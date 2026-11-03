/**
 * @file gpu_backend.cu
 * @brief The CUDA cCCE and gCCE paths for large clusters.
 */
#include "calculation/gpu/gpu_backend.h"
#include "calculation/gpu/compensated.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
using Z = cuDoubleComplex;
// an allocation that failed is GpuOutOfMemory (gpu_backend.h): not sticky, the context goes on
void check(cudaError_t s, const char* op){
    if (s == cudaErrorMemoryAllocation){ cudaGetLastError(); throw GpuOutOfMemory(std::string(op)+": "+cudaGetErrorString(s)); }
    if (s != cudaSuccess) throw std::runtime_error(std::string(op)+": "+cudaGetErrorString(s));
}
void check(cublasStatus_t s, const char* op){
    if (s == CUBLAS_STATUS_ALLOC_FAILED) throw GpuOutOfMemory(std::string(op)+": cuBLAS status "+std::to_string(s)+" (allocation failed)");
    if (s != CUBLAS_STATUS_SUCCESS) throw std::runtime_error(std::string(op)+": cuBLAS status "+std::to_string(s));
}
void check(cusolverStatus_t s, const char* op){
    if (s == CUSOLVER_STATUS_ALLOC_FAILED) throw GpuOutOfMemory(std::string(op)+": cuSOLVER status "+std::to_string(s)+" (allocation failed)");
    if (s != CUSOLVER_STATUS_SUCCESS) throw std::runtime_error(std::string(op)+": cuSOLVER status "+std::to_string(s));
}
// CCEX_TEST_GPU_OOM_BYTES (gpu_backend.h): the size from which every allocation fails as if the card were full
size_t testOutOfMemoryBytes(){
    static const size_t bytes = std::getenv("CCEX_TEST_GPU_OOM_BYTES") ?
        std::strtoull(std::getenv("CCEX_TEST_GPU_OOM_BYTES"),nullptr,10) : std::numeric_limits<size_t>::max();
    return bytes;
}
#define GPU_CHECK(op) check((op),#op)

void* allocateDevice(size_t bytes);
void releaseDevice(void* pointer) noexcept;

template<class T> struct Buffer {
    T* p=nullptr;
    size_t capacity=0;
    Buffer()=default;
    Buffer(const Buffer&)=delete;
    Buffer& operator=(const Buffer&)=delete;
    ~Buffer(){ if (p) releaseDevice(p); }
    void allocate(size_t n){
        if (p) throw std::runtime_error("GPU buffer already allocated");
        if (n > std::numeric_limits<size_t>::max()/sizeof(T))
            throw std::runtime_error("GPU allocation size overflow");
        if (n) p=static_cast<T*>(allocateDevice(n*sizeof(T)));
        capacity=n;
    }
    void reserve(size_t n){
        if(n<=capacity)return;
        if(n>std::numeric_limits<size_t>::max()/sizeof(T))throw std::runtime_error("GPU allocation size overflow");
        T* next=static_cast<T*>(allocateDevice(n*sizeof(T)));
        if(p)releaseDevice(p);p=next;capacity=n;
    }
};

struct Context {
    cudaStream_t stream=nullptr;
    cublasHandle_t blas=nullptr;
    cusolverDnHandle_t solver=nullptr;
    size_t cacheLimit=size_t(ccex::gpuCacheMiB)*1024*1024;
    bool recovering=false;              // gpuOutOfMemoryRecovery: released buffers are freed, not cached
    GpuMemoryStats memory{};
    std::unordered_map<void*,size_t> sizes;
    std::vector<void*> cached;
    void trim() noexcept {
        for(void* p:cached){ cudaFree(p); sizes.erase(p); }
        cached.clear(); memory.cachedBytes=0;
    }
    ~Context(){
        if (stream) cudaStreamSynchronize(stream);
        trim();
        if (solver) cusolverDnDestroy(solver);
        if (blas) cublasDestroy(blas);
        if (stream) cudaStreamDestroy(stream);
    }
};
std::unique_ptr<Context> context;

void* allocateDevice(size_t bytes){
    if (!context) throw std::runtime_error("CUDA backend was not initialized");
    auto& c=*context;
    size_t best=c.cached.size(), bestSize=std::numeric_limits<size_t>::max();
    for(size_t i=0;i<c.cached.size();++i){
        size_t size=c.sizes.at(c.cached[i]);
        // Do not pin a large solver workspace as a tiny long-lived state buffer.
        if(size>=bytes && size-bytes<=bytes && size<bestSize){best=i;bestSize=size;}
    }
    if(best<c.cached.size()){
        void* p=c.cached[best]; c.cached.erase(c.cached.begin()+best);
        c.memory.cachedBytes-=bestSize; ++c.memory.reuses; return p;
    }
    void* p=nullptr;
    cudaError_t status=bytes>=testOutOfMemoryBytes() ? cudaErrorMemoryAllocation : cudaMalloc(&p,bytes);
    if(status==cudaErrorMemoryAllocation && !c.cached.empty()){
        // Cached buffers must not prevent a larger, otherwise feasible request.
        cudaGetLastError(); c.trim(); status=bytes>=testOutOfMemoryBytes() ? cudaErrorMemoryAllocation : cudaMalloc(&p,bytes);
    }
    GPU_CHECK(status);
    try { c.sizes.emplace(p,bytes); } catch(...){cudaFree(p);throw;}
    ++c.memory.allocations;
    return p;
}

void releaseDevice(void* p) noexcept {
    if(!context){cudaFree(p);return;}
    auto& c=*context;
    auto found=c.sizes.find(p);
    if(found==c.sizes.end()){cudaFree(p);return;}
    const size_t bytes=found->second;
    // All users of these buffers run on this context's single stream. Reusing
    // an address is safe without synchronizing: subsequent accesses are ordered
    // after the previous owner's kernels on that same stream.
    if(!c.recovering && bytes<=c.cacheLimit && c.memory.cachedBytes<=c.cacheLimit-bytes){
        try {
            c.cached.push_back(p); c.memory.cachedBytes+=bytes;
            c.memory.peakCachedBytes=std::max(c.memory.peakCachedBytes,c.memory.cachedBytes);
            return;
        } catch(...) {}
    }
    cudaFree(p); c.sizes.erase(found);
}
// Owned by the first work object. Compare addresses on each call so changed
// membership and interleaved single-work evolution cannot leave stale pointers.
template<class T> struct PointerTable {
    std::vector<T> host;
    Buffer<T> device;
    void update(std::vector<T> next){
        if(next==host)return;
        device.reserve(next.size());
        GPU_CHECK(cudaMemcpyAsync(device.p,next.data(),next.size()*sizeof(T),cudaMemcpyHostToDevice,context->stream));
        host=std::move(next);++context->memory.pointerUploads;
    }
};
struct BatchWorkspace {
    PointerTable<Z*> pointers,density;
    PointerTable<const double*> eigenvalues;
    Buffer<Z> output;
};
const Z one=make_cuDoubleComplex(1,0), zero=make_cuDoubleComplex(0,0);

void gemm(int d, int cols, const Z* a, const Z* b, Z* out,
          cublasOperation_t transA=CUBLAS_OP_N, cublasOperation_t transB=CUBLAS_OP_N){
    GPU_CHECK(cublasZgemm(context->blas,transA,transB,d,cols,d,&one,a,d,b,d,&zero,out,d));
}

__global__ void resetStates(Z* out, const Z* initial, int d, int width, int total, bool vectorMode){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        int row=i%d, col=(i/d)%width;
        out[i] = vectorMode ? initial[row] : make_cuDoubleComplex(row==col ? 1.0 : 0.0,0);
    }
}

__device__ Z spectralPhase(double eigenvalue,double tau,double correction,bool refined){
    double sine,cosine;const double angle=eigenvalue*tau;
    sincos(angle,&sine,&cosine);Z phase=make_cuDoubleComplex(cosine,-sine);
    if(refined){
        double smallSine,smallCosine;
        sincos(__fma_rn(eigenvalue,tau,-angle)+correction*tau,&smallSine,&smallCosine);
        phase=cuCmul(phase,make_cuDoubleComplex(smallCosine,-smallSine));
    }
    return phase;
}

__global__ void phaseScale(Z* states, const double* evals, const double* times,
                           double fraction, int d, int width, int total, const double* corrections=nullptr){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        const double tau=times[i/(d*width)]*fraction;
        const Z phase=spectralPhase(evals[i%d],tau,corrections ? corrections[i%d] : 0,corrections!=nullptr);
        states[i]=cuCmul(states[i],phase);
    }
}

// Upulse acts on the qubit only. Apply it directly instead of allocating kron(P,I).
__global__ void applyPulse(const Z* in, Z* out, const Z* pulse,
                           int d, int q, int b, int total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        int row=i%d, bath=row%b, qubit=row/b;
        int base=(i/d)*d+bath;
        Z sum=make_cuDoubleComplex(0,0);
        for (int k=0; k<q; ++k) sum=cuCadd(sum,cuCmul(pulse[qubit+k*q],in[base+k*b]));
        out[i]=sum;
    }
}

__global__ void reduceState(const Z* states, Z* reduced, int d, int q, int b, int total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        int t=i/(q*q), row=i%q, col=(i/q)%q;
        const Z* psi=states+(size_t)t*d;
        Z sum=make_cuDoubleComplex(0,0);
        for (int k=0; k<b; ++k) sum=cuCadd(sum,cuCmul(psi[row*b+k],cuConj(psi[col*b+k])));
        reduced[i]=sum;
    }
}

__global__ void reduceDensity(const Z* density, Z* reduced, int d, int q, int b, int total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        int t=i/(q*q), row=i%q, col=(i/q)%q;
        const Z* rho=density+(size_t)t*d*d;
        Z sum=make_cuDoubleComplex(0,0);
        for (int k=0; k<b; ++k) sum=cuCadd(sum,rho[row*b+k+(size_t)(col*b+k)*d]);
        reduced[i]=sum;
    }
}

__global__ void fullDensity(const Z* states, Z* density, int d, size_t total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<total; i+=(size_t)blockDim.x*gridDim.x){
        size_t t=i/((size_t)d*d);
        int row=i%d, col=(i/d)%d;
        density[i]=cuCmul(states[t*d+row],cuConj(states[t*d+col]));
    }
}

int blocks(size_t n){ return (int)std::min<size_t>(4096,(n+255)/256); }

// Form the small eigensolver residual before projection, retaining cancellation
// beside the central-spin energy. One block owns each eigenvector.
__global__ void eigenResidual(const Z* h,const Z* vectors,const double* values,Z* residual,int d){
    const int col=blockIdx.x;
    for(int row=threadIdx.x;row<d;row+=blockDim.x){
        DoublePair real{0,0},imag{0,0};
        for(int k=0;k<d;++k){
            const Z v=vectors[k+(size_t)col*d];
            if(k==row){
                const auto diagonal=pairDifference({h[row+(size_t)row*d].x,0},{values[col],0});
                real=pairAdd(real,pairScale(diagonal,v.x));imag=pairAdd(imag,pairScale(diagonal,v.y));
            }else{
                const Z a=row>k ? h[row+(size_t)k*d] : cuConj(h[k+(size_t)row*d]);
                if(a.x==0 && a.y==0)continue;
                real=pairAdd(real,pairAdd(pairProduct(a.x,v.x),pairProduct(-a.y,v.y)));
                imag=pairAdd(imag,pairAdd(pairProduct(a.x,v.y),pairProduct(a.y,v.x)));
            }
        }
        residual[row+(size_t)col*d]=make_cuDoubleComplex(real.high+real.low,imag.high+imag.low);
    }
}

// The Jacobi refinement of eigen_refinement.h in a parallel (round-robin) order:
// a sweep is d-1 steps (one more for odd d), each rotating d/2 disjoint pairs at
// once, so E <- J^H E J is the column update followed by the row update. Every
// pair uses the CPU formula, skip rule and two-double energies; only the order
// of the rotations differs.
__device__ void atomicMaxNonNegative(unsigned long long* target,double value){
    atomicMax(target,(unsigned long long)__double_as_longlong(value));
}
__global__ void refineDiagonal(Z* e,double* values,double* corrections,int d,unsigned long long* largest){
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<d;i+=blockDim.x*gridDim.x){
        const auto energy=pairAdd({values[i],0},{e[i+(size_t)i*d].x,0});
        values[i]=energy.high;corrections[i]=energy.low;e[i+(size_t)i*d]=make_cuDoubleComplex(0,0);
        atomicMaxNonNegative(largest,fabs(energy.high));
    }
}
__global__ void refineSymmetrize(Z* e,int d){
    for(size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x;index<(size_t)d*d;index+=(size_t)blockDim.x*gridDim.x){
        const int i=index%d,j=index/d;
        if(i<=j)continue;
        const Z entry=cuCmul(make_cuDoubleComplex(.5,0),cuCadd(e[i+(size_t)j*d],cuConj(e[j+(size_t)i*d])));
        e[i+(size_t)j*d]=entry;e[j+(size_t)i*d]=cuConj(entry);
    }
}
// Pair i of step r: the circle method over m (= d rounded up to even) indices; an index d is the idle one.
__device__ bool jacobiPair(int i,int r,int m,int d,int& p,int& q){
    int a,b;
    if(i==0){a=r;b=m-1;}else{a=(r+i)%(m-1);b=(r-i+(m-1))%(m-1);}
    p=min(a,b);q=max(a,b);
    return q<d;
}
__global__ void jacobiRotations(const Z* e,double* values,double* corrections,int d,int m,int r,double tolerance,
        int* pairs,double* cosines,Z* sines){
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<m/2;i+=blockDim.x*gridDim.x){
        int p,q;pairs[2*i]=-1;
        if(!jacobiPair(i,r,m,d,p,q))continue;
        const Z off=e[p+(size_t)q*d];const double magnitude=cuCabs(off);
        if(magnitude<=tolerance)continue;
        const auto delta=pairDifference({values[q],corrections[q]},{values[p],corrections[p]});
        const double gap=delta.high+delta.low;
        const double t=gap==0 ? 1 : copysign(2*magnitude/(fabs(gap)+hypot(gap,2*magnitude)),gap);
        const double cosine=1/sqrt(1+t*t);
        const auto shift=pairProduct(t,magnitude);
        const auto a=pairDifference({values[p],corrections[p]},shift),b=pairAdd({values[q],corrections[q]},shift);
        values[p]=a.high;corrections[p]=a.low;values[q]=b.high;corrections[q]=b.low;
        pairs[2*i]=p;pairs[2*i+1]=q;cosines[i]=cosine;sines[i]=cuCmul(off,make_cuDoubleComplex(t*cosine/magnitude,0));
    }
}
// E <- E J and V <- V J: thread (k, pair) owns rows k of columns p and q.
__global__ void jacobiColumns(Z* e,Z* vectors,int d,int m,const int* pairs,const double* cosines,const Z* sines){
    const size_t total=(size_t)d*(m/2);
    for(size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x;index<total;index+=(size_t)blockDim.x*gridDim.x){
        const int k=index%d,i=index/d,p=pairs[2*i];
        if(p<0)continue;
        const int q=pairs[2*i+1];const Z c=make_cuDoubleComplex(cosines[i],0),s=sines[i];
        if(k!=p && k!=q){
            const Z a=e[k+(size_t)p*d],b=e[k+(size_t)q*d];
            e[k+(size_t)p*d]=cuCsub(cuCmul(c,a),cuCmul(cuConj(s),b));
            e[k+(size_t)q*d]=cuCadd(cuCmul(s,a),cuCmul(c,b));
        }
        const Z a=vectors[k+(size_t)p*d],b=vectors[k+(size_t)q*d];
        vectors[k+(size_t)p*d]=cuCsub(cuCmul(c,a),cuCmul(cuConj(s),b));
        vectors[k+(size_t)q*d]=cuCadd(cuCmul(s,a),cuCmul(c,b));
    }
}
// E <- J^H E: thread (k, pair) owns columns k of rows p and q; the pair's own 2x2 block is cleared
// (its energies were moved into values, as on the CPU).
__global__ void jacobiRows(Z* e,int d,int m,const int* pairs,const double* cosines,const Z* sines){
    const size_t total=(size_t)d*(m/2);
    for(size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x;index<total;index+=(size_t)blockDim.x*gridDim.x){
        const int k=index%d,i=index/d,p=pairs[2*i];
        if(p<0)continue;
        const int q=pairs[2*i+1];
        if(k==p || k==q){e[p+(size_t)k*d]=e[q+(size_t)k*d]=make_cuDoubleComplex(0,0);continue;}
        const Z c=make_cuDoubleComplex(cosines[i],0),s=sines[i];
        const Z a=e[p+(size_t)k*d],b=e[q+(size_t)k*d];
        e[p+(size_t)k*d]=cuCsub(cuCmul(c,a),cuCmul(s,b));
        e[q+(size_t)k*d]=cuCadd(cuCmul(cuConj(s),a),cuCmul(c,b));
    }
}
__global__ void largestEntry(const Z* e,size_t total,unsigned long long* largest){
    for(size_t index=(size_t)blockIdx.x*blockDim.x+threadIdx.x;index<total;index+=(size_t)blockDim.x*gridDim.x)
        atomicMaxNonNegative(largest,cuCabs(e[index]));
}

double readLargest(Buffer<unsigned long long>& largest){
    unsigned long long bits=0;
    GPU_CHECK(cudaMemcpyAsync(&bits,largest.p,sizeof bits,cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
    double value;std::memcpy(&value,&bits,sizeof value);return value;
}

// Returns false when maxSweeps sweeps end with an off-diagonal entry above the tolerance.
bool refineEigenbasis(Z* e,Z* vectors,double* values,double* corrections,int d,int maxSweeps){
    Buffer<unsigned long long> largest;largest.allocate(1);
    GPU_CHECK(cudaMemsetAsync(largest.p,0,sizeof(unsigned long long),context->stream));
    refineDiagonal<<<blocks(d),256,0,context->stream>>>(e,values,corrections,d,largest.p);
    GPU_CHECK(cudaGetLastError());
    refineSymmetrize<<<blocks((size_t)d*d),256,0,context->stream>>>(e,d);
    GPU_CHECK(cudaGetLastError());
    const double eps=std::numeric_limits<double>::epsilon();
    const double tolerance=64*eps*eps*readLargest(largest);
    const int m=d+(d&1),half=m/2;
    Buffer<int> pairs;pairs.allocate(2*(size_t)half);
    Buffer<double> cosines;cosines.allocate(half);
    Buffer<Z> sines;sines.allocate(half);
    const size_t work=(size_t)d*half;
    for(int sweep=0;sweep<maxSweeps;++sweep){
        for(int r=0;r<m-1;++r){
            jacobiRotations<<<blocks(half),256,0,context->stream>>>(e,values,corrections,d,m,r,tolerance,pairs.p,cosines.p,sines.p);
            GPU_CHECK(cudaGetLastError());
            jacobiColumns<<<blocks(work),256,0,context->stream>>>(e,vectors,d,m,pairs.p,cosines.p,sines.p);
            GPU_CHECK(cudaGetLastError());
            jacobiRows<<<blocks(work),256,0,context->stream>>>(e,d,m,pairs.p,cosines.p,sines.p);
            GPU_CHECK(cudaGetLastError());
        }
        GPU_CHECK(cudaMemsetAsync(largest.p,0,sizeof(unsigned long long),context->stream));
        largestEntry<<<blocks((size_t)d*d),256,0,context->stream>>>(e,(size_t)d*d,largest.p);
        GPU_CHECK(cudaGetLastError());
        if(readLargest(largest)<=tolerance)return true;
    }
    return false;
}

void diagonalize(int d, const void* hamiltonian, Buffer<Z>& vectors, Buffer<double>& values){
    vectors.allocate((size_t)d*d);
    values.allocate(d);
    GPU_CHECK(cudaMemcpyAsync(vectors.p,hamiltonian,(size_t)d*d*sizeof(Z),cudaMemcpyHostToDevice,context->stream));
    int workspaceSize=0;
    GPU_CHECK(cusolverDnZheevd_bufferSize(context->solver,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,
               d,vectors.p,d,values.p,&workspaceSize));
    Buffer<Z> workspace; workspace.allocate(workspaceSize);
    Buffer<int> info; info.allocate(1);
    GPU_CHECK(cusolverDnZheevd(context->solver,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,
               d,vectors.p,d,values.p,workspace.p,workspaceSize,info.p));
    int hostInfo=0;
    GPU_CHECK(cudaMemcpyAsync(&hostInfo,info.p,sizeof(int),cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
    if (hostInfo) throw std::runtime_error("Hermitian eigensolver failed: info="+std::to_string(hostInfo));
}

__global__ void resetCce(Z* a, Z* b, const Z* initialA, const Z* initialB,
                          const Z* overlap, int d, int width, int total, bool ensemble){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<(size_t)total; i+=(size_t)blockDim.x*gridDim.x){
        int row=i%d, col=(i/d)%width;
        // Both ensembles start from the same basis V_alpha: I in the alpha
        // eigenbasis, M^H in the beta eigenbasis, with M=V_alpha^H V_beta.
        a[i]=ensemble ? make_cuDoubleComplex(row==col ? 1.0 : 0.0,0) : initialA[row];
        b[i]=ensemble ? cuConj(overlap[col+(size_t)row*d]) : initialB[row];
    }
}

__global__ void overlapCce(const Z* a, const Z* b, Z* out, int size, int d,
                            bool ensemble, bool conjugate){
    __shared__ Z partial[256];
    const size_t base=(size_t)blockIdx.x*size;
    Z sum=make_cuDoubleComplex(0,0);
    for (size_t i=threadIdx.x; i<(size_t)size; i+=blockDim.x)
        sum=cuCadd(sum,cuCmul(cuConj(b[base+i]),a[base+i]));
    partial[threadIdx.x]=sum;
    __syncthreads();
    for (int stride=blockDim.x/2; stride; stride/=2){
        if (threadIdx.x<stride) partial[threadIdx.x]=cuCadd(partial[threadIdx.x],partial[threadIdx.x+stride]);
        __syncthreads();
    }
    if (!threadIdx.x){
        sum=partial[0];
        if (ensemble){ sum.x/=d; sum.y/=d; }
        if (conjugate) sum.y=-sum.y;
        out[blockIdx.x]=sum;
    }
}
}

namespace {
// Hahn echo, ensemble: L = (1/d) sum_ij |X_ij|^2 Uu_j conj(Uu_i) with X = M Dd M^H,
// M = V_alpha^H V_beta (the CPU kernel's formula; one d x d x d product per step).
// hahnScale writes Y_t = M diag(Dd(t)) for every step t of the batch.
__global__ void hahnScale(Z* y, const Z* overlap, const double* evalB, const double* times,
                          int d, size_t total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<total; i+=(size_t)blockDim.x*gridDim.x){
        const size_t t=i/((size_t)d*d); const int col=(i/d)%d, row=i%d;
        y[i]=cuCmul(overlap[row+(size_t)col*d],spectralPhase(evalB[col],times[t],0,false));
    }
}
__global__ void hahnPhases(Z* phases, const double* evalA, const double* times, int d, size_t total){
    for (size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x; i<total; i+=(size_t)blockDim.x*gridDim.x)
        phases[i]=spectralPhase(evalA[i%d],times[i/d],0,false);
}
// One block per step: sum over the d x d entries of X_t.
__global__ void hahnTrace(const Z* x, const Z* phases, Z* out, int d){
    __shared__ Z partial[256];
    const size_t base=(size_t)blockIdx.x*d*d;
    const Z* u=phases+(size_t)blockIdx.x*d;
    Z sum=make_cuDoubleComplex(0,0);
    for (size_t k=threadIdx.x; k<(size_t)d*d; k+=blockDim.x){
        const int i=k%d, j=k/d;
        const Z v=x[base+k];
        const double w=v.x*v.x+v.y*v.y;
        const Z p=cuCmul(u[j],cuConj(u[i]));
        sum=cuCadd(sum,make_cuDoubleComplex(w*p.x,w*p.y));
    }
    partial[threadIdx.x]=sum;
    __syncthreads();
    for (int stride=blockDim.x/2; stride; stride/=2){
        if (threadIdx.x<stride) partial[threadIdx.x]=cuCadd(partial[threadIdx.x],partial[threadIdx.x+stride]);
        __syncthreads();
    }
    if (!threadIdx.x){ sum=partial[0]; sum.x/=d; sum.y/=d; out[blockIdx.x]=sum; }
}
}

struct GpuGcce {
    int d,q,b,width,batch,npulse;
    bool vectorMode;
    Buffer<Z> eigenvectors, initial, pulses, states, scratch, densityScratch, reduced, full;
    Buffer<double> eigenvalues, eigenCorrections, times;
    BatchWorkspace batchWork;
};

bool gpuInitialize(int localRank, int localSize, char* description, int capacity){
    // Honor an explicit request to hide every GPU before entering the driver.
    // CUDA 12.8 / WSL driver 577.00 aborts inside cudaGetDeviceCount for an empty
    // CUDA_VISIBLE_DEVICES (also reproduced in a standalone cudart-only program).
    const char* visible=std::getenv("CUDA_VISIBLE_DEVICES");
    if (visible && (!*visible || !std::strcmp(visible,"-1"))){
        snprintf(description,capacity,"CPU (CUDA initialization failed: CUDA_VISIBLE_DEVICES hides all GPUs)");
        return false;
    }
    try {
        int count=0;
        GPU_CHECK(cudaGetDeviceCount(&count));
        if (!count) throw std::runtime_error("no visible CUDA devices");
        // Scheduler-isolated ranks commonly see one device each. Otherwise map the
        // node-local rank to the visible-device ordinal, never the global MPI rank.
        int device=localRank%count;
        GPU_CHECK(cudaSetDevice(device));
        cudaDeviceProp properties;
        GPU_CHECK(cudaGetDeviceProperties(&properties,device));
        auto next=std::make_unique<Context>();
        GPU_CHECK(cudaStreamCreateWithFlags(&next->stream,cudaStreamNonBlocking));
        GPU_CHECK(cublasCreate(&next->blas));
        GPU_CHECK(cublasSetStream(next->blas,next->stream));
        GPU_CHECK(cusolverDnCreate(&next->solver));
        GPU_CHECK(cusolverDnSetStream(next->solver,next->stream));
        context=std::move(next);
        snprintf(description,capacity,"CUDA device=%d name=%s local_rank=%d local_size=%d visible_devices=%d",
                 device,properties.name,localRank,localSize,count);
        return true;
    } catch (const std::exception& e){
        context.reset();
        snprintf(description,capacity,"CPU (CUDA initialization failed: %s)",e.what());
        return false;
    }
}

void gpuFinalize(){ context.reset(); }

// A68: CCEX_CUDA_ARCH is the -arch both builds pass (CMake: its CUDA architectures, "native" there leaves the list nvcc
// compiled for, __CUDA_ARCH_LIST__, e.g. 860)
#define CCEX_STRINGIFY_(...) #__VA_ARGS__
#define CCEX_STRINGIFY(...) CCEX_STRINGIFY_(__VA_ARGS__)
void gpuBuildDescription(char* out, int capacity){
#ifdef CCEX_CUDA_ARCH
    std::string arch = CCEX_STRINGIFY(CCEX_CUDA_ARCH);
#else
    std::string arch = "native";
#endif
#ifdef __CUDA_ARCH_LIST__
    const int compiled[] = {__CUDA_ARCH_LIST__};
    if (arch == "native" || arch.empty()) arch = std::to_string(compiled[0]/10);
#endif
    snprintf(out, capacity, "CUDA %d.%d, sm_%s", CUDART_VERSION/1000, CUDART_VERSION%1000/10, arch.c_str());
}
int gpuVisibleDevices(){
    const char* visible=std::getenv("CUDA_VISIBLE_DEVICES");   // the guard of gpuInitialize
    if (visible && (!*visible || !std::strcmp(visible,"-1"))) return 0;
    int count=0;
    return cudaGetDeviceCount(&count)==cudaSuccess ? count : 0;
}
std::size_t gpuDeviceMemorySmallest(){
    // the same guard as gpuInitialize: an empty CUDA_VISIBLE_DEVICES aborts inside cudaGetDeviceCount (review 5f7834 R4)
    const char* visible=std::getenv("CUDA_VISIBLE_DEVICES");
    if (visible && (!*visible || !std::strcmp(visible,"-1"))) return 0;
    int count=0;
    if (cudaGetDeviceCount(&count)!=cudaSuccess) return 0;
    std::size_t smallest=0;
    for (int d=0; d<count; ++d){
        cudaDeviceProp prop;
        if (cudaGetDeviceProperties(&prop,d)==cudaSuccess && (!smallest || prop.totalGlobalMem<smallest))
            smallest=prop.totalGlobalMem;
    }
    return smallest;
}
void gpuSetCacheLimit(size_t bytes){
    if (!context) return;
    context->cacheLimit=bytes;
    if(context->memory.cachedBytes>bytes) context->trim();
}
void gpuOutOfMemoryRecovery(bool on){
    if (!context) return;
    context->recovering=on;
    if (on) context->trim();
}
GpuMemoryStats gpuMemoryStats(){return context ? context->memory : GpuMemoryStats{};}

GpuGcce* gpuGcceCreate(int dimension, int qubitDimension, bool vectorMode,
                      int batchSteps, const void* hamiltonian, const void* initial,
                      int pulseCount, const void* pulses, bool refine, int maxSweeps){
    if (!context) throw std::runtime_error("CUDA backend was not initialized");
    if (dimension<1 || qubitDimension<1 || dimension%qubitDimension || batchSteps<1 || pulseCount<0 || maxSweeps<0)
        throw std::runtime_error("invalid gCCE GPU dimensions");
    auto w=std::make_unique<GpuGcce>();
    w->d=dimension; w->q=qubitDimension; w->b=dimension/qubitDimension;
    w->vectorMode=vectorMode; w->width=vectorMode ? 1 : dimension;
    w->batch=batchSteps; w->npulse=pulseCount;
    const size_t square=(size_t)dimension*dimension;
    const size_t columns=(size_t)w->width*batchSteps;
    if (square>INT_MAX || columns>INT_MAX/(size_t)dimension ||
        (size_t)qubitDimension*qubitDimension*batchSteps>INT_MAX)
        throw std::runtime_error("GPU batch exceeds 32-bit indexing, reduce gpu_batch_steps or cluster dimension");
    diagonalize(dimension,hamiltonian,w->eigenvectors,w->eigenvalues);
    w->eigenCorrections.allocate(dimension);
    GPU_CHECK(cudaMemsetAsync(w->eigenCorrections.p,0,size_t(dimension)*sizeof(double),context->stream));
    if(refine){
        Buffer<Z> original,residual,projected;
        original.allocate(square);residual.allocate(square);projected.allocate(square);
        GPU_CHECK(cudaMemcpyAsync(original.p,hamiltonian,square*sizeof(Z),cudaMemcpyHostToDevice,context->stream));
        eigenResidual<<<dimension,256,0,context->stream>>>(original.p,w->eigenvectors.p,w->eigenvalues.p,residual.p,dimension);
        GPU_CHECK(cudaGetLastError());
        gemm(dimension,dimension,w->eigenvectors.p,residual.p,projected.p,CUBLAS_OP_C);
        if(!refineEigenbasis(projected.p,w->eigenvectors.p,w->eigenvalues.p,w->eigenCorrections.p,dimension,maxSweeps))throw std::runtime_error("gCCE eigensystem refinement did not converge: dim="+
            std::to_string(dimension)+" gcce_refinement_max_sweeps="+std::to_string(maxSweeps));
    }

    size_t initialSize=vectorMode ? dimension : square;
    w->initial.allocate(initialSize);
    GPU_CHECK(cudaMemcpyAsync(w->initial.p,initial,initialSize*sizeof(Z),cudaMemcpyHostToDevice,context->stream));
    size_t pulseSize=(size_t)pulseCount*qubitDimension*qubitDimension;
    w->pulses.allocate(pulseSize);
    if (pulseSize) GPU_CHECK(cudaMemcpyAsync(w->pulses.p,pulses,pulseSize*sizeof(Z),cudaMemcpyHostToDevice,context->stream));
    w->states.allocate(dimension*columns);
    w->scratch.allocate(dimension*columns);
    if (!vectorMode) w->densityScratch.allocate(square);
    w->reduced.allocate((size_t)qubitDimension*qubitDimension*batchSteps);
    w->times.allocate(batchSteps);
    GPU_CHECK(cudaStreamSynchronize(context->stream));
    return w.release();
}

void gpuGcceEvolve(GpuGcce* w, int count, const double* times,
                   const double* delayFractions, void* reduced, void* full){
    if (!w || count<1 || count>w->batch) throw std::runtime_error("invalid GPU time batch");
    int cols=w->width*count, total=w->d*cols;
    GPU_CHECK(cudaMemcpyAsync(w->times.p,times,count*sizeof(double),cudaMemcpyHostToDevice,context->stream));
    resetStates<<<blocks(total),256,0,context->stream>>>(w->states.p,w->initial.p,w->d,w->width,total,w->vectorMode);
    GPU_CHECK(cudaGetLastError());
    for (int p=0; p<=w->npulse; ++p){
        // All time points share the eigenvectors: GEMM processes their states as
        // columns. Each column group receives its own phase on the device.
        gemm(w->d,cols,w->eigenvectors.p,w->states.p,w->scratch.p,CUBLAS_OP_C);
        phaseScale<<<blocks(total),256,0,context->stream>>>(w->scratch.p,w->eigenvalues.p,w->times.p,
                                                          delayFractions[p],w->d,w->width,total,w->eigenCorrections.p);
        GPU_CHECK(cudaGetLastError());
        gemm(w->d,cols,w->eigenvectors.p,w->scratch.p,w->states.p);
        if (p<w->npulse){
            applyPulse<<<blocks(total),256,0,context->stream>>>(w->states.p,w->scratch.p,
                w->pulses.p+(size_t)p*w->q*w->q,w->d,w->q,w->b,total);
            GPU_CHECK(cudaGetLastError());
            std::swap(w->states.p,w->scratch.p);
        }
    }
    int reducedCount=w->q*w->q*count;
    size_t fullCount=(size_t)count*w->d*w->d;
    if (w->vectorMode){
        reduceState<<<blocks(reducedCount),256,0,context->stream>>>(w->states.p,w->reduced.p,w->d,w->q,w->b,reducedCount);
        GPU_CHECK(cudaGetLastError());
        if (full){
            if (!w->full.p) w->full.allocate((size_t)w->batch*w->d*w->d);
            fullDensity<<<blocks(fullCount),256,0,context->stream>>>(w->states.p,w->full.p,w->d,fullCount);
            GPU_CHECK(cudaGetLastError());
            GPU_CHECK(cudaMemcpyAsync(full,w->full.p,fullCount*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
        }
    }else{
        for (int t=0; t<count; ++t){
            size_t offset=(size_t)t*w->d*w->d;
            gemm(w->d,w->d,w->states.p+offset,w->initial.p,w->densityScratch.p);
            gemm(w->d,w->d,w->densityScratch.p,w->states.p+offset,w->scratch.p+offset,CUBLAS_OP_N,CUBLAS_OP_C);
        }
        reduceDensity<<<blocks(reducedCount),256,0,context->stream>>>(w->scratch.p,w->reduced.p,w->d,w->q,w->b,reducedCount);
        GPU_CHECK(cudaGetLastError());
        if (full) GPU_CHECK(cudaMemcpyAsync(full,w->scratch.p,fullCount*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
    }
    GPU_CHECK(cudaMemcpyAsync(reduced,w->reduced.p,(size_t)reducedCount*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
}

void gpuGcceDestroy(GpuGcce* work){ delete work; }

struct GpuCce {
    int d,width,batch,npulse;
    bool ensemble;
    Buffer<Z> overlap, initialA, initialB, a, b, scratch, result, phases;
    Buffer<double> evalA, evalB, durations;
    BatchWorkspace batchWork;
};

GpuCce* gpuCceCreate(int d, bool ensemble, int batchSteps, int pulseCount,
                    const void* hAlpha, const void* hBeta, const void* initial){
    if (!context) throw std::runtime_error("CUDA backend was not initialized");
    if (d<1 || batchSteps<1 || pulseCount<0 || pulseCount==INT_MAX || (!ensemble && !initial))
        throw std::runtime_error("invalid cCCE GPU dimensions or initial state");
    auto w=std::make_unique<GpuCce>();
    w->d=d; w->width=ensemble ? d : 1; w->batch=batchSteps; w->npulse=pulseCount; w->ensemble=ensemble;
    const size_t square=(size_t)d*d, columns=(size_t)w->width*batchSteps;
    if (square>INT_MAX || columns>INT_MAX/(size_t)d)
        throw std::runtime_error("GPU batch exceeds 32-bit indexing, reduce gpu_batch_steps or cluster dimension");
    Buffer<Z> va,vb;
    diagonalize(d,hAlpha,va,w->evalA);
    diagonalize(d,hBeta,vb,w->evalB);
    w->overlap.allocate(square);
    gemm(d,d,va.p,vb.p,w->overlap.p,CUBLAS_OP_C);
    if (!ensemble){
        Buffer<Z> psi; psi.allocate(d);
        w->initialA.allocate(d); w->initialB.allocate(d);
        GPU_CHECK(cudaMemcpyAsync(psi.p,initial,(size_t)d*sizeof(Z),cudaMemcpyHostToDevice,context->stream));
        gemm(d,1,va.p,psi.p,w->initialA.p,CUBLAS_OP_C);
        gemm(d,1,vb.p,psi.p,w->initialB.p,CUBLAS_OP_C);
        GPU_CHECK(cudaStreamSynchronize(context->stream));
    }
    w->a.allocate(d*columns); w->b.allocate(d*columns); w->scratch.allocate(d*columns);
    w->result.allocate(batchSteps);
    w->durations.allocate((size_t)batchSteps*(pulseCount+1));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
    return w.release();
}

void gpuCceEvolve(GpuCce* w, int count, const double* durations, void* coherence){
    if (!w || count<1 || count>w->batch) throw std::runtime_error("invalid cCCE GPU time batch");
    const int cols=w->width*count, total=w->d*cols;
    GPU_CHECK(cudaMemcpyAsync(w->durations.p,durations,(size_t)count*(w->npulse+1)*sizeof(double),
                             cudaMemcpyHostToDevice,context->stream));
    // Hahn echo on the ensemble with equal intervals: one product per step (hahn* kernels).
    bool hahn = w->ensemble && w->npulse==1;
    for (int t=0; hahn && t<count; ++t) hahn = durations[t]==durations[count+t];
    if (hahn){
        const int d=w->d;
        const size_t square=(size_t)d*d, all=square*count;
        if (!w->phases.p) w->phases.allocate((size_t)d*w->batch);
        hahnScale<<<blocks(all),256,0,context->stream>>>(w->a.p,w->overlap.p,w->evalB.p,w->durations.p,d,all);
        GPU_CHECK(cudaGetLastError());
        hahnPhases<<<blocks((size_t)d*count),256,0,context->stream>>>(w->phases.p,w->evalA.p,w->durations.p,d,(size_t)d*count);
        GPU_CHECK(cudaGetLastError());
        GPU_CHECK(cublasZgemmStridedBatched(context->blas,CUBLAS_OP_N,CUBLAS_OP_C,d,d,d,&one,
                  w->a.p,d,(long long)square,w->overlap.p,d,0,&zero,w->scratch.p,d,(long long)square,count));
        hahnTrace<<<count,256,0,context->stream>>>(w->scratch.p,w->phases.p,w->result.p,d);
        GPU_CHECK(cudaGetLastError());
        GPU_CHECK(cudaMemcpyAsync(coherence,w->result.p,(size_t)count*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
        GPU_CHECK(cudaStreamSynchronize(context->stream));
        return;
    }
    resetCce<<<blocks(total),256,0,context->stream>>>(w->a.p,w->b.p,w->initialA.p,w->initialB.p,
                                                    w->overlap.p,w->d,w->width,total,w->ensemble);
    GPU_CHECK(cudaGetLastError());
    for (int p=0; p<=w->npulse; ++p){
        bool odd=p%2;
        if (p){
            // A switches alpha->beta on odd intervals; B switches beta->alpha.
            gemm(w->d,cols,w->overlap.p,w->a.p,w->scratch.p,odd ? CUBLAS_OP_C : CUBLAS_OP_N);
            std::swap(w->a.p,w->scratch.p);
            gemm(w->d,cols,w->overlap.p,w->b.p,w->scratch.p,odd ? CUBLAS_OP_N : CUBLAS_OP_C);
            std::swap(w->b.p,w->scratch.p);
        }
        const double* interval=w->durations.p+(size_t)p*count;
        phaseScale<<<blocks(total),256,0,context->stream>>>(w->a.p,odd ? w->evalB.p : w->evalA.p,
                                                          interval,1.0,w->d,w->width,total);
        GPU_CHECK(cudaGetLastError());
        phaseScale<<<blocks(total),256,0,context->stream>>>(w->b.p,odd ? w->evalA.p : w->evalB.p,
                                                          interval,1.0,w->d,w->width,total);
        GPU_CHECK(cudaGetLastError());
    }
    // Express A in B's final basis before the state/Frobenius inner product.
    gemm(w->d,cols,w->overlap.p,w->a.p,w->scratch.p,w->npulse%2 ? CUBLAS_OP_N : CUBLAS_OP_C);
    // No conjugate: the CPU ensemble used to return one for an even pulse count
    // and this matched it; both now return the mean of the single samples.
    overlapCce<<<count,256,0,context->stream>>>(w->scratch.p,w->b.p,w->result.p,w->d*w->width,
        w->d,w->ensemble,false);
    GPU_CHECK(cudaGetLastError());
    GPU_CHECK(cudaMemcpyAsync(coherence,w->result.p,(size_t)count*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
}

void gpuCceDestroy(GpuCce* work){ delete work; }

namespace {
__global__ void resetCceBatch(Z* const* a, Z* const* b, Z* const* initialA,
        Z* const* initialB, Z* const* overlap, int d, int width,
        size_t perCluster, size_t total, bool ensemble){
    for(size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;i<total;i+=(size_t)blockDim.x*gridDim.x){
        size_t c=i/perCluster, offset=i%perCluster;
        int row=offset%d,col=(offset/d)%width;
        a[c][offset]=ensemble ? make_cuDoubleComplex(row==col ? 1.0 : 0.0,0) : initialA[c][row];
        b[c][offset]=ensemble ? cuConj(overlap[c][col+(size_t)row*d]) : initialB[c][row];
    }
}
__global__ void phaseCceBatch(Z* const* states, const double* const* evals,
        const double* durations, int d, int width, size_t perCluster, size_t total,
        const double* const* corrections=nullptr, double fraction=1){
    for(size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;i<total;i+=(size_t)blockDim.x*gridDim.x){
        size_t c=i/perCluster, offset=i%perCluster;
        const double tau=durations[offset/(d*width)]*fraction;
        const Z phase=spectralPhase(evals[c][offset%d],tau,corrections ? corrections[c][offset%d] : 0,corrections!=nullptr);
        states[c][offset]=cuCmul(states[c][offset],phase);
    }
}
__global__ void overlapCceBatch(Z* const* a, Z* const* b, Z* out,
        int size,int d,int count,bool ensemble,bool conjugate){
    __shared__ Z partial[256];
    const int c=blockIdx.x/count,t=blockIdx.x%count;
    const size_t base=(size_t)t*size;
    Z sum=make_cuDoubleComplex(0,0);
    for(size_t i=threadIdx.x;i<(size_t)size;i+=blockDim.x)
        sum=cuCadd(sum,cuCmul(cuConj(b[c][base+i]),a[c][base+i]));
    partial[threadIdx.x]=sum; __syncthreads();
    for(int stride=blockDim.x/2;stride;stride/=2){
        if(threadIdx.x<stride) partial[threadIdx.x]=cuCadd(partial[threadIdx.x],partial[threadIdx.x+stride]);
        __syncthreads();
    }
    if(!threadIdx.x){
        sum=partial[0];
        if(ensemble){sum.x/=d;sum.y/=d;}
        if(conjugate)sum.y=-sum.y;
        out[blockIdx.x]=sum;
    }
}
}

void gpuCceEvolveBatch(GpuCce* const* works,int n,int count,const double* durations,void* coherence){
    if(n<1 || !works || !works[0] || count<1) throw std::runtime_error("invalid cCCE cluster batch");
    if(n==1){gpuCceEvolve(works[0],count,durations,coherence);return;}
    const GpuCce& first=*works[0];
    if((size_t)n*count>INT_MAX) throw std::runtime_error("cCCE cluster/time batch exceeds indexing limit");
    std::vector<Z*> hostPointers((size_t)6*n);
    std::vector<const double*> hostEvals((size_t)2*n);
    for(int c=0;c<n;++c){
        GpuCce* w=works[c];
        if(!w || w->d!=first.d || w->width!=first.width || w->ensemble!=first.ensemble ||
           w->npulse!=first.npulse || count>w->batch)
            throw std::runtime_error("cCCE batch requires matching dimensions, state mode and pulses");
        hostPointers[c]=w->a.p;hostPointers[n+c]=w->b.p;hostPointers[2*n+c]=w->scratch.p;
        hostPointers[3*n+c]=w->overlap.p;hostPointers[4*n+c]=w->initialA.p;hostPointers[5*n+c]=w->initialB.p;
        hostEvals[c]=w->evalA.p;hostEvals[n+c]=w->evalB.p;
    }
    auto& batchWork=works[0]->batchWork;
    batchWork.pointers.update(std::move(hostPointers));batchWork.eigenvalues.update(std::move(hostEvals));
    batchWork.output.reserve(size_t(n)*count);
    auto& pointers=batchWork.pointers.device;auto& evalPointers=batchWork.eigenvalues.device;
    auto& output=batchWork.output;
    GPU_CHECK(cudaMemcpyAsync(first.durations.p,durations,(size_t)count*(first.npulse+1)*sizeof(double),cudaMemcpyHostToDevice,context->stream));
    Z** a=pointers.p;Z** b=pointers.p+n;Z** scratch=pointers.p+2*n;Z** overlap=pointers.p+3*n;
    const int cols=first.width*count;
    const size_t perCluster=(size_t)first.d*cols,total=perCluster*n;
    resetCceBatch<<<blocks(total),256,0,context->stream>>>(a,b,pointers.p+4*n,pointers.p+5*n,overlap,
        first.d,first.width,perCluster,total,first.ensemble);
    GPU_CHECK(cudaGetLastError());
    auto multiply=[&](Z** in,Z** out,cublasOperation_t trans){
        GPU_CHECK(cublasZgemmBatched(context->blas,trans,CUBLAS_OP_N,first.d,cols,first.d,
            &one,overlap,first.d,in,first.d,&zero,out,first.d,n));
    };
    for(int p=0;p<=first.npulse;++p){
        const bool odd=p%2;
        if(p){
            multiply(a,scratch,odd ? CUBLAS_OP_C : CUBLAS_OP_N);std::swap(a,scratch);
            multiply(b,scratch,odd ? CUBLAS_OP_N : CUBLAS_OP_C);std::swap(b,scratch);
        }
        phaseCceBatch<<<blocks(total),256,0,context->stream>>>(a,evalPointers.p+(odd ? n : 0),
            first.durations.p+(size_t)p*count,first.d,first.width,perCluster,total);
        GPU_CHECK(cudaGetLastError());
        phaseCceBatch<<<blocks(total),256,0,context->stream>>>(b,evalPointers.p+(odd ? 0 : n),
            first.durations.p+(size_t)p*count,first.d,first.width,perCluster,total);
        GPU_CHECK(cudaGetLastError());
    }
    multiply(a,scratch,first.npulse%2 ? CUBLAS_OP_N : CUBLAS_OP_C);
    overlapCceBatch<<<n*count,256,0,context->stream>>>(scratch,b,output.p,first.d*first.width,
        first.d,count,first.ensemble,false);
    GPU_CHECK(cudaGetLastError());
    GPU_CHECK(cudaMemcpyAsync(coherence,output.p,(size_t)n*count*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
}

namespace {
__global__ void resetGcceBatch(Z* const* states,Z* const* initial,int d,int width,
        size_t perCluster,size_t total,bool vectorMode){
    for(size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;i<total;i+=(size_t)blockDim.x*gridDim.x){
        size_t c=i/perCluster,offset=i%perCluster;int row=offset%d,col=(offset/d)%width;
        states[c][offset]=vectorMode ? initial[c][row] : make_cuDoubleComplex(row==col ? 1.0 : 0.0,0);
    }
}
__global__ void pulseGcceBatch(Z* const* in,Z* const* out,Z* const* pulses,
        int p,int d,int q,int b,size_t perCluster,size_t total){
    for(size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;i<total;i+=(size_t)blockDim.x*gridDim.x){
        size_t c=i/perCluster,offset=i%perCluster;
        int row=offset%d,bath=row%b,qubit=row/b;
        size_t base=(offset/d)*d+bath;
        const Z* pulse=pulses[c]+(size_t)p*q*q;
        Z sum=make_cuDoubleComplex(0,0);
        for(int k=0;k<q;++k)sum=cuCadd(sum,cuCmul(pulse[qubit+k*q],in[c][base+k*b]));
        out[c][offset]=sum;
    }
}
__global__ void reduceGcceBatch(Z* const* states,Z* result,int d,int q,int b,
        int count,size_t total,bool vectorMode){
    const size_t perCluster=(size_t)count*q*q;
    for(size_t i=(size_t)blockIdx.x*blockDim.x+threadIdx.x;i<total;i+=(size_t)blockDim.x*gridDim.x){
        const size_t c=i/perCluster,offset=i%perCluster;
        const int t=offset/(q*q),row=offset%q,col=(offset/q)%q;
        const Z* state=states[c]+(size_t)t*d*(vectorMode ? 1 : d);
        Z sum=make_cuDoubleComplex(0,0);
        for(int k=0;k<b;++k)sum=cuCadd(sum,vectorMode ?
            cuCmul(state[row*b+k],cuConj(state[col*b+k])) : state[row*b+k+(size_t)(col*b+k)*d]);
        result[i]=sum;
    }
}
}

void gpuGcceEvolveBatch(GpuGcce* const* works,int n,int count,const double* times,
        const double* delays,void* reduced){
    if(n<1 || !works || !works[0] || count<1)throw std::runtime_error("invalid gCCE cluster batch");
    if(n==1){gpuGcceEvolve(works[0],count,times,delays,reduced,nullptr);return;}
    const GpuGcce& first=*works[0];
    const size_t resultCount=(size_t)n*count*first.q*first.q;
    if(resultCount>INT_MAX)throw std::runtime_error("gCCE cluster/time batch exceeds indexing limit");
    std::vector<Z*> hostPointers((size_t)6*n);
    std::vector<const double*> hostEvals((size_t)2*n);
    for(int c=0;c<n;++c){
        GpuGcce* w=works[c];
        if(!w || w->d!=first.d || w->q!=first.q || w->vectorMode!=first.vectorMode ||
           w->npulse!=first.npulse || count>w->batch)
            throw std::runtime_error("gCCE batch requires matching dimensions, state mode and pulses");
        if(!w->vectorMode && !w->full.p)w->full.allocate((size_t)w->batch*w->d*w->d);
        hostPointers[c]=w->states.p;hostPointers[n+c]=w->scratch.p;hostPointers[2*n+c]=w->eigenvectors.p;
        hostPointers[3*n+c]=w->initial.p;hostPointers[4*n+c]=w->pulses.p;hostPointers[5*n+c]=w->full.p;
        hostEvals[c]=w->eigenvalues.p;hostEvals[n+c]=w->eigenCorrections.p;
    }
    auto& batchWork=works[0]->batchWork;
    batchWork.pointers.update(std::move(hostPointers));batchWork.eigenvalues.update(std::move(hostEvals));
    batchWork.output.reserve(resultCount);
    auto& pointers=batchWork.pointers.device;auto& evals=batchWork.eigenvalues.device;
    auto& output=batchWork.output;
    GPU_CHECK(cudaMemcpyAsync(first.times.p,times,count*sizeof(double),cudaMemcpyHostToDevice,context->stream));
    Z** states=pointers.p;Z** scratch=pointers.p+n;Z** vectors=pointers.p+2*n;
    const int cols=first.width*count;
    const size_t perCluster=(size_t)first.d*cols,total=perCluster*n;
    resetGcceBatch<<<blocks(total),256,0,context->stream>>>(states,pointers.p+3*n,first.d,first.width,perCluster,total,first.vectorMode);
    GPU_CHECK(cudaGetLastError());
    auto multiply=[&](Z** in,Z** out,cublasOperation_t trans){
        GPU_CHECK(cublasZgemmBatched(context->blas,trans,CUBLAS_OP_N,first.d,cols,first.d,
            &one,vectors,first.d,in,first.d,&zero,out,first.d,n));
    };
    for(int p=0;p<=first.npulse;++p){
        multiply(states,scratch,CUBLAS_OP_C);
        phaseCceBatch<<<blocks(total),256,0,context->stream>>>(scratch,evals.p,first.times.p,
            first.d,first.width,perCluster,total,evals.p+n,delays[p]);
        GPU_CHECK(cudaGetLastError());
        multiply(scratch,states,CUBLAS_OP_N);
        if(p<first.npulse){
            pulseGcceBatch<<<blocks(total),256,0,context->stream>>>(states,scratch,pointers.p+4*n,p,
                first.d,first.q,first.b,perCluster,total);
            GPU_CHECK(cudaGetLastError());
            std::swap(states,scratch);
        }
    }
    Z** finalStates=states;
    auto& densityPointers=batchWork.density.device;
    if(!first.vectorMode){
        const int batches=n*count;
        std::vector<Z*> host((size_t)4*batches);
        // An odd pulse count swaps the physical state/scratch buffers.
        for(int c=0;c<n;++c)for(int t=0;t<count;++t){
            const size_t index=(size_t)c*count+t,offset=(size_t)t*first.d*first.d;
            host[index]=(first.npulse%2 ? works[c]->scratch.p : works[c]->states.p)+offset;
            host[batches+index]=works[c]->initial.p;
            host[2*batches+index]=(first.npulse%2 ? works[c]->states.p : works[c]->scratch.p)+offset;
            host[3*batches+index]=works[c]->full.p+offset;
        }
        batchWork.density.update(std::move(host));
        GPU_CHECK(cublasZgemmBatched(context->blas,CUBLAS_OP_N,CUBLAS_OP_N,first.d,first.d,first.d,
            &one,densityPointers.p,first.d,densityPointers.p+batches,first.d,&zero,densityPointers.p+2*batches,first.d,batches));
        GPU_CHECK(cublasZgemmBatched(context->blas,CUBLAS_OP_N,CUBLAS_OP_C,first.d,first.d,first.d,
            &one,densityPointers.p+2*batches,first.d,densityPointers.p,first.d,&zero,densityPointers.p+3*batches,first.d,batches));
        finalStates=pointers.p+5*n;
    }
    reduceGcceBatch<<<blocks(resultCount),256,0,context->stream>>>(finalStates,output.p,first.d,first.q,first.b,count,resultCount,first.vectorMode);
    GPU_CHECK(cudaGetLastError());
    GPU_CHECK(cudaMemcpyAsync(reduced,output.p,resultCount*sizeof(Z),cudaMemcpyDeviceToHost,context->stream));
    GPU_CHECK(cudaStreamSynchronize(context->stream));
}
