// Device allocation lifetime probe, including persistent batch pointer tables; and (X12, review 20261013b-gpu-perf I1)
// the cache after an out-of-memory that came after some allocations: emptied by gpuOutOfMemoryRecovery(true), buffers
// released afterwards freed, and caching back after gpuOutOfMemoryRecovery(false).
#include "calculation/gpu/gpu_backend.h"
#include <cuda_runtime.h>
#include <complex>
#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using Z=std::complex<double>;
static void cycle(){
    constexpr int d=6,q=3,n=3,count=3;
    std::vector<Z> h(d*d),hb(d*d),rho(d*d),psi(d),pulse(q*q);
    for(int i=0;i<d;++i){h[i+i*d]=.1*i;hb[i+i*d]=.17*i;rho[i+i*d]=1./d;}
    h[1]=h[d]=.01;hb[2]=hb[2*d]=.02;psi[0]=1.;
    for(int i=0;i<q;++i)pulse[i+i*q]=1.;
    double times[]={.1,.2,.3},delays[]={.5,.5},durations[]={.05,.1,.15,.05,.1,.15};
    for(bool vector:{false,true}){
        std::vector<std::unique_ptr<GpuGcce,void(*)(GpuGcce*)>> gowned;
        std::vector<std::unique_ptr<GpuCce,void(*)(GpuCce*)>> cowned;
        std::vector<GpuGcce*> gs;std::vector<GpuCce*> cs;
        for(int i=0;i<n;++i){
            gowned.emplace_back(gpuGcceCreate(d,q,vector,count,h.data(),vector ? psi.data() : rho.data(),1,pulse.data()),gpuGcceDestroy);
            cowned.emplace_back(gpuCceCreate(d,!vector,count,1,h.data(),hb.data(),vector ? psi.data() : nullptr),gpuCceDestroy);
            gs.push_back(gowned.back().get());cs.push_back(cowned.back().get());
        }
        std::vector<Z> gr(n*count*q*q),cr(n*count);
        for(int k=0;k<2;++k){
            gpuGcceEvolveBatch(gs.data(),n,count,times,delays,gr.data());
            gpuCceEvolveBatch(cs.data(),n,count,durations,cr.data());
        }
    }
}
static size_t freeBytes(){
    size_t available=0,total=0;
    if(cudaDeviceSynchronize()!=cudaSuccess || cudaMemGetInfo(&available,&total)!=cudaSuccess)
        throw std::runtime_error("cannot query device memory");
    return available;
}
// X12: a dimension-63 work of 16384 steps fails part way (CCEX_TEST_GPU_OOM_BYTES 8 MiB, set in main), a small work lives
// across it as a batch's works would
static void recovery(){
    constexpr int d=6,q=3,big=63;
    std::vector<Z> h(d*d),rho(d*d),pulse(q*q),hBig(big*big),rhoBig(big*big);
    for(int i=0;i<d;++i){h[i+i*d]=.1*i;rho[i+i*d]=1./d;}
    for(int i=0;i<big;++i){hBig[i+i*big]=.01*i;rhoBig[i+i*big]=1./big;}
    for(int i=0;i<q;++i)pulse[i+i*q]=1.;
    const auto small=[&]{return gpuGcceCreate(d,q,false,3,h.data(),rho.data(),1,pulse.data());};
    gpuSetCacheLimit(size_t(64)*1024*1024);
    GpuGcce* alive=small();
    bool thrown=false;
    try { gpuGcceDestroy(gpuGcceCreate(big,q,false,16384,hBig.data(),rhoBig.data(),1,pulse.data())); }
    catch(const GpuOutOfMemory&){ thrown=true; }
    const size_t left=gpuMemoryStats().cachedBytes;
    if(!thrown || !left) throw std::runtime_error("recovery: no out-of-memory after some allocations (thrown "+std::to_string(thrown)+", cached "+std::to_string(left)+")");
    gpuOutOfMemoryRecovery(true);
    const size_t trimmed=gpuMemoryStats().cachedBytes;
    gpuGcceDestroy(alive);
    const size_t afterRelease=gpuMemoryStats().cachedBytes;
    gpuOutOfMemoryRecovery(false);
    gpuGcceDestroy(small());
    const size_t again=gpuMemoryStats().cachedBytes;
    std::cout<<"recovery: cached after the failure "<<left<<", after recovery "<<trimmed<<", after releasing a live work "<<afterRelease
             <<", after a new work once recovered "<<again<<'\n';
    if(trimmed || afterRelease || !again) throw std::runtime_error("recovery: the cache was not emptied, or not used again");
}
int main(){
    setenv("CCEX_TEST_GPU_OOM_BYTES","8388608",1);   // before the first allocation reads it; the cycles stay far below
    char description[512];if(!gpuInitialize(0,1,description,sizeof(description))){std::cerr<<description<<'\n';return 77;}
    try {
        for(size_t cache:{size_t(0),size_t(64)*1024*1024}){
            gpuSetCacheLimit(cache);for(int i=0;i<5;++i)cycle();
            gpuSetCacheLimit(0);const auto before=freeBytes();gpuSetCacheLimit(cache);
            for(int i=0;i<60;++i)cycle();
            gpuSetCacheLimit(0);const auto after=freeBytes();
            std::cout<<"cache_bytes="<<cache<<" cycles=60 before="<<before<<" after="<<after
                     <<" drift="<<static_cast<long long>(after)-static_cast<long long>(before)<<'\n';
            if(before!=after)throw std::runtime_error("device memory drift; inspect other GPU activity before attributing a leak");
        }
        std::cout<<"PASS device allocation lifetime probe\n";
        recovery();
        std::cout<<"PASS out-of-memory recovery: cache emptied, later releases freed, caching again after\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';gpuFinalize();return 1;}
    gpuFinalize();
}
