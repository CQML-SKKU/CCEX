/**
 * @file noise_gpu.cu
 * @brief The GPU kernels of method noise: one thread per time point, each summing over a batch of the bath's modes.
 */
// In a noise run almost all the time goes into these sums, not into the clusters' Hamiltonians and eigensystems (1 ppm P1,
// 70553 pairs: 0.2 s of eigensystems beside 18 s of sums at 1001 time points): the host keeps the eigensystems and hands
// over (frequency, weight) lists. The arithmetic is that of noise.cpp. One thread takes one time point and one piece of
// the modes: with a thread per time point alone, 101 points fill one block and the card runs about 10 times slower per
// phase than at 1001 points (RTX 3050).
#include "calculation/gpu/noise_gpu.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(cudaError_t s, const char* op){
    if (s != cudaSuccess) throw std::runtime_error(std::string("noise GPU: ")+op+": "+cudaGetErrorString(s));
}
#define NOISE_CHECK(op) check((op),#op)

struct Device {
    double* p=nullptr; std::size_t n=0;
    void reserve(std::size_t want){
        if (want<=n) return;
        if (p) NOISE_CHECK(cudaFree(p));
        p=nullptr; n=0;
        NOISE_CHECK(cudaMalloc((void**)&p,want*sizeof(double)));
        n=want;
    }
    void copy(const double* from, std::size_t count){
        reserve(count);
        NOISE_CHECK(cudaMemcpy(p,from,count*sizeof(double),cudaMemcpyHostToDevice));
    }
    ~Device(){ if (p) cudaFree(p); }
};
const int kThreads=128;
const long kPiece=256;        // the modes of one thread

// |Y(w, T)|^2 by the intervals: sum_j s_j L_j sinc(w L_j / 2) exp(i w (centre of j)), noise.cpp's phaseIntegral
__global__ void phaseSum(const double* times, const double* fractions, int nfraction, const double* omega,
                         const double* weight, long n, double* out, int ntime, long total){
    const long id=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=total) return;
    const double T=times[id/(total/ntime)];
    const long first=(id%(total/ntime))*kPiece, last=first+kPiece<n ? first+kPiece : n;
    double sum=0;
    for (long k=first; k<last; ++k){
        const double w=omega[k];
        double re=0, im=0, start=0, sign=1;
        for (int j=0; j<nfraction; ++j){
            const double L=fractions[j]*T, x=0.5*w*L;
            const double sinc=fabs(x)<1e-4 ? 1-x*x/6 : sin(x)/x;
            double s, c;
            sincos(w*(start+0.5*L),&s,&c);
            re+=sign*L*sinc*c;
            im+=sign*L*sinc*s;
            start+=L;
            sign=-sign;
        }
        sum+=weight[k]*(re*re+im*im);
    }
    out[id]=sum;
}
__global__ void cosineSum(const double* times, const double* omega, const double* weightRe, const double* weightIm,
                          long n, double* out, int ntime, long total){
    const long id=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if (id>=total) return;
    const double time=times[id/(total/ntime)];
    const long first=(id%(total/ntime))*kPiece, last=first+kPiece<n ? first+kPiece : n;
    double sum=0;
    for (long k=first; k<last; ++k){
        double s, c;
        sincos(omega[k]*time,&s,&c);
        sum+=weightRe[k]*c-weightIm[k]*s;
    }
    out[id]=sum;
}
}

struct GpuNoise {
    int ntime, nfraction;
    Device times, fractions, omega, re, im, out;
    std::vector<double> host;
    long total=0;             // threads of the kernel: time points x pieces of the modes
    unsigned blocks(std::size_t modes){
        total=ntime*(long)((modes+kPiece-1)/kPiece);
        out.reserve(total);
        return (unsigned)((total+kThreads-1)/kThreads);
    }
};

GpuNoise* gpuNoiseCreate(int ntime, const double* times, int nfraction, const double* fractions){
    GpuNoise* work=new GpuNoise;
    work->ntime=ntime;
    work->nfraction=nfraction;
    try {
        work->times.copy(times,ntime);
        work->fractions.copy(fractions,nfraction);
    } catch (...) { delete work; throw; }
    return work;
}

void gpuNoiseDestroy(GpuNoise* work){ delete work; }

static void collect(GpuNoise* work, double* sum){
    NOISE_CHECK(cudaGetLastError());
    work->host.resize(work->total);
    NOISE_CHECK(cudaMemcpy(work->host.data(),work->out.p,work->total*sizeof(double),cudaMemcpyDeviceToHost));
    const long pieces=work->total/work->ntime;
    for (long id=0; id<work->total; ++id) sum[id/pieces]+=work->host[id];
}

void gpuNoisePhase(GpuNoise* work, std::size_t n, const double* omega, const double* weight, double* sum){
    if (!n) return;
    work->omega.copy(omega,n);
    work->re.copy(weight,n);
    const unsigned blocks=work->blocks(n);
    phaseSum<<<blocks,kThreads>>>(work->times.p,work->fractions.p,work->nfraction,work->omega.p,work->re.p,(long)n,
                                  work->out.p,work->ntime,work->total);
    collect(work,sum);
}

void gpuNoiseCosine(GpuNoise* work, std::size_t n, const double* omega, const double* re, const double* im, double* sum){
    if (!n) return;
    work->omega.copy(omega,n);
    work->re.copy(re,n);
    work->im.copy(im,n);
    const unsigned blocks=work->blocks(n);
    cosineSum<<<blocks,kThreads>>>(work->times.p,work->omega.p,work->re.p,work->im.p,(long)n,work->out.p,work->ntime,
                                   work->total);
    collect(work,sum);
}
