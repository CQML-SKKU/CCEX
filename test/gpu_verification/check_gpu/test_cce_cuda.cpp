// Independent cCCE oracle in the original bath basis, using general matrix
// exponentials. The GPU instead works in two conditional eigenbases.
#include "calculation/gpu/gpu_backend.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Z=std::complex<double>;
using M=Eigen::MatrixXcd;
using V=Eigen::VectorXcd;

static double runCase(int d, bool ensemble, int pulses, bool sameHamiltonian, int clusters=1, bool batchCreate=false){
    std::vector<M> alpha, beta;
    std::vector<V> initial;
    std::vector<std::unique_ptr<GpuCce,void(*)(GpuCce*)>> owned;
    std::vector<GpuCce*> works;
    for(int c=0;c<clusters;++c){
        M a=M::Random(d,d), b=M::Random(d,d);
        M ha=(a+a.adjoint()).eval()/double(d);
        M hb=(b+b.adjoint()).eval()/double(d)+0.23*M::Identity(d,d);
        if (sameHamiltonian) hb=ha;
        V psi=V::Random(d); psi.normalize();
        alpha.push_back(ha);beta.push_back(hb);initial.push_back(psi);
        if(batchCreate) continue;
        owned.emplace_back(gpuCceCreate(d,ensemble,4,pulses,ha.data(),hb.data(),
                                       ensemble ? nullptr : psi.data()),gpuCceDestroy);
        works.push_back(owned.back().get());
    }
    if(batchCreate){     // the eigensystems of all the clusters in one call (gpuCceCreateBatch)
        std::vector<const void*> hA,hB,psi;
        for(int c=0;c<clusters;++c){hA.push_back(alpha[c].data());hB.push_back(beta[c].data());psi.push_back(initial[c].data());}
        works=gpuCceCreateBatch(clusters,d,ensemble,4,pulses,hA.data(),hB.data(),ensemble ? nullptr : psi.data());
        for(GpuCce* w:works) owned.emplace_back(w,gpuCceDestroy);
    }
    double maximum=0;
    for (int iteration=0; iteration<3; ++iteration){
        int count=iteration==1 ? 1 : 4;
        const double times[]={0,0.37,1.131,2.419};
        std::vector<double> durations((pulses+1)*count);
        for (int p=0; p<=pulses; ++p) for (int t=0; t<count; ++t){
            double time=iteration==1 ? 0.715 : times[t];
            if (ensemble){
                // The half interval in double, as the CPU ensemble path has computed it since 2026-09-19
                // (simulator_cce.cpp: tau = time/(2.0*npulse)); it was a float here until 2026-10-08.
                const double tau=pulses ? time/(2.0*pulses) : 0;
                durations[p*count+t]=!pulses ? time : (p==0 || p==pulses ? tau : 2.0*tau);
            }else{
                // Unequal intervals exercise the manual-sequence single path.
                durations[p*count+t]=time*(p+1)/(0.5*(pulses+1)*(pulses+2));
            }
        }
        std::vector<Z> result((size_t)clusters*count);
        gpuCceEvolveBatch(works.data(),clusters,count,durations.data(),result.data());
        const auto before=gpuMemoryStats();
        std::vector<Z> repeated(result.size());
        gpuCceEvolveBatch(works.data(),clusters,count,durations.data(),repeated.data());
        const auto after=gpuMemoryStats();
        if(before.pointerUploads!=after.pointerUploads || before.allocations!=after.allocations || result!=repeated)
            throw std::runtime_error("cCCE unchanged batch did not reuse tables/output exactly");
        for(int c=0;c<clusters;++c) for (int t=0; t<count; ++t){
            const M& ha=alpha[c]; const M& hb=beta[c]; const V& psi=initial[c];
            M ua=M::Identity(d,d), ub=ua;
            for (int p=0; p<=pulses; ++p){
                const M& ah=p%2 ? hb : ha;
                const M& bh=p%2 ? ha : hb;
                ua=((Z(0,-durations[p*count+t])*ah).exp()*ua).eval();
                ub=((Z(0,-durations[p*count+t])*bh).exp()*ub).eval();
            }
            Z expected;
            if (ensemble){
                // The mean over a complete basis of what one state gives below. Until
                // 2026-09-19 the ensemble returned the conjugate for an even pulse
                // count and this test expected it; CPU and CUDA both dropped that
                // then (db94f3c23), checked against exact three-spin answers.
                expected=(ub.adjoint()*ua).trace()/double(d);
            }else expected=(psi.adjoint()*ub.adjoint()*ua*psi)(0,0);
            double error=std::abs(result[(size_t)c*count+t]-expected);
            if (!std::isfinite(error) || error>3e-11 || std::abs(result[(size_t)c*count+t])>1+3e-11)
                throw std::runtime_error("cCCE conditional propagation mismatch");
            if (sameHamiltonian && std::abs(result[(size_t)c*count+t]-Z(1,0))>3e-11)
                throw std::runtime_error("identical Hamiltonians must preserve coherence");
            maximum=std::max(maximum,error);
        }
        if(clusters>2){
            std::swap(works[1],works[2]);std::swap(alpha[1],alpha[2]);
            std::swap(beta[1],beta[2]);std::swap(initial[1],initial[2]);
        }
    }
    std::cout<<"PASS cCCE dim="<<d<<" ensemble="<<ensemble<<" pulses="<<pulses
             <<" clusters="<<clusters<<" same_H="<<sameHamiltonian<<" max_error="<<maximum<<'\n';
    return maximum;
}

int main(){
    char description[512]={};
    if (!gpuInitialize(0,1,description,sizeof(description))){
        std::cerr<<description<<'\n'; return 77;
    }
    std::cout<<description<<'\n'<<std::scientific;
    try {
        double maximum=0;
        for (bool ensemble : {true,false}){
            for (int p=0; p<=4; ++p) maximum=std::max(maximum,runCase(3,ensemble,p,false));
            maximum=std::max(maximum,runCase(1,ensemble,0,false));
            maximum=std::max(maximum,runCase(4,ensemble,3,true));
            maximum=std::max(maximum,runCase(16,ensemble,2,false));
            for(int p=0;p<=4;++p) maximum=std::max(maximum,runCase(3,ensemble,p,false,5));
            maximum=std::max(maximum,runCase(1,ensemble,4,true,3));
            maximum=std::max(maximum,runCase(16,ensemble,3,false,8));
            for(int p=0;p<=4;++p) maximum=std::max(maximum,runCase(3,ensemble,p,false,5,true));
            maximum=std::max(maximum,runCase(1,ensemble,4,true,3,true));
            maximum=std::max(maximum,runCase(16,ensemble,3,false,8,true));
            maximum=std::max(maximum,runCase(64,ensemble,2,false,4,true));
        }
        auto memory=gpuMemoryStats();
        if(!memory.reuses || memory.peakCachedBytes>64*1024*1024)
            throw std::runtime_error("device buffer cache did not reuse within budget");
        gpuSetCacheLimit(0);
        runCase(3,true,2,false,3);
        if(gpuMemoryStats().cachedBytes)throw std::runtime_error("disabled cache retained buffers");
        std::cout<<"cCCE overall max_error="<<maximum<<'\n';
    }catch (const std::exception& e){
        std::cerr<<"FAIL: "<<e.what()<<'\n'; gpuFinalize(); return 1;
    }
    gpuFinalize();
}
