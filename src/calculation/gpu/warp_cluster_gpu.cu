// A40: the warp kernel (warp_cluster.h says what a cluster computes). One cluster per warp: its matrices in shared memory,
// the eigensystem by a parallel cyclic Jacobi (round-robin ordering: every round rotates disjoint index pairs at once, the
// lanes sharing the column and row updates), then the lanes take different time points (point i on lane i % 32), each
// running the whole interval chain of its point with the eigenvectors read from shared memory. The A40 step-0 prototype
// of this scheme ((internal notes), job 698) ran P1 dimension-12 gCCE clusters 1.6 times as fast
// as 32 CPU cores. Its one-thread-per-cluster variant was wrong at -O3 on sm_120 (the NVVM optimizer on dynamically
// indexed local arrays, as A26 step 0); this kernel keeps the matrices in shared memory and every per-lane vector loop
// unrolled, and the standing check test_warp_cluster_cuda compares it with CPU-dumped clusters of every compiled size.
#include "calculation/gpu/warp_cluster_gpu.h"
#include <cstdio>
#include <stdexcept>
#include <string>
#include <type_traits>

using smallcluster::C;
using smallcluster::cmake;
using smallcluster::cadd;
using smallcluster::csub;
using smallcluster::cmul;
using smallcluster::cconj;
using smallcluster::cscale;
using smallcluster::cabs2;
using warpcluster::Plan;

namespace {
void check(cudaError_t s, const char* op){
    if (s != cudaSuccess) throw std::runtime_error(std::string("warp-cluster GPU: ")+op+": "+cudaGetErrorString(s));
}

__constant__ Plan P;

__device__ __forceinline__ C cjmul(C a, C b){ return cmake(a.re*b.re+a.im*b.im, a.re*b.im-a.im*b.re); }   // conj(a) b

// exp(-i E tau) as refinedPhases (eigen_refinement.h) with zero corrections: the gCCE CPU path's phase
__device__ __forceinline__ C gccePhase(double e, double tau){
    const double a = e*tau, tail = fma(e,tau,-a);
    double s, c; sincos(a,&s,&c);
    const C z = cmake(c,-s);
    if (!(fabs(tail) < 0x1p-27)){ double s2, c2; sincos(tail,&s2,&c2); return cmul(z, cmake(c2,-s2)); }
    return cmul(z, cmake(1.0,-tail));
}

// u = ph[f] without a run-time index into the per-lane arrays (they would go to local memory with dynamic offsets, the
// pattern of the sm_120 miscompile; small_cluster_math.h's chain() selects its slot the same way)
template<int NF, int D> __device__ __forceinline__ void select(const C (&ph)[NF][D], int f, C (&u)[D]){
    #pragma unroll
    for (int ff=0; ff<NF; ++ff) if (ff == f || NF == 1){
        #pragma unroll
        for (int k=0; k<D; ++k) u[k] = ph[ff][k];
    }
}

// The parallel Jacobi's per-round rotations (N/2 <= 16 pairs)
struct Rotations { double c[16]; C s[16]; int p[16], q[16], on[16]; };

// h (D x D, column-major, the lower triangle read) = V diag(E) V^dagger with A, V row-major in shared memory; E (shared,
// D values) the eigenvalues. The rotation is small_cluster_math.h's (that of the A26 prototype, checked against Eigen to
// 1e-13); the threshold, the sweep limit, the convergence test, and the scaling by a power of two (the largest component
// in [0.5, 1), so that no finite h overflows or underflows in the rotations; review 20261012-233058-452278 P1-01/P3-02)
// are its too. False: not converged, or h, an eigenvalue or an eigenvector component not finite.
template<int D> __device__ bool warpJacobi(const C* __restrict__ h, C* A, C* V, double* E, Rotations& R, int lane){
    constexpr int N = D + (D & 1), NP = N/2, ROUNDS = N-1;
    static_assert(NP <= 16, "dimension");
    double largest = 0; int finite = 1;
    for (int t=lane; t<D*D; t+=32){
        const int i = t/D, j = t%D;
        C a = i>=j ? h[i+j*D] : cconj(h[j+i*D]);
        if (i == j) a.im = 0;
        finite = finite && isfinite(a.re) && isfinite(a.im);
        largest = fmax(largest, fmax(fabs(a.re), fabs(a.im)));
        A[t] = a; V[t] = cmake(i==j ? 1.0 : 0.0, 0.0);
    }
    for (int o=16; o; o>>=1) largest = fmax(largest, __shfl_xor_sync(0xffffffffu, largest, o));
    if (!__all_sync(0xffffffffu, finite)) return false;
    int exponent = 0;
    if (largest > 0){
        frexp(largest, &exponent);
        for (int t=lane; t<D*D; t+=32) A[t] = cmake(ldexp(A[t].re, -exponent), ldexp(A[t].im, -exponent));
    }
    __syncwarp();
    for (int sweep=0; sweep<smallcluster::jacobiSweeps; ++sweep){
        int rotated = 0;
        for (int round=0; round<ROUNDS; ++round){
            if (lane < NP){
                // circle method: N-1 fixed, the others turn; index D is the dummy of an odd D
                const int a = lane==0 ? ROUNDS : (round+lane)%ROUNDS, b = lane==0 ? round : (round-lane+ROUNDS)%ROUNDS;
                const int p = a<b ? a : b, q = a<b ? b : a;
                int on = 0; double c = 1; C s = cmake(0,0);
                if (q < D){
                    const C off = A[p*D+q];
                    const double m = hypot(off.re, off.im);
                    const double ap = A[p*D+p].re, aq = A[q*D+q].re;
                    if (m > 1e-17*(fabs(ap)+fabs(aq)) && m > 0){
                        const double g = aq-ap;
                        const double t = g==0 ? 1.0 : copysign(2*m/(fabs(g)+hypot(g,2*m)), g);
                        c = 1/sqrt(1+t*t); s = cscale(off, t*c/m); on = 1;
                    }
                }
                R.c[lane] = c; R.s[lane] = s; R.p[lane] = p; R.q[lane] = q; R.on[lane] = on;
                rotated |= on;
            }
            __syncwarp();
            for (int t=lane; t<2*D*NP; t+=32){            // columns: A J and V J (row k, pair pr)
                const int m = t/(D*NP), k = (t%(D*NP))/NP, pr = t%NP;
                if (!R.on[pr]) continue;
                C* X = m ? V : A;
                const int p = R.p[pr], q = R.q[pr];
                const double c = R.c[pr]; const C s = R.s[pr], sc = cconj(s);
                const C a = X[k*D+p], b = X[k*D+q];
                X[k*D+p] = csub(cscale(a,c), cmul(sc,b)); X[k*D+q] = cadd(cmul(s,a), cscale(b,c));
            }
            __syncwarp();
            for (int t=lane; t<D*NP; t+=32){              // rows: J^dagger A (column k, pair pr)
                const int k = t/NP, pr = t%NP;
                if (!R.on[pr]) continue;
                const int p = R.p[pr], q = R.q[pr];
                const double c = R.c[pr]; const C s = R.s[pr], sc = cconj(s);
                const C a = A[p*D+k], b = A[q*D+k];
                A[p*D+k] = csub(cscale(a,c), cmul(s,b)); A[q*D+k] = cadd(cmul(sc,a), cscale(b,c));
            }
            __syncwarp();
            if (lane < NP && R.on[lane]){ A[R.p[lane]*D+R.q[lane]] = cmake(0,0); A[R.q[lane]*D+R.p[lane]] = cmake(0,0); }
            __syncwarp();
        }
        if (!__any_sync(0xffffffffu, rotated)) break;
    }
    // the test on the scaled matrix; explicit finite flags, as fmax drops a NaN
    double off = 0, diag = 0;
    for (int t=lane; t<D*D; t+=32){
        finite = finite && isfinite(V[t].re) && isfinite(V[t].im);
        if (t/D == t%D){ finite = finite && isfinite(A[t].re); diag = fmax(diag, fabs(A[t].re)); }
        else { const double a = hypot(A[t].re, A[t].im); finite = finite && isfinite(a); off = fmax(off, a); }
    }
    for (int o=16; o; o>>=1){
        off = fmax(off, __shfl_xor_sync(0xffffffffu, off, o));
        diag = fmax(diag, __shfl_xor_sync(0xffffffffu, diag, o));
    }
    if (lane < D){ E[lane] = ldexp(A[lane*D+lane].re, exponent); finite = finite && isfinite(E[lane]); }
    __syncwarp();
    return __all_sync(0xffffffffu, finite) && off <= 1e-12*diag;
}

// The end of a cluster: status ok, or - a value not finite on any lane - every output zeroed and statusNotFinite (the
// host computes the cluster on the CPU), as the packed kernel
__device__ __forceinline__ void finish(C* o, int values, int bad, unsigned char& status, int lane){
    if (__any_sync(0xffffffffu, bad)){
        __syncwarp();
        for (int t=lane; t<values; t+=32) o[t] = cmake(0,0);
        if (lane == 0) status = smallcluster::statusNotFinite;
        return;
    }
    if (lane == 0) status = smallcluster::statusOk;
}

// Warps per block: as many as fit in 48 KB of static shared memory, at most 4
template<class S> struct WarpsFor { static constexpr int value = sizeof(S) > 24000 ? 1 : sizeof(S) > 16000 ? 2 : sizeof(S) > 12000 ? 3 : 4; };

// ---------------------------------------------------------------- gCCE
template<int D> struct GcceShared { C A[D*D]; C V[D*D]; double E[D]; Rotations R; };
// X11 (review 20261013b-gpu-perf O1): Hahn (one pulse) keeps F = V^dagger (Uq (x) I) V, the pulse in the eigenbasis
template<int D> struct GcceHahnShared : GcceShared<D> { C F[D*D]; };
template<int D, bool HAHN> using GcceSharedFor = typename std::conditional<HAHN, GcceHahnShared<D>, GcceShared<D>>::type;

// HAHN (warpClusterLaunch: dimension 12, two intervals): each time point applies the pulse as F, 2 D^2 products instead
// of 3 D^2 + QD (the last digits differ from the general chain, within the GPU tolerance)
// The HAHN kernel is allowed 255 registers (minimum one block per SM): with the default 168 nvcc spilled it on sm_120 and
// P1 gCCE2 warp runs took 1.19 times as long as before X11 (job 17444); with 255, 0.85-0.87 (job 17447). 0 for the
// others leaves their code as it was (SASS identical on sm_120).
template<int QD, int BD, int NF, bool HAHN>
__global__ void __launch_bounds__(32*WarpsFor<GcceSharedFor<QD*BD,HAHN>>::value, HAHN ? 1 : 0) gcceKernel(const C* __restrict__ in, const WarpEntry* __restrict__ list,
                                                                               long n, C* __restrict__ out, unsigned char* __restrict__ status){
    constexpr int D = QD*BD, W = WarpsFor<GcceSharedFor<D,HAHN>>::value;
    __shared__ GcceSharedFor<D,HAHN> all[W];
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    const long c = (long)blockIdx.x*W + wid;
    if (c >= n) return;                           // whole warps leave together
    GcceSharedFor<D,HAHN>& w = all[wid];
    const WarpEntry e = list[c];
    const C* h = in + e.in;
    const int nstep = P.time.nstep, per = P.dm ? QD*QD : 1;
    C* o = out + e.out*(long)nstep*per;
    if (!warpJacobi<D>(h, w.A, w.V, w.E, w.R, lane)){
        for (int t=lane; t<nstep*per; t+=32) o[t] = cmake(0,0);
        if (lane == 0) status[e.out] = smallcluster::statusNotConverged;
        return;
    }
    if constexpr (HAHN){                          // F = V^dagger T, T = (Uq (x) I) V (T over A, no longer needed)
        for (int t=lane; t<D*D; t+=32){
            const int k = t/D, j = t%D, a = k/BD, b = k%BD;
            C s = cmake(0,0);
            for (int q=0; q<QD; ++q) s = cadd(s, cmul(P.pulse[0][a+q*QD], w.V[(q*BD+b)*D+j]));
            w.A[t] = s;
        }
        __syncwarp();
        for (int t=lane; t<D*D; t+=32){
            const int i = t/D, j = t%D;
            C s = cmake(0,0);
            for (int k=0; k<D; ++k) s = cadd(s, cjmul(w.V[k*D+i], w.A[k*D+j]));
            w.F[t] = s;
        }
        __syncwarp();
    }
    // the starts in the eigenbasis, S[m*D+j] (over A, no longer needed): V^dagger psi0, or V^dagger (qpsi0 (x) e_m)
    const int nstart = P.ensemble ? BD : 1;
    C* S = w.A;
    if (!P.ensemble){
        if (lane < D){
            C a = cmake(0,0);
            for (int k=0; k<D; ++k) a = cadd(a, cjmul(w.V[k*D+lane], h[D*D+k]));
            S[lane] = a;
        }
    } else {
        for (int t=lane; t<BD*D; t+=32){
            const int m = t/D, j = t%D;
            C a = cmake(0,0);
            for (int q=0; q<QD; ++q) a = cadd(a, cjmul(w.V[(q*BD+m)*D+j], P.qpsi0[q]));
            S[t] = a;
        }
    }
    __syncwarp();
    const int intervals = P.time.intervals;
    int bad = 0;                                  // a value not finite (the cluster then goes back to the CPU)
    double tfree = 0.0;                           // the time of point i summed as the CPU sums it (tfree += deltat)
    for (int k=0; k<lane; ++k) tfree += P.time.deltat;
    for (int i=lane; i<nstep; i+=32){
        C ph[NF][D];
        #pragma unroll
        for (int f=0; f<NF; ++f){
            if (f >= P.time.distinct) break;
            const double tau = tfree*P.time.fraction[f];
            #pragma unroll
            for (int k=0; k<D; ++k) ph[f][k] = gccePhase(w.E[k], tau);
        }
        C L = cmake(0,0);
        C rho[QD*QD];
        #pragma unroll
        for (int k=0; k<QD*QD; ++k) rho[k] = cmake(0,0);
        for (int m=0; m<nstart; ++m){
            C x[D], y[D], u[D];
            select<NF,D>(ph, P.time.same[0], u);
            #pragma unroll
            for (int j=0; j<D; ++j) x[j] = cmul(u[j], S[m*D+j]);
            if constexpr (HAHN){                  // y = V u1 F u0 S
                #pragma unroll
                for (int j=0; j<D; ++j){ C a = cmake(0,0);
                    #pragma unroll
                    for (int k=0; k<D; ++k) a = cadd(a, cmul(w.F[j*D+k], x[k])); y[j] = a; }
                select<NF,D>(ph, P.time.same[1], u);
                #pragma unroll
                for (int j=0; j<D; ++j) x[j] = cmul(u[j], y[j]);
            }
            #pragma unroll
            for (int k=0; k<D; ++k){ C a = cmake(0,0);
                #pragma unroll
                for (int j=0; j<D; ++j) a = cadd(a, cmul(w.V[k*D+j], x[j])); y[k] = a; }
            if constexpr (!HAHN) for (int p=1; p<intervals; ++p){
                const C* U = P.pulse[p-1];
                #pragma unroll
                for (int a=0; a<QD; ++a)
                    #pragma unroll
                    for (int b=0; b<BD; ++b){ C t = cmake(0,0);
                        #pragma unroll
                        for (int q=0; q<QD; ++q) t = cadd(t, cmul(U[a+q*QD], y[q*BD+b])); x[a*BD+b] = t; }
                select<NF,D>(ph, P.time.same[p], u);
                #pragma unroll
                for (int j=0; j<D; ++j){ C a = cmake(0,0);
                    #pragma unroll
                    for (int k=0; k<D; ++k) a = cadd(a, cjmul(w.V[k*D+j], x[k])); y[j] = cmul(a, u[j]); }
                #pragma unroll
                for (int k=0; k<D; ++k){ C a = cmake(0,0);
                    #pragma unroll
                    for (int j=0; j<D; ++j) a = cadd(a, cmul(w.V[k*D+j], y[j])); x[k] = a; }
                #pragma unroll
                for (int k=0; k<D; ++k) y[k] = x[k];
            }
            // Psi(i,b) = y[i*BD+b]
            if (P.dm){
                #pragma unroll
                for (int a=0; a<QD; ++a)
                    #pragma unroll
                    for (int b=0; b<QD; ++b){ C s = cmake(0,0);
                        #pragma unroll
                        for (int k=0; k<BD; ++k) s = cadd(s, cmul(y[a*BD+k], cconj(y[b*BD+k]))); rho[a+b*QD] = cadd(rho[a+b*QD], s); }
            } else {
                #pragma unroll
                for (int b=0; b<BD; ++b){
                    C a = cmake(0,0), bb = cmake(0,0);
                    #pragma unroll
                    for (int q=0; q<QD; ++q){ a = cadd(a, cjmul(P.psia[q], y[q*BD+b])); bb = cadd(bb, cjmul(P.psib[q], y[q*BD+b])); }
                    L = cadd(L, cmul(a, cconj(bb)));
                }
            }
        }
        const double scale = 1.0/nstart;
        if (P.dm){
            #pragma unroll
            for (int k=0; k<QD*QD; ++k){
                const C v = P.ensemble ? cscale(rho[k], scale) : rho[k];
                bad = bad || !isfinite(v.re) || !isfinite(v.im);
                o[(long)i*per+k] = v;
            }
        } else {
            const C v = P.ensemble ? cscale(L, scale) : L;
            bad = bad || !isfinite(v.re) || !isfinite(v.im);
            o[i] = v;
        }
        for (int k=0; k<32; ++k) tfree += P.time.deltat;
    }
    finish(o, nstep*per, bad, status[e.out], lane);
}

// ---------------------------------------------------------------- cCCE
template<int D> struct CceShared { C A[D*D]; C Vu[D*D]; C Vd[D*D]; double Eu[D], Ed[D]; C r0[D], c0[D]; Rotations R; };

template<int D, int NF>
__global__ void __launch_bounds__(32*WarpsFor<CceShared<D>>::value) cceKernel(const C* __restrict__ in, const WarpEntry* __restrict__ list,
                                                                         long n, C* __restrict__ out, unsigned char* __restrict__ status){
    constexpr int W = WarpsFor<CceShared<D>>::value;
    __shared__ CceShared<D> all[W];
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    const long c = (long)blockIdx.x*W + wid;
    if (c >= n) return;
    CceShared<D>& w = all[wid];
    const WarpEntry e = list[c];
    const C* h = in + e.in;
    const int nstep = P.time.nstep;
    C* o = out + e.out*(long)nstep;
    bool ok = warpJacobi<D>(h, w.A, w.Vu, w.Eu, w.R, lane);
    __syncwarp();
    ok = warpJacobi<D>(h+D*D, w.A, w.Vd, w.Ed, w.R, lane) && ok;
    if (!ok){
        for (int t=lane; t<nstep; t+=32) o[t] = cmake(0,0);
        if (lane == 0) status[e.out] = smallcluster::statusNotConverged;
        return;
    }
    __syncwarp();
    // M = Vu^dagger Vd over A
    for (int t=lane; t<D*D; t+=32){
        const int i = t/D, j = t%D;
        C s = cmake(0,0);
        for (int k=0; k<D; ++k) s = cadd(s, cmul(cconj(w.Vu[k*D+i]), w.Vd[k*D+j]));
        w.A[t] = s;
    }
    if (!P.ensemble && lane < D){                 // psi^dagger Vd and Vu^dagger psi
        const C* psi = h+2*D*D;
        C r = cmake(0,0), s = cmake(0,0);
        for (int k=0; k<D; ++k){ r = cadd(r, cmul(cconj(psi[k]), w.Vd[k*D+lane])); s = cadd(s, cmul(cconj(w.Vu[k*D+lane]), psi[k])); }
        w.r0[lane] = r; w.c0[lane] = s;
    }
    __syncwarp();
    const C (&M)[D][D] = *reinterpret_cast<const C(*)[D][D]>(w.A);
    const C (&r0)[D] = w.r0;
    const C (&c0)[D] = w.c0;
    const bool hahn = P.time.intervals == 2 && P.time.same[0] == P.time.same[1];
    int bad = 0;
    double tfree = 0.0;
    for (int k=0; k<lane; ++k) tfree += P.time.deltat;
    for (int i=lane; i<nstep; i+=32){
        C up[NF][D], down[NF][D];
        #pragma unroll
        for (int f=0; f<NF; ++f){
            if (f >= P.time.distinct) break;
            const double tau = tfree*P.time.fraction[f];
            #pragma unroll
            for (int k=0; k<D; ++k){ up[f][k] = smallcluster::phase(w.Eu[k]*tau); down[f][k] = smallcluster::phase(w.Ed[k]*tau); }
        }
        C L;
        if (!P.ensemble){
            L = smallcluster::chain<D,NF>(P.time, M, r0, c0, up, down);
        } else if (hahn){
            // the Hahn-ensemble identity of small_cluster_math.h: sum_ab |X_ab|^2 conj(Du_a) Du_b / D, X = M Dd M^dagger
            const C (&u)[D] = up[0];
            const C (&d)[D] = down[0];
            C sum = cmake(0,0);
            #pragma unroll
            for (int a=0; a<D; ++a){
                C rowSum = cmake(0,0);
                #pragma unroll
                for (int b=0; b<D; ++b){
                    C x = cmake(0,0);
                    #pragma unroll
                    for (int k=0; k<D; ++k) x = cadd(x, cmul(cmul(M[a][k], d[k]), cconj(M[b][k])));
                    rowSum = cadd(rowSum, cscale(u[b], cabs2(x)));
                }
                sum = cadd(sum, cmul(cconj(u[a]), rowSum));
            }
            L = cscale(sum, 1.0/D);
        } else {
            L = cmake(0,0);
            for (int k=0; k<D; ++k){
                C rk[D], ck[D];                   // e_k^dagger Vd and Vu^dagger e_k
                #pragma unroll
                for (int j=0; j<D; ++j){ rk[j] = w.Vd[k*D+j]; ck[j] = cconj(w.Vu[k*D+j]); }
                L = cadd(L, smallcluster::chain<D,NF>(P.time, M, rk, ck, up, down));
            }
            L = cscale(L, 1.0/D);
        }
        bad = bad || !isfinite(L.re) || !isfinite(L.im);
        o[i] = L;
        for (int k=0; k<32; ++k) tfree += P.time.deltat;
    }
    finish(o, nstep, bad, status[e.out], lane);
}

template<class K> void launch(K kernel, int warps, long n, cudaStream_t stream, const C* in, const WarpEntry* list, C* out, unsigned char* status){
    kernel<<<(unsigned)((n+warps-1)/warps), 32*warps, 0, stream>>>(in, list, n, out, status);
    check(cudaGetLastError(), "launch");
}
template<int QD, int BD, int NF, bool HAHN> void gcceOne(long n, cudaStream_t s, const C* in, const WarpEntry* list, C* out, unsigned char* st){
    launch(gcceKernel<QD,BD,NF,HAHN>, WarpsFor<GcceSharedFor<QD*BD,HAHN>>::value, n, s, in, list, out, st);
}
template<int QD, int BD, bool HAHN> void gcceSlots(int nf, long n, cudaStream_t s, const C* in, const WarpEntry* list, C* out, unsigned char* st){
    if (nf == 1) gcceOne<QD,BD,1,HAHN>(n,s,in,list,out,st);
    else if (nf == 2) gcceOne<QD,BD,2,HAHN>(n,s,in,list,out,st);
    else gcceOne<QD,BD,4,HAHN>(n,s,in,list,out,st);
}
template<int D> void cceSlots(int nf, long n, cudaStream_t s, const C* in, const WarpEntry* list, C* out, unsigned char* st){
    if (nf == 1) launch(cceKernel<D,1>, WarpsFor<CceShared<D>>::value, n, s, in, list, out, st);
    else if (nf == 2) launch(cceKernel<D,2>, WarpsFor<CceShared<D>>::value, n, s, in, list, out, st);
    else launch(cceKernel<D,4>, WarpsFor<CceShared<D>>::value, n, s, in, list, out, st);
}
}

void warpClusterUpload(const Plan& plan){
    check(cudaDeviceSynchronize(), "synchronize before the plan");
    check(cudaMemcpyToSymbol(P, &plan, sizeof(Plan)), "plan");
}

void warpClusterLaunch(const Plan& plan, int dim, const C* in, const WarpEntry* list, long n, C* out, unsigned char* status,
                       cudaStream_t stream){
    if (n <= 0) return;
    const int nf = smallcluster::slotsFor(plan.time.distinct);
    if (!nf || !warpcluster::compiled(plan.mode, plan.qdim, dim)) throw std::runtime_error("warp-cluster GPU: not compiled for this cluster");
    if (plan.mode == warpcluster::modeCce){
        if (dim == 8) cceSlots<8>(nf,n,stream,in,list,out,status);
        else cceSlots<16>(nf,n,stream,in,list,out,status);
        return;
    }
    switch (dim/3){                               // qubit dimension 3 (warpcluster::compiled)
        case 2: gcceSlots<3,2,false>(nf,n,stream,in,list,out,status); break;
        case 3: gcceSlots<3,3,false>(nf,n,stream,in,list,out,status); break;
        case 4:                                   // X11: the pulse in the eigenbasis for Hahn, at this dimension only so far
            if (plan.time.intervals == 2) gcceSlots<3,4,true>(nf,n,stream,in,list,out,status);
            else gcceSlots<3,4,false>(nf,n,stream,in,list,out,status);
            break;
        case 6: gcceSlots<3,6,false>(nf,n,stream,in,list,out,status); break;
        case 7: gcceSlots<3,7,false>(nf,n,stream,in,list,out,status); break;
        case 8: gcceSlots<3,8,false>(nf,n,stream,in,list,out,status); break;
        default: gcceSlots<3,9,false>(nf,n,stream,in,list,out,status); break;
    }
}
