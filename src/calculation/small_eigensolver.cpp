// A42: the fixed-sweep complex Jacobi for the small cCCE eigensystems (small_eigensolver.h). Compiled with
// -ffp-contract=off (do_compile.sh rule, CMake source property) and the pragma below for a build that forgets the flag:
// nothing here may be fused, so that one lane computes what any other lane computes. No Eigen and no std::complex
// arithmetic in this file (only loads and stores): every operation is a double +, -, *, /, sqrt or fabs.
#if defined(__clang__)
#  pragma clang fp contract(off)
#else
#  pragma GCC optimize("fp-contract=off")
#endif
#include "../../include/calculation/small_eigensolver.h"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <utility>

namespace {
typedef std::complex<double> cd;

// The input range the Jacobi is trusted with (review 20261013-221715 X1); outside it ok is false and the caller uses Eigen.
// rangeMax: the squares of 2 N^2 parts below 1e300 sum without overflow (the norm, r, the off-diagonal sum). rangeMin: an
// a_pq zeroed by the r <= 1e-150 branch is then below 1e-20 of the largest element, so that zeroing cannot pass for
// convergence (and the squares of the largest elements stay normal numbers). Hamiltonians in rad/us are 1e-3..1e9.
constexpr double rangeMax = 1e150, rangeMin = 1e-130;

// L matrices of size N, structure of arrays [row][col][lane]. Rotation (p,q) with a_pq = r e^{i phi}:
// G(:,p) = c e_p - s e^{-i phi} e_q, G(:,q) = s e_p + c e^{-i phi} e_q (t = tan, the root of the smaller angle);
// A <- G^dagger A G, V <- V G; the 2x2 block is then set exactly (a_pp - t r, a_qq + t r, 0).
template<int N, int L>
void jacobi(const cd* const* H, double* const* w, cd* const* V, bool* ok){
    double ar[N][N][L], ai[N][N][L], vr[N][N][L], vi[N][N][L], norm[L];
    for (int l=0; l<L; ++l) norm[l] = 0.0;
    for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) for (int l=0; l<L; ++l){
        ar[i][j][l] = H[l][j*N+i].real(); ai[i][j][l] = H[l][j*N+i].imag();
        vr[i][j][l] = i == j ? 1.0 : 0.0;  vi[i][j][l] = 0.0;
        norm[l] = norm[l] + (ar[i][j][l]*ar[i][j][l] + ai[i][j][l]*ai[i][j][l]);
    }
    // N 2: the first rotation sets the one off-diagonal element to exactly 0; the second sweep is the identity and settles
    // the signs of zeros, the same bits as four sweeps (review 20261013-221715 X4)
    for (int sw=0; sw<(N == 2 ? 2 : smalleig::sweeps); ++sw) for (int p=0; p<N-1; ++p) for (int q=p+1; q<N; ++q){
        double c[L], s[L], er[L], ei[L], t[L], rr[L], dp[L], dq[L];
        for (int l=0; l<L; ++l){
            dp[l] = ar[p][p][l]; dq[l] = ar[q][q][l];
            const double xr = ar[p][q][l], xi = ai[p][q][l];
            const double r = std::sqrt(xr*xr + xi*xi);
            // identity when a_pq is negligible: below 1e-18 of the diagonal pair, or so small that its square lost bits
            // to underflow (e^{i phi} = a_pq/r would then not have modulus 1)
            const bool z = !(r > 1e-150) || !(r > 1e-18*(std::fabs(dp[l]) + std::fabs(dq[l])));
            const double rs = z ? 1.0 : r;
            const double th = (dq[l] - dp[l]) / (2.0*rs);
            const double tt = (th >= 0.0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th*th + 1.0));
            t[l] = z ? 0.0 : tt;
            c[l] = 1.0/std::sqrt(t[l]*t[l] + 1.0); s[l] = t[l]*c[l];
            er[l] = z ? 1.0 : xr/rs; ei[l] = z ? 0.0 : xi/rs;
            rr[l] = z ? 0.0 : r;
        }
        for (int k=0; k<N; ++k) for (int l=0; l<L; ++l){          // columns of A and V
            const double yr = ar[k][q][l]*er[l] + ai[k][q][l]*ei[l], yi = ai[k][q][l]*er[l] - ar[k][q][l]*ei[l];
            const double xr = ar[k][p][l], xi = ai[k][p][l];
            ar[k][p][l] = c[l]*xr - s[l]*yr; ai[k][p][l] = c[l]*xi - s[l]*yi;
            ar[k][q][l] = s[l]*xr + c[l]*yr; ai[k][q][l] = s[l]*xi + c[l]*yi;
            const double ur = vr[k][q][l]*er[l] + vi[k][q][l]*ei[l], ui = vi[k][q][l]*er[l] - vr[k][q][l]*ei[l];
            const double wr = vr[k][p][l], wi = vi[k][p][l];
            vr[k][p][l] = c[l]*wr - s[l]*ur; vi[k][p][l] = c[l]*wi - s[l]*ui;
            vr[k][q][l] = s[l]*wr + c[l]*ur; vi[k][q][l] = s[l]*wi + c[l]*ui;
        }
        for (int k=0; k<N; ++k) for (int l=0; l<L; ++l){          // rows of A
            const double yr = ar[q][k][l]*er[l] - ai[q][k][l]*ei[l], yi = ai[q][k][l]*er[l] + ar[q][k][l]*ei[l];
            const double xr = ar[p][k][l], xi = ai[p][k][l];
            ar[p][k][l] = c[l]*xr - s[l]*yr; ai[p][k][l] = c[l]*xi - s[l]*yi;
            ar[q][k][l] = s[l]*xr + c[l]*yr; ai[q][k][l] = s[l]*xi + c[l]*yi;
        }
        for (int l=0; l<L; ++l){
            ar[p][p][l] = dp[l] - t[l]*rr[l]; ar[q][q][l] = dq[l] + t[l]*rr[l];
            ar[p][q][l] = ai[p][q][l] = ar[q][p][l] = ai[q][p][l] = 0.0;
            ai[p][p][l] = ai[q][q][l] = 0.0;
        }
    }
    for (int l=0; l<L; ++l){
        if (!w[l]) continue;                                       // a padded lane
        double off = 0.0;
        for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) if (i != j) off = off + (ar[i][j][l]*ar[i][j][l] + ai[i][j][l]*ai[i][j][l]);
        ok[l] = off <= 1e-28*norm[l];                              // |off|_F <= 1e-14 |H|_F; false for a NaN
        // ascending, Eigen's selection sort (first minimum, swap when it is not in place)
        double d[N]; int perm[N];
        for (int i=0; i<N; ++i){ d[i] = ar[i][i][l]; perm[i] = i; }
        for (int i=0; i<N-1; ++i){
            int k = i;
            for (int j=i+1; j<N; ++j) if (d[j] < d[k]) k = j;
            if (k != i){ std::swap(d[i],d[k]); std::swap(perm[i],perm[k]); }
        }
        for (int i=0; i<N; ++i){
            w[l][i] = d[i];
            for (int r=0; r<N; ++r) V[l][i*N+r] = cd(vr[r][perm[i]][l], vi[r][perm[i]][l]);
        }
    }
}

// The checks of review 20261013-221715 X1, on each matrix by itself (not part of the lanes' arithmetic): H in range - every
// part finite and at most rangeMax (false for a NaN), the largest at least rangeMin unless H is 0 - and w, V finite.
template<int N>
bool inRange(const cd* H){
    const double* x = reinterpret_cast<const double*>(H);
    double big = 0.0; bool in = true;
    for (int k=0; k<2*N*N; ++k){
        const double a = std::fabs(x[k]);
        in = in & (a <= rangeMax);
        big = a > big ? a : big;
    }
    return in & ((big == 0.0) | (big >= rangeMin));
}
template<int N>
bool finite(const double* w, const cd* V){
    const double* x = reinterpret_cast<const double*>(V);
    bool in = true;
    for (int k=0; k<N; ++k) in = in & (std::fabs(w[k]) <= DBL_MAX);
    for (int k=0; k<2*N*N; ++k) in = in & (std::fabs(x[k]) <= DBL_MAX);
    return in;
}

template<int N>
void solveN(int count, const cd* const* H, double* const* w, cd* const* V, bool* ok){
    if (count == 1) jacobi<N,1>(H,w,V,ok);
    else {
        const cd* h[smalleig::maxLanes]; double* ww[smalleig::maxLanes]; cd* vv[smalleig::maxLanes]; bool k[smalleig::maxLanes];
        for (int l=0; l<smalleig::maxLanes; ++l){                 // padded lanes repeat lane 0 and are not written
            h[l] = H[l < count ? l : 0]; ww[l] = l < count ? w[l] : nullptr; vv[l] = l < count ? V[l] : nullptr;
        }
        jacobi<N,smalleig::maxLanes>(h,ww,vv,k);
        for (int l=0; l<count; ++l) ok[l] = k[l];
    }
    for (int l=0; l<count; ++l) ok[l] = ok[l] && inRange<N>(H[l]) && finite<N>(w[l],V[l]);
}

long long solvedCount = 0, fallbackCount = 0;
}

namespace smalleig {
bool on(){
    static const bool value = !(getenv("CCEX_SMALL_EIG") && !strcmp(getenv("CCEX_SMALL_EIG"),"0"));
    return value;
}
void solve(int D, int count, const cd* const* H, double* const* w, cd* const* V, bool* ok){
    if (D == 2) solveN<2>(count,H,w,V,ok); else solveN<4>(count,H,w,V,ok);
}
void count(long long solved, long long fallbacks){ solvedCount += solved; fallbackCount += fallbacks; }
void takeCounts(long long& solved, long long& fallbacks){
    solved = solvedCount; fallbacks = fallbackCount; solvedCount = fallbackCount = 0;
}
}
