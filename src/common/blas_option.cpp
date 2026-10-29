/**
 * @file blas_option.cpp
 * @brief The optional external BLAS/LAPACK (Intel MKL; Apple's Accelerate on macOS) calls; empty without CCEX_BLAS_DIRECT.
 */
// A46: the library calls of the external BLAS/LAPACK option (include/common/blas_option.h). Empty unless the build
// defines CCEX_BLAS_DIRECT. A file of its own, so that the code around the call sites compiles as before.
#include "common/blas_option.h"

#ifdef CCEX_BLAS_DIRECT
#include <algorithm>

extern "C" {
void zgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*,
            const std::complex<double>*, const int*, const std::complex<double>*, const int*,
            const std::complex<double>*, std::complex<double>*, const int*);
void zheevd_(const char*, const char*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*,
             const int*, double*, const int*, int*, const int*, int*);
}

#ifndef CCEX_BLAS_ACCELERATE
#include <cstdlib>
extern "C" int MKL_CBWR_Set(int) __attribute__((weak));   // MKL's; null with another library behind the calls
namespace {
// MKL's results depend on where in memory the arrays start unless its reproducible mode is on: with MKL 2026.2 two runs
// that take the same path wrote different last digits (check_input_forms column_evolution C1-C3, 2026-10-28). AUTO keeps
// the CPU's fastest code path (the run times did not change), STRICT makes it independent of the alignment. Before any
// MKL call (static initialization); an MKL_CBWR the user set is left to MKL
const int cbwrSet = MKL_CBWR_Set && !std::getenv("MKL_CBWR") ? MKL_CBWR_Set(2 /* MKL_CBWR_AUTO */ | 0x10000 /* MKL_CBWR_STRICT */) : 0;
}
#endif

void ccexBlas::hahnProduct(const Eigen::MatrixXcd& M, const Eigen::VectorXcd& d, Eigen::MatrixXcd& scratch,
                           Eigen::MatrixXcd& X){
    const int n = int(M.rows());
    const std::complex<double> one(1,0), zero(0,0);
    scratch = M.array().rowwise() * d.transpose().array();
    X.resize(n,n);
    zgemm_("N","C",&n,&n,&n,&one,scratch.data(),&n,M.data(),&n,&zero,X.data(),&n);
}

void ccexBlas::product(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, Eigen::MatrixXcd& C){
    const int m = int(A.rows()), k = int(A.cols()), n = int(B.cols());
    const std::complex<double> one(1,0), zero(0,0);
    C.resize(m,n);
    zgemm_("N","N",&m,&n,&k,&one,A.data(),&m,B.data(),&k,&zero,C.data(),&m);
}

void CcexEigensolver::compute(const Eigen::MatrixXcd& H){
    const int n = int(H.rows());
    lapack_ = n >= ccexBlas::eigMin;
    if (!lapack_){ es_.compute(H); return; }
    V_ = H;
    w_.resize(n);
    int info = 0;
    if (n != wsDim_){                                   // workspace query once per dimension
        std::complex<double> wq; double rq; int iq; const int q = -1;
        zheevd_("V","L",&n,V_.data(),&n,w_.data(),&wq,&q,&rq,&q,&iq,&q,&info);
        work_.resize(std::max(1,int(wq.real()))); rwork_.resize(std::max(1,int(rq))); iwork_.resize(std::max(1,iq));
        wsDim_ = n;
    }
    const int lw = int(work_.size()), lrw = int(rwork_.size()), liw = int(iwork_.size());
    zheevd_("V","L",&n,V_.data(),&n,w_.data(),work_.data(),&lw,rwork_.data(),&lrw,iwork_.data(),&liw,&info);
    info_ = info == 0 ? Eigen::Success : Eigen::NumericalIssue;
}
#endif
