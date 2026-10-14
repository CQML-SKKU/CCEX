#pragma once
// A46: Intel MKL for the large clusters (do_compile.sh and CMake turn it on when they find MKL; CCEX_BLAS=off turns it off).
// Without CCEX_BLAS_DIRECT CcexEigensolver is Eigen's solver, inline as before (plus A42's setGiven).
// With it, the large clusters go to the library (src/common/blas_option.cpp): the cCCE ensemble Hahn step's
// X = (M diag(d)) M^dagger by zgemm from dimension ccexBlas::gemmMin on, and the cCCE/gCCE cluster eigensolve by zheevd
// from ccexBlas::eigMin on. The thresholds are where MKL was faster and the h-BN runs gained (A46 step 0, idm job 17337:
// zgemm pays from 12, zheevd from 27; 21 and 63 ran as fast and changed fewer outputs). Below them Eigen runs, as before.
#include <Eigen/Dense>

#include <complex>

#ifndef CCEX_BLAS_DIRECT
// Eigen's solver (inline, as before A46), and A42's eigensystems given from outside (setGiven: the small-cluster Jacobi,
// src/calculation/cce_setup_small.cpp) - eigenvectors()/eigenvalues() then return those until the next compute()
class CcexEigensolver {
public:
    void compute(const Eigen::MatrixXcd& H){ given_ = false; es_.compute(H); }
    void setGiven(int n, const double* w, const std::complex<double>* V){
        given_ = true; V_ = Eigen::Map<const Eigen::MatrixXcd>(V,n,n); w_ = Eigen::Map<const Eigen::VectorXd>(w,n);
    }
    const Eigen::MatrixXcd& eigenvectors() const { return given_ ? V_ : es_.eigenvectors(); }
    const Eigen::VectorXd& eigenvalues() const { return given_ ? w_ : es_.eigenvalues(); }
    Eigen::ComputationInfo info() const { return given_ ? Eigen::Success : es_.info(); }
private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es_;
    bool given_ = false;
    Eigen::MatrixXcd V_;
    Eigen::VectorXd w_;
};
#else
#include <complex>
#include <vector>

namespace ccexBlas {
constexpr int gemmMin = 21;
constexpr int eigMin = 63;
// X = (M diag(d)) M^dagger by zgemm; scratch holds M diag(d)
void hahnProduct(const Eigen::MatrixXcd& M, const Eigen::VectorXcd& d, Eigen::MatrixXcd& scratch, Eigen::MatrixXcd& X);
}

// Eigen's SelfAdjointEigenSolver below ccexBlas::eigMin, LAPACK zheevd (lower triangle, as Eigen reads) from it on
class CcexEigensolver {
public:
    void compute(const Eigen::MatrixXcd& H);
    // A42: an eigensystem given from outside (the small-cluster Jacobi), held where zheevd's result is held
    void setGiven(int n, const double* w, const std::complex<double>* V){
        lapack_ = true; info_ = Eigen::Success;
        V_ = Eigen::Map<const Eigen::MatrixXcd>(V,n,n); w_ = Eigen::Map<const Eigen::VectorXd>(w,n);
    }
    const Eigen::MatrixXcd& eigenvectors() const { return lapack_ ? V_ : es_.eigenvectors(); }
    const Eigen::VectorXd& eigenvalues() const { return lapack_ ? w_ : es_.eigenvalues(); }
    Eigen::ComputationInfo info() const { return lapack_ ? info_ : es_.info(); }
private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es_;
    bool lapack_ = false;
    Eigen::MatrixXcd V_;
    Eigen::VectorXd w_;
    Eigen::ComputationInfo info_ = Eigen::Success;
    int wsDim_ = -1;                                   // the dimension the workspace below was sized for
    std::vector<std::complex<double>> work_;
    std::vector<double> rwork_;
    std::vector<int> iwork_;
};
#endif
