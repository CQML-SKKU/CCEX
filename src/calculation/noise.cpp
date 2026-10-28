/**
 * @file noise.cpp
 * @brief method noise: the Gaussian coherence W(T) = exp(-chi(T)) of the qubit's phase noise (semi-cce, ported).
 *
 * The phase is phi = int_0^T f(t) dw(t) dt, f = +-1 flipped by every ideal pi pulse, dw the qubit frequency noise (rad/ms),
 * and chi = <phi^2> / 2. With Y(w) = int_0^T f(t) exp(i w t) dt:
 *  - a model spectrum S(f) of the frequency noise (kHz^2/kHz, two-sided): chi = 4 pi^2 int_0^inf S(f) |Y(2 pi f)|^2 df;
 *  - a cluster: dw = Delta m_s sum_i A_zz,i I_z,i, the bath under H = H_bath + p sum_i (A_z,i . I_i) with
 *    p = (|m_alpha| + |m_beta|) / 2 (semi-cce's projection; user decision 2026-10-28), infinite temperature:
 *    chi_C = (1 / 2D) sum_ab |<a|dw|b>|^2 |Y(E_a - E_b)|^2. For a Hahn echo this is semi-cce's mode_sum.
 * chi adds over the CCE terms; the cCCE product of the clusters' exp(-chi_C) is exp(-(that sum)), so the clusters, their
 * distribution and the product are the cCCE ones (simulator.cpp). A cluster hands over -chi_C, and the product takes
 * exp(-chi_C) in its scaled form: a FID chi_C passes 745 within microseconds, where exp underflows to 0 and 0^-1 is inf.
 */
#include "../../include/calculation/noise.h"
#include "../../include/calculation/hamiltonian.h"
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

namespace {
// the interval lengths as fractions of T; interval j has the sign (-1)^j (every pulse an X or Y pi flip,
// validatePulsesForMethod)
std::vector<double> intervalFractions(Pulse* pls){
    std::vector<double> frac(Pulse_getNpulse(pls)+1);
    double** sequence = Pulse_getSequence(pls);
    for (size_t j=0; j<frac.size(); ++j) frac[j] = sequence[j][2];
    return frac;
}

// Y(w) = sum_j s_j L_j sinc(w L_j / 2) exp(i w (centre of j)): no division, exact at w = 0
std::complex<double> phaseIntegral(double w, double T, const std::vector<double>& frac){
    std::complex<double> y = 0;
    double start = 0, sign = 1;
    for (double fj : frac){
        const double L = fj*T, x = 0.5*w*L;
        const double sinc = std::fabs(x) < 1e-4 ? 1 - x*x/6 : std::sin(x)/x;
        y += sign*L*sinc*std::polar(1.0, w*(start + 0.5*L));
        start += L;
        sign = -sign;
    }
    return y;
}

// OU, C(t) = (2 pi sigma)^2 exp(-|t| / tau): chi = 2 pi^2 sigma^2 sum_ij s_i s_j int_i int_j exp(-|t1 - t2| / tau), exact
double ouChi(double sigma, double tau, double T, const std::vector<double>& frac){
    const size_t n = frac.size();
    double sum = 0, start_i = 0;
    for (size_t i=0; i<n; ++i){
        const double Li = frac[i]*T, x = Li/tau;
        // int_i int_i = 2 tau^2 (x + expm1(-x)); its series below 1e-3, where the two terms cancel
        sum += 2*tau*tau*(x < 1e-3 ? x*x*(0.5 - x*(1.0/6 - x*(1.0/24 - x/120))) : x + std::expm1(-x));
        double start_j = start_i + Li;
        for (size_t j=i+1; j<n; ++j){
            const double Lj = frac[j]*T, sign = ((j-i) % 2) ? -1 : 1;
            sum += 2*sign*tau*tau*std::exp(-(start_j - start_i - Li)/tau)*std::expm1(-Li/tau)*std::expm1(-Lj/tau);
            start_j += Lj;
        }
        start_i += Li;
    }
    return 2*M_PI*M_PI*sigma*sigma*sum;
}

// S(f) = A (f / 1 kHz)^-beta on [f_min, f_max]: 8-point Gauss-Legendre on panels no wider than 5 % of f (the spectrum)
// and 1/(8T) (the oscillation of |Y|^2). shortcut: the panels number about 8 T f_max, slow for f_max T >> 1e6
double powerChi(double A, double beta, double fmin, double fmax, double T, const std::vector<double>& frac){
    if (T <= 0) return 0;
    static const double x[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static const double w[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    double sum = 0;
    for (double a = fmin; a < fmax;){
        const double b = std::min(fmax, a + std::min(0.05*a, 0.125/T)), mid = 0.5*(a+b), half = 0.5*(b-a);
        for (int k=0; k<4; ++k)
            for (double f : {mid - half*x[k], mid + half*x[k]})
                sum += w[k]*half*A*std::pow(f,-beta)*std::norm(phaseIntegral(2*M_PI*f,T,frac));
        a = b;
    }
    return 4*M_PI*M_PI*sum;
}
}

double noiseModelChi(Config* cnf, Pulse* pls, double T){
    const std::vector<double> frac = intervalFractions(pls);
    if (!strcmp(cnf->noise_model,"OU")) return ouChi(cnf->noise_sigma,cnf->noise_tau_c,T,frac);
    return powerChi(cnf->noise_amplitude,cnf->noise_beta,cnf->noise_f_min,cnf->noise_f_max,T,frac);
}

void noiseModelRun(Config* cnf, Pulse* pls, Output* op){
    if (rank != 0) return;
    const int nstep = Config_getNstep(cnf);
    const double deltat = Config_getDeltat(cnf);
    std::vector<MatrixXcd> result(nstep);
    for (int i=0; i<nstep; ++i) result[i] = MatrixXcd::Constant(1,1,std::exp(-noiseModelChi(cnf,pls,i*deltat)));
    Output_save(op,result.data(),result.data(),nstep,deltat,0);
}

void calCoherenceNoiseInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result){
    const int nspin = BathArray_getNspin(ba), bdim = BathArray_dim(ba);
    auto cleanup = [nspin](MatrixXcd** p){ if(p){ for(int i=0;i<nspin;++i) delete[] p[i]; delete[] p; } };
    std::unique_ptr<MatrixXcd*,decltype(cleanup)> sigmas(BathArray_PauliOperators(ba),cleanup);

    const double msa = findZbasisSubLevel(qa->psia), msb = findZbasisSubLevel(qa->psib);
    const double p = 0.5*(std::fabs(msa) + std::fabs(msb));
    MatrixXcd H = HamilBath(ba,sigmas.get(),cnf);
    MatrixXcd dw = MatrixXcd::Zero(bdim,bdim);
    for (int ib=0; ib<nspin; ++ib){
        const MatrixXcd tensor = BathArray_getBath_i_hypf_j(ba,ib,0);
        MatrixXcd v = MatrixXcd::Zero(1,3);
        for (int c=0; c<3; ++c) v(0,c) = p*tensor(2,c).real();
        addExpandedHamiltonian(H,sigmas.get(),calHamiltonianSingleInt(v,sigmas.get()[ib]),nspin,ib);
        v.setZero();
        v(0,2) = (msa - msb)*tensor(2,2).real();
        addExpandedHamiltonian(dw,sigmas.get(),calHamiltonianSingleInt(v,sigmas.get()[ib]),nspin,ib);
    }
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
    const MatrixXcd& U = es.eigenvectors();
    const Eigen::VectorXd& E = es.eigenvalues();
    const Eigen::MatrixXd M2 = (U.adjoint()*dw*U).cwiseAbs2();

    const std::vector<double> frac = intervalFractions(pls);
    const int nstep = Config_getNstep(cnf);
    const double deltat = Config_getDeltat(cnf);
    for (int istep=0; istep<nstep; ++istep){
        const double T = istep*deltat;
        double chi = 0.5*M2.diagonal().sum()*std::norm(phaseIntegral(0,T,frac));   // the static part (FID only)
        for (int a=0; a<bdim; ++a)
            for (int b=a+1; b<bdim; ++b) chi += M2(a,b)*std::norm(phaseIntegral(E(a)-E(b),T,frac));
        result[istep] = MatrixXcd::Constant(1,1,-chi/bdim);   // -chi_C, not exp(-chi_C): see noiseFactor (simulator.cpp)
    }
}
