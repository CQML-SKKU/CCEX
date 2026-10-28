/**
 * @file noise.cpp
 * @brief method noise: the Gaussian coherence W(T) = exp(-chi(T)) of the qubit's phase noise (semi-cce, ported).
 *
 * The phase is phi = int_0^T f(t) dw(t) dt, f = +-1 flipped by every ideal pi pulse, dw the qubit frequency noise (rad/ms),
 * and chi = <phi^2> / 2. With Y(w) = int_0^T f(t) exp(i w t) dt:
 *  - a model spectrum S(f) of the frequency noise (kHz^2/kHz, two-sided): chi = 4 pi^2 int_0^inf S(f) |Y(2 pi f)|^2 df;
 *  - a cluster: dw = Delta m_s sum_i A_zz,i I_z,i, the secular hyperfine term (semi-cce's), with
 *    hamiltonian.noise.pseudo_secular also A_zx,i I_x,i + A_zy,i I_y,i; the bath under H = H_bath + p sum_i (A_z,i . I_i) with
 *    p = (|m_alpha| + |m_beta|) / 2 (semi-cce's projection; user decision 2026-10-28), infinite temperature:
 *    chi_C = (1 / 2D) sum_ab |<a|dw|b>|^2 |Y(E_a - E_b)|^2. For a Hahn echo this is semi-cce's mode_sum.
 * hamiltonian.noise.compute chooses what a cluster gives (noiseCluster) and what the run writes (noiseFinish):
 *  - "coherence": -chi_C at the run's times; the coherence exp(-chi).
 *  - "spectrum" (semi-cce's psd): the correlation C_C(t) of B = sum_i A_zz,i I_z,i (kHz) at infinite temperature on the
 *    grid of hamiltonian.noise.fft, C = (1/D) sum_ab |<a|B|b>|^2 cos(w_ab t); then S(f) by FFT of C (semi-cce's
 *    psd_from_corr, its window, pad_factor, subtract_mean) and chi = Delta m_s^2 2 int S(f) 2 pi^2 |Y(2 pi f)|^2 df by the
 *    trapezoid rule (coherence_from_psd, clamp_negative): C(t), S(f) and the coherence are written.
 *  - "correlation" (PyCCE's CCENoise): C_C(t) = Re Tr[rho {B(t), B} / 2] at the run's times, rho the cluster's pure
 *    initial state with configuration.s_range (and the mean field of the other spins, as PyCCE's external_spins_field)
 *    or 1/D: C(t) alone is written.
 * Each adds over the CCE terms: a cluster's values times its power are summed (simulator.cpp), not multiplied, and the
 * sum is turned into the outputs once, so no exp(-chi_C) of one cluster can underflow (a FID chi_C passes 745 within
 * microseconds).
 */
#include "../../include/calculation/noise.h"
#include "../../include/calculation/hamiltonian.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>
#include <unsupported/Eigen/FFT>

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

namespace {
// the cluster's Hamiltonian with p = (|m_alpha| + |m_beta|) / 2 (HamilBath has the mean field unless it is off), and
// B = sum_i A_z,i . I_i (rad/ms) with A_zz alone unless pseudo_secular
void clusterOperators(QubitArray* qa, BathArray* ba, Config* cnf, MatrixXcd& H, MatrixXcd& B){
    const int nspin = BathArray_getNspin(ba), bdim = BathArray_dim(ba);
    auto cleanup = [nspin](MatrixXcd** p){ if(p){ for(int i=0;i<nspin;++i) delete[] p[i]; delete[] p; } };
    std::unique_ptr<MatrixXcd*,decltype(cleanup)> sigmas(BathArray_PauliOperators(ba),cleanup);
    const double msa = findZbasisSubLevel(qa->psia), msb = findZbasisSubLevel(qa->psib);
    const double p = 0.5*(std::fabs(msa) + std::fabs(msb));
    H = HamilBath(ba,sigmas.get(),cnf);
    B = MatrixXcd::Zero(bdim,bdim);
    for (int ib=0; ib<nspin; ++ib){
        const MatrixXcd tensor = BathArray_getBath_i_hypf_j(ba,ib,0);
        MatrixXcd v = MatrixXcd::Zero(1,3);
        for (int c=0; c<3; ++c) v(0,c) = p*tensor(2,c).real();
        addExpandedHamiltonian(H,sigmas.get(),calHamiltonianSingleInt(v,sigmas.get()[ib]),nspin,ib);
        v.setZero();
        for (int c = cnf->noise_pseudo_secular ? 0 : 2; c<3; ++c) v(0,c) = tensor(2,c).real();
        addExpandedHamiltonian(B,sigmas.get(),calHamiltonianSingleInt(v,sigmas.get()[ib]),nspin,ib);
    }
}

double deltaMs(QubitArray* qa){ return findZbasisSubLevel(qa->psia) - findZbasisSubLevel(qa->psib); }

std::string outputName(Output* op, int istate, const char* suffix){
    std::string name = Output_getOutfile(op);
    if (istate) name += "_state" + std::to_string(istate);
    return name + suffix;
}

void writeColumns(const std::string& name, const char* header, const std::vector<double>& x, const std::vector<double>& y){
    FILE* fp = fopen(name.c_str(),"w");
    if (!fp){ fprintf(stderr,"Error: cannot open file %s\n",name.c_str()); ccexFatal(1); }
    fprintf(fp,"# %s\n",header);
    for (size_t i=0; i<x.size(); ++i) fprintf(fp,"%.12e\t%+.12e\n",x[i],y[i]);
    if (fclose(fp)){ fprintf(stderr,"Error: cannot write file %s\n",name.c_str()); ccexFatal(1); }
}

// semi-cce psd_from_corr: S(f) = dt Re FFT of C made even (length M = 2N - 1, or pad_factor times that, odd), the
// padding at C's tail mean (0 with subtract_mean); f = k / (M dt), k = 0 .. M/2
void spectrumOf(std::vector<double> C, double dt, Config* cnf, std::vector<double>& f, std::vector<double>& S){
    const int N = int(C.size());
    const int tail = std::max(1, N/5);
    double Cinf = 0;
    for (int n=N-tail; n<N; ++n) Cinf += C[n];
    Cinf /= tail;
    if (cnf->noise_fft_subtract_mean) for (double& c : C) c -= Cinf;
    if (!strcmp(cnf->noise_fft_window,"hann"))
        for (int n=0; n<N; ++n) C[n] *= 0.5*(1 + std::cos(M_PI*n/(N-1)));
    else if (!strcmp(cnf->noise_fft_window,"tukey")){   // alpha 0.5: the last half tapered
        const int ts = int(N*0.5);
        for (int n=ts; n<N; ++n) C[n] *= 0.5*(1 + std::cos(M_PI*(n-ts)/(N-ts)));
    }else if (!strcmp(cnf->noise_fft_window,"exp"))
        for (int n=0; n<N; ++n) C[n] *= std::exp(-3.0*n/(N-1));
    const int Mmin = 2*N - 1;
    int Mpad = Mmin*cnf->noise_fft_pad_factor;
    if (Mpad % 2 == 0) ++Mpad;
    const int M = std::max(Mmin,Mpad);
    std::vector<double> x(M, cnf->noise_fft_subtract_mean ? 0 : Cinf);
    x[0] = C[0];
    for (int n=1; n<N; ++n) x[n] = x[M-n] = C[n];
    Eigen::FFT<double> fft;
    std::vector<std::complex<double>> X;
    fft.fwd(X,x);
    f.resize(M/2+1);
    S.resize(M/2+1);
    for (int k=0; k<=M/2; ++k){ f[k] = k/(M*dt); S[k] = dt*X[k].real(); }
}
}

int noiseLength(Config* cnf){
    return !strcmp(cnf->noise_compute,"spectrum") ? cnf->noise_fft_nstep + 1 : Config_getNstep(cnf);
}

void noiseCluster(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, std::vector<double>& value){
    const int bdim = BathArray_dim(ba);
    MatrixXcd H, B;
    clusterOperators(qa,ba,cnf,H,B);
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
    const MatrixXcd& U = es.eigenvectors();
    const Eigen::VectorXd& E = es.eigenvalues();
    const MatrixXcd M = U.adjoint()*B*U;
    value.assign(noiseLength(cnf),0.0);

    if (!strcmp(cnf->noise_compute,"coherence")){
        const Eigen::MatrixXd M2 = M.cwiseAbs2()*(deltaMs(qa)*deltaMs(qa));
        const std::vector<double> frac = intervalFractions(pls);
        const double deltat = Config_getDeltat(cnf);
        for (size_t i=0; i<value.size(); ++i){
            const double T = i*deltat;
            double chi = 0.5*M2.diagonal().sum()*std::norm(phaseIntegral(0,T,frac));   // the static part (FID only)
            for (int a=0; a<bdim; ++a)
                for (int b=a+1; b<bdim; ++b) chi += M2(a,b)*std::norm(phaseIntegral(E(a)-E(b),T,frac));
            value[i] = -chi/bdim;
        }
        return;
    }
    // C(t) = Re sum_ab G_ab exp(i w_ab t), G_ab = M_ab ((M R)_ba + (R M)_ba) / 2 in kHz^2 (B / 2 pi), R the initial state
    const bool pure = !strcmp(cnf->noise_compute,"correlation") && Config_getNstate(cnf) > 0;
    MatrixXcd G;
    if (pure){
        const MatrixXcd psi = U.adjoint()*BathArray_Psi0(ba);
        const MatrixXcd R = psi*psi.adjoint();
        const MatrixXcd MR = M*R, RM = R*M;
        G = M.cwiseProduct((MR + RM).transpose())*(0.5/(4*M_PI*M_PI));
    }else
        G = M.cwiseAbs2().cast<std::complex<double>>()*(1.0/(4*M_PI*M_PI*bdim));
    const double dt = !strcmp(cnf->noise_compute,"spectrum") ? cnf->noise_fft_time_max/cnf->noise_fft_nstep : Config_getDeltat(cnf);
    for (size_t i=0; i<value.size(); ++i){
        const double t = i*dt;
        double c = 0;
        for (int a=0; a<bdim; ++a)
            for (int b=0; b<bdim; ++b) c += (G(a,b)*std::polar(1.0,(E(a)-E(b))*t)).real();
        value[i] = c;
    }
}

bool noiseFinish(QubitArray* qa, Config* cnf, Pulse* pls, Output* op, int istate,
                 std::vector<double>& with, std::vector<double>& without, MatrixXcd* wD, MatrixXcd* nD){
    const int n = int(with.size());
    std::vector<double> sumWith(n), sumWithout(n);
    MPI_Reduce(with.data(),sumWith.data(),n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(without.data(),sumWithout.data(),n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank != 0) return true;
    const int nstep = Config_getNstep(cnf);
    const double deltat = Config_getDeltat(cnf);
    if (!strcmp(cnf->noise_compute,"coherence")){
        for (int i=0; i<nstep; ++i){
            wD[i] = MatrixXcd::Constant(1,1,std::exp(sumWith[i]));
            nD[i] = MatrixXcd::Constant(1,1,std::exp(sumWithout[i]));
        }
        return true;
    }
    const bool spectrum = !strcmp(cnf->noise_compute,"spectrum");
    const double dt = spectrum ? cnf->noise_fft_time_max/cnf->noise_fft_nstep : deltat;
    std::vector<double> t(n);
    for (int i=0; i<n; ++i) t[i] = i*dt;
    writeColumns(outputName(op,istate,"_noise_correlation"),
                 "t (ms)  C(t) (kHz^2): the correlation of sum_i A_zz,i I_z,i (with pseudo_secular A_z,i . I_i), Re <{B(t), B}>/2",
                 t,sumWith);
    if (!spectrum) return false;
    std::vector<double> f, S;
    spectrumOf(sumWith,dt,cnf,f,S);
    writeColumns(outputName(op,istate,"_noise_spectrum"),"f (kHz)  S(f) (kHz^2/kHz, two-sided): dt Re FFT of C(t)",f,S);
    if (cnf->noise_fft_clamp_negative) for (double& s : S) s = std::max(s,0.0);
    const std::vector<double> frac = intervalFractions(pls);
    const double dm2 = deltaMs(qa)*deltaMs(qa);
    for (int i=0; i<nstep; ++i){
        const double T = i*deltat;
        double chi = 0, prev = 0;
        for (size_t k=0; k<f.size(); ++k){   // trapezoid of S(f) F(2 pi f T) / f^2, F / f^2 = 2 pi^2 |Y(2 pi f)|^2
            const double g = S[k]*2*M_PI*M_PI*std::norm(phaseIntegral(2*M_PI*f[k],T,frac));
            if (k) chi += 0.5*(g + prev)*(f[k] - f[k-1]);
            prev = g;
        }
        wD[i] = nD[i] = MatrixXcd::Constant(1,1,std::exp(-2*dm2*chi));
    }
    return true;
}
