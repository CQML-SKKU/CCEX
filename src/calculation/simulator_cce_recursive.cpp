#include "../../include/calculation/simulator.h"
#include <cmath>

// A18: the ensemble cCCE (nstate 0) with evolution_update = recursive. The two ensemble branches of calCoherenceCceInto
// - CPMG / Hahn / Ramsey (tau = t/(2 npulse), Hahn t/2, Ramsey t; the sequence is not read) and the Manual sequence
// (interval by interval, tau = t f) - with the same products, but the phases exp(-i E tau) of each tau computed by sin
// and cos at steps 0, 16, 32, ... only, and in between the previous step's times exp(-i E dtau) (dtau = how much tau
// grows per step: deltat/2, deltat, deltat/(2 npulse) and twice that, or f deltat), as evolveSingleState does (A16 b).
// A file of its own: written into calCoherenceCceInto (as a template of both modes) it changed how g++ compiled the
// exact branches there (the sign of a zero imaginary part in check_numerics hbn_cce2_hf3_qd0).
namespace {
constexpr int phaseAnchor = 16;
// exp(-i E tau) as the Hahn and Manual branches form it (sin and cos of the real argument)
void realPhases(const Eigen::VectorXcd& E, double tau, Eigen::VectorXcd& out){
    out.resize(E.size());
    for (int k=0; k<E.size(); k++){
        const double a = E(k).real()*tau;
        out(k) = doublec(std::cos(a), -std::sin(a));
    }
}
// exp(-i x), x = E tau, as the Ramsey and CPMG branches form it
void complexPhases(const Eigen::VectorXcd& x, Eigen::VectorXcd& out){
    out = x.array().cos() + doublec(0,-1)*x.array().sin();
}
// The phases of one tau at step i: by sin and cos every phaseAnchor steps (and the factors for one step at step 0),
// otherwise the previous step's times those factors
struct Carried {
    Eigen::VectorXcd up, down, stepUp, stepDown;
    void advance(int i, double tau, double dtau, bool complexForm, const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown){
        auto form = [&](double t, Eigen::VectorXcd& u, Eigen::VectorXcd& d){
            if (complexForm){ complexPhases(valUp*t,u); complexPhases(valDown*t,d); }
            else            { realPhases(valUp,t,u);    realPhases(valDown,t,d); }
        };
        if (i % phaseAnchor == 0){
            form(tau,up,down);
            if (i == 0) form(dtau,stepUp,stepDown);
        } else {
            up.array() *= stepUp.array(); down.array() *= stepDown.array();
        }
    }
};
}

// X8: the ensemble Ramsey reads only Tr[M Dd M^dagger Du^*] = sum_ij |M_ij|^2 Dd_j conj(Du_i) = Du^dagger W Dd, W = |M|^2
// (once per cluster): D^2 per step instead of a D^3 product
namespace {
void ramseyWeights(const MatrixXcd& M, MatrixXcd& W){ W = M.cwiseAbs2().cast<doublec>(); }
doublec ramseyTrace(const Eigen::VectorXcd& up, const MatrixXcd& W, const Eigen::VectorXcd& down){
    return (up.adjoint() * W * down)(0,0);
}
}

// X8: the exact ensemble Ramsey of calCoherenceCceInto (phases by cos and sin at every step), with the trace above; as
// there, the mean over the bath, conjugated (an even number of pulses, 2026-09-19)
void evolveEnsembleRamseyCce(const MatrixXcd& M, const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown,
                             int nstep, double deltat, MatrixXcd* result){
    const int bdim = int(M.rows());
    MatrixXcd W;
    ramseyWeights(M,W);
    Eigen::VectorXcd U_up, U_down;
    double tfree = 0.0;
    for (int i = 0; i < nstep; i++){
        complexPhases(valUp*tfree,U_up);
        complexPhases(valDown*tfree,U_down);
        doublec overlap = ramseyTrace(U_up,W,U_down);
        overlap /= doublec(bdim*1.0,0.0);
        result[i] = MatrixXcd::Zero(1,1);
        result[i](0,0) = std::conj(overlap);
        tfree += deltat;
    }
}

void evolveEnsembleRecursiveCce(const MatrixXcd& M, const MatrixXcd& M_up, const MatrixXcd& M_down,
                                const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown, Pulse* pulse,
                                bool manualIntervals, int nstep, double deltat, MatrixXcd* result){

    const int bdim = int(M.rows());
    const int npulse = Pulse_getNpulse(pulse);
    double** sequence = Pulse_getSequence(pulse);
    int* sequence_indices = Pulse_getSequenceIndices(pulse);
    double tfree = 0.0;

    if (!manualIntervals){
        Carried one, two;                    // tau and, for CPMG, 2 tau
        MatrixXcd X(bdim,bdim), L_mat, sum_L_a, sum_L_b;
        MatrixXcd L_a_even, L_b_even, L_a_odd, L_b_odd, L_a2_even, L_b2_even, L_a2_odd, L_b2_odd;
        Eigen::MatrixXd X2(bdim,bdim);
        MatrixXcd ramseyW;                   // X8: W = |M|^2 for Ramsey (ramseyTrace)
        if (npulse == 0) ramseyWeights(M,ramseyW);
        doublec overlap;
        for (int i = 0; i < nstep; i++){
            if (npulse == 1){
                // Hahn echo: Tr[X Du X^dagger Du^*] with X = M Dd M^dagger
                one.advance(i, tfree/2.0, deltat/2.0, false, valUp, valDown);
#ifdef CCEX_BLAS_DIRECT
                if (bdim >= ccexBlas::gemmMin) ccexBlas::hahnProduct(M, one.down, L_mat, X);   // A46: zgemm
                else
#endif
                X.noalias() = (M.array().rowwise() * one.down.transpose().array()).matrix() * M.adjoint();
                X2 = X.cwiseAbs2();
                overlap = (one.up.adjoint() * X2.cast<doublec>() * one.up)(0,0);
            } else if (npulse == 0){
                // Ramsey: Tr[M U_down M^dagger U_up^dagger] = U_up^dagger W U_down
                one.advance(i, tfree, deltat, true, valUp, valDown);
                overlap = ramseyTrace(one.up,ramseyW,one.down);
            } else {
                // CPMG: tau at both ends, 2 tau in between
                const double tau = tfree/(2*npulse), dtau = deltat/(2*npulse);
                one.advance(i, tau, dtau, true, valUp, valDown);
                two.advance(i, 2*tau, 2*dtau, true, valUp, valDown);
                const Eigen::VectorXcd &U_up = one.up, &U_down = one.down, &U_up2 = two.up, &U_down2 = two.down;
                L_a_even = M.array().rowwise()           * U_down.transpose().array();   // M * U_down
                L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array();       // M^dagger * U_up^dagger
                if (npulse % 2 != 0){
                    L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().array();  // M^dagger * U_up
                    L_b_odd = M.array().rowwise()           * U_down.adjoint().array();  // M * U_down^dagger
                }
                if (npulse >= 3){
                    L_a2_even = M.array().rowwise()           * U_down2.transpose().array();
                    L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();
                }
                L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().array();
                L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array();

                sum_L_a = L_a_even;
                sum_L_b = L_b_even;
                int loopN = 1;
                for (int nn = 1; nn <= npulse; nn++){
                    if (nn == npulse){
                        if (nn % 2 == 0){ sum_L_a = sum_L_a * L_a_even; sum_L_b = L_b_even * sum_L_b; }
                        else            { sum_L_a = sum_L_a * L_a_odd;  sum_L_b = L_b_odd * sum_L_b;  }
                    } else {
                        if (loopN % 2 != 0){ sum_L_a = sum_L_a * L_a2_odd;  sum_L_b = L_b2_odd * sum_L_b;  }
                        else               { sum_L_a = sum_L_a * L_a2_even; sum_L_b = L_b2_even * sum_L_b; }
                        loopN++;
                    }
                }
                overlap = (sum_L_a.array() * sum_L_b.transpose().array()).sum();
            }
            // as calCoherenceCceInto: the mean over the bath, conjugated for an even number of pulses (2026-09-19)
            overlap /= doublec(bdim*1.0,0.0);
            if (npulse%2==0) overlap = std::conj(overlap);
            result[i] = MatrixXcd::Zero(1,1);
            result[i](0,0) = overlap;
            tfree += deltat;
        }
        return;
    }

    // Manual sequence: the interval-by-interval chain of calCoherenceCceInto, Tr[UdownDagger Uup]/d
    std::vector<Carried> phases(npulse+1);
    MatrixXcd UdownDagger(bdim,bdim), Uup(bdim,bdim), Uscratch(bdim,bdim);
    std::vector<MatrixXcd> Mdagger_UupTauDagger(npulse+1), M_UdownTauDagger(npulse+1), Mdagger_UupTau(npulse+1), M_UdownTau(npulse+1);
    for (int i = 0; i < nstep; i++){
        for (int ipulse=0; ipulse<npulse+1; ipulse++){
            const int same = sequence_indices[ipulse];
            if (same == ipulse){
                Carried& p = phases[ipulse];
                p.advance(i, tfree*sequence[ipulse][2], deltat*sequence[ipulse][2], false, valUp, valDown);
                Mdagger_UupTauDagger[ipulse] = (M.adjoint().array().rowwise() * p.up.adjoint().array()).matrix();
                M_UdownTauDagger[ipulse]     = (M.array().rowwise()           * p.down.adjoint().array()).matrix();
                Mdagger_UupTau[ipulse]       = (M.adjoint().array().rowwise() * p.up.transpose().array()).matrix();
                M_UdownTau[ipulse]           = (M.array().rowwise()           * p.down.transpose().array()).matrix();
            } else {
                Mdagger_UupTauDagger[ipulse] = Mdagger_UupTauDagger[same];
                M_UdownTauDagger[ipulse]     = M_UdownTauDagger[same];
                Mdagger_UupTau[ipulse]       = Mdagger_UupTau[same];
                M_UdownTau[ipulse]           = M_UdownTau[same];
            }
            if (ipulse == 0){
                const Carried& p = phases[0];
                UdownDagger = (M_down.array().rowwise() * p.down.adjoint().array()).matrix();       // M_down U_down^dagger
                Uup.noalias() = (M.adjoint().array().rowwise() * p.up.transpose().array()).matrix() * M_up.adjoint();
            } else if (ipulse%2 == 1){
                Uscratch.noalias() = UdownDagger * Mdagger_UupTauDagger[ipulse]; UdownDagger.swap(Uscratch);
                Uscratch.noalias() = M_UdownTau[ipulse] * Uup;                   Uup.swap(Uscratch);
            } else {
                Uscratch.noalias() = UdownDagger * M_UdownTauDagger[ipulse];     UdownDagger.swap(Uscratch);
                Uscratch.noalias() = Mdagger_UupTau[ipulse] * Uup;               Uup.swap(Uscratch);
            }
        }
        const doublec overlap = (UdownDagger.array() * Uup.transpose().array()).sum();
        result[i] = MatrixXcd::Zero(1,1);
        result[i](0,0) = overlap / doublec(bdim*1.0,0.0);
        tfree += deltat;
    }
}
