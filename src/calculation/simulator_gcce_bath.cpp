#include "../../include/calculation/simulator_gcce_bath.h"
#include "../../include/calculation/hamiltonian.h"
#include "../../include/physical/pulse_bath.h"
#include <cmath>
#include <string>
#include <vector>

namespace {
// exp(-i theta sigma / 2) in the {balphams, bbetams} subspace of spin bspin (adding_MECCE cal_Upulse_when_pulseiter_isFalse)
MatrixXcd bathRotation(const Pulse* pls, int row) {
    MatrixXcd* s = getGeneralPauliOperators(getSpinor(pls->bspin, pls->balphams), getSpinor(pls->bspin, pls->bbetams));
    const char axis = pls->tot_bpulse_axes[row];
    const int k = axis == 'X' ? 1 : axis == 'Y' ? 2 : axis == 'Z' ? 3 : 0;   // sigma[0]=I, [1]=X, [2]=Y, [3]=Z
    const double theta = pls->tot_bpulse_angles[row] * M_PI / 180.0;
    MatrixXcd U = (doublec(0.0, -1.0) * s[k] * theta / 2.0).exp();
    delete[] s;
    return U;
}
// the rule of adding_MECCE src/simulator_gcce.cpp:546-553: the name, and any pump tone within the window of the detuning
bool named(const BathSpin* bs, const Pulse* pls) { return !strcasecmp(bs->name, pls->bpulse_defect); }
bool inWindow(const BathSpin* bs, const Pulse* pls) {
    for (int k = 0; k < pls->n_bpulse_shift; ++k)
        if (std::abs(pls->bpulse_energy_shifts[k] - bs->detuning) < kResonanceWindowRadkHz) return true;
    return false;
}
bool driven(const BathSpin* bs, const Pulse* pls) { return named(bs, pls) && inWindow(bs, pls); }
}

MatrixXcd bathPulseOperator(BathArray* ba, const Pulse* pls, int row) {
    const int nspin = ba ? BathArray_getNspin(ba) : 0;
    MatrixXcd Ub = MatrixXcd::Identity(1, 1);
    for (int j = 0; j < nspin; ++j) {
        const float S = BathArray_getBath_i_spin(ba, j);
        MatrixXcd Uj;
        if (driven(ba->bath[j], pls)) {
            if (std::fabs(S - pls->bspin) > 1e-6) {   // validateBathPulseTargets stops this before any cluster
                fprintf(stderr, "Error: bath pulse on %s: pulse.bath.spin %g differs from the spin %g of that bath spin\n",
                        ba->bath[j]->name, pls->bspin, (double)S);
                ccexFatal(1);
            }
            Uj = bathRotation(pls, row);
        } else {
            Uj = MatrixXcd::Identity((int)(2 * S + 1), (int)(2 * S + 1));
        }
        Ub = kron(Ub, Uj);
    }
    return Ub;
}

void applyBathPulsesGcce(GcceWork* w, MatrixXcd* Upulses, BathArray* ba, const Pulse* pls) {
    if (!Pulse_hasBathPulses(pls) || !ba || BathArray_getNspin(ba) == 0) return;
    bool changed = false;
    for (int i = 0; i < pls->npulse; ++i) {
        if (pls->tot_bpulse_axes[i] == 'I' || std::fabs(pls->tot_bpulse_angles[i]) < 1e-12) continue;
        MatrixXcd Ub = bathPulseOperator(ba, pls, i);
        if (Ub.isIdentity(0.0)) continue;   // nothing driven in this cluster: keep kron(U, 1) and the block path
        w->UpulsesExpanded[i] = kron(Upulses[i], Ub);
        changed = true;
    }
    if (changed) w->pulseBdim = 0;   // applyPulseGcce's block path assumes kron(U, 1)
}

void validateBathPulseTargets(BathArray* ba, DefectArray* dfa, const Pulse* pls) {
    if (!Pulse_hasBathPulses(pls)) return;
    auto refuse = [](const char* msg) { if (rank == 0) fprintf(stderr, "Error: %s\n", msg); MPI_Abort(MPI_COMM_WORLD, 1); };
    // every spin a cluster can hold: the bath file's, and the Defect sub-spins createBathArray adds to them
    std::vector<const BathSpin*> spins;
    const int nspin = BathArray_getNspin(ba);
    for (int j = 0; j < nspin; ++j) spins.push_back(ba->bath[j]);
    if (DefectArray_getNdefect(dfa) > 0)
        for (int j = 0; j < nspin; ++j)
            for (int k = 0; k < DefectArray_getNaddspins_i(dfa, j); ++k) spins.push_back(DefectArray_getSubbath_i_isp(dfa, j, k));
    char msg[512];
    int nnamed = 0, nhit = 0;
    for (const BathSpin* bs : spins) {
        if (!named(bs, pls)) continue;
        ++nnamed; nhit += inWindow(bs, pls);
        if (std::fabs(bs->spin - pls->bspin) > 1e-6) {   // the rotation would not fit the spin's space (kron size)
            snprintf(msg, sizeof msg, "pulse.bath.spin %g differs from the spin %g of the bath spins named %s "
                     "(pulse.bath.defect)", pls->bspin, (double)bs->spin, pls->bpulse_defect);
            refuse(msg);
        }
    }
    if (nnamed == 0) {   // a typo: the run would be a plain echo
        snprintf(msg, sizeof msg, "pulse.bath.defect \"%s\" names no bath spin of this run (the names of the bath file "
                 "and of the Defect sub-spins; the comparison ignores case)", pls->bpulse_defect);
        refuse(msg);
    }
    // the two levels the pulse rotates between (getSpinor): |m| <= S, S - m whole, and not the same level
    const double S = pls->bspin;
    auto level = [S](double m) { return std::fabs(m) <= S + 1e-6 && std::fabs(std::remainder(S - m, 1.0)) < 1e-6; };
    if (std::fabs(pls->balphams - pls->bbetams) < 1e-6 || !level(pls->balphams) || !level(pls->bbetams)) {
        snprintf(msg, sizeof msg, "pulse.bath.alphams %g and betams %g must be two different levels of spin %g "
                 "(m = -S, -S+1, ..., S)", pls->balphams, pls->bbetams, S);
        refuse(msg);
    }
    if (nhit == 0 && rank == 0) {   // spins of that name, none in reach of a tone: the pulses act as identity
        std::string tones;
        for (int k = 0; k < pls->n_bpulse_shift; ++k) {
            char tone[64];
            snprintf(tone, sizeof tone, "%s%.10g", k ? ", " : "", RADKHZ_TO_MHZ(pls->bpulse_energy_shifts[k]));
            tones += tone;
        }
        fflush(stdout);   // the warning on a line of its own in a log that holds both
        fprintf(stderr, "Warning: pulse.bath drives no spin: no spin named %s has a detuning within %g rad*kHz of the "
                "tones [%s] MHz; the bath pulses act as identity\n", pls->bpulse_defect, kResonanceWindowRadkHz,
                tones.c_str());
    }
}
