#ifndef __CCEX_SIMULATOR_GCCE_BATH_H_
#define __CCEX_SIMULATOR_GCCE_BATH_H_
// DEER bath pulses (pulse.bath) on the gCCE CPU paths. A file of its own: simulator_gcce.cpp only calls
// applyBathPulsesGcce, so its own code compiles as before (simulator.h, the note above firstFreeGcce).
#include "calculation/simulator.h"

// U_b for merged row `row`: exp(-i theta sigma_axis / 2) on every bath spin whose name is bpulse_defect and whose
// detuning is within kResonanceWindowRadkHz of a pump tone, identity on the others (kron over the cluster's spins).
// ba == NULL or no spin: 1x1 identity.
constexpr double kResonanceWindowRadkHz = 10.0;   // adding_MECCE src/simulator_gcce.cpp:546-553 (arbitrary, kept)
MatrixXcd bathPulseOperator(BathArray* ba, const Pulse* pls, int row);
// After prepareGccePulses: rows with a bath pulse get UpulsesExpanded[row] = kron(Upulses[row], U_b); if any row
// changed, pulseBdim = 0 so applyPulseGcce takes the dense product (its block path assumes kron(U, 1)).
void applyBathPulsesGcce(GcceWork* w, MatrixXcd* Upulses, BathArray* ba, const Pulse* pls);
// Once the bath and the Defect detunings are set (main.cpp), before any cluster: refused (rank 0 says why, ccexRefuse: every rank ends with exit 1, A55)
// when no bath spin (bath file or Defect sub-spin) has the name pulse.bath.defect, when such a spin's S is not
// pulse.bath.spin, or when pulse.bath.alphams/betams are not two different levels of that S. When none of those spins
// is within the window of a tone, rank 0 warns once and the run goes on (an off-resonant pump is a control).
void validateBathPulseTargets(BathArray* ba, DefectArray* dfa, const Pulse* pls);
#endif
