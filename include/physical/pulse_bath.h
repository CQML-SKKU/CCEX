#ifndef __CCEX_PULSE_BATH_H_
#define __CCEX_PULSE_BATH_H_
#include <vector>
#include "../../include/physical/pulse.h"
struct SequenceRow { double frac; char axis; double angle; };      // interval fraction, then the pulse at its end
struct MergedRow   { double frac; char qAxis; double qAngle; char bAxis; double bAngle; };
// Rows of each input are intervals whose pulse sits at the interval's end; the last row's pulse is never applied.
// Pulses with axis 'I' or angle 0 are dropped. Edges closer than kMergeTol are one edge (both pulses there).
// Every input pulse lands on exactly one output edge. Output: the merged intervals, last row 'I'/0 for both.
// Two pulses of the same input within kMergeTol share an edge and the later one wins (so do two bath pulses each
// within kMergeTol of one qubit edge): the caller refuses that input.
constexpr double kMergeTol = 1e-9;
std::vector<MergedRow> mergePulseSequences(const std::vector<SequenceRow>& qubit, const std::vector<SequenceRow>& bath);

// bnpulse > 0: replaces npulse, sequence, pulse_axes, pulse_angles, sequence_indices with the merged sequence,
// fills tot_bpulse_axes/angles, and sets pulsename to "Manual". bnpulse == 0: does nothing.
// Refused (rank 0 says why, ccexRefuse: every rank ends with exit 1, A55): a last bath row that is not identity, and two pulses of one channel that the
// merge would put on one edge (one of them would be lost). Call once, after every option source is applied.
void Pulse_resolveBathSequence(Pulse* pls);
bool Pulse_hasBathPulses(const Pulse* pls);   // bnpulse > 0
#endif
