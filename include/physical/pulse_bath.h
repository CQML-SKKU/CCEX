#ifndef __CCEX_PULSE_BATH_H_
#define __CCEX_PULSE_BATH_H_
#include <vector>
struct SequenceRow { double frac; char axis; double angle; };      // interval fraction, then the pulse at its end
struct MergedRow   { double frac; char qAxis; double qAngle; char bAxis; double bAngle; };
// Rows of each input are intervals whose pulse sits at the interval's end; the last row's pulse must be identity.
// Pulses with axis 'I' or angle 0 are dropped. Edges closer than kMergeTol are one edge (both pulses there).
// Every input pulse lands on exactly one output edge. Output: the merged intervals, last row 'I'/0 for both.
// Two pulses of the same input within kMergeTol share an edge and the later one wins: the caller refuses that input.
constexpr double kMergeTol = 1e-9;
std::vector<MergedRow> mergePulseSequences(const std::vector<SequenceRow>& qubit, const std::vector<SequenceRow>& bath);
#endif
