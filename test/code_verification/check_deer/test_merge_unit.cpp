#include "physical/pulse_bath.h"
#include <cmath>
#include <cstdio>
// Defined in main.cpp for the real binary; the test provides its own (as test_rotation_unit.cpp).
bool verbosity = false;
int  rank = 0;
int  nprocess = 1;
static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { std::printf("  FAIL  %s\n", msg); ++fails; } else std::printf("  ok    %s\n", msg); } while (0)
static bool near(double a, double b) { return std::fabs(a - b) < 1e-12; }
static double sum(const std::vector<MergedRow>& m) { double s = 0; for (auto& r : m) s += r.frac; return s; }

int main() {
    const std::vector<SequenceRow> hahn = {{0.5, 'X', 180.0}, {0.5, 'I', 0.0}};
    {   // bath pi exactly on the qubit pi: one edge, both pulses
        auto m = mergePulseSequences(hahn, {{0.5, 'X', 180.0}, {0.5, 'I', 0.0}});
        CHECK(m.size() == 2, "coincident: 2 rows");
        CHECK(near(m[0].frac, 0.5) && m[0].qAxis == 'X' && m[0].bAxis == 'X', "coincident: X on both at 0.5");
        CHECK(m[1].qAxis == 'I' && m[1].bAxis == 'I', "coincident: last row identity");
    }
    {   // bath pi at 0.505: two edges, each pulse once
        auto m = mergePulseSequences(hahn, {{0.505, 'X', 180.0}, {0.495, 'I', 0.0}});
        CHECK(m.size() == 3, "near: 3 rows");
        CHECK(near(m[0].frac, 0.5) && m[0].qAxis == 'X' && m[0].bAxis == 'I', "near: qubit only at 0.5");
        CHECK(near(m[1].frac, 0.005) && m[1].qAxis == 'I' && m[1].bAxis == 'X', "near: bath only at 0.505");
        CHECK(near(sum(m), 1.0), "near: fractions sum to 1");
    }
    {   // within kMergeTol: one edge (the branch would apply the bath pulse on two edges here)
        auto m = mergePulseSequences(hahn, {{0.5 + 5e-10, 'X', 180.0}, {0.5 - 5e-10, 'I', 0.0}});
        int nb = 0; for (auto& r : m) nb += (r.bAxis != 'I');
        CHECK(m.size() == 2 && nb == 1, "5e-10 apart: one edge, bath pulse applied once");
    }
    {   // 0.009 apart must stay two edges (da3fbe017: the old 0.01 window merged these and applied pi twice)
        auto m = mergePulseSequences(hahn, {{0.509, 'X', 180.0}, {0.491, 'I', 0.0}});
        int nb = 0; for (auto& r : m) nb += (r.bAxis != 'I');
        CHECK(m.size() == 3 && nb == 1, "0.009 apart: two edges, bath pulse once");
    }
    {   // Ramsey qubit + bath pi at 0.5
        auto m = mergePulseSequences({{1.0, 'I', 0.0}}, {{0.5, 'X', 180.0}, {0.5, 'I', 0.0}});
        CHECK(m.size() == 2 && m[0].qAxis == 'I' && m[0].bAxis == 'X', "ramsey: bath-only edge");
    }
    {   // CPMG-2 qubit + no bath pulse: unchanged
        auto m = mergePulseSequences({{0.25, 'X', 180.0}, {0.5, 'X', 180.0}, {0.25, 'I', 0.0}}, {{1.0, 'I', 0.0}});
        CHECK(m.size() == 3 && near(m[1].frac, 0.5) && m[2].qAxis == 'I', "cpmg2 alone unchanged");
    }
    {   // R6: two pulses of one channel within kMergeTol share an edge and one is lost, so Pulse_resolveBathSequence
        // refuses such an input (Task A2); the merge itself keeps its contract for every other input
        auto m = mergePulseSequences({{0.5, 'X', 180.0}, {5e-10, 'Y', 90.0}, {0.5 - 5e-10, 'I', 0.0}}, {{1.0, 'I', 0.0}});
        int nq = 0; for (auto& r : m) nq += (r.qAxis != 'I');
        CHECK(m.size() == 2 && nq == 1, "same-channel pair within kMergeTol: one edge, one pulse survives (caller refuses)");
    }
    {   // pulses of different channels 2e-9 apart (beyond kMergeTol) stay on two edges, each once
        auto m = mergePulseSequences(hahn, {{0.5 + 2e-9, 'X', 180.0}, {0.5 - 2e-9, 'I', 0.0}});
        int nq = 0, nb = 0; for (auto& r : m) { nq += (r.qAxis != 'I'); nb += (r.bAxis != 'I'); }
        CHECK(m.size() == 3 && nq == 1 && nb == 1, "2e-9 apart: two edges, each pulse once");
    }
    {   // a bath pulse at 1 + 5e-10 (its fractions sum to 1 within the reader's 1e-9): the last interval is 0, not
        // negative (a negative fraction would run that interval backwards in time)
        auto m = mergePulseSequences(hahn, {{1.0 + 5e-10, 'X', 180.0}, {0.0, 'I', 0.0}});
        CHECK(m.size() == 3 && m[1].bAxis == 'X' && m[2].frac >= 0.0, "edge past 1 within 1e-9: last fraction clamped to 0");
    }
    std::printf(" merge unit : FAIL=%d\n", fails);
    return fails ? 1 : 0;
}
