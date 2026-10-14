#include "../../include/physical/pulse_bath.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace {
struct Edge { double at; char qAxis = 'I'; double qAngle = 0.0; char bAxis = 'I'; double bAngle = 0.0; };
bool isIdentity(const SequenceRow& r) { return r.axis == 'I' || std::fabs(r.angle) < 1e-12; }
// the edge within kMergeTol of `at`, or a new one
Edge& edgeAt(std::vector<Edge>& edges, double at) {
    for (auto& e : edges) if (std::fabs(e.at - at) < kMergeTol) return e;
    edges.push_back(Edge{at});
    return edges.back();
}
}

std::vector<MergedRow> mergePulseSequences(const std::vector<SequenceRow>& qubit, const std::vector<SequenceRow>& bath) {
    std::vector<Edge> edges;
    double at = 0.0;   // qubit edges first: a bath edge within kMergeTol takes the qubit edge's position
    for (size_t i = 0; i + 1 < qubit.size(); ++i) {
        at += qubit[i].frac;
        if (isIdentity(qubit[i])) continue;
        Edge& e = edgeAt(edges, at); e.qAxis = qubit[i].axis; e.qAngle = qubit[i].angle;
    }
    at = 0.0;
    for (size_t i = 0; i + 1 < bath.size(); ++i) {
        at += bath[i].frac;
        if (isIdentity(bath[i])) continue;
        Edge& e = edgeAt(edges, at); e.bAxis = bath[i].axis; e.bAngle = bath[i].angle;
    }
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) { return a.at < b.at; });
    std::vector<MergedRow> out;
    double prev = 0.0;
    for (const auto& e : edges) {
        out.push_back(MergedRow{e.at - prev, e.qAxis, e.qAngle, e.bAxis, e.bAngle});
        prev = e.at;
    }
    // an edge may lie up to 1e-9 past 1 (the readers' tolerance on the sum): the last interval is then 0, not negative
    out.push_back(MergedRow{std::max(0.0, 1.0 - prev), 'I', 0.0, 'I', 0.0});
    return out;
}

// ---- Part 2: the merged sequence installed into the Pulse (Task A2) ----

bool Pulse_hasBathPulses(const Pulse* pls) { return pls->bnpulse > 0; }

namespace {
void refuseSequence(const char* msg) {
    if (rank == 0) fprintf(stderr, "Error: %s\n", msg);
    ccexRefuse();
}
// the pulses the merge places: every row's but the last's, identity ones dropped
int placedPulses(const std::vector<SequenceRow>& rows) {
    int n = 0;
    for (size_t i = 0; i + 1 < rows.size(); ++i) n += !isIdentity(rows[i]);
    return n;
}
}

void Pulse_resolveBathSequence(Pulse* pls) {
    if (!Pulse_hasBathPulses(pls)) return;
    std::vector<SequenceRow> q, b;
    for (int i = 0; i <= pls->npulse; ++i)  q.push_back({pls->sequence[i][2], pls->pulse_axes[i], pls->pulse_angles[i]});
    for (int i = 0; i <= pls->bnpulse; ++i) b.push_back({pls->bsequence_fracs[i], pls->bpulse_axes[i], pls->bpulse_angles[i]});
    // a row's pulse sits at the end of its interval: the last row's would come after the run (the qubit's last row is
    // ignored as everywhere else in CCEX - the generated CPMG writes X 180 there)
    if (!isIdentity(b.back()))
        refuseSequence("the last row of pulse.bath.sequence must be identity (axis I or angle 0): its pulse would come "
                       "after the last interval and never be applied");
    const std::vector<MergedRow> m = mergePulseSequences(q, b);
    // each pulse exactly once (R6): two of one channel on one edge keep the later one only
    int mq = 0, mb = 0;
    for (const auto& r : m) { mq += !isIdentity({r.frac, r.qAxis, r.qAngle}); mb += !isIdentity({r.frac, r.bAxis, r.bAngle}); }
    if (mq < placedPulses(q))
        refuseSequence("pulse.sequence has two qubit pulses at the same time; bath pulses need them apart "
                       "(more than 1e-9 of the sequence)");
    if (mb < placedPulses(b))
        refuseSequence("pulse.bath.sequence has two bath pulses at the same time; give each its own time (more than 1e-9 "
                       "of the sequence apart, and not both within 1e-9 of one qubit pulse)");
    // the qubit arrays are re-made at the merged length, with the same allocators as Pulse_alloc* (pulse.cpp)
    Pulse_freeSequence(pls); Pulse_freeAxes(pls); Pulse_freeAngles(pls); Pulse_freeSequenceIndices(pls);
    pls->npulse = (int)m.size() - 1;
    Pulse_allocSequence(pls); Pulse_allocAxes(pls); Pulse_allocAngles(pls); Pulse_allocSequenceIndices(pls);
    pls->tot_bpulse_axes   = (char*)calloc(m.size() + 1, sizeof(char));
    pls->tot_bpulse_angles = (double*)calloc(m.size(), sizeof(double));
    double at = 0.0;
    for (size_t i = 0; i < m.size(); ++i) {
        pls->sequence[i][0] = at; at += m[i].frac; pls->sequence[i][1] = at; pls->sequence[i][2] = m[i].frac;
        pls->pulse_axes[i] = m[i].qAxis;      pls->pulse_angles[i] = m[i].qAngle;
        pls->tot_bpulse_axes[i] = m[i].bAxis; pls->tot_bpulse_angles[i] = m[i].bAngle;
    }
    assign_sequence_indices(pls);
    Pulse_setPulsename(pls, (char*)"Manual");
}
