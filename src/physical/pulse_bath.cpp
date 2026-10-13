#include "../../include/physical/pulse_bath.h"
#include <algorithm>
#include <cmath>

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
    out.push_back(MergedRow{1.0 - prev, 'I', 0.0, 'I', 0.0});
    return out;
}
