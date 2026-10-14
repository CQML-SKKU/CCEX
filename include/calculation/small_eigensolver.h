#ifndef CCEX_SMALL_EIGENSOLVER_H
#define CCEX_SMALL_EIGENSOLVER_H
// A42: the eigensystems of the small cCCE cluster Hamiltonians (bath dimension 2 and 4) by a cyclic complex Jacobi with a
// FIXED number of sweeps, several matrices at once (one per lane). Every lane does the same operations - no branch on the
// data - and the file is compiled without contraction (-ffp-contract=off in both builds and a pragma), so a matrix gets the
// same bits alone (count 1) as in any lane of a group: the batched setup of the A30 window and the scalar path of
// calCoherenceCceInto agree bit for bit. Not Eigen's arithmetic: the eigensystems differ from Eigen's in the last digits
// (P1 coherences by <= 7.1e-10, (internal notes)). A lane that has not converged after the
// sweeps is reported (ok false) and its caller solves that matrix with Eigen; CCEX_SMALL_EIG=0 turns the whole thing off
// (Eigen for every matrix, the code before A42). gCCE and the GPU paths do not use it.
#include <complex>

namespace smalleig {
constexpr int maxLanes = 8;
constexpr int sweeps = 4;                                         // for D 4; D 2 needs and does 2
// false when CCEX_SMALL_EIG=0 (read once)
bool on();
// Solves count (1..maxLanes) Hermitian D x D matrices, D 2 or 4: H[l] column major (both triangles read); the eigenvalues
// ascending into w[l][0..D), the eigenvectors, column major in the same order, into V[l][0..D*D). ok[l] is false when lane
// l's off-diagonal part after the sweeps is not below 1e-14 of the matrix (Frobenius norms), when H[l] has a part that is
// not finite or above 1e150 or its largest part is below 1e-130 (and not 0), or when w or V is not finite; its w and V
// are then not to be used.
void solve(int D, int count, const std::complex<double>* const* H, double* const* w, std::complex<double>* const* V, bool* ok);
// [profile] counts of this rank (solved, fell back to Eigen); counted by the callers
void count(long long solved, long long fallbacks);
void takeCounts(long long& solved, long long& fallbacks);
}
#endif
