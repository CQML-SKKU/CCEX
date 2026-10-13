#include "calculation/superoperator.h"
#include <cstdio>
// Defined in main.cpp for the real binary; the test provides its own (as test_merge_unit.cpp).
bool verbosity = false;
int  rank = 0;
int  nprocess = 1;
static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { std::printf("  FAIL  %s\n", msg); ++fails; } else std::printf("  ok    %s\n", msg); } while (0)
int main() {
    srand(7);
    const int N = 3;
    MatrixXcd A = MatrixXcd::Random(N, N), B = MatrixXcd::Random(N, N), R = MatrixXcd::Random(N, N);
    CHECK((vec_to_mat(mat_to_vec(R)) - R).norm() < 1e-14, "vec/mat round trip");
    CHECK((op_to_supop(A, B) * mat_to_vec(R) - mat_to_vec(A * R * B)).norm() < 1e-12, "vec(A R B) = (A x B^T) vec(R)");
    MatrixXcd H = MatrixXcd::Random(N, N); H = (H + H.adjoint()).eval();
    MatrixXcd comm = doublec(0, -1) * (H * R - R * H);
    CHECK((coherent_superoperator(H) * mat_to_vec(R) - mat_to_vec(comm)).norm() < 1e-12, "coherent = -i[H, .]");
    MatrixXcd C = get_jump_operator_matrix("-", 3, 2.0);
    MatrixXcd D = C * R * C.adjoint() - 0.5 * (C.adjoint() * C * R + R * C.adjoint() * C);
    CHECK((collapse_superoperator_single(C) * mat_to_vec(R) - mat_to_vec(D)).norm() < 1e-12, "dissipator = C R C^dag - 1/2{C^dag C, R}");
    Eigen::RowVectorXcd tr = mat_to_vec(MatrixXcd::Identity(N, N)).transpose();   // Tr(X) = vec(I)^T vec(X), row-major
    CHECK((tr * collapse_superoperator_single(C)).norm() < 1e-12, "dissipator preserves the trace");
    // S = 1/2: sqrt(g) S_z = sqrt(g) diag(1/2, -1/2); S = 1: S_- has sqrt(2) below the diagonal
    MatrixXcd z = get_jump_operator_matrix("z", 2, 4.0);
    CHECK(std::abs(z(0, 0) - doublec(1.0, 0)) < 1e-14 && std::abs(z(1, 1) - doublec(-1.0, 0)) < 1e-14, "z, S=1/2, rate 4");
    MatrixXcd m = get_jump_operator_matrix("-", 3, 1.0);
    CHECK(std::abs(m(1, 0) - doublec(std::sqrt(2.0), 0)) < 1e-14 && std::abs(m(2, 1) - doublec(std::sqrt(2.0), 0)) < 1e-14, "S_-, S=1");
    CHECK((get_jump_operator_matrix("+", 3, 1.0) - m.adjoint()).norm() < 1e-14, "S_+ = S_-^dag");
    CHECK(isJumpOperatorType("Z") && isJumpOperatorType("-") && !isJumpOperatorType("+-") && !isJumpOperatorType("-+"), "operator names");
    std::printf(" superoperator unit : FAIL=%d\n", fails);
    return fails ? 1 : 0;
}
