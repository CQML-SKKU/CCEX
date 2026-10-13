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
    // complex jump operators ("y" is imaginary, "x" real-symmetric): a C vs C* or C^T vs C^dag mix-up shows here
    for (const char* t : {"x", "y", "+", "z"}) {
        for (int dim : {2, 3}) {
            MatrixXcd Cc = get_jump_operator_matrix(t, dim, 1.7);
            MatrixXcd Rc = MatrixXcd::Random(dim, dim);
            MatrixXcd Dc = Cc * Rc * Cc.adjoint() - 0.5 * (Cc.adjoint() * Cc * Rc + Rc * Cc.adjoint() * Cc);
            char msg[96]; std::snprintf(msg, sizeof msg, "dissipator, complex R, op %s, dim %d", t, dim);
            CHECK((collapse_superoperator_single(Cc) * mat_to_vec(Rc) - mat_to_vec(Dc)).norm() < 1e-12, msg);
            Eigen::RowVectorXcd trc = mat_to_vec(MatrixXcd::Identity(dim, dim)).transpose();
            std::snprintf(msg, sizeof msg, "trace preserved, op %s, dim %d", t, dim);
            CHECK((trc * collapse_superoperator_single(Cc)).norm() < 1e-12, msg);
        }
    }
    CHECK(get_jump_operator_matrix("y", 2, 1.0).imag().norm() > 0.1, "y is genuinely complex");
    // closed-form pure dephasing: S=1/2, C = sqrt(g) S_z, rho0 = |+x><+x|: rho_01(t) = 1/2 exp(-g t / 2), Tr = 1
    {
        const double g = 2.0;
        MatrixXcd Lz = collapse_superoperator_single(get_jump_operator_matrix("z", 2, g));
        MatrixXcd rho0(2, 2); rho0 << 0.5, 0.5, 0.5, 0.5;
        bool ok = true;
        for (double t : {0.1, 0.7, 1.9, 4.0}) {
            MatrixXcd rt = vec_to_mat(simple_incoherent_propagator(t, Lz) * mat_to_vec(rho0));
            ok = ok && std::abs(rt(0, 1) - doublec(0.5 * std::exp(-g * t / 2), 0)) < 1e-12
                    && std::abs(rt(1, 0) - doublec(0.5 * std::exp(-g * t / 2), 0)) < 1e-12
                    && std::abs(rt.trace() - doublec(1, 0)) < 1e-12;
        }
        CHECK(ok, "pure dephasing: rho_01 = 1/2 exp(-g t/2), trace 1");
    }
    // S = 3/2 ladder: sqrt(3), 2, sqrt(3)
    {
        MatrixXcd l = get_jump_operator_matrix("-", 4, 1.0);
        CHECK(std::abs(l(1, 0) - doublec(std::sqrt(3.0), 0)) < 1e-14 && std::abs(l(2, 1) - doublec(2.0, 0)) < 1e-14
              && std::abs(l(3, 2) - doublec(std::sqrt(3.0), 0)) < 1e-14 && std::abs(l.sum() - doublec(2 * std::sqrt(3.0) + 2, 0)) < 1e-13, "S_-, S=3/2");
    }
    // expand_operator and incoherent_superoperator on a 2-spin toy bath: "e" (S=1/2) and "n" (S=1)
    {
        MatrixXcd sig0[1] = {MatrixXcd::Identity(2, 2)}, sig1[1] = {MatrixXcd::Identity(3, 3)};
        MatrixXcd* bs[2] = {sig0, sig1};
        MatrixXcd I2 = MatrixXcd::Identity(2, 2), I3 = MatrixXcd::Identity(3, 3);
        MatrixXcd o2 = get_jump_operator_matrix("y", 2, 1.0), o3 = get_jump_operator_matrix("x", 3, 1.0);
        CHECK((expand_operator(o2, bs, 2, 0) - kron(o2, I3)).norm() < 1e-14, "expand: operator on spin 0 = kron(op, I)");
        CHECK((expand_operator(o3, bs, 2, 1) - kron(I2, o3)).norm() < 1e-14, "expand: operator on spin 1 = kron(I, op)");
        BathSpin* sp[2] = {new BathSpin(), new BathSpin()};
        std::strcpy(sp[0]->name, "e"); sp[0]->spin = 0.5f;
        std::strcpy(sp[1]->name, "n"); sp[1]->spin = 1.0f;
        BathArray ba; ba.nspin = 2; ba.bath = sp; ba.prop_nspecies = 0; ba.prop_names = NULL; ba.prop_gyros = NULL; ba.prop_spins = NULL;
        JumpOperatorArray* joa = JumpOperatorArray_init();
        JumpOperatorArray_alloc(joa, 2);
        JumpOperatorArray_setOp(joa, 0, "E", "y", 3.0);        // case-insensitive name match, rate in rad*kHz as given
        JumpOperatorArray_setOp(joa, 1, "qubit", "z", 5.0);    // matches no bath spin
        MatrixXcd L = incoherent_superoperator(&ba, bs, joa, 6);
        MatrixXcd Le = collapse_superoperator_single(kron(get_jump_operator_matrix("y", 2, 3.0), I3));
        CHECK(L.rows() == 36 && (L - Le).norm() < 1e-12, "incoherent: jump on 'e' only (name match, qubit adds nothing)");
        JumpOperatorArray_setOp(joa, 1, "n", "x", 0.5);
        MatrixXcd L2 = incoherent_superoperator(&ba, bs, joa, 6);
        MatrixXcd Ln = collapse_superoperator_single(kron(I2, get_jump_operator_matrix("x", 3, 0.5)));
        CHECK((L2 - Le - Ln).norm() < 1e-12, "incoherent: 'e' and 'n' jumps add");
        CHECK(incoherent_superoperator(&ba, bs, NULL, 6).norm() == 0.0, "incoherent: no jumps -> 0");
        JumpOperatorArray_freeAll(joa);
        delete sp[0]; delete sp[1];
    }
    std::printf(" superoperator unit : FAIL=%d\n", fails);
    return fails ? 1 : 0;
}
