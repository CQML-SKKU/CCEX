#include "../../include/calculation/superoperator.h"
#include "../../include/calculation/hamiltonian.h"
#include "../../include/common/memory.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include <cmath>

/* ================================================================
 * JumpOperatorArray management
 * ================================================================ */

JumpOperatorArray* JumpOperatorArray_init(){
    JumpOperatorArray* joa = (JumpOperatorArray*)allocArray1d(1, sizeof(JumpOperatorArray));
    joa->njump = 0;
    joa->ops = NULL;
    return joa;
}

void JumpOperatorArray_alloc(JumpOperatorArray* joa, int njump){
    joa->njump = njump;
    joa->ops = (JumpOperator**)allocArray2d(njump, 1, sizeof(JumpOperator));
    for (int i = 0; i < njump; i++){
        joa->ops[i]->bath_name[0] = '\0';
        joa->ops[i]->op_type[0] = '\0';
        joa->ops[i]->rate = 0.0;
    }
}

void JumpOperatorArray_freeAll(JumpOperatorArray* joa){
    if (joa->ops != NULL){
        freeArray2d((void***)&(joa->ops), joa->njump);
    }
    freeArray1d((void**)&joa);
}

void JumpOperatorArray_setOp(JumpOperatorArray* joa, int idx, const char* bath_name, const char* op_type, double rate_radkhz){
    strncpy(joa->ops[idx]->bath_name, bath_name, MAX_CHARARRAY_LENGTH - 1);
    joa->ops[idx]->bath_name[MAX_CHARARRAY_LENGTH - 1] = '\0';
    strncpy(joa->ops[idx]->op_type, op_type, sizeof(joa->ops[idx]->op_type) - 1);
    joa->ops[idx]->op_type[sizeof(joa->ops[idx]->op_type) - 1] = '\0';
    joa->ops[idx]->rate = rate_radkhz;
}

void JumpOperatorArray_report(const JumpOperatorArray* joa){
    if (rank == 0){
        printTitle("Jump Operators (ME-CCE)");
        printf("  Number of jump operators: %d\n", joa->njump);
        for (int i = 0; i < joa->njump; i++){
            printf("  [%d] bath_name: %s, operator: %s, rate: %.6e rad*kHz\n",
                   i, joa->ops[i]->bath_name, joa->ops[i]->op_type, joa->ops[i]->rate);
        }
        printLine();
    }
}

/* ================================================================
 * Superoperator utility functions
 * ================================================================ */

MatrixXcd op_to_supop(const MatrixXcd& left, const MatrixXcd& right){
    // S[rho] = left * rho * right
    // In vectorized form (row-major): S = kron(left, right^T)
    return kroneckerProduct(left, right.transpose()).eval();
}

VectorXcd mat_to_vec(const MatrixXcd& mat){
    // Row-major vectorization: stack rows
    int N = mat.rows();
    VectorXcd vec(N * N);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            vec(i * N + j) = mat(i, j);
        }
    }
    return vec;
}

MatrixXcd vec_to_mat(const VectorXcd& vec){
    int N2 = vec.size();
    int N = (int)std::round(std::sqrt((double)N2));
    MatrixXcd mat(N, N);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            mat(i, j) = vec(i * N + j);
        }
    }
    return mat;
}

MatrixXcd coherent_superoperator(const MatrixXcd& H){
    int N = H.rows();
    MatrixXcd eye = MatrixXcd::Identity(N, N);
    return doublec(0.0, -1.0) * (op_to_supop(H, eye) - op_to_supop(eye, H));
}

MatrixXcd projected_coherent_superoperator(const MatrixXcd& Ha, const MatrixXcd& Hb){
    int N = Ha.rows();
    MatrixXcd eye = MatrixXcd::Identity(N, N);
    return doublec(0.0, -1.0) * (op_to_supop(Ha, eye) - op_to_supop(eye, Hb));
}

MatrixXcd collapse_superoperator_single(const MatrixXcd& C){
    int N = C.rows();
    MatrixXcd eye = MatrixXcd::Identity(N, N);
    MatrixXcd CdagC = C.adjoint() * C;

    // D = C⊗C* - 1/2(C†C⊗I) - 1/2(I⊗(C†C)^T)
    // Note: C* = C.conjugate(), and in op_to_supop the right is transposed,
    //       so op_to_supop(C, C.adjoint()) = kron(C, (C†)^T) = kron(C, C*)
    MatrixXcd cn_rho_cndag = op_to_supop(C, C.adjoint());
    MatrixXcd cndag_cn_rho = op_to_supop(CdagC, eye);
    MatrixXcd rho_cndag_cn = op_to_supop(eye, CdagC);

    return cn_rho_cndag - 0.5 * cndag_cn_rho - 0.5 * rho_cndag_cn;
}

bool isJumpOperatorType(const char* t){
    return !strcasecmp(t, "-") || !strcasecmp(t, "+") || !strcasecmp(t, "z") || !strcasecmp(t, "x") || !strcasecmp(t, "y");
}

MatrixXcd get_jump_operator_matrix(const char* op_type, int dim, double rate){
    // Get the operator matrix for the given type
    // dim = 2S+1

    MatrixXcd op = MatrixXcd::Zero(dim, dim);
    double sqrt_rate = std::sqrt(rate);

    if (strcasecmp(op_type, "-") == 0){
        // Lowering operator: sigma_-
        // For spin-1/2: |0><1| (maps |up> to |down>)
        // General: S_- with matrix elements <m-1|S_-|m> = sqrt(S(S+1) - m(m-1))
        double S = (dim - 1.0) / 2.0;
        for (int i = 1; i < dim; i++){
            double m = S - i + 1; // m value of the bra state (row i-1 has m=S, row i has m=S-1, ...)
            // Actually for standard convention: row index r corresponds to m = S - r
            // S_- |m> = sqrt(S(S+1)-m(m-1)) |m-1>
            // <m-1|S_-|m> at matrix element (r+1, r) where r corresponds to m=S-r
            double m_val = S - (i - 1); // m value
            op(i, i - 1) = sqrt_rate * std::sqrt(S * (S + 1) - m_val * (m_val - 1));
        }
    }
    else if (strcasecmp(op_type, "+") == 0){
        // Raising operator: sigma_+ = (sigma_-)†
        double S = (dim - 1.0) / 2.0;
        MatrixXcd Sminus = MatrixXcd::Zero(dim, dim);
        for (int i = 1; i < dim; i++){
            double m_val = S - (i - 1);
            Sminus(i, i - 1) = std::sqrt(S * (S + 1) - m_val * (m_val - 1));
        }
        op = sqrt_rate * Sminus.adjoint();
    }
    else if (strcasecmp(op_type, "z") == 0){
        // sigma_z / S_z
        double S = (dim - 1.0) / 2.0;
        for (int i = 0; i < dim; i++){
            op(i, i) = sqrt_rate * (S - i);
        }
    }
    else if (strcasecmp(op_type, "x") == 0){
        // S_x = (S_+ + S_-) / 2
        double S = (dim - 1.0) / 2.0;
        MatrixXcd Sminus = MatrixXcd::Zero(dim, dim);
        for (int i = 1; i < dim; i++){
            double m_val = S - (i - 1);
            Sminus(i, i - 1) = std::sqrt(S * (S + 1) - m_val * (m_val - 1));
        }
        MatrixXcd Splus = Sminus.adjoint();
        op = sqrt_rate * (Splus + Sminus) / 2.0;
    }
    else if (strcasecmp(op_type, "y") == 0){
        // S_y = (S_+ - S_-) / (2i)
        double S = (dim - 1.0) / 2.0;
        MatrixXcd Sminus = MatrixXcd::Zero(dim, dim);
        for (int i = 1; i < dim; i++){
            double m_val = S - (i - 1);
            Sminus(i, i - 1) = std::sqrt(S * (S + 1) - m_val * (m_val - 1));
        }
        MatrixXcd Splus = Sminus.adjoint();
        op = sqrt_rate * (Splus - Sminus) / doublec(0.0, 2.0);
    }
    else{
        fprintf(stderr, "Error: Unknown jump operator type: %s\n", op_type);
        ccexFatal(1);
    }

    return op;
}

MatrixXcd expand_operator(const MatrixXcd& single_op, MatrixXcd** bsigmas, int nspin, int ispin){
    // Expand single-spin operator to full Hilbert space using tensor products
    // I ⊗ I ⊗ ... ⊗ single_op ⊗ ... ⊗ I

    if (nspin == 1){
        return single_op;
    }

    MatrixXcd result;
    bool initialized = false;

    for (int i = 0; i < nspin; i++){
        int dim_i = bsigmas[i][0].rows(); // Identity matrix dimension
        MatrixXcd op_i;
        if (i == ispin){
            op_i = single_op;
        } else {
            op_i = MatrixXcd::Identity(dim_i, dim_i);
        }

        if (!initialized){
            result = op_i;
            initialized = true;
        } else {
            result = kron(result, op_i);
        }
    }

    return result;
}

MatrixXcd incoherent_superoperator(BathArray* ba, MatrixXcd** bsigmas, const JumpOperatorArray* joa, int bdim){

    int N2 = bdim * bdim;
    MatrixXcd L_incoh = MatrixXcd::Zero(N2, N2);

    if (joa == NULL || joa->njump == 0){
        return L_incoh;
    }

    int nspin = BathArray_getNspin(ba);

    for (int ij = 0; ij < joa->njump; ij++){
        const char* target_name = joa->ops[ij]->bath_name;
        const char* op_type = joa->ops[ij]->op_type;
        double rate = joa->ops[ij]->rate;

        // Apply this jump operator to all bath spins with matching name
        for (int ib = 0; ib < nspin; ib++){
            char* spin_name = BathArray_getBath_i_name(ba, ib);

            if (strcasecmp(spin_name, target_name) == 0){
                // Get single-spin dimension
                int dim_i = BathArray_dimBath_i(ba, ib);

                // Build the single-spin jump operator matrix
                MatrixXcd C_single = get_jump_operator_matrix(op_type, dim_i, rate);

                // Expand to full Hilbert space
                MatrixXcd C_full = expand_operator(C_single, bsigmas, nspin, ib);

                // Add the collapse superoperator
                L_incoh += collapse_superoperator_single(C_full);
            }
        }
    }

    return L_incoh;
}

MatrixXcd simple_incoherent_propagator(double t, const MatrixXcd& L){
    // U(t) = exp(L * t)
    // Note: CCEX Hamiltonians are already in radkHz (angular frequency),
    //       so no 2*pi factor is needed (unlike PyCCE which uses kHz).
    MatrixXcd exponent = L * t;
    return exponent.exp();
}

MatrixXcd propagate_superpropagators(const MatrixXcd& u_before_pi, const MatrixXcd& u_after_pi, int number){
    // Combine propagators for CPMG sequence with `number` pi-pulses
    // v_he = u_after * u_before  (one Hahn echo unit)

    MatrixXcd v_he = u_after_pi * u_before_pi;

    if (number == 1){
        return v_he;
    }

    MatrixXcd v_he_reversed = u_before_pi * u_after_pi;
    MatrixXcd v_cp = v_he_reversed * v_he; // One CPMG unit (2 pi-pulses)

    if (number == 2){
        return v_cp;
    }

    // Matrix power for repeated CPMG units
    MatrixXcd nonunitary = MatrixXcd::Identity(v_cp.rows(), v_cp.cols());
    int full_cycles = number / 2;

    // Simple matrix power (could optimize with binary exponentiation for large N)
    for (int i = 0; i < full_cycles; i++){
        nonunitary = v_cp * nonunitary;
    }

    if (number % 2 == 1){
        nonunitary = v_he * nonunitary;
    }

    return nonunitary;
}
