#ifndef __CCEX_SUPEROPERATOR_H_
#define __CCEX_SUPEROPERATOR_H_

#include "common/utilities.h"
#include "physical/bath.h"

/**
 * @brief Jump operator definition for ME-CCE / ME-gCCE
 * @details Defines a collapse operator C = sqrt(rate) S_op acting on the bath spins named bath_name
 *          (or on the qubit, named "qubit" or by its name, in ME-gCCE).
 */
typedef struct JumpOperator {
    char bath_name[MAX_CHARARRAY_LENGTH]; /**< Bath spin name (or "qubit") to apply this jump operator */
    char op_type[4];                      /**< Operator type: "-" (lowering), "+" (raising), "z" (dephasing), "x", "y" */
    double rate;                          /**< Rate of the jump process (Unit: rad*kHz) */
} JumpOperator;

/**
 * @brief Collection of jump operators for ME-CCE
 */
typedef struct JumpOperatorArray {
    int njump;           /**< Number of jump operators */
    JumpOperator** ops;  /**< Array of jump operator pointers */
} JumpOperatorArray;

// JumpOperatorArray management
JumpOperatorArray* JumpOperatorArray_init();
void JumpOperatorArray_alloc(JumpOperatorArray* joa, int njump);
/** @brief Stores rate_radkhz as given; the MHz -> rad*kHz conversion is done once, by the input parser. */
void JumpOperatorArray_setOp(JumpOperatorArray* joa, int i, const char* bath_name, const char* op_type, double rate_radkhz);
void JumpOperatorArray_freeAll(JumpOperatorArray* joa);
void JumpOperatorArray_report(const JumpOperatorArray* joa);

/** @brief True for "-", "+", "z", "x", "y" (case-insensitive). */
bool isJumpOperatorType(const char* op_type);

// Superoperator utility functions

/**
 * @brief Convert operator pair to superoperator: S[rho] = A * rho * B  -->  supop = kron(A, B^T)
 * @note Uses row-major vectorization convention (matching PyCCE)
 */
MatrixXcd op_to_supop(const MatrixXcd& A, const MatrixXcd& B);

/**
 * @brief Vectorize density matrix: N x N -> N^2 x 1 (row-major order)
 */
Eigen::VectorXcd mat_to_vec(const MatrixXcd& M);

/**
 * @brief Un-vectorize: N^2 x 1 -> N x N (row-major order)
 */
MatrixXcd vec_to_mat(const Eigen::VectorXcd& v);

/**
 * @brief Coherent superoperator: L_coh = -i(H x I - I x H^T)
 */
MatrixXcd coherent_superoperator(const MatrixXcd& H);

/**
 * @brief Projected coherent superoperator: L_proj = -i(H_a x I - I x H_b^T)
 * @details Used in conventional CCE where Hamiltonians are projected on qubit alpha/beta states
 */
MatrixXcd projected_coherent_superoperator(const MatrixXcd& Ha, const MatrixXcd& Hb);

/**
 * @brief Build collapse superoperator for a single jump operator C
 * @details D = C x C* - 1/2 (C^dag C) x I - 1/2 I x (C^dag C)^T
 */
MatrixXcd collapse_superoperator_single(const MatrixXcd& C);

/**
 * @brief Get sqrt(rate) S_op for a spin of Hilbert space dimension dim = 2S+1
 * @param op_type Operator type string: "-", "+", "z", "x", "y"
 * @param dim Hilbert space dimension of the spin (2S+1)
 * @param rate Square root of the rate is multiplied to the operator
 */
MatrixXcd get_jump_operator_matrix(const char* op_type, int dim, double rate);

/**
 * @brief Expand a single-spin operator to the full Hilbert space
 * @param single_op Operator acting on spin at index ispin
 * @param bsigmas Pauli operators for all spins (used for dimensions)
 * @param nspin Total number of spins
 * @param ispin Index of the spin this operator acts on
 */
MatrixXcd expand_operator(const MatrixXcd& single_op, MatrixXcd** bsigmas, int nspin, int ispin);

/**
 * @brief Build incoherent superoperator for all applicable bath spins in the cluster
 * @param ba BathArray of the cluster
 * @param bsigmas Pauli operators for each bath spin
 * @param joa JumpOperatorArray containing all jump operator definitions
 * @param bdim Total bath Hilbert space dimension
 */
MatrixXcd incoherent_superoperator(BathArray* ba, MatrixXcd** bsigmas, const JumpOperatorArray* joa, int bdim);

/**
 * @brief Compute super-propagator: U(t) = exp(L * t)
 */
MatrixXcd simple_incoherent_propagator(double t, const MatrixXcd& L);

/**
 * @brief Propagate CPMG super-propagators
 * @param Ub Superoperator before pi-pulse
 * @param Ua Superoperator after pi-pulse (swapped alpha/beta)
 * @param n Number of pi-pulses
 */
MatrixXcd propagate_superpropagators(const MatrixXcd& Ub, const MatrixXcd& Ua, int n);

#endif // __CCEX_SUPEROPERATOR_H_
