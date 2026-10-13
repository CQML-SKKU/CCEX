#include "../../include/calculation/hamiltonian.h"
#include "../../include/calculation/simulator.h"
#include "../../include/calculation/superoperator.h"
#include "../../include/calculation/small_hamiltonian.h"

// ME-gCCE: the gCCE (qubit and bath spins in the cluster space, the full Hamiltonian) with Lindblad jump operators,
// propagated in Liouville space (row-major vectorization, exp(L tau) per interval). Bath jumps act as 1_q (x) C_B,
// qubit jumps (bath_name "qubit" or the qubit's name) as C_q (x) 1_B in every cluster, the 0th cluster included, which
// the coherence is divided by as in the gCCE (PyCCE 1.1.1 LindbladgCCE has the same structure).
// Ported from adding_MECCE src/simulator_megcce.cpp as of 6022d7998, without the bath-pulse arguments (the DEER bath
// pulses of a megCCE run are refused at input) and without its full density-matrix output (savemode
// cluster_density_matrix and quantity dm are refused at input, main.cpp validateMasterEquation).

// L_incoh of the jump operators in the qubit x bath space: C_q (x) 1_B for a qubit jump, 1_q (x) C_B for each bath spin of
// the cluster with the jump's name. ba == NULL (the 0th cluster): the qubit jumps alone.
static MatrixXcd incoherent_superoperator_megcce(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas,
                                                 const JumpOperatorArray* joa, int bdim, int qdim){
    int totdim = qdim * bdim;
    int N2 = totdim * totdim;
    MatrixXcd L_incoh = MatrixXcd::Zero(N2, N2);

    if (joa == NULL || joa->njump == 0){
        return L_incoh;
    }

    int nspin = BathArray_getNspin(ba);
    MatrixXcd Iq = MatrixXcd::Identity(qdim, qdim);
    MatrixXcd Ib = MatrixXcd::Identity(bdim, bdim);

    for (int ij = 0; ij < joa->njump; ij++){
        const char* target_name = joa->ops[ij]->bath_name;
        const char* op_type = joa->ops[ij]->op_type;
        double rate = joa->ops[ij]->rate;

        // qubit jump operator: bath_name "qubit" or the name of the qubit (the test of validateMasterEquation)
        bool isQubit = (strcasecmp(target_name, "qubit") == 0 ||
                        QubitArray_getQubitIdx_fromName(qa, target_name) >= 0);
        if (isQubit){
            MatrixXcd C_q    = get_jump_operator_matrix(op_type, qdim, rate);
            MatrixXcd C_full = kron(C_q, Ib);
            L_incoh += collapse_superoperator_single(C_full);
            continue;
        }

        // bath spin jump operators: C_full = 1_q (x) C_bath
        for (int ib = 0; ib < nspin; ib++){
            char* spin_name = BathArray_getBath_i_name(ba, ib);

            if (strcasecmp(spin_name, target_name) == 0){
                int dim_i = BathArray_dimBath_i(ba, ib);

                // the single-spin jump operator, expanded to the bath space, then to qubit x bath
                MatrixXcd C_single = get_jump_operator_matrix(op_type, dim_i, rate);
                MatrixXcd C_bath = expand_operator(C_single, bsigmas, nspin, ib);
                MatrixXcd C_full = kron(Iq, C_bath);

                L_incoh += collapse_superoperator_single(C_full);
            }
        }
    }

    return L_incoh;
}

// The super-propagator of one time point: S_total = S_free(tau_{N+1}) S_pulse_N ... S_pulse_1 S_free(tau_1), with
// S_free = exp(L tau) and S_pulse = U (x) U^* for U = U_q (x) 1_B; the intervals follow the pulse sequence (its fractions,
// a Manual sequence included), as calPropagatorGcce.
static MatrixXcd calSuperPropagatorMegcce(BathArray* ba, const MatrixXcd& L, Pulse* pulse, double tfree,
                                          MatrixXcd* Upulses){

    int bdim = BathArray_dim(ba);
    int npulse = Pulse_getNpulse(pulse);
    double** sequence = Pulse_getSequence(pulse);
    int* sequence_indices = Pulse_getSequenceIndices(pulse);

    int N2 = L.rows(); // superoperator dimension (totdim^2)

    MatrixXcd Stotal;
    MatrixXcd* Sfrees = new MatrixXcd[npulse + 1];

    for (int ipulse = 0; ipulse < npulse + 1; ipulse++){

        double tau = tfree * sequence[ipulse][2];
        int sameTauIndex = sequence_indices[ipulse];

        // pulse superoperator
        MatrixXcd Spulse;
        if (ipulse < npulse){
            // the qubit pulse over the bath (DEER bath pulses: Phase C)
            MatrixXcd Upulse_full = kron(Upulses[ipulse], MatrixXcd::Identity(bdim, bdim));

            // S[rho] = U rho U^dagger
            Spulse = op_to_supop(Upulse_full, Upulse_full.adjoint());
        } else {
            // identity superoperator for the last interval
            Spulse = MatrixXcd::Identity(N2, N2);
        }

        // free evolution superoperator, exp(L tau), once per distinct interval
        MatrixXcd Sfree;
        if (sameTauIndex == ipulse){
            Sfree = (L * tau).exp();
        } else {
            Sfree = Sfrees[sameTauIndex];
        }

        // total super-propagator
        if (ipulse == 0 && ipulse != npulse){
            Stotal = Spulse * Sfree;
        } else if (ipulse == 0 && ipulse == npulse){
            Stotal = Sfree;
        } else if (ipulse == npulse){
            Stotal = Sfree * Stotal;
        } else {
            Stotal = Spulse * Sfree * Stotal;
        }

        Sfrees[ipulse] = Sfree;
    }

    delete[] Sfrees;
    return Stotal;
}

// op: the run's Output, as calCoherenceGcce takes it; nothing is written to it here (the full density matrices of savemode
// cluster_density_matrix are not ported, and the coherence goes back to calculate)
MatrixXcd* calCoherenceMegcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, Output* op){

    JumpOperatorArray* joa = cnf->jump_operators;

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);

    int qdim = QubitArray_dim(qa);
    int bdim = BathArray_dim(ba);
    int totdim = qdim * bdim;

    ////////////////////////////////
    // Hamiltonian (the gCCE's)
    ////////////////////////////////
    MatrixXcd Htot = MatrixXcd::Zero(totdim, totdim);

    // Pauli matrices of the qubit
    MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

    // qubit Hamiltonian
    MatrixXcd Hq = HamilQubit(qa, ba, qsigmas, cnf);
    MatrixXcd Hq_expand = kron(Hq, MatrixXcd::Identity(bdim, bdim));

    // qubit pulses (pulse.per_qubit and pulse.detuning_factor != 1 are refused at input)
    MatrixXcd* Upulses = new MatrixXcd[Pulse_getNpulse(pulse)];
    cachedUpulses(Upulses, qa, pulse);

    MatrixXcd** bsigmas = NULL;

    if (nspin > 0){
        // bath Hamiltonian
        bsigmas = BathArray_PauliOperators(ba);
        MatrixXcd Hb = HamilBath(ba, bsigmas, cnf);
        MatrixXcd Hb_expand = kron(MatrixXcd::Identity(qdim, qdim), Hb);

        // qubit - bath interaction (full, not secular)
        MatrixXcd Hqb = HamilQubitBath(qa, ba, qsigmas, bsigmas, cnf);

        Htot = Hq_expand + Hb_expand + Hqb;
    } else {
        Htot = Hq;
    }

    ////////////////////////////////
    // Lindbladian superoperator
    ////////////////////////////////

    // coherent part: -i (H x 1 - 1 x H^T)
    MatrixXcd L_coh = coherent_superoperator(Htot);

    // incoherent part: the dissipators of the qubit and bath jump operators
    MatrixXcd L_incoh = incoherent_superoperator_megcce(qa, ba, bsigmas, joa, bdim, qdim);

    MatrixXcd L = L_coh + L_incoh;

    ////////////////////////////////
    // Density matrix
    ////////////////////////////////
    bool isEnsemble = (Config_getNstate(cnf) == 0);

    MatrixXcd qrho0 = QubitArray_Rho0(qa);
    MatrixXcd brho0 = BathArray_Rho0(ba, isEnsemble);
    MatrixXcd rho0 = kron(qrho0, brho0);
    Eigen::VectorXcd rho0_vec = mat_to_vec(rho0);

    ////////////////////////////////
    // Propagation
    ////////////////////////////////
    // the coherence only (quantity dm is refused at input)
    if (strcasecmp(Config_getQuantity(cnf), "coherence") != 0){
        fprintf(stderr, "Error : calCoherenceMegcce : quantity %s is not coherence\n", Config_getQuantity(cnf));
        ccexFatal(1);
    }
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    int nstep = Config_getNstep(cnf);
    double deltat = Config_getDeltat(cnf);
    double tfree = 0.0;
    MatrixXcd* result = new MatrixXcd[nstep];

    for (int i = 0; i < nstep; i++){

        // super-propagator of this time point
        MatrixXcd Sprop = calSuperPropagatorMegcce(ba, L, pulse, tfree, Upulses);

        // density matrix at t
        Eigen::VectorXcd rhot_vec = Sprop * rho0_vec;
        MatrixXcd rhot = vec_to_mat(rhot_vec);

        // trace over the bath
        MatrixXcd reducedRhot = rhot;
        for (int ib = nspin - 1; ib >= 0; ib--){
            int bdim_i = BathArray_dimBath_i(ba, ib);
            reducedRhot = partialtrace(reducedRhot, bdim_i, bdim_i);
        }

        // coherence: <alpha| rho_q |beta> (the gCCE's element)
        result[i] = psia.adjoint() * reducedRhot * psib;

        tfree += deltat;
    }

    // free
    if (bsigmas != NULL){
        for (int i = 0; i < nspin; i++){
            delete[] bsigmas[i];
        }
        delete[] bsigmas;
    }
    for (int i = 0; i < nqubit; i++){
        delete[] qsigmas[i];
    }
    delete[] qsigmas;
    delete[] Upulses;

    return result;
}
