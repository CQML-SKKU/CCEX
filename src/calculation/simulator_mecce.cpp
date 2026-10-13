#include "../../include/calculation/hamiltonian.h"
#include "../../include/calculation/simulator.h"
#include "../../include/calculation/superoperator.h"
#include "../../include/calculation/small_hamiltonian.h"

// ME-CCE (PRL 132, 250401 Eq. 2): the bath density matrix projected on the qubit states, <alpha|rho|beta>, under
// -i H_alpha rho + i rho H_beta and Lindblad jumps on the bath spins, propagated in Liouville space (row-major
// vectorization, exp(L t) per step). The qubit is not in the cluster space and L_0 = 1, as in the cCCE; qubit jump
// operators are refused at input (main.cpp validateMasterEquation; megCCE takes them).
// Ported from adding_MECCE src/simulator_mecce.cpp as of 74a27780b (before its qubit jump operators).
MatrixXcd* calCoherenceMecce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls){

    JumpOperatorArray* joa = cnf->jump_operators;

    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);

    int nqubit = QubitArray_getNqubit(qa);
    if (nqubit > 1){
        fprintf(stderr,"Error : calCoherenceMecce : nqubit > 1, current nqubit = %d\n",nqubit);
        ccexFatal(1);
    }

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // H_alpha, H_beta of the cCCE (bath + ms_a,b * secular hyperfine, with the knight and hfmedi terms of
    // hamiltonian.cCCE when they are on): the function calCoherenceCceInto uses
    MatrixXcd Halpha, Hbeta;
    cceClusterHamiltonians(qa,ba,cnf,Halpha,Hbeta);

    // spin operators of the cluster, for the jump operators
    MatrixXcd** bsigmas = BathArray_PauliOperators(ba);

    ////////////////////////////////
    // Lindbladian superoperator
    ////////////////////////////////

    // coherent part, projected: -i (H_alpha x 1 - 1 x H_beta^T)
    MatrixXcd L_coh = projected_coherent_superoperator(Halpha, Hbeta);

    // incoherent part: the dissipators of the jump operators on the bath spins (by name; a qubit name matches none)
    MatrixXcd L_incoh = incoherent_superoperator(ba, bsigmas, joa, bdim);

    MatrixXcd L = L_coh + L_incoh;

    // after an odd number of pi pulses alpha and beta are exchanged
    MatrixXcd L_flipped = projected_coherent_superoperator(Hbeta, Halpha) + L_incoh;

    ////////////////////////////////
    // Initial state
    ////////////////////////////////
    const bool isEnsemble = Config_getNstate(cnf) == 0;

    MatrixXcd rho0;
    if (isEnsemble){
        rho0 = MatrixXcd::Identity(bdim, bdim) / doublec(bdim, 0.0);
    }else{
        MatrixXcd psi0 = BathArray_Psi0(ba);
        rho0 = psi0 * psi0.adjoint();
    }
    Eigen::VectorXcd rho0_vec = mat_to_vec(rho0);

    ////////////////////////////////
    // Propagation
    ////////////////////////////////

    int nstep = Config_getNstep(cnf);
    double deltat = Config_getDeltat(cnf);
    MatrixXcd* result = new MatrixXcd[nstep];

    int npulse = Pulse_getNpulse(pls);
    double tfree = 0.0;

    for (int i = 0; i < nstep; i++){

        MatrixXcd super_prop;

        if (npulse == 0){
            // Ramsey
            super_prop = simple_incoherent_propagator(tfree, L);
        }else{
            // Hahn echo / CPMG: tau before the first pulse and after the last, 2 tau between
            double tau = tfree / (2.0 * npulse);

            MatrixXcd U_before = simple_incoherent_propagator(tau, L);
            MatrixXcd U_after  = simple_incoherent_propagator(tau, L_flipped);

            super_prop = propagate_superpropagators(U_before, U_after, npulse);
        }

        Eigen::VectorXcd rho_t_vec = super_prop * rho0_vec;
        MatrixXcd rho_t = vec_to_mat(rho_t_vec);

        // coherence: Tr(rho_t) (the ensemble's rho_0 = 1/bdim carries the normalization)
        result[i] = MatrixXcd::Zero(1, 1);
        result[i](0, 0) = rho_t.trace();

        tfree += deltat;
    }

    // free
    for (int i = 0; i < nspin; i++){
        delete[] bsigmas[i];
    }
    delete[] bsigmas;

    return result;
}
