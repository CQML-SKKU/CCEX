#include "../../include/calculation/simulator.h"
#include "../../include/calculation/hamiltonian.h"
#include "../../include/common/memory.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/gpu/eigen_refinement.h"
#include "../../include/calculation/small_hamiltonian.h"
#include "../../include/calculation/sparse_propagator.h"
#include "../../include/calculation/simulator_gcce_bath.h"
#ifdef CCEX_USE_CUDA
#include "../../include/calculation/gpu/gpu_backend.h"
#include <memory>

struct GcceBatch {
    std::vector<std::unique_ptr<GpuGcce,void(*)(GpuGcce*)>> owned;
    std::vector<GpuGcce*> works;
    std::vector<std::unique_ptr<MatrixXcd[]>> results;
    std::vector<MatrixXcd> alpha,beta;
    std::vector<double> build,setup;
    bool outOfMemory=false;          // a cluster's setup ran out of GPU memory (backend auto): the batch is given up
};
#endif

#define sPI 3.14159265358979323846

// exp(-i*Htot*tau) applied to w->psi, in place, without ever forming the propagator:
//     exp(-i*H*tau)|v> = M ( exp(-i*lambda*tau) .* (M^dag |v>) )
// Two matrix-vector products and one elementwise scaling, reading the decomposition
// prepareGcceWork already computed. `phase` is exp(-i*lambda*tau) for this tau.
// evolution=vector only.
static inline void applyUfreeGcce(GcceWork* w, const Eigen::VectorXcd& phase){
    w->psiscratch.noalias() = w->Madj * w->psi;
    w->psiscratch.array()  *= phase.array();
    w->psi.noalias()        = w->M * w->psiscratch;
}


static void validateGcceHamiltonian(const MatrixXcd& Htot){
    // The eigendecomposition below stands in for exp(-i*Htot*tau) only if Htot really
    // is Hermitian. It is a Hamiltonian, so it must be -- check rather than assume,
    // because SelfAdjointEigenSolver reads one triangle and would silently accept a
    // non-Hermitian matrix by ignoring the other half.
    if (!gcceHamiltonianFinite(Htot)){
        fprintf(stderr,"Error: gCCE Hamiltonian contains non-finite entries\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    // A25: the same test without the temporaries (small_hamiltonian.cpp)
    double asym = 0.0, scale = 0.0;
    if (!gcceHamiltonianHermitian(Htot, asym, scale)){
        fprintf(stderr,"Error : prepareGcceWork : Htot is not Hermitian "
                       "(max|H-H^dag| = %.3e, max|H| = %.3e)\n", asym, scale);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
}

void prepareGcceWork(GcceWork* w, const MatrixXcd& Htot, bool useEigen, bool refine, int maxSweeps){
    w->useEigen = useEigen;
    if (!useEigen) return;
    validateGcceHamiltonian(Htot);

    // A25: one solver kept from cluster to cluster (no allocation once the dimension has been seen; the constructor called
    // the same compute())
    static CcexEigensolver es;   // A46: Eigen's solver unless built with CCEX_BLAS_DIRECT
    es.compute(Htot);
    if (es.info() != Eigen::Success){
        fprintf(stderr,"Error : prepareGcceWork : eigendecomposition of Htot failed\n");
        ccexFatal(1);   // one rank's cluster: the run stops as a whole (D38)
    }
    w->evals = es.eigenvalues();
    w->M     = es.eigenvectors();
    w->evalCorrections=Eigen::VectorXd::Zero(w->evals.size());
    try {if(refine)refineEigensystem(Htot,w->M,w->evals,w->evalCorrections,maxSweeps);}
    catch(const std::exception& error){
        fprintf(stderr,"Error: gCCE eigensystem dim=%ld gcce_refinement.max_sweeps=%d: %s\n",
            static_cast<long>(Htot.rows()),maxSweeps,error.what());
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    w->Madj  = w->M.adjoint();
}

Matrix3cd Sx_spin1() {
    Matrix3cd Sx = Matrix3cd::Zero();
    double factor = 1.0 / sqrt(2.0);
    //double factor = 1.0;

    Sx(0,1) = factor;
    Sx(1,0) = factor;
    Sx(1,2) = factor;
    Sx(2,1) = factor;

    return Sx;
}

MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, Output* op, GcceBatch* batchCollector){
    return calCoherenceGcceInto(qa, ba, cnf, pulse, op, nullptr, nullptr, batchCollector);
}

// A14: with into, the results go there (1 x 1 or qdim x qdim each, reused when already so) and alphaRho holds psia^dagger rho;
// without, a new array as before
MatrixXcd* calCoherenceGcceInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, Output* op, MatrixXcd* into,
                                Eigen::RowVectorXcd* alphaRhoInto, GcceBatch* batchCollector){
    const double profileStart = cnf->profile ? MPI_Wtime() : 0.0;

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int ntotspin = nqubit + nspin;

    int qdim = QubitArray_dim(qa);
    int bdim = BathArray_dim(ba);
    int totdim = qdim*bdim;

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // A25: one qubit with bath spins: Htot and Hq from buffers kept from one cluster to the next (small_hamiltonian.cpp; the
    // same bits as the general functions below, which stay for a cluster without bath spins and for more than one qubit)
    static MatrixXcd HtotKept, HqKept;
    MatrixXcd HtotOwn, HqOwn;
    MatrixXcd** qsigmas = nullptr;
    const bool kept = nspin > 0 && nqubit == 1;
    // A2: propagator = sparse from sparse_min_dim on (evolution = vector; main.cpp refuses the rest): the Hamiltonian in
    // compressed rows from its local terms, no dense Htot (HtotKept is not written and not read on this path)
    const bool sparse = sparsePropagatorFor(cnf,totdim);
    GcceSparseWork sw;
    if (kept && sparse){
        LocalHamiltonian terms;
        gcceClusterTerms(qa,ba,cnf,terms,HqKept);
        csrFromLocal(terms,sw.H);
    }else if (kept){
        gcceClusterHamiltonian(qa,ba,cnf,HtotKept,HqKept);
    }else{
        // Declare the total Hamiltonian
        HtotOwn = MatrixXcd::Zero(totdim,totdim);

        // Pauli matrices of Qubit
        qsigmas = QubitArray_PauliOperators(qa);

        // Qubit Hamiltonian
        HqOwn = HamilQubit(qa,ba,qsigmas,cnf);
        MatrixXcd Hq_expand = kron(HqOwn,MatrixXcd::Identity(bdim,bdim));

        // Bath Hamiltonian
        if (nspin > 0){
            MatrixXcd** bsigmas = BathArray_PauliOperators(ba);
            MatrixXcd Hb = HamilBath(ba,bsigmas,cnf);
            MatrixXcd Hb_expand = kron(MatrixXcd::Identity(qdim,qdim),Hb);

            // Qubit - Bath Interaction Hamiltonian
            MatrixXcd Hqb = HamilQubitBath(qa,ba,qsigmas,bsigmas,cnf);

            // Total Hamiltonian
            HtotOwn = Hq_expand + Hb_expand + Hqb;

            for (int i=0; i<nspin; i++){
                delete[] bsigmas[i];
            }
            delete[] bsigmas;

        }else{
            HtotOwn = HqOwn;
        }
    }
    MatrixXcd& Htot = kept ? HtotKept : HtotOwn;
    MatrixXcd& Hq = kept ? HqKept : HqOwn;
    if (sparse && !kept) csrFromDense(HtotOwn,sw.H);   // no bath spins: Htot is the qubit's

    ////////////////////////////////
    // Density matrix
    ////////////////////////////////
    const double profileBuildEnd = cnf->profile ? MPI_Wtime() : 0.0;
    // Check if this is the ensemble calculation
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = false;
    if (nstate == 0){isEnsemble = true;}

    // evolution=vector propagates |psi> instead of rho, which needs rho0 to BE a pure
    // state. QubitArray_Rho0 always is (psi0*psi0^dag) and BathArray_Rho0 is too except
    // for the ensemble average, so the requirement is nstate > 0. Refuse rather than
    // silently falling back -- a run that quietly used a different algorithm than the
    // one it was asked for is worse than a run that stops.
    const bool useVector = (strcasecmp(Config_getEvolution(cnf),"vector")==0);
    const bool recursivePhases = (strcasecmp(cnf->evolution_update,"recursive")==0);   // A16 b, A18; propagator=eigen (main.cpp)
    if (useVector){
        if (isEnsemble){
            fprintf(stderr,"Error : evolution=vector needs a pure rho0, but nstate=0 "
                           "(ensemble) makes the bath density matrix mixed. "
                           "Use nstate > 0, or evolution=matrix.\n");
            exit(1);
        }
        if (strcasecmp(Config_getPropagator(cnf),"expm")==0){
            fprintf(stderr,"Error : evolution=vector requires propagator=eigen "
                           "(with expm the per-step matrix exponential dominates and "
                           "there is nothing to gain).\n");
            exit(1);
        }
    }

    MatrixXcd rho0;
    Eigen::VectorXcd psi0;
    if (useVector){
        // BathArray_Psi0 leaves its result empty when nspin==0; the bath factor is the
        // 1x1 identity there, so the qubit state is already the whole state.
        psi0 = (nspin > 0) ? Eigen::VectorXcd(kron(QubitArray_getPsi0(qa),BathArray_Psi0(ba)))
                           : Eigen::VectorXcd(QubitArray_getPsi0(qa));
    }else{
        MatrixXcd qrho0 = QubitArray_Rho0(qa);
        MatrixXcd brho0 = BathArray_Rho0(ba, isEnsemble);
        rho0 = kron(qrho0,brho0);
    }

    ////////////////////////////////
    // Propagation
    ////////////////////////////////
    // Quantity
    char* quantity = Config_getQuantity(cnf);

    // alpha, beta
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    // Alloc the result variable
    int nstep     = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);
    double tfree  = 0.0;
    MatrixXcd* result = into ? into : new MatrixXcd[nstep];
    Eigen::RowVectorXcd alphaRhoOwn;
    Eigen::RowVectorXcd& alphaRho = alphaRhoInto ? *alphaRhoInto : alphaRhoOwn;

    MatrixXcd* result_tot = nullptr; // Density matrix before trace out (total density matrix)
    if (strcasecmp(op->savemode,"cluster_density_matrix")==0){
        result_tot = new MatrixXcd[nstep]; //!
    }

    MatrixXcd* Upulses = new MatrixXcd[pulse->npulse];
    if (pulse->detuning_factor == 1.0){
        // A25: without pulseiter they depend on psia, psib, the angles and the axes only - formed again when one changes
        if (!Pulse_getPulseiter(pulse)) cachedUpulses(Upulses, qa, pulse);
        else                            calUpulses(Upulses, qa, pulse);
        //calTDUpulses(Upulses, qa, pulse, Hq);
    } else{
        calTDUpulses(Upulses, qa, pulse, Hq);
    }

    // Everything fixed for this cluster: the eigendecomposition that replaces a
    // matrix exponential on every step.
    GcceWork gw;
    bool useCuda = computeUseCuda(cnf,totdim);
#ifdef CCEX_USE_CUDA
    std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> gpu(nullptr,gpuGcceDestroy);
    if (useCuda){
        validateGcceHamiltonian(Htot);
        std::vector<doublec> pulseData((size_t)pulse->npulse*qdim*qdim);
        for (int p=0; p<pulse->npulse; ++p)
            std::copy(Upulses[p].data(),Upulses[p].data()+qdim*qdim,
                      pulseData.data()+(size_t)p*qdim*qdim);
        try {
            gpu.reset(gpuGcceCreate(totdim,qdim,useVector,std::min(cnf->gpu_batch_steps,std::max(1,nstep)),
                Htot.data(),useVector ? psi0.data() : rho0.data(),pulse->npulse,pulseData.data(),cnf->gcce_refinement,cnf->gcce_refinement_max_sweeps));
        } catch (const GpuOutOfMemory& e){
            // backend auto: this cluster on the CPU (one by one), or, in a cluster batch, the batch given up
            if (!computeCudaOutOfMemory(cnf,totdim,e.what())){
                fprintf(stderr,"Error: CUDA gCCE setup: %s\n",e.what());
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            if (batchCollector) batchCollector->outOfMemory=true;
            else useCuda=false;
        } catch (const std::exception& e){
            fprintf(stderr,"Error: CUDA gCCE setup: %s\n",e.what());
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
#endif
    // the CPU's work for the cluster: on the CPU from the start, or after the GPU ran out of memory (backend auto)
    const auto prepareOnCpu=[&](){
    if (sparse){
        if (!useVector){
            fprintf(stderr,"Error : propagator=sparse needs evolution=vector\n");
            ccexFatal();
        }
        double asym, scale;
        if (!csrFiniteHermitian(sw.H,asym,scale)){
            fprintf(stderr,"Error : propagator=sparse : the gCCE Hamiltonian (dim %d) is not finite and Hermitian "
                           "(max|H-H^dag| = %.3e, max|H| = %.3e)\n",totdim,asym,scale);
            ccexFatal();
        }
        sw.qdim = qdim; sw.bdim = bdim;
        sw.levels = Eigen::SelfAdjointEigenSolver<MatrixXcd>(Hq).eigenvectors();   // the parts of GcceSparseWork::evolve
    }else{
        // eigen (also below sparse_min_dim with propagator = sparse) or expm
        prepareGcceWork(&gw,Htot,strcasecmp(Config_getPropagator(cnf),"expm")!=0,
                        cnf->gcce_refinement,cnf->gcce_refinement_max_sweeps);
        prepareGccePulses(&gw,Upulses,pulse->npulse,bdim);
        applyBathPulsesGcce(&gw,Upulses,ba,pulse);           // DEER: kron(U_q, U_b) where a bath spin is driven
    }
    };
    if (!useCuda) prepareOnCpu();
    const double profileSetupEnd = cnf->profile ? MPI_Wtime() : 0.0;

#ifdef CCEX_USE_CUDA
    if(batchCollector){
        if(!useCuda || result_tot){
            fprintf(stderr,"Error: gCCE cluster batching requires CUDA and reduced output\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        batchCollector->owned.emplace_back(gpu.release(),gpuGcceDestroy);
        batchCollector->works.push_back(batchCollector->owned.back().get());
        batchCollector->results.emplace_back(result);
        batchCollector->alpha.push_back(psia);batchCollector->beta.push_back(psib);
        batchCollector->build.push_back(profileBuildEnd-profileStart);
        batchCollector->setup.push_back(profileSetupEnd-profileBuildEnd);
        if(qsigmas){for(int i=0;i<nqubit;++i)delete[] qsigmas[i];delete[] qsigmas;}
        delete[] Upulses;
        return result;
    }
    std::vector<MatrixXcd> gpuReduced;
    if (useCuda){
        gpuReduced.resize(nstep);
        const int batch = std::min(cnf->gpu_batch_steps,std::max(1,nstep));
        std::vector<doublec> reduced((size_t)batch*qdim*qdim);
        std::vector<doublec> full(result_tot ? (size_t)batch*totdim*totdim : 0);
        std::vector<double> times(batch), delays(pulse->npulse+1);
        // Use the same representative interval as the CPU's delay cache.
        for (int p=0; p<=pulse->npulse; ++p)
            delays[p]=pulse->sequence[pulse->sequence_indices[p]][2];
        double time=0.0;
        try {
            for (int start=0; start<nstep; start+=batch){
                int count=std::min(batch,nstep-start);
                for (int t=0; t<count; ++t){ times[t]=time; time+=deltat; }
                gpuGcceEvolve(gpu.get(),count,times.data(),delays.data(),reduced.data(),
                              result_tot ? full.data() : nullptr);
                for (int t=0; t<count; ++t){
                    gpuReduced[start+t]=Eigen::Map<const MatrixXcd>(reduced.data()+(size_t)t*qdim*qdim,qdim,qdim);
                    if (result_tot) result_tot[start+t]=Eigen::Map<const MatrixXcd>(full.data()+(size_t)t*totdim*totdim,totdim,totdim);
                }
            }
        } catch (const GpuOutOfMemory& e){
            // backend auto: the whole cluster again on the CPU (what the card had is freed)
            if (!computeCudaOutOfMemory(cnf,totdim,e.what())){
                fprintf(stderr,"Error: CUDA gCCE evolution: %s\n",e.what());
                MPI_Abort(MPI_COMM_WORLD,1);
            }
            gpu.reset();
            useCuda=false;
            prepareOnCpu();
        } catch (const std::exception& e){
            fprintf(stderr,"Error: CUDA gCCE evolution: %s\n",e.what());
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
#endif

    MatrixXcd reducedRhot;
    for (int i=0; i<nstep; i++){

#ifdef CCEX_USE_CUDA
        if (useCuda){
            reducedRhot=gpuReduced[i];
        }else
#endif
        if (useVector){
            // |psi(t)> = Utot |psi0>, without ever forming Utot.
            // evolution_update = recursive: the step number, so the phases between every 16th step are carried over (A16 b);
            // exact: the call without it, as before (one call with a variable step changed the bits of a dm result, g++)
            const Eigen::VectorXcd& psit = sparse ? propagateStateGcceSparse(sw, psi0, pulse, Upulses, tfree, cnf->sparse_tol)
                                         : recursivePhases ? propagateStateGcce(psi0, pulse, tfree, &gw, i, deltat)
                                                           : propagateStateGcce(psi0, pulse, tfree, &gw);

            // Tracing the bath out of |psi><psi| is a reshape. kron() lays the state out
            // as psi(iq*bdim + ib), so reading it as a qdim x bdim matrix Psi gives
            //   (Tr_bath |psi><psi|)_{ij} = sum_b psi(i,b) conj(psi(j,b)) = (Psi Psi^dag)_{ij}
            // which is exactly what the partialtrace chain below computes term by term.
            Eigen::Map<const Eigen::Matrix<doublec,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
                Psi(psit.data(), qdim, bdim);
            reducedRhot.noalias() = Psi * Psi.adjoint();

            // cluster_density_matrix wants the full qubit+bath density matrix; from a pure state that is
            // an outer product, O(n^2), cheaper than the O(n^3) the matrix path pays.
            if (strcasecmp(op->savemode,"cluster_density_matrix")==0){
                result_tot[i] = psit * psit.adjoint();
            }
        }else{
            // evolution_update = recursive (A18): the phases carried over as for vector; exact: as before
            MatrixXcd Utot = recursivePhases ? calPropagatorGcceRecursive(pulse, tfree, &gw, i, deltat)
                                             : calPropagatorGcce(qa, Htot, pulse, tfree, Upulses, &gw);
            // Density matrix for time
            MatrixXcd rhot = Utot * rho0 * Utot.adjoint();

            // Save the density matrix
            if (strcasecmp(op->savemode,"cluster_density_matrix")==0){
                result_tot[i] = rhot;
            }

            // Trace for tha bath state
            reducedRhot = rhot;
            for (int ib=nspin-1; ib>=0; ib--){
                int bdim_i = BathArray_dimBath_i(ba, ib);
                reducedRhot = partialtrace(reducedRhot, bdim_i, bdim_i);
            }
        }

        // Get the phase for qubit' two states
        if (strcasecmp(quantity,"coherence")==0){ // Coherence
            // psia^dagger rho, then times psib: the two products Eigen made in temporaries, written into reused storage (A14)
            alphaRho.noalias() = psia.adjoint() * reducedRhot;
            result[i].resize(1,1);
            result[i].noalias() = alphaRho * psib;

        }else if (strcasecmp(quantity,"dm")==0){ // Density matrix
            result[i] = reducedRhot;

        }else{
            fprintf(stderr,"Error : Quantity is neither coherence or dm\n");
            exit(1);
        }
        // Update tFree
        tfree += deltat;
    }

    if (cnf->profile) computeProfile(cnf,"gcce",useCuda ? "cuda" : sparse ? "sparse" : "cpu",totdim,
        profileBuildEnd-profileStart,profileSetupEnd-profileBuildEnd,MPI_Wtime()-profileSetupEnd);
    if (cnf->profile && sparse)
        printf("[sparse] rank=%d method=gcce dim=%d nnz=%ld matvecs=%lld substeps=%lld rejected=%lld\n",
               rank, totdim, long(sw.H.val.size()), sw.stats.matvecs, sw.stats.substeps, sw.stats.rejected);

    // Save the density matrix
    // to op's outfile as it is: the caller names it after the cluster and the state (Output_fullOutput). Labelling it here
    // with 0..nspin+nqubit-1 gave clusters of the same size one file name (2026-10-13).
    if (strcasecmp(op->savemode,"cluster_density_matrix")==0){
        Output_save(op,result_tot,result_tot,nstep,deltat,0);
    }

    // free
    if (qsigmas){
        for (int i=0; i<nqubit; i++){
            delete[] qsigmas[i];
        }
        delete[] qsigmas;
    }
    delete[] Upulses;
    if (result_tot != nullptr){
        delete[] result_tot;
    }

    return result;
}

MatrixXcd getPauliOperator(MatrixXcd* sigmas, char axis) {
    switch (axis) {
        case 'I': return sigmas[0];
        case 'X': return sigmas[1];
        case 'Y': return sigmas[2];
        case 'Z': return sigmas[3];
        default:
            throw std::invalid_argument("Invalid axis character");
    }
}

double getPhase(char axis) {
    double PI = 3.14159265358979323846;
    switch (axis) {
        case 'X': return 0.0;
        case 'Y': return PI/2;
        default:
            throw std::invalid_argument("Invalid axis character (For Time-dependent pulse operator)");
    }
}

double getAngle(double degree) {
    return degree * (M_PI / 180.0);  // Radian
} 
/**
 * kron(Upulses[i], I_bdim) does not depend on the step, so both evolution paths
 * expand it once per cluster. The matrix path used to rebuild it, and an identity
 * for the last interval that nothing read, on every pulse of every step.
 */
void prepareGccePulses(GcceWork* w, MatrixXcd* Upulses, int npulse, int bdim){

    w->UpulsesExpanded.resize(npulse);
    for (int ipulse=0; ipulse<npulse; ipulse++){
        w->UpulsesExpanded[ipulse] = kron(Upulses[ipulse],MatrixXcd::Identity(bdim,bdim));
    }
    w->pulseBdim = bdim;
    w->startEigen.resize(0);   // a new cluster: firstFreeGcce forms M^dag |psi0> again
}

/**
 * evolution=vector. The same operator as calPropagatorGcce, applied to a state instead
 * of assembled. That function builds, for npulse pulses,
 *     Utotal = Ufree_N (Upulse_{N-1} Ufree_{N-1}) ... (Upulse_0 Ufree_0)
 * so applying it to a state runs right to left: Ufree_0 first, then alternating the
 * pulse and the next Ufree. Every application is a matrix-VECTOR product, O(n^2), where
 * assembling Utotal costs a chain of O(n^3) matrix-matrix products.
 *
 * Ufree(tau) is never formed either -- see applyUfreeGcce. Requires propagator=eigen;
 * Config_setEvolution documents why, and calCoherenceGcce enforces it.
 */
const Eigen::VectorXcd& propagateStateGcce(const Eigen::VectorXcd& psi0, Pulse* pulse, double tfree, GcceWork* w,
                                           int istep, double deltat){

    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    int     npulse           = Pulse_getNpulse(pulse);

    // A pulse sequence repeats its delays -- CPMG runs npulse-1 identical inner taus --
    // and the matrix path has always exploited that through sequence_indices. Do the
    // same here: exp(-i*lambda*tau) costs one sincos per eigenvalue, so recomputing it
    // for a tau already seen is the most expensive thing to get wrong. An identical tau
    // gives an identical phase vector, so reusing it is bit-for-bit the same as
    // recomputing it.
    // A16 (b), evolution_update = recursive: with istep >= 0 (the istep-th step of a uniform grid deltat), the phases are computed every
    // gccePhaseAnchor steps and in between are the previous step's times exp(-i lambda f deltat)
    const int gccePhaseAnchor = 16;
    if ((int)w->phases.size() < npulse+1){ w->phases.resize(npulse+1); w->phaseSteps.resize(npulse+1); }
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        if (sequence_indices[ipulse] == ipulse){
            if (istep < 0 || istep % gccePhaseAnchor == 0){
                refinedPhases(w->evals,w->evalCorrections,tfree*sequence[ipulse][2],w->phases[ipulse]);
                if (istep == 0)
                    refinedPhases(w->evals,w->evalCorrections,deltat*sequence[ipulse][2],w->phaseSteps[ipulse]);
            } else {
                w->phases[ipulse].array() *= w->phaseSteps[ipulse].array();
            }
        }
    }

    firstFreeGcce(w, psi0, w->phases[sequence_indices[0]]);   // A16, simulator_gcce_vector.cpp

    for (int ipulse=1; ipulse<npulse+1; ipulse++){
        applyPulseGcce(w, ipulse-1, w->psi, w->psiscratch);
        w->psi.swap(w->psiscratch);
        applyUfreeGcce(w, w->phases[sequence_indices[ipulse]]);
    }

    return w->psi;
}


MatrixXcd calPropagatorGcce(QubitArray* qa, const MatrixXcd& Htot, Pulse* pulse, double tfree, MatrixXcd* Upulses, GcceWork* w){

    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    int     npulse           = Pulse_getNpulse(pulse);

    /////////////////////////////////////////////////////////////
    // Propagator for total
    MatrixXcd Utotal;
    // Reused across the nstep loop instead of new/delete on every call.
    if ((int)w->Ufrees.size() < npulse+1){ w->Ufrees.resize(npulse+1); }

    // Propagator(total) =  U(tauN+1) (U_pulseN) ... (U_pulse1) U(tau1)
    // Pulse Index  0    1    2    3    4    5    6    7    8    |
    // Pulse (8#)   |____|____|____|____|____|____|____|____|____|
    //              |tau1|tau2|           .......           |tau9|
    //              |seq1|seq2|           .......           |seq9|
    // Pulse delay  0                                          tfree
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        
        // Propagator for free evolution
        double       tau = tfree * sequence[ipulse][2];
        int sameTauIndex = sequence_indices[ipulse];

        // Calculate free evolution operators (Ufree)
        if (sameTauIndex == ipulse){ // if Ufree haven't been calculated
            if (w->useEigen){
                // U = exp(-iHtau), from the decomposition prepareGcceWork already did:
                //   exp(-i*H*tau) = M * diag(exp(-i*lambda_k*tau)) * M^dag
                // Htot does not change with the step, so one decomposition serves every
                // tau instead of a matrix exponential per step.
                refinedPhases(w->evals,w->evalCorrections,tau,w->phase);
                w->scratch.noalias()        = w->M * w->phase.asDiagonal();
                w->Ufrees[ipulse].noalias() = w->scratch * w->Madj;
            }else{
                // propagator=expm : Eigen's general matrix exponential, as CCEX has
                // always done it. Kept for reproducing numbers generated before the
                // eigendecomposition path existed.
                w->Ufrees[ipulse] = ((-1.0) * doublec(0.0,1.0) * Htot * tau).exp();
            }
        }
        else{ // if Ufree have been calculated
            // get previously calculated Ufree for the same tau
            w->Ufrees[ipulse] = w->Ufrees[sameTauIndex];
        }
        const MatrixXcd& Ufree = w->Ufrees[ipulse];

        // =================================== //
        //     Calculate U_total operator      //
        // =================================== //
        if (ipulse==0 && ipulse != npulse){
            Utotal = w->UpulsesExpanded[ipulse] * Ufree;
        }
        else if (ipulse==0 && ipulse == npulse){
            Utotal = Ufree;
        }
        else if (ipulse==npulse){
            Utotal = Ufree * Utotal;
        }
        else {
            Utotal = w->UpulsesExpanded[ipulse] * Ufree * Utotal;
        }
    }
    

    return Utotal;
}

MatrixXcd cal_Upulse_when_pulseiter_isFalse(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma){
    MatrixXcd exponent = ((-1.0) * doublec(0.0,1.0) * sigma * angle / 2.0); 
    MatrixXcd Upulse = exponent.exp();
    return Upulse;
}

MatrixXcd cal_TDUpulse_when_pulseiter_isTrue(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, double angle, int ipulse) {

    int nqubit             = QubitArray_getNqubit(qa);
    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);
    char axis              = pulse->pulse_axes[ipulse];
    double Ediff           = calQubit_Energy_difference(qa, Hq);
    std::cout << "Ediff : " << Ediff << std::endl;
    std::cout << "Qubit Hamiltonian : \n" << Hq  << std::endl;

    MatrixXcd U_total = MatrixXcd::Identity(1, 1); // 초기 유니타리

    // 전체 유니타리 구성
    for (int iq = 0; iq < nqubit; iq++) {
        MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa, iq);
        MatrixXcd beta  = QubitArray_getQubit_i_beta(qa, iq);
        std::cout << "alpha: \n" << alpha << std::endl;
        std::cout << "beta : \n" << beta << std::endl;
        exit(1);
        if (alpha.size() == 0 || beta.size() == 0) {
            fprintf(stderr, "[Error] alpha or beta not set for qubit %d\n", iq);
            exit(EXIT_FAILURE);
        }

        MatrixXcd* sigmas = getGeneralPauliOperators(alpha, beta);
        MatrixXcd ox = getPauliOperator(sigmas, 'X');
        MatrixXcd oy = getPauliOperator(sigmas, 'Y');
        MatrixXcd tempTDUpulse = cal_TDUpulse(qa, Hq, pulse, ox, oy, axis, angle, Ediff, detuning_factor, pulse->pulse_time, 10000);
        U_total = kron(U_total, tempTDUpulse);

        delete[] sigmas;
    }

    return U_total;
}


MatrixXcd cal_Upulse_when_pulseiter_isTrue(QubitArray* qa, Pulse* pulse, double angle, int ipulse) {
    int nqubit = QubitArray_getNqubit(qa);
    MatrixXcd U_total = MatrixXcd::Identity(1, 1); // 초기 유니타리 (1x1 단위행렬)

    for (int iq = 0; iq < nqubit; iq++) {
        MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa, iq);
        MatrixXcd beta  = QubitArray_getQubit_i_beta(qa, iq);

        if (alpha.rows() == 0 || alpha.cols() == 0 || beta.rows() == 0 || beta.cols() == 0) {
            fprintf(stderr, "[Error] Function - cal_Upulse_when_pulseiter_isTrue: alpha or beta not set for qubit %d\n", iq);
            exit(EXIT_FAILURE);
        }

        MatrixXcd* sigmas = getGeneralPauliOperators(alpha, beta); // sigma[0]=I, [1]=X, [2]=Y, [3]=Z
        MatrixXcd  sigma = getPauliOperator(sigmas, (pulse->pulse_axes[ipulse]));

        // Expand (tensor product)
        MatrixXcd sigma_expanded = MatrixXcd::Identity(1, 1);
        for (int jq = 0; jq < nqubit; jq++) {
            int dim = QubitArray_dimQubit_i(qa, jq);
            MatrixXcd part;
            if (jq == iq) {
                part = sigma;
            } else {
                part = MatrixXcd::Identity(dim, dim);
            }
            sigma_expanded = kron(sigma_expanded, part);
        }

        MatrixXcd exponent = ((-1.0) * doublec(0.0,1.0) * angle * sigma_expanded / 2.0);
        MatrixXcd U_iq = exponent.exp();

        U_total = U_iq * U_total;

        delete[] sigmas;
    }

    return U_total;
}

MatrixXcd cal_Upulse_when_pulseiter_isTrue_whenOnlyXPulse(QubitArray* qa, Pulse* pulse, double angle){

    int nqubit = QubitArray_getNqubit(qa);
    int  qdim  = QubitArray_dim(qa);
    MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);

    for (int iq=0; iq<nqubit; iq++){
    
        MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa,iq);
        MatrixXcd beta = QubitArray_getQubit_i_beta(qa,iq);
    
        if (alpha.rows() == 0 && alpha.cols() == 0){
            fprintf(stderr,"Error : calPropagatorGcce : alpha or bata is not set\n");
            fprintf(stderr,"If pulse.per_qubit is true, you have to give alpha,beta\n");
            exit(1);
        }
    
        /////////////////////////////////////////////////////////////
        // Generate the general Pauli operators for i-th qubit
        MatrixXcd sigmaExpanded_x;
        MatrixXcd* sigmas = getGeneralPauliOperators(alpha,beta);
        /////////////////////////////////////////////////////////////
        MatrixXcd sigma;
        for (int jq=0; jq<iq; jq++){
            if (iq==jq){
                sigma = sigmas[1];
            }else{
                int jqdim = QubitArray_dimQubit_i(qa,jq);
                sigma = MatrixXcd::Identity(jqdim,jqdim);
            }
    
            if (jq==0){
                sigmaExpanded_x = sigma;
            }else{
                sigmaExpanded_x = kron(sigmaExpanded_x,sigma);
            }
        }
        /////////////////////////////////////////////////////////////
        // -i*sigma_x_i
        // MatrixXcd Upulse_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x); // Before 25.08.08
        MatrixXcd exponent_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x * M_PI / 2.0); // Changed 25.08.08
        MatrixXcd Upulse_iq = exponent_iq.exp(); // Changed 25.08.08
        Upulse = Upulse * Upulse_iq;
    }
    return Upulse;
}

void calUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse){

    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
    int     npulse        = Pulse_getNpulse(pulse);
    double* pulse_angles  = Pulse_getPulseAngles(pulse);
    char*   pulse_axes    = Pulse_getPulseAxes(pulse);
    bool    pulseiter     = Pulse_getPulseiter(pulse);
    int     qdim          = QubitArray_dim(qa);

    /////////////////////////////////////////////////////////////
    // Propagator for pulse
    //  Upulse = exp(-i*sigma_x*pi/2) -> Pi-Pulse
    for (int ipulse=0; ipulse<npulse; ipulse++){
        MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);
        double angle = getAngle(pulse_angles[ipulse]);  // Radian
        if (!pulseiter){
            MatrixXcd* sigmas = QubitArray_PauliOperator_fromPsiaPsib(qa);
            MatrixXcd  sigma = getPauliOperator(sigmas, (pulse_axes[ipulse]));
            Upulse = cal_Upulse_when_pulseiter_isFalse(qa, pulse, angle, sigma);
            delete[] sigmas;

        } else{
            if ((strcasecmp(pulse->pulsename, "CPMG") == 0) || (strcasecmp(pulse->pulsename, "HahnEcho") == 0) || (strcasecmp(pulse->pulsename, "Ramsey") == 0)){
                Upulse = cal_Upulse_when_pulseiter_isTrue_whenOnlyXPulse(qa, pulse, angle);
            }else{
                Upulse = cal_Upulse_when_pulseiter_isTrue(qa, pulse, angle, ipulse);
            }
        }
        Upulses[ipulse] = Upulse;
    }
    //return Upulses;
}

//MatrixXcd cal_Upulse_withDetuning(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma, MatrixXcd sigma_z){
//
//    printf("Function: cal_Upulse_withDetuning \n");
//    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
//    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);
//    double ptime = (pulse->pulse_time) / 1e6;
//    //std::cout << "Detun/Omega: " << detuning / angle * ptime << std::endl;
//
//    MatrixXcd factor = (detuning * ptime * sigma_z) + (angle * sigma);
//    MatrixXcd U_detuning = ((-1.0) * doublec(0.0,1.0) * factor / 2.0).exp(); 
//    return U_detuning;
//}

double calQubit_Energy_difference(QubitArray* qa, MatrixXcd Hq){
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    double Ediff = 0.0;
    int alphaIdx;
    int betaIdx;
    //std::cout << "psia: " << psia << std::endl;
    //std::cout << "psib: " << psib << std::endl;
    for (int i = 0; i < psia.rows(); ++i) {
        if (psia(i, 0) != std::complex<double>(0.0, 0.0)) {
            alphaIdx = i;
            //std::cout << "non-zero element idx for pisa: " << i << std::endl;
        }
        if (psib(i, 0) != std::complex<double>(0.0, 0.0)) {
            betaIdx = i;
            //std::cout << "non-zero element idx for psib: " << i <<std::endl;
        }
    }
    Ediff = std::abs(Hq(alphaIdx, alphaIdx).real() - Hq(betaIdx, betaIdx).real());
    //std::cout << "Ediff: " << Ediff << std::endl;
    return Ediff;
}

void calTDUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse, MatrixXcd Hq){

    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
    int     npulse        = Pulse_getNpulse(pulse);
    double* pulse_angles  = Pulse_getPulseAngles(pulse);
    char*   pulse_axes    = Pulse_getPulseAxes(pulse);
    bool    pulseiter     = Pulse_getPulseiter(pulse);
    int     qdim          = QubitArray_dim(qa);
    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);

    /////////////////////////////////////////////////////////////
    // Calculate time-dependent pulse operator. 
    //printf("Calculating pulse energy level difference! \n");
    double Ediff = calQubit_Energy_difference(qa, Hq);

    for (int ipulse=0; ipulse<npulse; ipulse++){
        MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);
        double angle = getAngle(pulse_angles[ipulse]);  // Radian
        if (!pulseiter){
            MatrixXcd* sigmas = QubitArray_PauliOperator_fromPsiaPsib(qa);
            MatrixXcd  ox = getPauliOperator(sigmas, 'X');
            MatrixXcd  oy = getPauliOperator(sigmas, 'Y');
            char axis     = pulse_axes[ipulse];
            // printf(" || Calculating Time-dependent pulse operator || \n");
            Upulse = cal_TDUpulse(qa, Hq, pulse, ox, oy, axis, angle, Ediff, detuning_factor, pulse->pulse_time, 10000);
            delete[] sigmas;

        } else{
            printf("There's code developed for pulse iter time-dependent pulse operator!\n");
            Upulse = cal_Upulse_when_pulseiter_isTrue(qa, pulse, angle, ipulse);
            exit(1);
        }
        Upulses[ipulse] = Upulse;
    }
}

MatrixXcd cal_TDUpulse(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, MatrixXcd ox, MatrixXcd oy, char axis, double pulse_angle, double Ediff, double detuning_factor, double pulse_time, int Nstep) {

    //printf("Function: cal_TDUpulse\n");

    // Omega(Rabi_frequency) = Angle * pulse_duration_time //
    //double Omega = pulse_angle / pulse->pulse_time * 1e6;
    
    double ptime = pulse_time / 1e6;
    double Omega_Rabi = pulse_angle / ptime;

    // detuning = Ediff - omega // 
    double omega    = Ediff * detuning_factor;
    double detuning = Ediff - omega; // unit : rad*kHz
    double phase    = getPhase(axis);
    double Omega    = std::sqrt(Omega_Rabi*Omega_Rabi + detuning*detuning);

    const int qdim = ox.rows();
    double dt = ptime / static_cast<double>(Nstep);

     printf("Qubit Hamiltonian : \n");
     std::cout << Hq << std::endl;
     printf("\n");
     std::cout << "Pulse angle: " << pulse_angle << "  rad"  << std::endl;
     std::cout << "Pulse Time : " << pulse_time << "  ns" << std::endl;
     //std::cout << "dt         : " << dt*1e6 << "  ns" << std::endl;
     std::cout << "w0         : " << Ediff*0.001  << "  rad*MHz" << std::endl;
     std::cout << "Rabi       : " << Omega_Rabi*0.001  << "  rad*MHz" << std::endl;
     std::cout << "w          : " << omega     *0.001  << "  rad*MHz" << std::endl;
     std::cout << "Omega      : " << Omega     *0.001  << "  rad*MHz" << std::endl;
     std::cout << "detuning   : " << detuning  *0.001  << "  rad*MHz" << std::endl;
     //std::cout << "phase      : " << phase                      << std::endl;

    MatrixXcd TDUpulse = MatrixXcd::Identity(qdim, qdim);
    //MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);
    for (int n = 0; n < Nstep; ++n) {
        double t = n * dt;
        //MatrixXcd H = Hq + Omega * cos(omega * t) * ox / 2.0;
        MatrixXcd H = Hq + Omega * ( (cos(omega*t + phase)*ox) + (sin(omega*t + phase)*oy) ) / 2.0;
        MatrixXcd U = ((-1.0) * doublec(0.0, 1.0) * H * dt).exp();
        TDUpulse = U * TDUpulse;  
    }
    std::cout << "TDUpuslse : \n" << TDUpulse << std::endl;
    //exit(1);

    return TDUpulse;
}

#ifdef CCEX_USE_CUDA
std::vector<MatrixXcd*> calCoherenceGcceCudaBatch(QubitArray* qa,
        const std::vector<BathArray*>& baths,Config* cnf,Pulse* pulse,Output* output){
    try {
        GcceBatch work;
        const int n=int(baths.size()),q=QubitArray_dim(qa),nstep=cnf->nstep;
        const int batch=std::min(cnf->gpu_batch_steps,std::max(1,nstep));
        for(BathArray* bath:baths){
            calCoherenceGcce(qa,bath,cnf,pulse,output,&work);
            if(work.outOfMemory)return {};   // backend auto: the caller computes the clusters one by one (on the CPU)
        }
        std::vector<doublec> reduced((size_t)n*batch*q*q);
        std::vector<double> times(batch),delays(pulse->npulse+1);
        for(int p=0;p<=pulse->npulse;++p)delays[p]=pulse->sequence[pulse->sequence_indices[p]][2];
        const double start=cnf->profile ? MPI_Wtime() : 0;
        double time=0;
        for(int first=0;first<nstep;first+=batch){
            const int count=std::min(batch,nstep-first);
            for(int t=0;t<count;++t){times[t]=time;time+=(double)cnf->deltat;}
            gpuGcceEvolveBatch(work.works.data(),n,count,times.data(),delays.data(),reduced.data());
            for(int c=0;c<n;++c)for(int t=0;t<count;++t){
                const auto rho=Eigen::Map<const MatrixXcd>(reduced.data()+((size_t)c*count+t)*q*q,q,q);
                if(!strcasecmp(cnf->quantity,"coherence"))
                    work.results[c][first+t]=work.alpha[c].adjoint()*rho*work.beta[c];
                else work.results[c][first+t]=rho;
            }
        }
        const double evolution=(cnf->profile ? MPI_Wtime() : 0)-start;
        const int dimension=q*BathArray_dim(baths[0]);
        std::vector<MatrixXcd*> results;
        for(int c=0;c<n;++c){
            computeProfile(cnf,"gcce","cuda",dimension,work.build[c],work.setup[c],evolution/n);
            results.push_back(work.results[c].release());
        }
        computeClusterBatch(cnf,dimension,n);
        return results;
    }catch(const GpuOutOfMemory& error){
        if(computeCudaOutOfMemory(cnf,QubitArray_dim(qa)*BathArray_dim(baths[0]),error.what()))return {};
        fprintf(stderr,"Error: CUDA gCCE cluster batch: %s\n",error.what());
        MPI_Abort(MPI_COMM_WORLD,1);return {};
    }catch(const std::exception& error){
        fprintf(stderr,"Error: CUDA gCCE cluster batch: %s\n",error.what());
        MPI_Abort(MPI_COMM_WORLD,1);return {};
    }
}
#endif
