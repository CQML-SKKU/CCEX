#include "../../include/calculation/simulator.h"
#include "../../include/calculation/hamiltonian.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/hamiltonian_check.h"
#include "../../include/calculation/small_hamiltonian.h"
#include "../../include/calculation/sparse_propagator.h"
namespace {
std::pair<MatrixXcd,MatrixXcd> cceHamiltonians(QubitArray* qa,BathArray* ba,Config* cnf) {
    const int nspin=BathArray_getNspin(ba);
    auto cleanup=[nspin](MatrixXcd** p){if(p){for(int i=0;i<nspin;++i)delete[] p[i];delete[] p;}};
    std::unique_ptr<MatrixXcd*,decltype(cleanup)> sigmas(BathArray_PauliOperators(ba),cleanup);
    MatrixXcd bath=HamilBath(ba,sigmas.get(),cnf);
    std::unique_ptr<MatrixXcd[]> interaction(HamilQubitBathSecularApp(qa,ba,sigmas.get(),cnf));
    interaction[0]=bath+interaction[0];   // in place: two bdim x bdim matrices fewer
    interaction[1]=bath+interaction[1];
    return {std::move(interaction[0]),std::move(interaction[1])};
}
#ifdef CCEX_USE_CUDA
void validateCceCudaHamiltonians(const MatrixXcd& a,const MatrixXcd& b) {
    for(const MatrixXcd* h:{&a,&b})
        if(!finiteAndHermitian(*h))
            throw std::runtime_error("cCCE Hamiltonian must be finite and Hermitian");
}
void cceCudaDurations(const Pulse* pulse,bool ensemble,int count,double dt,
                      double& time,std::vector<double>& durations) {
    const int npulse=pulse->npulse;
    for(int t=0;t<count;++t) {
        for(int p=0;p<=npulse;++p) {
            double duration;
            if(ensemble) {
                if(!npulse)duration=time;
                else {
                    const double tau=time/(2.0*npulse);
                    duration=(p==0 || p==npulse) ? tau : 2.0*tau;
                }
            } else duration=time*pulse->sequence[pulse->sequence_indices[p]][2];
            durations[size_t(p)*count+t]=duration;
        }
        time+=dt;
    }
}
#endif
}

#ifdef CCEX_USE_CUDA
#include "../../include/calculation/gpu/gpu_backend.h"
#include <memory>

static MatrixXcd* calCoherenceCceCuda(const MatrixXcd& ha, const MatrixXcd& hb,
    BathArray* ba, Config* cnf, Pulse* pulse, double start, double buildEnd){
    try {
        validateCceCudaHamiltonians(ha,hb);
        if (strcasecmp(Config_getQuantity(cnf),"coherence"))
            throw std::runtime_error("cCCE requires quantity=coherence");
        const bool ensemble=Config_getNstate(cnf)==0;
        const int d=ha.rows(), nstep=Config_getNstep(cnf), npulse=Pulse_getNpulse(pulse);
        const int batch=std::min(cnf->gpu_batch_steps,std::max(1,nstep));
        MatrixXcd psi;
        if (!ensemble) psi=BathArray_Psi0(ba);
        std::unique_ptr<GpuCce,void(*)(GpuCce*)> work(
            gpuCceCreate(d,ensemble,batch,npulse,ha.data(),hb.data(),ensemble ? nullptr : psi.data()),gpuCceDestroy);
        const double setupEnd=cnf->profile ? MPI_Wtime() : 0.0;
        std::unique_ptr<MatrixXcd[]> result(new MatrixXcd[nstep]);
        std::vector<doublec> coherence(batch);
        std::vector<double> durations((size_t)batch*(npulse+1));
        double time=0.0, deltat=(double)Config_getDeltat(cnf);
        for (int first=0; first<nstep; first+=batch){
            const int count=std::min(batch,nstep-first);
            cceCudaDurations(pulse,ensemble,count,deltat,time,durations);
            gpuCceEvolve(work.get(),count,durations.data(),coherence.data());
            for (int t=0; t<count; ++t) result[first+t]=MatrixXcd::Constant(1,1,coherence[t]);
        }
        if (cnf->profile) computeProfile(cnf,"cce","cuda",d,buildEnd-start,setupEnd-buildEnd,MPI_Wtime()-setupEnd);
        return result.release();
    } catch (const GpuOutOfMemory& e){
        // backend auto: nullptr, and the caller computes this cluster on the CPU (what the card had is freed)
        if (computeCudaOutOfMemory(cnf,int(ha.rows()),e.what())) return nullptr;
        fprintf(stderr,"Error: CUDA cCCE: %s\n",e.what());
        MPI_Abort(MPI_COMM_WORLD,1);
        return nullptr;
    } catch (const std::exception& e){
        fprintf(stderr,"Error: CUDA cCCE: %s\n",e.what());
        MPI_Abort(MPI_COMM_WORLD,1);
        return nullptr;
    }
}
#endif

namespace {
// A16 (b), evolution_update = recursive: the phases of the time steps between two of these come from the previous step's
// (phaseAnchor - 1 products at most)
constexpr int phaseAnchor = 16;
struct SingleStateRun {
    const Eigen::VectorXcd& valUp; const Eigen::VectorXcd& valDown;
    double** sequence; const int* sequenceIndices; int npulse, nstep; double deltat;
    MatrixXcd* result; CceWorkspace& ws;
};
// One state: only psi^dagger (UdownDagger Uup) psi is needed, and both chains are products of M or M^dagger with
// diagonal phases. Applied to the state from its two ends they are matrix-vector products - D^2 per factor instead of
// D^3 - with the same factors as the matrix branch of calCoherenceCceInto:
//   UdownDagger = M_down D_d0^* [M^dagger D_u^* | M D_d^*] ...
//   Uup         = ... [M D_d | M^dagger D_u] M^dagger D_u0 M_up^dagger
// Mat/Row/Col are MatrixXcd/RowVectorXcd/VectorXcd or their fixed-size forms (A16). The phases of an interval with an
// earlier equal one (sequenceIndices) are read from that one instead of copied (A16; the same values).
// out = in * A. A fixed size makes Eigen use its unrolled product, which sums in another order than the general kernel
// (measured with icpx and AVX-512 on idm; not with g++ and AVX2), so for a fixed size the general kernel is called
// directly with the arguments Eigen's dispatch gives it: the same bits, without the dispatch. Dynamic: as before.
template<class Mat,class Row>
inline void rowTimes(const Row& in,const Mat& A,Row& out){
    if constexpr (Mat::RowsAtCompileTime == Eigen::Dynamic){
        out.noalias() = in * A;
    } else {
        typedef Eigen::internal::const_blas_data_mapper<doublec,Eigen::Index,Eigen::RowMajor> Lhs;
        typedef Eigen::internal::const_blas_data_mapper<doublec,Eigen::Index,Eigen::ColMajor> Rhs;
        const Eigen::Index n = A.rows();
        out.setZero();
        Eigen::internal::general_matrix_vector_product<Eigen::Index,doublec,Lhs,Eigen::RowMajor,false,doublec,Rhs,false>::run(
            n, n, Lhs(A.data(),n), Rhs(in.data(),1), out.data(), 1, doublec(1.0,0.0));
    }
}
// out = A * in for evolution_update = exact, the same way (P1-04 of review 20261012-233058-452278): Eigen's unrolled
// fixed-size product sums and signs zeros otherwise than the general kernel (dynamic (1, -0), fixed (1, +0); last bits
// without -march=native), so a fixed size calls the column-major kernel with the arguments of Eigen's dispatch - out starts
// from the dispatch's setZero. recursive keeps the unrolled product (faster; its fixed size differs from the dynamic path
// in the last bits anyway: g++ -march=native vectorizes the phase step for a run-time size only).
template<class Mat,class Col>
inline void colTimes(const Mat& A,const Col& in,Col& out){
    if constexpr (Mat::RowsAtCompileTime == Eigen::Dynamic){
        out.noalias() = A * in;
    } else {
        typedef Eigen::internal::const_blas_data_mapper<doublec,Eigen::Index,Eigen::ColMajor> Lhs;
        typedef Eigen::internal::const_blas_data_mapper<doublec,Eigen::Index,Eigen::RowMajor> Rhs;
        const Eigen::Index n = A.rows();
        out.setZero();
        Eigen::internal::general_matrix_vector_product<Eigen::Index,doublec,Lhs,Eigen::ColMajor,false,doublec,Rhs,false>::run(
            n, n, Lhs(A.data(),n), Rhs(in.data(),1), out.data(), 1, doublec(1.0,0.0));
    }
}
// row * col for exact: for a fixed size the run-time-sized reduction of the dynamic path through maps (same order)
template<class Row,class Col>
inline doublec rowDotCol(const Row& row,const Col& col){
    if constexpr (Row::SizeAtCompileTime == Eigen::Dynamic){
        return (row * col)(0,0);
    } else {
        const Eigen::Map<const Eigen::RowVectorXcd> r(row.data(),row.size());
        const Eigen::Map<const Eigen::VectorXcd> c(col.data(),col.size());
        return (r * c)(0,0);
    }
}
// Recursive: evolution_update = recursive (A16 b); exact: sin and cos at every step. A template argument, not a run-time
// flag: with the flag g++ compiled both modes differently (last bits of diverged CCE4/pCCE points for recursive)
template<bool Recursive,class Mat,class Row,class Col>
void evolveSingleState(const SingleStateRun& r,const Mat& M,const Mat& Madj,const Row& rowStart,const Col& colStart,
                       Row& row,Row& rowNext,Col& col,Col& colNext,Col& scaled){
    typedef Eigen::Map<const Col> Phase;
    const int bdim=int(M.rows()), npulse=r.npulse;
    // per interval: bdim up, bdim down, and with Recursive their factors for one step (A16 b)
    constexpr int slots = Recursive ? 4 : 2;
    std::vector<doublec>& buf=r.ws.phases;
    if (buf.size() < size_t(slots*bdim*(npulse+1))) buf.resize(size_t(slots*bdim*(npulse+1)));
    double tfree=0.0;
    for (int i = 0; i < r.nstep; i++) {
        for (int ipulse=0; ipulse<npulse+1; ipulse++){
            const int same = r.sequenceIndices[ipulse];
            doublec* up = buf.data()+size_t(slots*bdim*same);
            doublec* down = up+bdim;
            doublec* stepUp = Recursive ? up+2*bdim : nullptr;     // exact holds 2 slots: no pointer past them
            doublec* stepDown = Recursive ? up+3*bdim : nullptr;
            if (same == ipulse){
                // exact: exp(-i E tau) from sin and cos at every step. recursive (A16 b): so every phaseAnchor steps, and
                // in between the previous step's value times exp(-i E f deltat) - tau grows by f deltat per step (f = this
                // interval's fraction)
                if (!Recursive || i % phaseAnchor == 0){
                    const double tau = tfree * r.sequence[ipulse][2];
                    for (int k=0; k<bdim; k++){
                        const double a_up = r.valUp(k).real()*tau, a_down = r.valDown(k).real()*tau;
                        up[k]   = doublec(std::cos(a_up),   -std::sin(a_up));
                        down[k] = doublec(std::cos(a_down), -std::sin(a_down));
                    }
                    if (Recursive && i == 0){
                        const double dtau = r.deltat * r.sequence[ipulse][2];
                        for (int k=0; k<bdim; k++){
                            const double a_up = r.valUp(k).real()*dtau, a_down = r.valDown(k).real()*dtau;
                            stepUp[k]   = doublec(std::cos(a_up),   -std::sin(a_up));
                            stepDown[k] = doublec(std::cos(a_down), -std::sin(a_down));
                        }
                    }
                } else {
                    for (int k=0; k<bdim; k++){ up[k] *= stepUp[k]; down[k] *= stepDown[k]; }
                }
            }
            const Phase Uu(up,bdim), Ud(down,bdim);
            if (ipulse==0){
                row = rowStart.array() * Ud.adjoint().array();              // psi^dagger M_down D_d^*
                scaled = (Uu.array() * colStart.array()).matrix();
                if constexpr (Recursive) col.noalias() = Madj * scaled;           // M^dagger D_u M_up^dagger psi
                else colTimes(Madj,scaled,col);
            } else if (ipulse%2==1){
                rowTimes(row,Madj,rowNext); row = rowNext.array() * Uu.adjoint().array();   // * M^dagger D_u^*
                scaled = (Ud.array() * col.array()).matrix();
                if constexpr (Recursive) colNext.noalias() = M * scaled;          // M D_d *
                else colTimes(M,scaled,colNext);
                col.swap(colNext);
            } else {
                rowTimes(row,M,rowNext);    row = rowNext.array() * Ud.adjoint().array();   // * M D_d^*
                scaled = (Uu.array() * col.array()).matrix();
                if constexpr (Recursive) colNext.noalias() = Madj * scaled;       // M^dagger D_u *
                else colTimes(Madj,scaled,colNext);
                col.swap(colNext);
            }
        }
        r.result[i].resize(1,1);                         // no allocation when the caller's buffer already holds 1 x 1
        if constexpr (Recursive) r.result[i](0,0) = (row * col)(0,0);
        else r.result[i](0,0) = rowDotCol(row,col);
        tfree += r.deltat;
    }
}
}

MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse){
    MatrixXcd* result = new MatrixXcd[Config_getNstep(cnf)];
    CceWorkspace ws;
    calCoherenceCceInto(qa, ba, cnf, pulse, result, ws);
    return result;
}

void calCoherenceCceInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, MatrixXcd* result, CceWorkspace& ws){
    
    double stime = MPI_Wtime();
    double etime = 0.0;

    int bdim = BathArray_dim(ba);

    int nqubit = QubitArray_getNqubit(qa);
    if (nqubit > 1){
        fprintf(stderr,"Error : calCoherenceCce : nqubit > 1, current nqubit = %d\n",nqubit);
        exit(1);
    }

    // A2: propagator = sparse from sparse_min_dim on (one state; main.cpp refuses the rest) - no dense Hamiltonian
    if (sparsePropagatorFor(cnf,bdim)){
        calCoherenceCceSparse(qa,ba,cnf,pulse,result,stime);
        return;
    }

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // A25: the two Hamiltonians in ws, without the general functions' allocations (small_hamiltonian.cpp; the same bits)
    cceClusterHamiltonians(qa,ba,cnf,ws.Halpha,ws.Hbeta);

#ifdef CCEX_USE_CUDA
    if (computeUseCuda(cnf,bdim)){
        const double buildEnd=cnf->profile ? MPI_Wtime() : 0.0;
        MatrixXcd* cuda=calCoherenceCceCuda(ws.Halpha,ws.Hbeta,ba,cnf,pulse,stime,buildEnd);
        if (cuda){                     // nullptr: out of GPU memory under backend auto, on the CPU below
            for (int i=0; i<Config_getNstep(cnf); i++) result[i] = cuda[i];
            delete[] cuda;
            return;
        }
    }
#endif

    const double profileBuildEnd = cnf->profile ? MPI_Wtime() : 0.0;

    // Check if this is the ensemble calculation; one state: the initial state of the bath
    const bool isEnsemble = Config_getNstate(cnf) == 0;
    MatrixXcd psi0;
    if (!isEnsemble){
        psi0 = BathArray_Psi0(ba);
    }
    if (!cceSetupSmallAlone(ws, isEnsemble ? nullptr : &psi0)) cceSetup(ws, isEnsemble ? nullptr : &psi0);   // A42

    // Quantity
    char* quantity = Config_getQuantity(cnf);
    if (strcasecmp(quantity,"coherence")!=0){
        fprintf(stderr,"Error : calCoherenceCce : quantity is not coherence\n");
        exit(1);
    }
    const double profileSetupEnd = cnf->profile ? MPI_Wtime() : 0.0;

    cceEvolve(cnf, pulse, result, ws);

    // free
    if (cnf->profile) computeProfile(cnf,"cce","cpu",bdim,
        profileBuildEnd-stime,profileSetupEnd-profileBuildEnd,MPI_Wtime()-profileSetupEnd);
}

void cceSetup(CceWorkspace& ws, const MatrixXcd* psi){
    ////////////////////////////////
    // Diagonalize the Hamiltonian
    ////////////////////////////////
    // the solvers, M and the eigenvalues live in ws (A14: no allocation once the dimension has been seen)
    ws.up.compute(ws.Halpha);
    ws.down.compute(ws.Hbeta);
    const MatrixXcd& M_up   = ws.up.eigenvectors();
    const MatrixXcd& M_down = ws.down.eigenvectors();
    // M_up is the eigenvector matrix of a Hermitian matrix, so its inverse is its
    // adjoint: no LU inverse (7-24 % of this setup at dim 256-2048, 2026-09-30)
    MatrixXcd& M = ws.M;
    M.noalias() = M_up.adjoint()*M_down;
    ws.valUp = ws.up.eigenvalues();
    ws.valDown = ws.down.eigenvalues();
    if (psi){
        // One state: the two ends of the chain (evolveSingleState)
        const MatrixXcd& psi0 = *psi;
        ws.rowStart.noalias() = psi0.adjoint() * M_down;                        // psi^dagger M_down
        ws.colStart.noalias() = M_up.adjoint() * psi0;                          // M_up^dagger psi
        ws.Madj = M.adjoint();
    }
}

void cceEvolve(Config* cnf, Pulse* pulse, MatrixXcd* result, CceWorkspace& ws){
    const int bdim = int(ws.M.rows());
    const MatrixXcd& M_up   = ws.up.eigenvectors();
    const MatrixXcd& M_down = ws.down.eigenvectors();
    const MatrixXcd& M = ws.M;
    const Eigen::VectorXcd& es_up_Val = ws.valUp;
    const Eigen::VectorXcd& es_down_Val = ws.valDown;
    const bool isEnsemble = Config_getNstate(cnf) == 0;

    ////////////////////////////////
    // Propagation
    ////////////////////////////////

    // Npulse
    int npulse = Pulse_getNpulse(pulse);
    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);

    // Alloc the result variable
    int nstep = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);

    double tfree = 0.0;

    // The ensemble block below is written for CPMG: it takes tau = t/(2*npulse) and
    // never reads the sequence. A Manual sequence has its own intervals, so it goes
    // through the interval-by-interval block instead, and ends in a trace over the
    // bath rather than an expectation value in one state.
    const bool manualIntervals = strcasecmp(pulse->pulsename,"Manual")==0;

    if (isEnsemble && strcasecmp(cnf->evolution_update,"recursive")==0){
        // A18: both ensemble branches below with the phases carried over from step to step (a file of its own, so that
        // the exact branches here compile as before)
        evolveEnsembleRecursiveCce(M,M_up,M_down,es_up_Val,es_down_Val,pulse,manualIntervals,nstep,deltat,result);

    }else if (isEnsemble && !manualIntervals && npulse==0){
        // X8: Ramsey reads only Tr[M U_down M^dagger U_up^dagger] = U_up^dagger |M|^2 U_down (D^2 per step instead of a
        // D^3 product); in simulator_cce_recursive.cpp, so that the branches here compile as before (written here, or
        // with the old Ramsey branch below removed, the single-state branch moved in the last bits: p1_pcce2_sk2)
        evolveEnsembleRamseyCce(M,es_up_Val,es_down_Val,nstep,deltat,result);

    }else if (isEnsemble && !manualIntervals){

        MatrixXcd U_up; MatrixXcd U_down;
        MatrixXcd U_up2; MatrixXcd U_down2;

        MatrixXcd L_mat; 
        MatrixXcd L_a_even; MatrixXcd L_b_even;
        MatrixXcd L_a_odd; MatrixXcd L_b_odd;
        MatrixXcd L_a2_even; MatrixXcd L_b2_even;
        MatrixXcd L_a2_odd; MatrixXcd L_b2_odd;
        MatrixXcd sum_L_a; MatrixXcd sum_L_b;

        // Hahn echo: the chain is Tr[X Du X^dagger Du^*] with X = M Dd M^dagger (the
        // second half, M Dd^* M^dagger, is X^dagger exactly), which is
        // sum_ij |X_ij|^2 Du_j conj(Du_i): one D^3 product per step instead of three.
        MatrixXcd X(bdim,bdim);
        Eigen::MatrixXd X2(bdim,bdim);
        doublec overlap;

        for (int i = 0; i < nstep; i++) {
            //Note that the free evolution time between the pulses
            // is t_free/2.0 for the Hahn-echo sequence.

            if(npulse==1){
                const double tau = tfree/2.0;
                Eigen::VectorXcd Uu(bdim), Ud(bdim);
                for (int k=0; k<bdim; k++){
                    const double a_up = es_up_Val(k).real()*tau, a_down = es_down_Val(k).real()*tau;
                    Uu(k) = doublec(std::cos(a_up),   -std::sin(a_up));
                    Ud(k) = doublec(std::cos(a_down), -std::sin(a_down));
                }
                X.noalias() = (M.array().rowwise() * Ud.transpose().array()).matrix() * M.adjoint();
                X2 = X.cwiseAbs2();
                overlap = (Uu.adjoint() * X2.cast<doublec>() * Uu)(0,0);
            }
            //T2_star (Ramsey): no longer reached (X8, evolveEnsembleRamseyCce above), but kept: without it g++ compiles the
            // single-state branch of this function differently (check_numerics p1_pcce2_sk2 moved in the last bits)
            else if(npulse==0){

                Eigen::VectorXcd es_up_Val_tau      =es_up_Val*tfree;
                Eigen::VectorXcd es_down_Val_tau    =es_down_Val*tfree;
                Eigen::VectorXcd U_up               =es_up_Val_tau.array().cos() + doublec(0,-1)*es_up_Val_tau.array().sin();
                Eigen::VectorXcd U_down             =es_down_Val_tau.array().cos() + doublec(0,-1)*es_down_Val_tau.array().sin();

                // Coherence function
                //L_mat     =  M * U_down * M.adjoint()* U_up.adjoint();
                L_mat     =  (M.array().rowwise() * U_down.transpose().eval().array() ).matrix() 
                            * (M.adjoint().array().rowwise() * U_up.adjoint().array()).matrix();



            }else{
            //T2 (CPMG sequence)
            //if you want to use the Hahn-echo, set the N=1
                double tau = tfree/(2*npulse);

                //define diagonal matrix to vector 
                Eigen::VectorXcd es_up_Val_tau      = es_up_Val*tau;
                Eigen::VectorXcd es_down_Val_tau    = es_down_Val*tau;
                Eigen::VectorXcd U_up               = es_up_Val_tau.array().cos() + doublec(0,-1)*es_up_Val_tau.array().sin();
                Eigen::VectorXcd U_down             = es_down_Val_tau.array().cos() + doublec(0,-1)*es_down_Val_tau.array().sin();

                Eigen::VectorXcd es_up_Val_tau2     = es_up_Val*2*tau;
                Eigen::VectorXcd es_down_Val_tau2   = es_down_Val*2*tau;
                Eigen::VectorXcd U_up2              = es_up_Val_tau2.array().cos() + doublec(0,-1)*es_up_Val_tau2.array().sin();
                Eigen::VectorXcd U_down2            = es_down_Val_tau2.array().cos() + doublec(0,-1)*es_down_Val_tau2.array().sin();
                

                //pulse multiple
                L_a_even = M.array().rowwise()           * U_down.transpose().eval().array();   //M * U_down
                L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array();       //M.adjoint * U_up.adjoint

                if(npulse % 2 != 0){
                    L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().eval().array();  //M.adjoint() * U_up
                    L_b_odd = M.array().rowwise()           * U_down.adjoint().array();  //M * U_down.adjoint
                }

                switch (npulse){
                default :
                case 3:
                    L_a2_even = M.array().rowwise()           * U_down2.transpose().eval().array();  //M * U_down2
                    L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();      //M.adjoint * U_up2.adjoint
                case 2:
                    L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().eval().array();     //M.adjoint * U_up2
                    L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array();     //M * U_down2.adjoint
                    break;
                case 1:
                    break;
                }

                sum_L_a = L_a_even; // M * U_down
                sum_L_b = L_b_even; // M.adjoint * U_up.adjoint


                int loopN=1;
                for (int nn = 1; nn <= npulse; nn ++){
                    if ( nn == npulse){
                        if ( nn % 2==0){  //even number pulse
                            //L_a = L_a * M * U_down;
                            //L_b = M.adjoint() * U_up.adjoint() * L_b;

                            //L_a_even = M.array().rowwise()           * U_down.transpose().eval().array(); //M * U_down
                            //L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array(); //M.adjoint * U_up.adjoint

                            sum_L_a = sum_L_a * L_a_even;
                            sum_L_b = L_b_even * sum_L_b;

                        }else{ //odd number pulse
                            //L_a = L_a * M.adjoint() * U_up;
                            //L_b = M * U_down.adjoint() * L_b;

                            //L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().eval().array(); //M.adjoint() * U_up
                            //L_b_odd = M.array().rowwise()           * U_down.adjoint().array(); //M * U_down.adjoint
                            
                            sum_L_a = sum_L_a * L_a_odd;
                            sum_L_b = L_b_odd * sum_L_b;

                        }
                    }else{
                        if (loopN % 2 != 0){ //odd number term
                            //L_a = L_a * M.adjoint()* U_up2;
                            //L_b = M * U_down2.adjoint() * L_b;

                            //L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().eval().array(); //M.adjoint * U_up2
                            //L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array(); //M * U_down2.adjoint

                            sum_L_a = sum_L_a * L_a2_odd;
                            sum_L_b = L_b2_odd * sum_L_b;

                        }else{ //even number term 
                            //L_a = L_a * M * U_down2;
                            //L_b = M.adjoint() * U_up2.adjoint() * L_b;

                            //L_a2_even = M.array().rowwise()           * U_down2.transpose().eval().array(); //M * U_down2
                            //L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();     //M.adjoint * U_up2.adjoint

                            sum_L_a = sum_L_a * L_a2_even;
                            sum_L_b = L_b2_even * sum_L_b;

                        }
                        loopN++;
                    }
                }
                // only the trace of sum_L_a * sum_L_b is used
                overlap = (sum_L_a.array() * sum_L_b.transpose().array()).sum();
            }

            //////////////////////////////////////////////////////////////
            // The products above come out as Tr[U_d U_u^dagger ...] when the pulse
            // count is even and Tr[U_d^dagger U_u ...] when it is odd - one is the
            // conjugate of the other. A single sample is always the second, and an
            // ensemble has to be the mean of its samples: measured against the exact
            // answer for three spins, Ramsey and CPMG-2 had the imaginary part with
            // the opposite sign (2026-09-19).
            if (npulse==0) overlap = L_mat.trace();
            overlap /= doublec(bdim*1.0,0.0);
            if (npulse%2==0) overlap = std::conj(overlap);
            result[i] = MatrixXcd::Zero(1,1);
            result[i](0,0) = overlap;
            tfree += deltat;
        }
    //////////////////////////////////////////////////////////////
    //  Single Sample approch -> Huijin's code
    //////////////////////////////////////////////////////////////
    }else if (!isEnsemble){

        // One state: matrix-vector products from both ends (evolveSingleState); the start vectors and M^dagger from cceSetup
        const Eigen::RowVectorXcd& rowStart = ws.rowStart;
        const Eigen::VectorXcd&    colStart = ws.colStart;
        const MatrixXcd& Madj = ws.Madj;
        const SingleStateRun run{es_up_Val,es_down_Val,sequence,sequence_indices,npulse,nstep,deltat,result,ws};
        const bool recursive = strcasecmp(cnf->evolution_update,"recursive")==0;
        // A16: bath dimension 2 or 4 (every P1 cluster) with fixed-size types. exact: the bits of the dynamic path, zeros'
        // signs included - the products and the final sum through Eigen's general kernels (rowTimes, colTimes, rowDotCol);
        // recursive: the unrolled column products, last bits may differ from the dynamic path (P1-04 of review
        // 20261012-233058-452278). Any other dimension as before.
        if (bdim==4){
            const Eigen::Matrix4cd M4=M, Madj4=Madj;
            Eigen::RowVector4cd r0=rowStart, r, rn; Eigen::Vector4cd c0=colStart, c, cn, sc;
            if (recursive) evolveSingleState<true>(run,M4,Madj4,r0,c0,r,rn,c,cn,sc);
            else           evolveSingleState<false>(run,M4,Madj4,r0,c0,r,rn,c,cn,sc);
        } else if (bdim==2){
            const Eigen::Matrix2cd M2=M, Madj2=Madj;
            Eigen::RowVector2cd r0=rowStart, r, rn; Eigen::Vector2cd c0=colStart, c, cn, sc;
            if (recursive) evolveSingleState<true>(run,M2,Madj2,r0,c0,r,rn,c,cn,sc);
            else           evolveSingleState<false>(run,M2,Madj2,r0,c0,r,rn,c,cn,sc);
        } else {
            ws.row.resize(bdim); ws.rowNext.resize(bdim); ws.col.resize(bdim); ws.colNext.resize(bdim); ws.scaled.resize(bdim);
            if (recursive) evolveSingleState<true>(run,M,Madj,rowStart,colStart,ws.row,ws.rowNext,ws.col,ws.colNext,ws.scaled);
            else           evolveSingleState<false>(run,M,Madj,rowStart,colStart,ws.row,ws.rowNext,ws.col,ws.colNext,ws.scaled);
        }

    }else{

        // Hoisted out of the nstep loop. Every one of these is fully overwritten before
        // it is read on each step (ipulse==0 always assigns UupTau/UdownTau/UdownDagger
        // /Uup, and each [ipulse] slot is assigned in both branches), so lifting the
        // declarations changes no arithmetic -- it only stops the allocator from
        // rebuilding the same bdim-sized buffers nstep times per cluster. Eigen's
        // operator= reuses the existing buffer when the shape already matches, which is
        // exactly the case here because bdim is fixed for one calCoherenceCce call.
        MatrixXcd UdownDagger(bdim,bdim);
        MatrixXcd Uup(bdim,bdim);
        // Destination for the two products that would otherwise read their own left-hand
        // side. Eigen answers "A = A * B" by computing into a hidden temporary and
        // copying back, which is a fresh allocation per pulse per step; writing into a
        // reused buffer with noalias() and swapping the data pointers runs the same gemm
        // on the same operands in the same order, and swap() only exchanges pointers.
        MatrixXcd Uscratch(bdim,bdim);

        MatrixXcd* Mdagger_UupTauDagger = new MatrixXcd[npulse+1];
        MatrixXcd* M_UdownTauDagger = new MatrixXcd[npulse+1];
        MatrixXcd* Mdagger_UupTau = new MatrixXcd[npulse+1];
        MatrixXcd* M_UdownTau = new MatrixXcd[npulse+1];

        Eigen::VectorXcd UupTau(bdim);
        Eigen::VectorXcd UdownTau(bdim);

         //#pragma omp parallel for
         for (int i = 0; i < nstep; i++) {
            //double stime_a = MPI_Wtime();

            // Assign the U operators for each sequence
            // and calculate total U operators
            for (int ipulse=0; ipulse<npulse+1; ipulse++){

                //double       tau = tfree * tau_fractions[ipulse];
                double tau       = tfree * sequence[ipulse][2];
                int SameTauIndex = sequence_indices[ipulse];

                //printf("[%d] tau = %.10lf\n",ipulse,tau);

                // Calculate operators 
                if (SameTauIndex == ipulse){
                    // U = exp(-iHt)
                    // U = M(Ud)Mdagger
                    // Ud = exp(-i *Hd * t)
                    // ev_up_tau = Hd_up * tau
                    // ev_down_tau = Hd_down * tau
                    // UupTau = Ud_up
                    // UdownTau = Ud_down
                    //
                    // es_up_Val / es_down_Val are the eigenvalues of Hermitian matrices,
                    // so their imaginary parts are exactly zero and es_*_Val*tau is
                    // exactly real. For a real argument the complex cosine and sine
                    // collapse onto the real ones:
                    //     cos(x+0i) = cos(x)*cosh(0) - i*sin(x)*sinh(0)
                    // with cosh(0) == 1.0 and sinh(0) == 0.0 exactly, so what is stored
                    // here is bit for bit what the complex path produced. It just skips
                    // libm's ccosh/csinh, which the profile put at ~19% of the run.
                    for (int k=0; k<bdim; k++){
                        double a_up   = es_up_Val(k).real()   * tau;
                        double a_down = es_down_Val(k).real() * tau;
                        UupTau(k)   = doublec(std::cos(a_up),   -std::sin(a_up));
                        UdownTau(k) = doublec(std::cos(a_down), -std::sin(a_down));
                    }
    
                    //std::cout << "Eigenvalues\n" << std::endl;
                    //std::cout << es_up_tau << std::endl;
                    //std::cout << es_down_tau << std::endl;
                    //std::cout << "Uoperators\n" << std::endl;
                    //std::cout << UupTau << std::endl;
                    //std::cout << UdownTau << std::endl;

                    //Mdagger_UupTauDagger = M_dagger * UupTauDagger
                    //M_UdownTauDagger = M * UdownTauipulse,Dagger
                    //Mdagger_UupTau = M_dagger * UupTau
                    //M_UdownTau = M * UdownTau
                    Mdagger_UupTauDagger[ipulse] = (M.adjoint().array().rowwise() * UupTau.adjoint().array()).matrix();
                    M_UdownTauDagger[ipulse]     = (M.array().rowwise()           * UdownTau.adjoint().array()).matrix();
                    // .eval() materialised the transposed phase vector into a fresh
                    // buffer on every step. transpose() is a pure view and nothing here
                    // aliases it, so dropping the eval() leaves each coefficient product
                    // identical and removes one allocation per line per step.
                    Mdagger_UupTau[ipulse]       = (M.adjoint().array().rowwise() * UupTau.transpose().array()).matrix();
                    M_UdownTau[ipulse]           = (M.array().rowwise()           * UdownTau.transpose().array()).matrix();
                }
                else{
                    Mdagger_UupTauDagger[ipulse] = Mdagger_UupTauDagger[SameTauIndex];
                    M_UdownTauDagger[ipulse]     = M_UdownTauDagger[SameTauIndex]    ;
                    Mdagger_UupTau[ipulse]       = Mdagger_UupTau[SameTauIndex]      ;
                    M_UdownTau[ipulse]           = M_UdownTau[SameTauIndex]          ;
                }

                if (ipulse==0){
                    // M_down * U_down_dagger 
                    UdownDagger = (M_down.array().rowwise() 
                                * UdownTau.adjoint().array()).matrix();
                    // M_down_dagger * M_up * U_up * M_up_dagger 
                    Uup.noalias() = (M.adjoint().array().rowwise() * UupTau.transpose().array()).matrix()
                                * (M_up.adjoint());

                    //std::cout << "M_down * U_down_dagger" << std::endl;
                    //std::cout << UdownDagger << std::endl;
                    //std::cout << "M * U_up * M_up_dagger" << std::endl;
                    //std::cout << Uup << std::endl;
                }
                else if (ipulse!=0 && ipulse%2==1){
                    // odd
                    Uscratch.noalias() = UdownDagger * Mdagger_UupTauDagger[ipulse];
                    UdownDagger.swap(Uscratch);
                    Uscratch.noalias() = M_UdownTau[ipulse] * Uup;
                    Uup.swap(Uscratch);
                    //std::cout << "M_dagger * U_up_dagger" << std::endl;
                    //std::cout << Mdagger_UupTauDagger[ipulse] << std::endl;
                    //std::cout << "M * U_down" << std::endl;
                    //std::cout << M_UdownTau[ipulse] << std::endl;
                }
                else if (ipulse!=0 && ipulse%2==0){
                    //even
                    Uscratch.noalias() = UdownDagger * M_UdownTauDagger[ipulse];
                    UdownDagger.swap(Uscratch);
                    Uscratch.noalias() = Mdagger_UupTau[ipulse] * Uup;
                    Uup.swap(Uscratch);
                }
                else{printf("Error calsingle.cpp\n");exit(1);}
            }

            //std::cout << "UdownDagger" << UdownDagger << std::endl;
            //std::cout << "Uup" << Uup << std::endl;
            //std::cout << psi0 << std::endl;
            // (an ensemble with a Manual sequence: one state goes through the branch above)
            // Tr[UdownDagger*Uup]/d, summed without forming the product: the mean over
            // a complete basis of psi^dagger UdownDagger Uup psi.
            doublec overlap = (UdownDagger.array() * Uup.transpose().array()).sum();
            result[i] = MatrixXcd::Zero(1,1);
            result[i](0,0) = overlap / doublec(bdim*1.0,0.0);

            //printf("%.10lf  %.10lf %+.10lf j \n",tfree,result[i](0,0).real(),result[i](0,0).imag());
            tfree += deltat;                      

            //etime = MPI_Wtime();
            //printf("          Wall time step4-%d# = %.5f s\n",i, stime_a - etime);
        }

        delete[] Mdagger_UupTauDagger;
        delete[] M_UdownTauDagger;
        delete[] Mdagger_UupTau;
        delete[] M_UdownTau;
    }

    //etime = MPI_Wtime();
    //printf("          Wall time step done = %.5f s\n",stime - etime);
}

#ifdef CCEX_USE_CUDA
std::vector<MatrixXcd*> calCoherenceCceCudaBatch(QubitArray* qa,
        const std::vector<BathArray*>& baths, Config* cnf, Pulse* pulse){
    try {
        if(baths.empty() || QubitArray_getNqubit(qa)>1 || strcasecmp(cnf->quantity,"coherence"))
            throw std::runtime_error("invalid cCCE cluster batch");
        const int n=int(baths.size()),d=BathArray_dim(baths[0]),nstep=cnf->nstep,npulse=pulse->npulse;
        const int batch=std::min(cnf->gpu_batch_steps,std::max(1,nstep));
        const bool ensemble=cnf->nstate==0;
        std::vector<std::unique_ptr<GpuCce,void(*)(GpuCce*)>> owned;
        std::vector<GpuCce*> works;
        std::vector<std::unique_ptr<MatrixXcd[]>> results;
        std::vector<double> build(n),setup(n);
        for(int c=0;c<n;++c){
            BathArray* ba=baths[c];
            if(BathArray_dim(ba)!=d)throw std::runtime_error("mixed cCCE batch dimensions");
            const double start=cnf->profile ? MPI_Wtime() : 0;
            auto conditional=cceHamiltonians(qa,ba,cnf);
            const MatrixXcd& ha=conditional.first;const MatrixXcd& hbeta=conditional.second;
            const double ready=cnf->profile ? MPI_Wtime() : 0;
            validateCceCudaHamiltonians(ha,hbeta);
            MatrixXcd psi;
            if(!ensemble)psi=BathArray_Psi0(ba);
            owned.emplace_back(gpuCceCreate(d,ensemble,batch,npulse,ha.data(),hbeta.data(),
                                            ensemble ? nullptr : psi.data()),gpuCceDestroy);
            works.push_back(owned.back().get());
            results.emplace_back(new MatrixXcd[nstep]);
            build[c]=ready-start;setup[c]=(cnf->profile ? MPI_Wtime() : 0)-ready;
        }
        const double evolutionStart=cnf->profile ? MPI_Wtime() : 0;
        std::vector<double> durations((size_t)batch*(npulse+1));
        std::vector<doublec> coherence((size_t)n*batch);
        double time=0;
        for(int first=0;first<nstep;first+=batch){
            const int count=std::min(batch,nstep-first);
            cceCudaDurations(pulse,ensemble,count,double(cnf->deltat),time,durations);
            gpuCceEvolveBatch(works.data(),n,count,durations.data(),coherence.data());
            for(int c=0;c<n;++c)for(int t=0;t<count;++t)
                results[c][first+t]=MatrixXcd::Constant(1,1,coherence[(size_t)c*count+t]);
        }
        const double evolution=(cnf->profile ? MPI_Wtime() : 0)-evolutionStart;
        std::vector<MatrixXcd*> out;
        for(int c=0;c<n;++c){
            computeProfile(cnf,"cce","cuda",d,build[c],setup[c],evolution/n);
            out.push_back(results[c].release());
        }
        computeClusterBatch(cnf,d,n);
        return out;
    }catch(const GpuOutOfMemory& error){
        // backend auto: nothing (what the card had is freed), and the caller computes the clusters one by one, here on the CPU
        if(computeCudaOutOfMemory(cnf,BathArray_dim(baths[0]),error.what()))return {};
        fprintf(stderr,"Error: CUDA cCCE cluster batch: %s\n",error.what());
        MPI_Abort(MPI_COMM_WORLD,1);return {};
    }catch(const std::exception& error){
        fprintf(stderr,"Error: CUDA cCCE cluster batch: %s\n",error.what());
        MPI_Abort(MPI_COMM_WORLD,1);return {};
    }
}
#endif
