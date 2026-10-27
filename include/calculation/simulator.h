/**
 * @file simulator.h
 * @brief The calculation: cluster Hamiltonians, time evolution and the coherence of each cluster and of the whole bath.
 */
#include "calculation/gpu/accelerator_defaults.h"
#ifndef __CCEX_SIMULATOR_H_
#define __CCEX_SIMULATOR_H_

#include <unsupported/Eigen/MatrixFunctions>
#include <vector>   // GcceWork/CceWork members; g++ gets it transitively, icpc does not
#include "common/utilities.h"
#include "physical/qubit.h"
#include "physical/bath.h"
#include "physical/defect.h"
#include "common/output.h"
#include "physical/pulse.h"
#include "clustering/cluster.h"
#include "input/general.h"
#include "common/blas_option.h"   // A46: CcexEigensolver (Eigen's solver unless built with CCEX_BLAS_DIRECT)

//MatrixXcd cal_Upulse_withDetuning(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma, MatrixXcd sigma_z);
void calUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse);
void calTDUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse, MatrixXcd Hq);
void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cc, Pulse* pulse, Cluster* cluster, Output* output, int*** localClusters);

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit);

struct GcceBatch;
MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op, GcceBatch* batch=nullptr);
MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);
MatrixXcd* calCoherenceMecce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);   // bath jumps only; no 0th cluster (L_0 = 1, as CCE)
MatrixXcd* calCoherenceCceDriven(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, const int* idx);   // pulse.bath_drive (simulator_cce_drive.cpp); idx: the cluster spins in ba of calculate()
void BathDrive_check(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);           // once, on the whole bath: refusals, warnings, log
void BathDrive_prepareMeanField(BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls);   // s_range: the samples and their fields
MatrixXcd* calCoherenceMegcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op);   // ba == NULL: 0th cluster
/**
 * Storage calCoherenceCceInto reuses from one cluster to the next: the eigensolvers, the products and the
 * vectors of the time evolution. Owned by the caller (calculate); nothing in it carries a value from one call to the next - every
 * member is written before it is read.
 */
struct CceWorkspace {
    MatrixXcd Halpha, Hbeta;                         // A25: the conditional Hamiltonians
    CcexEigensolver up, down;
    MatrixXcd M, Madj;
    Eigen::VectorXcd valUp, valDown;
    Eigen::RowVectorXcd rowStart, row, rowNext;
    Eigen::VectorXcd colStart, col, colNext, scaled;
    std::vector<doublec> phases;                     // single state: the phases of each distinct interval (A16)
    // A62 (cceEvolveColumns): one column block's phases per distinct interval, the row side (columns x D) and the column
    // side (D x columns) of the chain and their next values
    std::vector<MatrixXcd> blockUp, blockDown;
    MatrixXcd blockRow, blockRowNext, blockCol, blockColNext;
};
/** calCoherenceGcce writing into into[0..nstep) when given (else a new array), reusing *alphaRho for psia^dagger rho. */
MatrixXcd* calCoherenceGcceInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op, MatrixXcd* into,
                                Eigen::RowVectorXcd* alphaRho, GcceBatch* batch=nullptr);
/** calCoherenceCce writing into result[0..nstep) (1 x 1 each, reused when already so) with the storage of ws; the same arithmetic. */
void calCoherenceCceInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result, CceWorkspace& ws);
/** The two halves of calCoherenceCceInto after the Hamiltonians: cceSetup takes ws.Halpha and ws.Hbeta to the
 *  eigensystems, M, the eigenvalues and, for one state (psi, d x 1, else nullptr), rowStart, colStart and Madj; cceEvolve
 *  is the time evolution from these into result[0..nstep). The batched small-cluster path (cce_batch.h) shares cceSetup. */
void cceSetup(CceWorkspace& ws, const MatrixXcd* psi);
void cceEvolve(Config* cnf, Pulse* pls, MatrixXcd* result, CceWorkspace& ws);
/** (simulator_cce_columns.cpp) The one-state time evolution of cceEvolve with
 *  the time points of a block of `columns` steps as matrix columns, so that each interval is one matrix-matrix product
 *  (D x D times D x columns) instead of `columns` matrix-vector products. The phases are those of evolveSingleState
 *  (exact: sin and cos at every step; recursive: so at steps 0, 16, 32, ... and the previous step's times one step's
 *  factor in between, carried across blocks). The same products in another summation order: the last bits differ.
 *  cceEvolve calls it for bath dimension >= column_evolution.min_dim when column_evolution is on (the column_evolution
 *  input block; on by default with min_dim 256, columns 64). */
void cceEvolveColumns(const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown, double** sequence,
                      const int* sequenceIndices, int npulse, int nstep, double deltat, bool recursive, int columns,
                      MatrixXcd* result, CceWorkspace& ws);
/** (cce_setup_small.cpp) cceSetup with the eigensystems of ws.Halpha (s 0) and ws.Hbeta (s 1) given - the small-cluster
 *  Jacobi's values w[s][0..d) and vectors V[s] (column major), set into ws.up/ws.down; a null w[s] has that one solved by
 *  Eigen. The batch window (cce_batch_window.h) solves its clusters' matrices in groups (small_eigensolver.h) and calls this. */
void cceSetupSmall(CceWorkspace& ws, const MatrixXcd* psi, const double* const* w, const doublec* const* V);
/** For bath dimension 2 and 4 (and CCEX_SMALL_EIG not 0) solves the two Hamiltonians one at a time with the Jacobi
 *  and calls cceSetupSmall - the bits a cluster gets in a group of the window - and returns true; otherwise false (the
 *  caller then runs cceSetup, Eigen). */
bool cceSetupSmallAlone(CceWorkspace& ws, const MatrixXcd* psi);
/** The ensemble (nstate 0) time evolution of calCoherenceCceInto with evolution_update = recursive
 *  (simulator_cce_recursive.cpp): M = M_up^dagger M_down, the eigenvalues of the two conditional Hamiltonians. */
void evolveEnsembleRecursiveCce(const MatrixXcd& M, const MatrixXcd& M_up, const MatrixXcd& M_down,
                                const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown, Pulse* pls,
                                bool manualIntervals, int nstep, double deltat, MatrixXcd* result);
/** X8: the ensemble Ramsey (npulse 0, not Manual) of calCoherenceCceInto with evolution_update = exact, by
 *  Tr[M Dd M^dagger Du^*] = Du^dagger |M|^2 Dd (simulator_cce_recursive.cpp; the recursive branch uses the same trace). */
void evolveEnsembleRamseyCce(const MatrixXcd& M, const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown,
                             int nstep, double deltat, MatrixXcd* result);
#ifdef CCEX_USE_CUDA
std::vector<MatrixXcd*> calCoherenceCceCudaBatch(QubitArray* qa,
    const std::vector<BathArray*>& baths, Config* cnf, Pulse* pulse);
std::vector<MatrixXcd*> calCoherenceGcceCudaBatch(QubitArray* qa,
    const std::vector<BathArray*>& baths, Config* cnf, Pulse* pulse, Output* output);
#endif

/**
 * Scratch for the gCCE propagator, owned by the caller and reused across the nstep
 * loop. Htot does not depend on the time step -- calCoherenceGcce builds it once and
 * the loop only varies tau -- so one eigendecomposition covers every tau:
 *     exp(-i*Htot*tau) = M * diag(exp(-i*lambda_k*tau)) * M^dag
 * That replaces a matrix exponential per step with a handful of scalar exponentials.
 * This is what simulator_cce.cpp has always done; gCCE was calling Eigen's general
 * .exp(), which runs a complex Schur decomposition and does not know Htot is Hermitian.
 */
struct GcceWork {
    bool useEigen;                       // false => legacy per-step matrix exponential
    MatrixXcd M, Madj;                   // eigenvectors of Htot, and their adjoint
    Eigen::VectorXd  evals;              // eigenvalues (real -- Htot is Hermitian)
    Eigen::VectorXd  evalCorrections;    // small Rayleigh corrections kept separately
    Eigen::VectorXd  evals_tau;          // scratch
    Eigen::VectorXcd phase;              // scratch
    MatrixXcd scratch;                   // scratch
    std::vector<MatrixXcd> Ufrees;       // per-pulse free-evolution propagators

    // evolution=vector only (see Config_setEvolution): |psi(t)>, the eigenbasis
    // intermediate of one Ufree application, the phase factor of each distinct tau, and
    // the pulse propagators already expanded over the bath.
    Eigen::VectorXcd psi, psiscratch;
    std::vector<Eigen::VectorXcd> phases, phaseSteps;   // phaseSteps: exp(-i lambda f deltat) per interval (A16 b)
    std::vector<MatrixXcd> UpulsesExpanded;
    // A16, evolution=vector: M^dag |psi0> (the same every step: formed by firstFreeGcce on its first call for the
    // cluster, emptied by prepareGccePulses), and the bath dimension of the pulses, for kron(U, 1) block by block
    Eigen::VectorXcd startEigen;
    int pulseBdim = 0;
};
// A16 (simulator_gcce_vector.cpp, a file of their own: added to simulator_gcce.cpp they changed how g++ compiled the
// rest of that file and the last bit of some gCCE results). w->psi = exp(-i H tau0) psi0 with phase = exp(-i lambda tau0),
// keeping M^dag psi0 for the cluster; and out = UpulsesExpanded[ipulse] * in.
void firstFreeGcce(GcceWork* w, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& phase);
void applyPulseGcce(GcceWork* w, int ipulse, const Eigen::VectorXcd& in, Eigen::VectorXcd& out);

/** Builds everything in GcceWork that is fixed for the cluster. Call once per cluster. */
void prepareGcceWork(GcceWork* w, const MatrixXcd& Htot, bool useEigen, bool refine=true,
                     int maxSweeps=ccex::refinementSweeps);
MatrixXcd calPropagatorGcce(QubitArray* qa, const MatrixXcd& Htot, Pulse* pls, double tfree, MatrixXcd* Upulses, GcceWork* w);
/** calPropagatorGcce (propagator=eigen) for evolution_update = recursive: the calls come for istep = 0, 1, 2, ... with
 *  tfree = istep deltat (accumulated), and the phases between every 16th step are carried over by one step's factor. */
MatrixXcd calPropagatorGcceRecursive(Pulse* pls, double tfree, GcceWork* w, int istep, double deltat);

/** Caches kron(Upulses[i], I_bdim), the pulses expanded to the qubit+bath space, in w->UpulsesExpanded once per
 *  cluster (it does not depend on the step), and resets the kept M^dag psi0 for the new cluster. */
void prepareGccePulses(GcceWork* w, MatrixXcd* Upulses, int npulse, int bdim);
/** evolution=vector: applies the same pulse sequence as calPropagatorGcce, but to the
 *  state vector psi0 instead of assembling the propagator. Result in w->psi.
 *  prepareGccePulses must have been called for this cluster first, and psi0 is the same in every call for the cluster
 *  (M^dag psi0 is formed on the first call and kept). istep >= 0 (evolution_update = recursive): the calls come for istep = 0, 1, 2, ... with
 *  tfree = istep deltat (accumulated), and the phases between every 16th step are carried over by one step's factor. */
const Eigen::VectorXcd& propagateStateGcce(const Eigen::VectorXcd& psi0, Pulse* pls, double tfree, GcceWork* w,
                                           int istep=-1, double deltat=0.0);

MatrixXcd HamilQubit(QubitArray* qa, BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilBath(BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilQubitBath(QubitArray* qa, BathArray* ba, MatrixXcd** qsigmas, MatrixXcd** bsigmas, Config* cnf);
MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf);
double calQubit_Energy_difference(QubitArray* qa, MatrixXcd Hq);
MatrixXcd cal_TDUpulse(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, MatrixXcd ox, MatrixXcd oy, char axis, double pulse_angle, double Ediff, double detuning_factor, double pulse_time, int Nstep);


MatrixXcd** gatherSigmas(MatrixXcd** qsigmas, MatrixXcd** bsigmas, int nqubit, int nspin);
MatrixXcd** gatherSigmas_singlequbit(MatrixXcd* qsigma, MatrixXcd** bsigmas, int nqubit, int nspin);
int close_state_index(MatrixXcd state, MatrixXcd eigenVectors);

#endif // __CCEX_SIMULATOR_H_
