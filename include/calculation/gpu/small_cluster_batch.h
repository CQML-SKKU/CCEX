#ifndef CCEX_SMALL_CLUSTER_BATCH_H
#define CCEX_SMALL_CLUSTER_BATCH_H
// A26 step 1: small cCCE clusters (bath dimension 2 and 4) computed on the GPU in large batches, one thread per cluster
// (small_cluster_gpu.cu): the packed kernel, "gpu_options": {"packed_max_dim": 2 or 4} (0: off, and nothing of it runs; the
// default is 4 in the GPU build with backend auto or cuda, method CCE and quantity coherence, else 0 -
// Config_resolvePackedMaxDim). The rest of the run is unchanged.
//
// calculate() asks get(ic) for each local cluster. The first time a window of clusters is reached, the rank builds their
// Hamiltonians (cceClusterHamiltonians, the CPU path's own function) and packs them into its part of a node shared-memory
// window; before waiting for that window it packs the next one, so that the card computes one window while the rank builds
// the next and assembles the one before. One process per node drives the card: node rank 0 runs a feeder thread that
// collects the packed windows of all ranks of the node, computes them in combined batches and writes the coherences back
// into the shared memory. Every rank keeps its clusters and does the rest of its work (Hamiltonians, the product, the MPI
// reduction) as before; only the eigensystems and the time evolution move to the card. One CUDA context per node instead of
// one per rank: 32 ranks starting CUDA at once took 3.5 s and their teardown about 2 s more on idm4-gpu (A26 step 1, job 683).
// A cluster the batch does not take (another dimension, a Hamiltonian that is not finite and Hermitian, a Jacobi that did
// not converge, and the clusters met while the card is still starting) gets nullptr, and calculate()
// computes it on the CPU as before. A card that cannot be opened turns the batch off under backend auto (every cluster on
// the CPU, as the auto contract says) and stops the run under cuda or cpu, also when the ranks had finished before.
//
// A40 (gpu_options.warp_max_dim, on by default in the GPU build with backend auto or cuda): the same windows and feeder also carry the clusters of the warp
// kernel (warp_cluster.h) - gCCE clusters of total dimension 6 to 27 with a qubit of dimension 3, cCCE clusters of bath
// dimension 8 and 16 - of a total dimension up to warp_max_dim. Their inputs go into an area of the window of their own
// (a window closes early when it is full), their outputs after the packed kernel's.
#include "calculation/simulator.h"
#include "calculation/gpu/small_cluster_math.h"
#include "calculation/gpu/warp_cluster.h"
#include <complex>

class SmallClusterBatch {
public:
    // Collective (every rank calls it at the same place). Active with gpu_options.packed_max_dim 2 or 4, method CCE with one
    // qubit, and a pulse plan of at most 4 distinct intervals - or with gpu_options.warp_max_dim and a run the warp kernel
    // takes (warpClusterRun). Asked for in a build without CUDA or with another method the run stops in smallClusterWarmUp.
    // localclusters, order: this rank's clusters (localclusters[n]) - the window is no larger than the largest number of
    // one order, and the batch is off (no memory, no card, no feeder) when no rank has a cluster of a dimension it takes,
    // or - with neither packed_max_dim nor warp_max_dim given - when the cluster evaluations (clusters x states) it would
    // take over all ranks are fewer than 20000 per rank (starting the card costs more than it saves; A40 follow-up,
    // 2026-10-13). A48: with packed_max_dim not given and a run the CPU batch takes (ccebatchRunTakes), the packed kernel
    // only for one state with at least 40000 clusters of bath dimension 2 and 4 per rank per state (off for the ensemble);
    // the CPU batch computes the rest. The warp input region
    // is there only when some rank has a cluster of a warp dimension.
    // savemode: the run's (the warp kernel is off for gCCE cluster_density_matrix).
    SmallClusterBatch(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE, int*** localclusters,
                      int order, const char* savemode = nullptr);
    ~SmallClusterBatch();
    SmallClusterBatch(const SmallClusterBatch&) = delete;
    SmallClusterBatch& operator=(const SmallClusterBatch&) = delete;
    bool on() const { return impl != nullptr; }
    // Collective, once per state before its clusters (after psia, psib and the qubit's psi0 are set): the gCCE warp plan
    // (pulses and qubit states) of the state; the run stops if it differs from the state before (it cannot: those are
    // set once per run). Nothing otherwise.
    void beginState(QubitArray* qa, Pulse* pls);
    // whether a cluster of this total dimension is the batch's (packed: 2 and 4 up to packed_max_dim; warp: the compiled ones up to warp_max_dim):
    // the others keep the ordinary CUDA cluster batching (review 20261012-233058-452278 P1-05)
    bool packs(int dim) const;
    // A48: where a window of the CPU batch (CceBatchWindow) that starts at cluster ic - one get() just left to the CPU -
    // has to stop, so that it takes none of the clusters this batch computes: the end of the window the CPU computes while
    // the card starts, ic+1 for one the card gave back or does not take, no bound (INT_MAX) when the packed kernel is off
    int cpuEnd(int ic) const;
    // the local clusters of order n (localclusters[n]: [0][0]-1 clusters, each [power, spins...]) of the current state
    void beginOrder(int n, int** clusters, int ncluster);
    // cluster ic's values when the GPU computed it, else nullptr: nstep points of per values each (1, the coherence; the
    // qdim x qdim reduced density matrix, column-major, for gCCE with quantity dm), and its bath dimension
    const std::complex<double>* get(int ic, int& dim, int& per);
    // collective, after the last state: stops the feeder, releases the shared memory and prints one [gpu-small] line by
    // rank 0 (counts and seconds over the ranks); the batch is off afterwards
    void report();
private:
    struct Impl;
    Impl* impl = nullptr;
};

// Collective, early in the run (computeBackendInit, before the bath is read and the clusters are made): when the packed
// or the warp kernel is asked for, every rank allocates its part of the node shared memory and node rank 0 starts the CUDA
// context (and page-locks that memory) in a thread of its own, so that the feeder finds it ready (about 0.2 s on
// idm4-gpu, otherwise in the first state). smallClusterWarmDown (computeBackendFinish) joins it if no batch took it;
// SmallClusterBatch calls it at once when no rank has a cluster of a dimension the batch takes (GPU 다음 셋 (1)).
// With qa and pls (review 20261013-221715 X7) nothing starts for a kernel the input already rules out (the constructor's
// checks that need no cluster: warpClusterRun, savemode, the packed plan, A48's ensemble rule).
void smallClusterWarmUp(const Config* cnf, QubitArray* qa = nullptr, Pulse* pls = nullptr, const char* savemode = nullptr);
// Collective, once the bath and the defects are read (A55): for a one-state standard cCCE run whose packed kernel the
// constructor decides by its threshold, smallClusterWarmUp starts nothing; this counts the clusters of bath dimension 2
// and 4 from the bath (the spins and the pairs closer than rdip) and starts the card only when they reach that threshold.
// Nothing otherwise.
void smallClusterWarmUpAfterBath(const Config* cnf, BathArray* ba, DefectArray* dfa);
void smallClusterWarmDown();

// The time plan of the run as the batch uses it (interval fractions, distinct slots); false if it has more than
// 4 distinct intervals or more than maxIntervals. Both builds.
bool smallClusterPlan(Config* cnf, Pulse* pls, smallcluster::Plan& plan, bool& ensemble, bool& recursive);
// Test data for the GPU check (test_small_cluster_cuda): with CCEX_SMALL_CLUSTER_DUMP=<file>, rank 0 appends its first
// CCEX_SMALL_CLUSTER_DUMP_N (default 400) cCCE or pCCE clusters of dimension 2 and 4 - plan, Halpha, Hbeta, psi and the CPU
// coherence just computed. Any build; nothing happens without the variable.
void smallClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result);

// A40: whether the warp kernel can take this run's clusters (any build): method CCE (quantity coherence) or gCCE (quantity
// coherence or dm, one qubit, propagator eigen, no gcce_refinement, savemode not cluster_density_matrix, pulses without pulseiter or
// detuning, evolution matrix for an ensemble), at most 4 distinct pulse intervals. Fills the time part of plan (and mode,
// ensemble, dm, qdim); false with why (a phrase) otherwise.
bool warpClusterRun(QubitArray* qa, Config* cnf, Pulse* pls, warpcluster::Plan& plan, const char*& why);
// The gCCE part of the plan: the pulses (cachedUpulses), psia, psib and the qubit psi0; psia and psib must be set.
void warpClusterQubit(QubitArray* qa, Pulse* pls, warpcluster::Plan& plan);
// Test data for the GPU check (test_warp_cluster_cuda): with CCEX_WARP_CLUSTER_DUMP=<file>, rank 0 appends its first
// CCEX_WARP_CLUSTER_DUMP_N (default 400) clusters of each dimension the warp kernel has (any build, warp_max_dim not needed) -
// the plan, the input and the CPU result just computed. Nothing happens without the variable.
void warpClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result);
#endif
