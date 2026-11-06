/**
 * @file compute_backend.h
 * @brief The choice of compute backend (cpu, cuda, auto) per cluster and the GPU start-up.
 */
#ifndef CCEX_COMPUTE_BACKEND_H
#define CCEX_COMPUTE_BACKEND_H

#include "../../input/general.h"
#include "../../physical/pulse.h"
#include "../../physical/qubit.h"
#include "../../physical/bath.h"
#include "../../physical/defect.h"
#include <set>
#include <vector>

// Called after all JSON/CLI overrides, before reading bath/tensor files. qa, pls, savemode (optional): for the early start
// of the small-cluster batch (smallClusterWarmUp, review 20261013-221715 X7).
void computeBackendInit(const Config* config, QubitArray* qa = nullptr, Pulse* pls = nullptr, const char* savemode = nullptr);
bool computeUseCuda(const Config* config, int dimension);
// computeUseCuda's answer without its side effects (it starts the card at the first cluster it sends there): whether a
// cluster of this dimension would be sent to the GPU, for the memory check before the calculation (A67)
bool computeMayUseCuda(const Config* config, int dimension);
// The methods whose clusters may go to the GPU: cCCE (cluster_mode standard or pcce) and gCCE with propagator=eigen.
// computeBackendInit and the evolution_update default (Config_validateEvolutionUpdate, A36) ask the same.
bool computeBackendMethodSupported(const Config* config);
// the GPU cluster batches take consecutive clusters of one dimension: this rank's clusters are then ordered by size, when
// cluster_batch_mb is above its default 64. Not on the CPU ranks of gpu_options.ranks_per_gpu, which do not batch, nor with
// the default, whose batches hold one large cluster: the ranks of a card or of the CPU then all run their large clusters at
// once, which cost more (GH200, ranks_per_gpu 4: CPU ranks ordered 386 s against 328 s for the 1 atm case; the defaults
// ordered 198 s against 165 s for O-rich; batches of 32 in 8192 MiB unordered 223 s against 137 s)
bool computeBatchByDimension(const Config* config);
// A cluster of this dimension ran out of GPU memory (GpuOutOfMemory, gpu_backend.h; what = its message). backend auto
// (user 2026-10-13, A37 job 699): true - the caller computes the cluster on the CPU, and from now on this rank sends no
// cluster of this dimension or larger to the GPU (computeUseCuda says no); one line on stdout. Otherwise (cuda) false: the
// caller stops the run as for any other GPU error.
bool computeCudaOutOfMemory(const Config* config, int dimension, const char* what);
// keepCard (A64): another configuration of configuration.k_range follows - the CUDA context stays open for it (the next
// computeBackendInit takes it as it is); the small-cluster batch and the profile are finished either way
void computeBackendFinish(bool keepCard = false);
// the memory of the node's smallest GPU when this run may put clusters on the GPU (backend cuda or auto, a CUDA build
// and a supported method), else 0; for the probe memory budget (A1). One counter serves the node, and the ranks of a
// node may all land on one GPU, so the smallest device bounds what may run at once (review 5f7834 R3).
std::size_t computeDeviceMemory(const Config* config);
// The measured per-cluster times, written out and forgotten. A separate process
// per calculation flushed these on exit; the in-process refinement loop has to
// flush them at the end of each step, or the cost table reads an empty log.
void computeProfileFlush();
void computeProfile(const Config* config, const char* method, const char* backend,
                    int dimension, double build, double setup, double evolution);
void computeClusterBatch(const Config* config, int dimension, int count);

// A98: execution.gpu_options.ranks_per_gpu = k > 0 (backend auto, cCCE). computeGpuSplitSetup, collective, before the
// clusters are made: on every node the first k x (visible GPUs) ranks are the GPU ranks, and the others never open the
// GPU (computeUseCuda says no), so a node without MPS runs every core and only k contexts per GPU. rankClusterRows then
// gives a cluster of gpu_options.min_dim and above to a GPU rank, weighted by the time its dimension took a GPU rank over
// the time on the CPU, both as loaded as in the run (a made-up cluster of that dimension, computeGpuSplitCalibrate), or
// to a CPU rank when that one would finish it first; a dimension as slow or slower on a GPU rank stays on the CPU of
// every rank.
void computeGpuSplitSetup(Config* config, QubitArray* qa, BathArray* ba, DefectArray* dfa, Pulse* pls);
// Every rank, after computeBackendInit and before computeGpuSplitSetup: the gpu_options the input leaves out, for cCCE on
// the GPU (backend auto or cuda, a GPU visible): ranks_per_gpu 4 when a node has more than 4 ranks a GPU (backend auto,
// cluster.mode standard), min_dim 256 with ranks_per_gpu (its timing sends a dimension to the CPU where the GPU is
// slower), evolution_batch_size 32, cluster_batch_mb two thirds of the card over the ranks that open it. Measured on a GH200
// node (288 ranks, 4 cards) and an RTX PRO 6000 (32 ranks, 1 card): 4 ranks a card with batches of 32 the fastest on both
// (2026-11-02, jobs 1200264 and idm 18162). Nothing in a CPU-only build.
void computeGpuAutoSettings(Config* config);
bool computeGpuSplitOn();
const std::vector<int>& computeGpuSplitRanks();   // the GPU ranks, ascending
// the dimension of a cluster of n spins (row[1..n]) when it goes to the GPU and is not measured yet, else 0
int computeGpuSplitPending(int n, const int* row);
// collective, the same dimensions on every rank: the pending ones, timed on the GPU ranks
void computeGpuSplitCalibrate(const std::set<int>& dimensions);
// GPU-rank time / CPU time of the cluster's dimension when a GPU rank may take it (that time < 1), else 0
double computeGpuSplitWeight(int n, const int* row);

#endif
