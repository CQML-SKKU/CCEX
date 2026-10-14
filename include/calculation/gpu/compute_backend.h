#ifndef CCEX_COMPUTE_BACKEND_H
#define CCEX_COMPUTE_BACKEND_H

#include "../../input/general.h"
#include "../../physical/pulse.h"
#include "../../physical/qubit.h"

// Called after all JSON/CLI overrides, before reading bath/tensor files. qa, pls, savemode (optional): for the early start
// of the small-cluster batch (smallClusterWarmUp, review 20261013-221715 X7).
void computeBackendInit(const Config* config, QubitArray* qa = nullptr, Pulse* pls = nullptr, const char* savemode = nullptr);
bool computeUseCuda(const Config* config, int dimension);
// The methods whose clusters may go to the GPU: cCCE (cluster_mode standard or pcce) and gCCE with propagator=eigen.
// computeBackendInit and the evolution_update default (Config_validateEvolutionUpdate, A36) ask the same.
bool computeBackendMethodSupported(const Config* config);
// A cluster of this dimension ran out of GPU memory (GpuOutOfMemory, gpu_backend.h; what = its message). backend auto
// (user 2026-10-13, A37 job 699): true - the caller computes the cluster on the CPU, and from now on this rank sends no
// cluster of this dimension or larger to the GPU (computeUseCuda says no); one line on stdout. Otherwise (cuda) false: the
// caller stops the run as for any other GPU error.
bool computeCudaOutOfMemory(const Config* config, int dimension, const char* what);
void computeBackendFinish();
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

#endif
