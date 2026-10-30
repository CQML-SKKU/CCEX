/**
 * @file gpu_backend.h
 * @brief The interface of the CUDA cCCE and gCCE paths for large clusters (no Eigen or MPI types cross it).
 */
#ifndef CCEX_GPU_BACKEND_H
#define CCEX_GPU_BACKEND_H
#include <cstddef>
#include <stdexcept>
#include "accelerator_defaults.h"

// No Eigen/MPI types cross the CUDA translation-unit boundary. Matrices are
// column-major complex doubles, represented as interleaved real/imag doubles.
// Host code and nvcc may use different SIMD alignment settings.
struct GpuGcce;
// Thrown by every call below when the card has no memory left for it (cudaMalloc, cuBLAS or cuSOLVER allocation
// failure). What the call had allocated is freed and the context stays usable, so with backend auto the caller can
// compute the cluster on the CPU instead (compute_backend.h, computeCudaOutOfMemory). Test hook:
// CCEX_TEST_GPU_OOM_BYTES=<n> makes every device allocation of n bytes or more fail this way.
struct GpuOutOfMemory : std::runtime_error { using std::runtime_error::runtime_error; };
bool gpuInitialize(int localRank, int localSize, char* description, int capacity);
void gpuFinalize();
// A68: "CUDA 12.8, sm_86" for the banner's Build line: the CUDA runtime this binary was compiled against (CUDART_VERSION)
// and the architecture the kernels were built for (CCEX_CUDA_ARCH, passed by both builds)
void gpuBuildDescription(char* out, int capacity);
// the global memory of the smallest visible device, 0 without one (or with CUDA_VISIBLE_DEVICES hiding them all)
std::size_t gpuDeviceMemorySmallest();
// the number of visible devices, 0 without one (A98)
int gpuVisibleDevices();
void gpuSetCacheLimit(std::size_t bytes);
// X12 (review 20261013b-gpu-perf I1): after an out-of-memory, true empties the reuse cache and frees every buffer
// released from then on (the partial work and the batch objects are destroyed after the exception was caught), so that
// none stays held on the card; false (the next cluster) caches again. The cache limit is not changed.
void gpuOutOfMemoryRecovery(bool on);
struct GpuMemoryStats {
    std::size_t allocations, reuses, cachedBytes, peakCachedBytes, pointerUploads;
};
GpuMemoryStats gpuMemoryStats();
GpuGcce* gpuGcceCreate(int dimension, int qubitDimension, bool vectorMode,
                      int batchSteps, const void* hamiltonian, const void* initial,
                      int pulseCount, const void* pulses, bool refine=true,
                      int maxSweeps=ccex::refinementSweeps);
// times uses the exact host time grid, including its repeated-addition rounding.
// reduced: count consecutive qdim x qdim matrices; full may be null.
void gpuGcceEvolve(GpuGcce* work, int count, const double* times,
                   const double* delayFractions, void* reduced, void* full);
void gpuGcceEvolveBatch(GpuGcce* const* works, int clusterCount, int count,
    const double* times, const double* delayFractions, void* reduced);
void gpuGcceDestroy(GpuGcce* work);

struct GpuCce;
// ensemble=true uses the normalized bath trace; otherwise initial is a pure
// bath state. Hamiltonians, initial and coherence are complex-double arrays.
GpuCce* gpuCceCreate(int dimension, bool ensemble, int batchSteps, int pulseCount,
                    const void* hAlpha, const void* hBeta, const void* initial);
// durations[p*count+t] is the exact host duration for interval p and time t.
// The ensemble path preserves legacy cCCE's even-pulse conjugation convention.
void gpuCceEvolve(GpuCce* work, int count, const double* durations, void* coherence);
// Works must share dimension, state mode, pulse count and time schedule.
// Output is cluster-major: clusterCount blocks of count complex scalars.
void gpuCceEvolveBatch(GpuCce* const* works, int clusterCount, int count,
                       const double* durations, void* coherence);
void gpuCceDestroy(GpuCce* work);
#endif
