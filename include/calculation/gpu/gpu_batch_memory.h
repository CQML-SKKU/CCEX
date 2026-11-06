/**
 * @file gpu_batch_memory.h
 * @brief The memory estimate that decides how many clusters and time points go to the GPU in one batch.
 */
#ifndef CCEX_GPU_BATCH_MEMORY_H
#define CCEX_GPU_BATCH_MEMORY_H
#include <stdexcept>
// Conservative admission estimate, not a cudaMemGetInfo guarantee. The D^2
// setup allowance remains for eigenvectors/solver temporaries in vector mode.
inline long double gpuClusterMemoryEstimate(int dimension,int qubitDimension,
        int timeBatch,int pulses,bool gcce,bool vectorMode) {
    if(dimension<1 || qubitDimension<1 || timeBatch<1 || pulses<0)
        throw std::invalid_argument("invalid GPU batch estimate dimensions");
    const long double d=dimension,q=qubitDimension,t=timeBatch,p=pulses;
    const long double width=vectorMode ? 1 : d;
    const long double complexEntries=3*t*d*width+8*d*d+4*d+t*(gcce ? q*q : 1)+p*q*q;
    const long double realEntries=2*d+t*(gcce ? 1 : p+1);
    return 16*complexEntries+8*realEntries+sizeof(void*)*(8+4*t);
}
// The workspace of the batched eigensolver for one cCCE cluster (its two Hamiltonians): cusolverDnXsyevBatched asked for
// 7 complex numbers per matrix element from D = 256 to 1296 (56 MiB for 8 at 256, 3377 MiB for 32 at 1296).
inline long double gpuCceBatchWorkspaceEstimate(int dimension) {
    const long double d=dimension;
    return 16*14*d*d;
}
#endif
