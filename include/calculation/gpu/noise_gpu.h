/**
 * @file noise_gpu.h
 * @brief The card side of method noise: sums over the bath's modes at every time point (no Eigen or MPI types cross it).
 */
#ifndef CCEX_NOISE_GPU_H
#define CCEX_NOISE_GPU_H
#include <cstddef>

struct GpuNoise;
// the time points and the pulse sequence (the interval lengths as fractions of T, interval j of sign (-1)^j); on the
// card of this process (computeBackendInit opened it). Throws std::runtime_error, as the two calls below.
GpuNoise* gpuNoiseCreate(int ntime, const double* times, int nfraction, const double* fractions);
void gpuNoiseDestroy(GpuNoise* work);
// sum[t] += sum_k weight_k |Y(omega_k, times[t])|^2, Y(w, T) = int_0^T f(t) exp(i w t) dt
void gpuNoisePhase(GpuNoise* work, std::size_t n, const double* omega, const double* weight, double* sum);
// sum[t] += sum_k (re_k cos(omega_k times[t]) - im_k sin(omega_k times[t]))
void gpuNoiseCosine(GpuNoise* work, std::size_t n, const double* omega, const double* re, const double* im, double* sum);
#endif
