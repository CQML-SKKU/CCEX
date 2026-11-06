/**
 * @file noise.h
 * @brief method noise: the coherence exp(-chi) of Gaussian phase noise, from a model spectrum or from the bath clusters.
 */
#ifndef __CCEX_NOISE_H_
#define __CCEX_NOISE_H_
#include "../../include/calculation/simulator.h"
#include <vector>

/** chi(T) of the model spectrum of cnf (noise_model "OU" or "power") at the total time T (ms) */
double noiseModelChi(Config* cnf, Pulse* pls, double T);
/** The model spectrum's coherence exp(-chi) at the run's time points, written as the run's coherence file (rank 0) */
void noiseModelRun(Config* cnf, Pulse* pls, Output* op);
/** The number of values a cluster gives: the run's time points, or those of hamiltonian.noise.fft (compute spectrum) */
int noiseLength(Config* cnf);
/** One cluster's values (hamiltonian.noise.compute: -chi_C at the run's times, or its correlation C_C(t)) times its power,
 *  added to with and, for a positive power, to without. On the GPU (execution.backend cuda) its modes wait in a batch
 *  until noiseFinish, or until the batch is full. */
void noiseCluster(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, int power, std::vector<double>& with,
                  std::vector<double>& without);
/** The sums of every rank's clusters (values times powers; without: the positive powers alone) to the outputs: the
 *  coherence into wD and nD on rank 0, the correlation and spectrum files. Returns false where no coherence is written
 *  (compute correlation). Collective. */
bool noiseFinish(QubitArray* qa, Config* cnf, Pulse* pls, Output* op, int istate,
                 std::vector<double>& with, std::vector<double>& without, MatrixXcd* wD, MatrixXcd* nD);

#endif
