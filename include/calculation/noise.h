/**
 * @file noise.h
 * @brief method noise: the coherence exp(-chi) of Gaussian phase noise, from a model spectrum or from the bath clusters.
 */
#ifndef __CCEX_NOISE_H_
#define __CCEX_NOISE_H_
#include "../../include/calculation/simulator.h"

/** chi(T) of the model spectrum of cnf (noise_model "OU" or "power") at the total time T (ms) */
double noiseModelChi(Config* cnf, Pulse* pls, double T);
/** The model spectrum's coherence exp(-chi) at the run's time points, written as the run's coherence file (rank 0) */
void noiseModelRun(Config* cnf, Pulse* pls, Output* op);
/** One cluster: -chi_C at every time point (1 x 1 each), chi_C the mode sum of the cluster's eigenspectrum; the cCCE
 *  product takes exp(-chi_C) (noiseFactor, simulator.cpp) and is then exp(-sum of the CCE terms of chi) */
void calCoherenceNoiseInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result);

#endif
