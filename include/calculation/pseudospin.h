/**
 * @file pseudospin.h
 * @brief method pseudospin: the cCCE with its pair clusters in the pseudospin model of nuclear flip-flops.
 */
#ifndef __CCEX_PSEUDOSPIN_H_
#define __CCEX_PSEUDOSPIN_H_
#include "../../include/calculation/simulator.h"

/** A pair cluster: P L_1 L_2 at every time point (1 x 1 each), P the pair's pseudospin term and L_k the full single-spin
 *  coherences, so that the cCCE assembly gives prod_k L_k prod_pairs P. Writes the pair's channels to the channel file. */
void calCoherencePseudospinInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result);
/** Every rank: open its part of "{outfile}_pseudospin_channels" (before calculate) */
void pseudospinChannelsBegin(Output* op);
/** Every rank: close its part; rank 0 then joins the parts into the one file (after calculate) */
void pseudospinChannelsEnd();

#endif
