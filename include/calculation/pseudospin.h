/**
 * @file pseudospin.h
 * @brief method pseudospin: the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops, alone.
 */
#ifndef __CCEX_PSEUDOSPIN_H_
#define __CCEX_PSEUDOSPIN_H_
#include "../../include/calculation/simulator.h"

/** A pair cluster: P at every time point (1 x 1 each), the pair's pseudospin term; the single spins' clusters are 1, so
 *  that the cCCE assembly gives prod_pairs P. Writes the pair's channels to the channel file. */
void calCoherencePseudospinInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result);
/** Every rank: open its part of "{outfile}_pseudospin_channels" (before calculate) */
void pseudospinChannelsBegin(Output* op);
/** Every rank: close its part; rank 0 then joins the parts into the one file (after calculate) */
void pseudospinChannelsEnd();

#endif
