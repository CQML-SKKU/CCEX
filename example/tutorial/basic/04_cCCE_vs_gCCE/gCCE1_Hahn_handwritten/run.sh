#!/bin/sh
# basic/04_cCCE_vs_gCCE/gCCE1_Hahn_handwritten: gCCE1 density matrix of a hand-written one-spin system
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png "gCCE1 (dm)=results/ccein_state1_wiDiv"
