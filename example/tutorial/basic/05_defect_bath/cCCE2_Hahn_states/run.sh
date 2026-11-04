#!/bin/sh
# basic/05_defect_bath/cCCE2_Hahn_states: cCCE2 Hahn echo of the NV in a 1 ppm P1 bath, 2 configurations x 10 sampled states
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png --unit us "mean of 20=results/ccein_k*_state*_wiDiv" "k1 state 1=results/ccein_k1_state1_wiDiv"
