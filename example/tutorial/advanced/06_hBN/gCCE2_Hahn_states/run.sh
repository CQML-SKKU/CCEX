#!/bin/sh
# advanced/06_hBN/gCCE2_Hahn_states: gCCE2 density matrix of V_B- for one sampled bath state
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png --unit us "gCCE2 (dm)=results/ccein_state1_wiDiv"
