#!/bin/sh
# basic/06_defect_bath/cCCE2_Hahn_states: cCCE2 Hahn echo of the NV in 1 ppm P1, 10 configurations x 10 sampled states
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --unit us --members "mean of 100=results/ccein_k*_state*_wiDiv"
