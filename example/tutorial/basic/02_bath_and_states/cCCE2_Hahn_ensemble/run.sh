#!/bin/sh
# basic/02_bath_and_states/cCCE2_Hahn_ensemble: the cCCE2 Hahn echo of the NV in ten 13C bath configurations (k_range) and their mean
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 10 results/bath 110 5 > results/bath.log   # results/bath/config_1..10, no 13C within 5 A
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --members "mean of 10=results/ccein_k*_wiDiv"
