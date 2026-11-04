#!/bin/sh
# basic/02_bath_and_states/cCCE2_Hahn_ensemble: the cCCE2 Hahn echo of the NV in four 13C bath configurations (k_range) and their mean
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
sh $ROOT/example/tutorial/make_natab_baths.sh 4 results/bath     # results/bath/config_1..4 (about 5 s the first time)
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png "mean of 4=results/ccein_k*_wiDiv" \
  "k1=results/ccein_k1_wiDiv" "k2=results/ccein_k2_wiDiv" "k3=results/ccein_k3_wiDiv" "k4=results/ccein_k4_wiDiv"
