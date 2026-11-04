#!/bin/sh
# basic/04_cCCE_vs_gCCE/gCCE2_Hahn_levelCrossing: the NV ms -1 / 0 at 1024 G, near the level crossing, by cCCE2 and gCCE2, ten configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 10 results/bath 110 5 > results/bath.log   # results/bath/config_1..10, no 13C within 5 A
for name in cCCE gCCE; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png "cCCE2=results/ccein_cCCE_k*_wiDiv" "gCCE2=results/ccein_gCCE_k*_wiDiv"
