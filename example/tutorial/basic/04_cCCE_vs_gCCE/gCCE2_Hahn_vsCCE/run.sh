#!/bin/sh
# basic/04_cCCE_vs_gCCE/gCCE2_Hahn_vsCCE: the NV in 13C by cCCE2 and by gCCE2
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in cCCE gCCE; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "cCCE2=results/ccein_cCCE_wiDiv" "gCCE2=results/ccein_gCCE_wiDiv"
