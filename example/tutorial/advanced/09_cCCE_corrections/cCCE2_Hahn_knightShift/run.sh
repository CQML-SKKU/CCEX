#!/bin/sh
# advanced/09_cCCE_corrections/cCCE2_Hahn_knightShift: cCCE2 in a tilted field with and without the Knight shift, and with delta_A as well
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in gcce wi_knight wi_knight_deltaA wo_knight; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --dashed "gCCE (reference)=results/ccein_gcce_wiDiv" \
    "cCCE with knight=results/ccein_wi_knight_wiDiv" "cCCE with knight and delta_A=results/ccein_wi_knight_deltaA_wiDiv" \
    "cCCE without=results/ccein_wo_knight_wiDiv"
