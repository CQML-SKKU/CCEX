#!/bin/sh
# advanced/10_cCCE_corrections/cCCE2_Hahn_knightShift: cCCE2 in a tilted field with and without the Knight shift, and with delta_A as well
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in wi_knight wo_knight wi_knight_deltaA; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "with Knight shift=results/ccein_wi_knight_wiDiv" "without=results/ccein_wo_knight_wiDiv" "Knight and delta_A=results/ccein_wi_knight_deltaA_wiDiv"
