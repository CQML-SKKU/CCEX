#!/bin/sh
# 04_hBN/gCCE2_Hahn_Overhauser: gCCE2 with and without the Overhauser field of the bath on the qubit
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in with_Overhauser without_Overhauser; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit ns "with the Overhauser field=results/ccein_with_Overhauser_state1_wiDiv" "without=results/ccein_without_Overhauser_state1_wiDiv"
