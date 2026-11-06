#!/bin/sh
# advanced/11_dissipation/megCCE2_Hahn_qubitRelaxation: gCCE2 against megCCE2 with a relaxation jump on the qubit
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in gCCE megCCE; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "gCCE=results/ccein_gCCE_wiDiv" "megCCE, qubit jump S- at 1 rad/ms=results/ccein_megCCE_wiDiv"
