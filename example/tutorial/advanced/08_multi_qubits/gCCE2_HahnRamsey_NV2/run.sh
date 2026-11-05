#!/bin/sh
# advanced/08_multi_qubits/gCCE2_HahnRamsey_NV2: 2 NV centres and a 13C pair, gCCE2 density matrix, Hahn echo and Ramsey
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in hahn ramsey; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "Hahn=results/ccein_hahn_state1_wiDiv" "Ramsey=results/ccein_ramsey_state1_wiDiv"
