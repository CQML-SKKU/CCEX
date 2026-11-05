#!/bin/sh
# advanced/08_multi_qubits/gCCE2_Hahn_oneVsTwoNV: gCCE2 Hahn echo of one NV and of two NVs 25 A apart in the same 13C bath
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in single multi2; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "one NV=results/ccein_single_wiDiv" "two NVs=results/ccein_multi2_wiDiv"
