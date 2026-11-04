#!/bin/sh
# basic/03_pulses/cCCE2_CPMG_ensemble: Hahn echo, CPMG-2 and CPMG-8 of the NV in 13C (pulse.count 1, 2, 8)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
for n in 1 2; do       # the same input with another pulse.count, written by set_input_keys.py as results/count<n>.json
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/count$n.json pulse.count=$n
  mpirun -n $NP $CCEX -f results/count$n.json -o results > results/count$n.log
done
python3 $PLOT results/coherence.png "Hahn (count 1)=results/count1_wiDiv" "CPMG-2=results/count2_wiDiv" "CPMG-8=results/ccein_wiDiv"
