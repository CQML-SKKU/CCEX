#!/bin/sh
# 11_noise/noiseModel_OU: the coherence under Ornstein-Uhlenbeck frequency noise, Ramsey, Hahn echo and CPMG-4
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
mkdir -p results
mpirun -n 1 $CCEX -f ccein.jsonc -o results > results/ccein.log
for n in 0 4; do       # the same input with another pulse.count, written by set_input_keys.py as results/count<n>.json
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/count$n.json pulse.count=$n
  mpirun -n 1 $CCEX -f results/count$n.json -o results > results/count$n.log
done
python3 $PLOT results/coherence.png --unit us "Ramsey=results/count0_wiDiv" "Hahn echo=results/ccein_wiDiv" \
    "CPMG-4=results/count4_wiDiv"
