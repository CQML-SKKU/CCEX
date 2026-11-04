#!/bin/sh
# basic/03_pulses/cCCE2_Manual_ensemble: a CPMG-2 written as a Manual sequence against the built-in CPMG-2
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/builtin.json 'pulse={"count": 2}'
mpirun -n $NP $CCEX -f results/builtin.json -o results > results/builtin.log
python3 $PLOT results/coherence.png "built-in CPMG-2=results/builtin_wiDiv" "Manual=results/ccein_wiDiv"
