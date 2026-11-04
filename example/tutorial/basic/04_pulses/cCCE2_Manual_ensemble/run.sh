#!/bin/sh
# basic/04_pulses/cCCE2_Manual_ensemble: a CPMG-2 written as a Manual sequence against the built-in CPMG-2, ten configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 10 results/bath 110 5 > results/bath.log   # results/bath/config_1..10, no 13C within 5 A
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/builtin.json 'pulse={"count": 2}'
mpirun -n $NP $CCEX -f results/builtin.json -o results > results/builtin.log
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png "built-in CPMG-2=results/builtin_k*_wiDiv" "Manual=results/ccein_k*_wiDiv"
