#!/bin/sh
# basic/02_bath_and_states/cCCE2_Hahn_states: 20 sampled bath states (s_range [1, 20]) against the ensemble (no s_range) of the same bath
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results/states
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
# the same input without s_range (the ensemble), written by set_input_keys.py as results/ensemble.json
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/ensemble.json configuration.s_range=null
mpirun -n $NP $CCEX -f results/ensemble.json -o results > results/ensemble.log
# the same 20 states from files: make_state_files.py writes the files ccein_files.jsonc names (results/states/state{s})
python3 $ROOT/example/tools/StateFiles/make_state_files.py ccein_files.jsonc --force > results/states.log
mpirun -n $NP $CCEX -f ccein_files.jsonc -o results > results/ccein_files.log
python3 $PLOT results/coherence_files.png "ensemble=results/ensemble_wiDiv" \
  "mean of 20 states from files=results/ccein_files_state*_wiDiv"
python3 $PLOT results/coherence.png "ensemble=results/ensemble_wiDiv" "mean of 20 states=results/ccein_state*_wiDiv" \
  "state 1=results/ccein_state1_wiDiv"
