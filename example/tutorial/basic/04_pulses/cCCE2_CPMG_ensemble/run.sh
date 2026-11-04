#!/bin/sh
# basic/04_pulses/cCCE2_CPMG_ensemble: Hahn echo, CPMG-2 and CPMG-4 of the NV in 13C, each at its 13C revival times,
# averaged over ten bath configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 10 results/bath 110 5 > results/bath.log   # results/bath/config_1..10, no 13C within 5 A
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
# the same input with N = 1 and 2 pulses, deltat = 2N / nu_L and the same 4 ms window
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/count1.json pulse.count=1 deltat=0.003735387 nstep=1070
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/count2.json pulse.count=2 deltat=0.007470773 nstep=535
mpirun -n $NP $CCEX -f results/count1.json -o results > results/count1.log
mpirun -n $NP $CCEX -f results/count2.json -o results > results/count2.log
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png \
  "Hahn=results/count1_k*_wiDiv" "CPMG-2=results/count2_k*_wiDiv" "CPMG-4=results/ccein_k*_wiDiv"
