#!/bin/sh
# basic/03_convergence/cCCE2_Hahn_timeAxis: cCCE2 Hahn echo with a 1 us step, a 10 us step, the 13C revivals and a 0.3 ms window
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
K=$ROOT/example/tools/InputKeys/set_input_keys.py
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 20 ../results/bath 150 5 > results/bath.log   # shared, made once
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $K ccein.jsonc results/step10us.json deltat=0.01 nstep=200
mpirun -n $NP $CCEX -f results/step10us.json -o results > results/step10us.log
python3 $K ccein.jsonc results/revivals.json deltat=0.003735387 nstep=536
mpirun -n $NP $CCEX -f results/revivals.json -o results > results/revivals.log
python3 $K ccein.jsonc results/short.json nstep=300
mpirun -n $NP $CCEX -f results/short.json -o results > results/short.log
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --bold-first "revivals=results/revivals_k*_wiDiv" "1 us step=results/ccein_k*_wiDiv" "10 us step=results/step10us_k*_wiDiv"
python3 $ROOT/example/tutorial/plot_coherence.py results/short.png "0.3 ms window=results/short_k*_wiDiv"
