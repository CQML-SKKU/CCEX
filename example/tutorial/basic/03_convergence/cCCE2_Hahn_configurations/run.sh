#!/bin/sh
# basic/03_convergence/cCCE2_Hahn_configurations: cCCE2 Hahn echo, the mean of the first 1 to 20 configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 20 ../results/bath 150 5 > results/bath.log   # shared, made once
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $ROOT/example/tutorial/plot_convergence.py results/coherence.png --name configurations \
  --first "results/ccein_k*_wiDiv" 1 2 3 5 7 10 15 20
python3 $ROOT/example/tutorial/plot_coherence.py results/spread.png --spread "20 configurations=results/ccein_k*_wiDiv"
