#!/bin/sh
# basic/03_convergence/cCCE2_Hahn_rdip: cCCE2 Hahn echo, rdip 4 to 14 A, ten configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 20 ../results/bath 150 5 > results/bath.log   # shared, made once
for v in 4 5 6 7 8 9 10 11 12 13 14; do
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/rdip$v.json rdip=$v
  mpirun -n $NP $CCEX -f results/rdip$v.json -o results > results/rdip$v.log
done
python3 $ROOT/example/tutorial/plot_convergence.py results/coherence.png --name rdip \
  --each "results/rdip{v}_k*_wiDiv" 4 5 6 7 8 9 10 11 12 13 14
