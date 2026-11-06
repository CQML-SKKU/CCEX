#!/bin/sh
# basic/03_convergence/cCCE2_Hahn_rbath: cCCE2 Hahn echo, rbath 20 to 70 A, ten configurations
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 20 ../results/bath 150 5 > results/bath.log   # shared, made once
for v in 20 25 30 35 40 45 50 55 60 65 70; do
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/rbath$v.json rbath=$v
  mpirun -n $NP $CCEX -f results/rbath$v.json -o results > results/rbath$v.log
done
python3 $ROOT/example/tutorial/plot_convergence.py results/coherence.png --name rbath \
  --each "results/rbath{v}_k*_wiDiv" 20 25 30 35 40 45 50 55 60 65 70
