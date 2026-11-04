#!/bin/sh
# basic/03_convergence/cCCE_Hahn_order: cCCE Hahn echo at order 1 to 4, configuration 1 (order 4 takes about 85 s)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 20 ../results/bath 150 5 > results/bath.log   # shared, made once
for v in 1 2 3 4; do
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/order$v.json order=$v "configuration.k_range=[1, 1]"
  mpirun -n $NP $CCEX -f results/order$v.json -o results > results/order$v.log
done
python3 $ROOT/example/tutorial/plot_convergence.py results/coherence.png --name order \
  --each "results/order{v}_k*_wiDiv" 1 2 3 4
