#!/bin/sh
# advanced/05_P1/cCCE4_Hahn_states: the P1 bath at order 2 and at order 4, to check that order 2 is enough
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-2}
mkdir -p results
for name in order2 order4; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us --dashed "order 2=results/ccein_order2_k*_state*_wiDiv" \
    "order 4=results/ccein_order4_k*_state*_wiDiv"
