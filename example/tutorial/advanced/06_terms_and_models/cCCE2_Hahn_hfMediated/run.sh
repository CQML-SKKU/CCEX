#!/bin/sh
# advanced/06_terms_and_models/cCCE2_Hahn_hfMediated: cCCE2 with and without the hyperfine-mediated term
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in wi_hfmedi wo_hfmedi; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "with hfmedi=results/ccein_wi_hfmedi_wiDiv" "without=results/ccein_wo_hfmedi_wiDiv"
