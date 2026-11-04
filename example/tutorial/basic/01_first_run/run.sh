#!/bin/sh
# basic/01_first_run: the Hahn echo of an NV in natural-abundance diamond (cCCE2, ensemble); a few seconds on 4 ranks
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log      # writes results/ccein_wiDiv, ccein_noDiv
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png "cCCE2 Hahn=results/ccein_wiDiv"
