#!/bin/sh
# advanced/06_terms_and_models/cCCE2_Hahn_dsj: cCCE2 with cluster.mode dsj against the standard clusters
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein_standard.jsonc -o results > results/ccein_standard.log
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png "standard=results/ccein_standard_wiDiv" "dsj=results/ccein_wiDiv"
