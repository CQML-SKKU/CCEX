#!/bin/sh
# advanced/08_analysis/cCCE2_Hahn_clusterReport: what the run reports about its clusters - a contribution table, the
# cluster plan alone, and the coherence of every cluster of a small bath
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in contributions plan curves; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "rbath 50 A=results/ccein_contributions_wiDiv" "rbath 15 A=results/ccein_curves_wiDiv"
