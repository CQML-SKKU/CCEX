#!/bin/sh
# advanced/11_cluster_shapes/cCCE2_Hahn_adaptive: adaptive refinement of the bundles (merge policy)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
# exit status 0: converged; 2: "no feasible improving merge" (the refinement stopped there; the last run is the result)
mpirun -n $NP $CCEX -f ccein.jsonc > results/ccein.log || echo "adaptive refinement ended with status $?"
python3 -c 'import json; s = json.load(open("results/adaptive/summary.json")); print("status:", s["status"], "| rounds:", len(s["rounds"]), "| final bundles:", len(s["final_bundles"]))'
FINAL=$(python3 -c 'import json; print(json.load(open("results/adaptive/summary.json"))["final_directory"])')
python3 $PLOT results/coherence.png "first run=results/adaptive/0001-run/coherence_wiDiv" "final run=$FINAL/coherence_wiDiv"
