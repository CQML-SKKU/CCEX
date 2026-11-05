#!/bin/sh
# advanced/10_cluster_modes/cCCE2_Hahn_bundle: cCCE2 with ordinary clusters and with bundles
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in standard bundle; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
grep -h '^\[adaptive\]' results/ccein_bundle.log || true     # how the spins were bundled
python3 $PLOT results/coherence.png "standard=results/ccein_standard_wiDiv" "bundle=results/ccein_bundle_wiDiv"
