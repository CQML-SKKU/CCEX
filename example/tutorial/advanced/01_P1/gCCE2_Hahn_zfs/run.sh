#!/bin/sh
# advanced/01_P1/gCCE2_Hahn_zfs: gCCE2 Hahn echo at zero field with and without the NV zero-field splitting
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in with_ZFS without_ZFS; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "with ZFS=results/ccein_with_ZFS_state1_wiDiv" "without ZFS=results/ccein_without_ZFS_state1_wiDiv"
