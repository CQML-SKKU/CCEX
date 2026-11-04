#!/bin/sh
# basic/05_defect_bath/cCCE2_Hahn_meanField: cCCE2 in 1 ppm P1 with the static field of the other bath spins off, z (default) and xyz
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in z off xyz; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "z (default)=results/ccein_z_state*_wiDiv" \
    "off=results/ccein_off_state*_wiDiv" "xyz=results/ccein_xyz_state*_wiDiv"
