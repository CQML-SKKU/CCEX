#!/bin/sh
# 03_P1/cCCE2_Hahn_VpBath: cCCE2 Hahn echo of the NV in a Vp bath at 10, 12 and 50 ppm
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for ppm in 10 12 50; do
  mpirun -n $NP $CCEX -f ccein_${ppm}ppm.jsonc -o results > results/ccein_${ppm}ppm.log
done
python3 $PLOT results/coherence.png --unit us --log-time "10 ppm=results/ccein_10ppm_state1_wiDiv" "12 ppm=results/ccein_12ppm_state1_wiDiv" \
  "50 ppm=results/ccein_50ppm_state1_wiDiv"
