#!/bin/sh
# advanced/05_P1/gCCE2_Hahn_DEER: the NV Hahn echo in 10 ppm P1 (gCCE2) with and without a DEER pump on the central P1 line
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in nopump pump; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "Hahn echo=results/ccein_nopump_state*_wiDiv" "with the P1 pump (DEER)=results/ccein_pump_state*_wiDiv"
