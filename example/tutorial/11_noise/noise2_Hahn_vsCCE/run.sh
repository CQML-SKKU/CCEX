#!/bin/sh
# 11_noise/noise2_Hahn_vsCCE: the Hahn echo of the NV in 1 ppm P1 by method noise and by cCCE2, 10 sampled states
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in noise cce; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "cCCE2, mean of 10=results/ccein_cce_state*_wiDiv" \
    "noise, mean of 10=results/ccein_noise_state*_wiDiv"
