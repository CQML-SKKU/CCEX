#!/bin/sh
# advanced/06_hBN/cCCE2_Hahn_tensors: V_B- in h-BN, cCCE2, five configurations, without the DFT tensors, with each and with both
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in hf0_qd0 hf3_qd2 hf0_qd2 hf3_qd0; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "no DFT tensor=results/ccein_hf0_qd0_k*_wiDiv" \
  "hyperfine + quadrupole=results/ccein_hf3_qd2_k*_wiDiv" \
  "quadrupole only=results/ccein_hf0_qd2_k*_wiDiv" "hyperfine only=results/ccein_hf3_qd0_k*_wiDiv"
