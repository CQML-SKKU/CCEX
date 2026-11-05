#!/bin/sh
# advanced/05_P1/cCCE2_Hahn_onsite14N: cCCE2 of the same P1 bath with the 14N as a detuning (2e) and as a spin (2e2n)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in 2e 2e2n; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us "2e (14N a detuning)=results/ccein_2e_k*_state*_wiDiv" \
    "2e2n (14N a spin)=results/ccein_2e2n_k*_state*_wiDiv"
