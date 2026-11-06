#!/bin/sh
# advanced/09_hamiltonian_terms/gCCE2_Hahn_meanField: the two mean-field keys of hamiltonian.common on the P1 bath of
# advanced/05_P1 (gCCE2, ten sampled states): both at their defaults, the bath one off, the qubit one on
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in default bath_off qubit_on; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us --dashed "defaults (bath z, qubit off)=results/ccein_default_state*_wiDiv" \
    "bath_bath_mean_field off=results/ccein_bath_off_state*_wiDiv" "qubit_bath_mean_field on=results/ccein_qubit_on_state*_wiDiv"
