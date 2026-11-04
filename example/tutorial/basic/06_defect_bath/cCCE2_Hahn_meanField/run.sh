#!/bin/sh
# basic/06_defect_bath/cCCE2_Hahn_meanField: cCCE2 in 1 ppm P1 with the static field of the other bath spins z (default), off and xyz
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
for mf in off xyz; do
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/$mf.json hamiltonian.common.bath_bath_mean_field=$mf
  mpirun -n $NP $CCEX -f results/$mf.json -o results > results/$mf.log
done
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --unit us --log-time --dashed \
  "z (default)=results/ccein_k*_state*_wiDiv" "off=results/off_k*_state*_wiDiv" "xyz=results/xyz_k*_state*_wiDiv"
