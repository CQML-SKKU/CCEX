#!/bin/sh
# advanced/08_coordinate_rotation/cCCE2_Hahn_rotationDefect: a P1 bath with Defect tensors, as generated and turned into
# the NV axis along [111], with the Defect species declared "frame": "qubit" (right here) and "file"
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in aligned cubic_qubit cubic_file; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit us --dashed "aligned (reference)=results/ccein_aligned_k1_state*_wiDiv" \
    "[111] bath, Defect frame qubit (right)=results/ccein_cubic_qubit_k1_state*_wiDiv" "[111] bath, Defect frame file (wrong)=results/ccein_cubic_file_k1_state*_wiDiv"
