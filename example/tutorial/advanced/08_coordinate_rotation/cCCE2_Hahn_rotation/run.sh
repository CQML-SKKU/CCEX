#!/bin/sh
# advanced/08_coordinate_rotation/cCCE2_Hahn_rotation: the bath as generated, with the NV axis along [111], with the rotation on, and without it
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in aligned cubic_rotated cubic_norotation; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --dashed "aligned=results/ccein_aligned_wiDiv" "[111] bath, rotation on=results/ccein_cubic_rotated_wiDiv" \
  "[111] bath, rotation off=results/ccein_cubic_norotation_wiDiv"
