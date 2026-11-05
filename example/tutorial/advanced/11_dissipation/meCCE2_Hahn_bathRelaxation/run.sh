#!/bin/sh
# advanced/11_dissipation/meCCE2_Hahn_bathRelaxation: the NV Hahn echo by cCCE2 and by meCCE2 with 13C relaxation (jump operators)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in cCCE meCCE; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png "cCCE=results/ccein_cCCE_wiDiv" "meCCE, 13C relaxation=results/ccein_meCCE_wiDiv"
