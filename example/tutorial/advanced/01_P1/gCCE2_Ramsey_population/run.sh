#!/bin/sh
# advanced/01_P1/gCCE2_Ramsey_population: the ms = 0 population of the NV in a 10 ppm P1 bath (gCCE2, density matrix)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png --free-y "ms = 0 population=results/ccein_wiDiv"
