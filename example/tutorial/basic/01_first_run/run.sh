#!/bin/sh
# basic/01_first_run: the Hahn echo of an NV in natural-abundance diamond (cCCE2, ensemble) in three steps: make a bath,
# run CCEX, plot the coherence with a T2 fit; a few seconds on 4 ranks
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results
# 1. a bath: one natural-abundance 13C configuration around an NV (tools/BathGenerator, seed 1, a 110 A cell), to see
#    what a bath file is. The run below reads the bath kept in bath/, so the results do not depend on this step.
NP=$NP sh $ROOT/example/tutorial/make_natab_baths.sh 1 results/bath > results/bath.log 2>&1 &&
  python3 -c 'import sys, numpy as np
nv = np.loadtxt(sys.argv[2], max_rows=1)[:3]
x = np.loadtxt(sys.argv[1], skiprows=1, usecols=(0, 1, 2)); d = np.linalg.norm(x - nv, axis=1)
print(f"1. results/bath/config_1: {len(x)} 13C, {(d < 50).sum()} within 50 A of the NV, the nearest at {d.min():.2f} A")' \
  results/bath/config_1 results/bath/config_defect ||
  echo "1. no bath made (tools/BathGenerator needs python3 with mpi4py, numpy and psutil; see results/bath.log)"
# 2. the calculation: writes results/ccein_wiDiv and ccein_noDiv
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
echo "2. results/ccein_wiDiv written"
# 3. the plot, with exp(-(t/T2)^n) fitted to the upper envelope of |L| (the fitted points are drawn)
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --fit "cCCE2 Hahn=results/ccein_wiDiv"
