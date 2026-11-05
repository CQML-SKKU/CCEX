#!/bin/sh
# advanced/05_P1/cCCE2_Hahn_makeBath: a P1 bath made with example/tools/BathGenerator (P1 on the carbon sites), then
# the cCCE2 Hahn echo in it
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results/bath
# POSCAR_pure_P: the diamond cell with its carbon sites named P; POSCAR_NV_P: the same with the NV (N and the vacancy
# V), inserted at the centre. prob_P1_100ppm: P1 (spin 1/2) on 0.01 % of the P sites, no spin on the N.
# The random numbers are seeded (1), so the same bath comes out every time. The generator may abort at interpreter exit
# after writing ("double free"), so the file is checked instead of its exit status.
python3 -c 'import sys, random, runpy, numpy as np; np.random.seed(1); random.seed(1); sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name="__main__")' \
    $ROOT/example/tools/BathGenerator/BathGenerator.py -o results/bath/bath_P1_100ppm -i POSCAR_pure_P -I POSCAR_NV_P \
    -p prob_P1_100ppm --size 300 --Defect V > results/bath/generator.log 2>&1 || true
[ -s results/bath/bath_P1_100ppm ] || { echo "the bath was not made: see results/bath/generator.log"; exit 1; }
echo "P1 in the cube: $(head -1 results/bath/bath_P1_100ppm | awk '{printf "%d", $1}')"
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
python3 $PLOT results/coherence.png --unit us "mean of 10=results/ccein_state*_wiDiv"
