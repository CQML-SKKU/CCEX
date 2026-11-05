#!/bin/sh
# advanced/06_hBN/hBN_makeBath: an h-BN bath around V_B- made with example/tools/BathGenerator, the relaxed V_B- cell of
# the VASP run of makeTensors_hBN inserted at its centre
set -e
ROOT=$(git rev-parse --show-toplevel)
mkdir -p results
# POSCAR_VB_flat: the relaxed cell with the vacancy (VB), the CONTCAR of the VASP run without its velocity block.
# POSCAR_pure_flat: the same lattice without the defect, the sites on the ideal grid and the vacancy filled with B
# (the generator needs the two cells to have the same lattice). ../bath/prob_hBN: the natural abundances of
# 10B, 11B, 14N and 15N. The random numbers are seeded (1), so the same bath comes out every time. The generator may
# abort at interpreter exit after writing ("double free"), so the file is checked instead of its exit status.
( python3 -c 'import sys, random, runpy, numpy as np; np.random.seed(1); random.seed(1); sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name="__main__")' \
    $ROOT/example/tools/BathGenerator/BathGenerator.py -o results/bath_hBN -i POSCAR_pure_flat -I POSCAR_VB_flat \
    -p ../bath/prob_hBN --size 60 --Defect VB --opt2D || true ) > results/generator.log 2>&1
[ -s results/bath_hBN ] || { echo "the bath was not made: see results/generator.log"; exit 1; }
echo "spins in the bath: $(head -1 results/bath_hBN | awk '{printf "%d", $1}')"
awk 'NR > 1 {print $4}' results/bath_hBN | sort | uniq -c
echo "the defect position (x y z): $(head -1 results/bath_hBN_defect)"
