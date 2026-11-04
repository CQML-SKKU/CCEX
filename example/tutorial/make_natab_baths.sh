#!/bin/sh
# make_natab_baths.sh N OUTDIR [SIZE]: N natural-abundance 13C bath configurations around an NV, OUTDIR/config_1 .. config_N
# and OUTDIR/config_defect (the NV position), made with tools/BathGenerator from basic/01_first_run/bath/genbath/ (diamond
# cells, 1.1 % 13C) in a SIZE A cell (default 110: rbath up to 53 A). Configuration k uses the random seed k, so every
# call makes the same files; files already there are kept. The configurations are shared over NP MPI ranks (default 4,
# make_natab_baths.py; without mpirun they are made one after another): about 5 s per configuration (SIZE 110), 40 in
# 51 s with NP=4 (i7-14700K, WSL). Needs python3 with mpi4py, numpy and psutil.
set -e
N=$1; OUT=$2; SIZE=${3:-110}
[ -n "$OUT" ] || { echo "usage: make_natab_baths.sh N OUTDIR [SIZE]"; exit 2; }
ROOT=$(git rev-parse --show-toplevel)
GEN=$ROOT/example/tutorial/basic/01_first_run/bath/genbath
mkdir -p $OUT
python3 -c 'import importlib.util, sys; sys.exit(importlib.util.find_spec("mpi4py") is None)' ||
  { echo "make_natab_baths.sh: BathGenerator.py needs python3 with mpi4py (and numpy, psutil)"; exit 1; }
if command -v mpirun > /dev/null; then
  for k in $(seq 1 $N); do
    [ -s $OUT/config_$k ] && continue
    mpirun -n ${NP:-4} python3 $ROOT/example/tutorial/make_natab_baths.py $N $OUT $SIZE
    exit
  done
  exit 0
fi
echo "make_natab_baths.sh: no mpirun, the configurations are made one after another"
for k in $(seq 1 $N); do
  [ -s $OUT/config_$k ] && continue
  # the generator may abort at interpreter exit after writing ("double free"); the file is checked instead
  python3 -c 'import sys, random, runpy, numpy as np; k = int(sys.argv[1]); np.random.seed(k); random.seed(k); sys.argv = sys.argv[2:]; runpy.run_path(sys.argv[0], run_name="__main__")' \
    $k $ROOT/example/tools/BathGenerator/BathGenerator.py -o $OUT/config_$k -i $GEN/POSCAR_pure -I $GEN/POSCAR_NV -p $GEN/prob_natab \
    --size $SIZE --Defect V --Defect_center 0.5 0.5 0.5 > $OUT/process_$k 2>&1 || true
  [ -s $OUT/config_$k ] || { echo "bath $k was not made: see $OUT/process_$k"; exit 1; }
  mv $OUT/config_${k}_defect $OUT/config_defect
done
