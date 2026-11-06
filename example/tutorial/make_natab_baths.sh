#!/bin/sh
# make_natab_baths.sh N OUTDIR [SIZE [RMIN]]: N natural-abundance 13C bath configurations around an NV, OUTDIR/config_1 .. config_N
# and OUTDIR/config_defect (the NV position), made with tools/BathGenerator from basic/01_first_run/bath/genbath/ (diamond
# cells, 1.1 % 13C) in a SIZE A cell (default 110: rbath up to 53 A). Configuration k uses the random seed k, so every
# call makes the same files; files already there are kept. The configurations are shared over NP MPI ranks (default 4,
# make_natab_baths.py; without mpirun they are made one after another): about 5 s per configuration (SIZE 110), 40 in
# 51 s with NP=4 on a desktop CPU. RMIN > 0 (A, default 0) leaves out the spins closer than RMIN to the NV, whose
# hyperfine coupling is too strong for a weak bath spin. Needs python3 with mpi4py, numpy and psutil.
set -e
N=$1; OUT=$2; SIZE=${3:-110}; RMIN=${4:-0}
[ -n "$OUT" ] || { echo "usage: make_natab_baths.sh N OUTDIR [SIZE [RMIN]]"; exit 2; }
ROOT=$(git rev-parse --show-toplevel)
GEN=$ROOT/example/tutorial/basic/01_first_run/bath/genbath
mkdir -p $OUT
python3 -c 'import importlib.util, sys; sys.exit(importlib.util.find_spec("mpi4py") is None)' ||
  { echo "make_natab_baths.sh: BathGenerator.py needs python3 with mpi4py (and numpy, psutil)"; exit 1; }
if command -v mpirun > /dev/null; then
  for k in $(seq 1 $N); do
    [ -s $OUT/config_$k ] && continue
    mpirun -n ${NP:-4} python3 $ROOT/example/tutorial/make_natab_baths.py $N $OUT $SIZE
    break
  done
else
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
fi
# RMIN: the spins closer than RMIN to the NV out of every configuration (the count of the header lowered to match);
# a second call removes nothing more
[ "$RMIN" = 0 ] || python3 - "$OUT" "$N" "$RMIN" <<'EOF'
import math, sys
out, n, rmin = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
nv = [float(x) for x in open(f"{out}/config_defect").readline().split()[:3]]
for k in range(1, n + 1):
    lines = open(f"{out}/config_{k}").read().splitlines()
    keep = [l for l in lines[1:] if math.dist([float(x) for x in l.split()[:3]], nv) >= rmin]
    first = lines[0].split()[0]
    head = lines[0].replace(first, "%.5f" % (float(first) - (len(lines) - 1 - len(keep))), 1)
    open(f"{out}/config_{k}", "w").write(head + "\n" + "\n".join(keep) + "\n")
EOF
