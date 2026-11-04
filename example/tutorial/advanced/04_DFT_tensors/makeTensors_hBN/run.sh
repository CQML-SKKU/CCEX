#!/bin/sh
# advanced/04_DFT_tensors/makeTensors_hBN: the hyperfine (Afile) and quadrupole (Qfile) tensor files of V_B- in h-BN made
# from a VASP run with example/tools/SpinTensorReader, and checked against the files kept with that example
set -e
ROOT=$(git rev-parse --show-toplevel)
STR=$ROOT/example/tools/SpinTensorReader
EX=$STR/examples/02.hBN_VB_VASP
mkdir -p results && cd results
# VASP CONTCAR (the relaxed cell, the defect marked by its name VB), then the OUTCAR for --A (hyperfine) and --Q (EFG);
# the tool writes Afile and Qfile in the working directory. Needs numpy, pandas, psutil.
python3 $STR/main.py VASP $EX/VASP_CONTCAR_flat VB --A $EX/VASP_OUTCAR_flat --Q $EX/VASP_OUTCAR_flat > spintensorreader.log
cmp Afile $EX/Afile_v2 && cmp Qfile $EX/Qfile_v2 && echo "Afile and Qfile are the files of the example (Afile_v2, Qfile_v2)"
sed -n 2p Afile                               # the provenance line CCEX reads: tensor format v3, MHz, S = 1
awk 'NR > 22 {print $1}' Afile | sort | uniq -c       # one row per atom of the cell: 53 B and 54 N
tail -n 3 Afile | cut -c1-60                  # the three 14N next to the vacancy: x y z (Angstrom from the defect), iso (MHz)
