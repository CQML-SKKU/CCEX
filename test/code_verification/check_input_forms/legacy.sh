#!/usr/bin/env bash
#
# The legacy-form inputs every other stage builds on.
# ---------------------------------------------------------------------------
# Writes _work/<system>/legacy.json for each system: the physics written the way
# v1.1.0 accepted it, with no key this branch added. axis2.sh, reject.sh and
# errcorr.sh derive their variants from these files.
#
# Until 2026-09-24 this was the generating half of axis1.sh, which then ran each
# legacy input under the v1.1.0 binary and this one and required byte identity.
# That axis is retired: since the fixes of 2026-09-19 (deltat as a double, the
# even-pulse conjugate) the numbers differ from v1.1.0 on purpose, and
# ../regression/ measures that difference against v1.1.0 and v1.2.2 instead.
#
# Usage:  bash test/code_verification/check_input_forms/legacy.sh [system ...]
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"


# --- ccein generation -------------------------------------------------------
# Each system writes ONE legacy ccein; the two runs differ only in the binary and
# in "outfile", which is patched per run so the outputs sit side by side.
gen(){  # gen <sysdir> <python body on stdin>
    local sys="$1"; mkdir -p "$WORK/$sys"; python3 - "$WORK/$sys" "$REPO" "$FIX" "$REG" "$EX"
}

# regression case -> absolute-path ccein, config 1 only (byte identity does not
# need the whole ensemble)
gen_reg(){  # gen_reg <sysname> <case> <tag>
    local sys="$1" case="$2" tag="$3"
    mkdir -p "$WORK/$sys"
    python3 - "$WORK/$sys" "$REG/$case" "$tag" <<'PY'
import json, sys, os, glob
work, casedir, tag = sys.argv[1:4]
d = json.load(open(f"{casedir}/inputfiles/ccein_{tag}.json"))
# the file keys (A47: Qubit.file, bath_state_file, Defect.state_file / axis_file)
for block, k in ((d.get("Qubit"),"file"),(d,"gyrofile"),(d,"bath_state_file"),(d.get("Defect"),"state_file"),
                 (d.get("Defect"),"axis_file")):
    if isinstance(block, dict) and isinstance(block.get(k), str):
        block[k] = os.path.abspath(os.path.join(casedir, block[k]))
baths = sorted(glob.glob(f"{casedir}/inputfiles/bath_*_1"))
d["bathfile"] = [os.path.abspath(baths[0])]
d["outfile"] = "./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    {os.path.basename(casedir)}/{tag}: bath={os.path.basename(baths[0])} "
      f"order={d.get('order')} rbath={d.get('rbath')} nstep={d.get('nstep')} nstate={d.get('nstate')}")
PY
}


# --- the systems ------------------------------------------------------------
declare -A SYSTEMS=(
  [S1a]="NV natab, CCE-2, ensemble bath"
  [S1b]="NV natab, CCE-2, sampled bath"
  [S2a]="NV natab, gCCE-2"
  [S2b]="NV + P1 10ppm, CCE-2, Defect n_nuclear_spins=0"
  [S2c]="NV + P1 10ppm, gCCE-2, Defect n_nuclear_spins=0"
  [S3]="P1 bath, Defect n_nuclear_spins=1 (14N), n_configurations=4, axis+state files"
  [S4]="hBN VB, hyperfine dft+contact + quadrupole plain (DFT tensor files, vertex boundary)"
  [S4b]="hBN VB, the same but with the MinDif/MaxDif axis-aligned boundary"
  [S5]="2 NV, gCCE, explicit intmap self/ZFS tensors"
)
ORDER=(S1a S1b S2a S2b S2c S3 S4 S4b S5)

# a generator that fails has written no legacy.json; without this the stage still
# ended "FAIL=0", and only the systems axis2.sh runs were caught later
not_written(){ echo "  FAIL  $sys : legacy.json not written"; FAIL=$((FAIL+1)); }
WANT=("${@:-${ORDER[@]}}")

for sys in "${WANT[@]}"; do
    echo "=== $sys : ${SYSTEMS[$sys]:-?}"
    case "$sys" in
      S1a) gen_reg S1a Diamond_NV_natab_cce   nstate0 || not_written ;;
      S1b) gen_reg S1b Diamond_NV_natab_cce   nstate1 || not_written ;;
      S2a) gen_reg S2a Diamond_NV_natab_gcce  nstate0 || not_written ;;
      S2b) gen_reg S2b Diamond_NV_P1_10ppm_cce  nstate0 || not_written ;;
      S2c) gen_reg S2c Diamond_NV_P1_10ppm_gcce nstate0 || not_written ;;

      S3)  mkdir -p "$WORK/S3"
           python3 - "$WORK/S3" "$REPO" "$FIX" <<'PY' || not_written
import json, re, sys
work, repo, fix = sys.argv[1:4]
src = open(f"{repo}/test/code_verification/bath_data/FullP1SpinBath/cce.json").read()
src = re.sub(r'!\S*','',src); src = re.sub(r'#[^\n]*','',src)   # CCEX json dialect
d = json.loads(src)
F = f"{fix}/single_data"
d["Qubit"]["file"]=f"{F}/bath_DiaP1_1ppm_Defect"; d["gyrofile"]=f"{F}/DiaP1_gyro"
d["bathfile"]=[f"{F}/bath_DiaP1_1ppm_1"];     d["Defect"]["axis_file"]=f"{F}/bathJT_DiaP1_1ppm_1"
d["bath_state_file"]=f"{F}/state_DiaP1_1ppm_"; d["Defect"]["state_file"]=f"{F}/stateEx_DiaP1_1ppm_"
# Reduced from the original (rbath 850 / rdip 620 / nstep 300, 469 spins, sized for
# 24 ranks). Equivalence needs reproducibility, not convergence; rbath 360 keeps
# 45 spins and 518 pairs, which still exercises every Defect code path.
d["rbath"]=360; d["rdip"]=360; d["nstep"]=20; d["nstate"]=1
d["outfile"]="./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
sp = d['Defect']['species'][0]
print(f"    P1 Defect: rbath={d['rbath']} nstep={d['nstep']} n_configurations={sp['n_configurations']} n_nuclear_spins={sp['n_nuclear_spins']}")
PY
           # bath_state_file passed BOTH in the ccein and on the command line, on purpose:
           # a missing -s falls back to randomized states that still look plausible.
           ;;

      S4|S4b)
           # S4  : Afile_vertex / Qfile_vertex -> CheckBD_vertex
           # S4b : Afile        / Qfile        -> CheckBD_Range, the axis-aligned
           #       MinDif/MaxDif box. Two different boundary tests over the same
           #       system, and only one of them is reachable per file.
           SUF=""; [ "$sys" = "S4" ] && SUF="_vertex"
           mkdir -p "$WORK/$sys"
           python3 - "$WORK/$sys" "$FIX" "$SUF" <<'PY' || not_written
import json, sys
work, fix, suf = sys.argv[1:4]
H=f"{fix}/hexagonal_data"
# rbath 20 -> 6 : 3897 spins -> 101, and 121 s -> 0.4 s on one core. The A/Q files
# cover this region comfortably, so both readers stay on their file path.
d = {"method":"cCCE","quantity":"coherence",
     "Qubit":{"file":f"{H}/defect","alphams":1,"betams":0},"gyrofile":f"{fix}/ensemble_data/Gyro_h_BN",
     "bathfile":[f"{H}/bath_1"],
     "order":2,"bfield":30000,"rbath":6,"rdip":4,"deltat":0.002,"nstep":20,
     "nstate":0,"pulse":{"count":1},
     "hamiltonian":{"common":{
     # missing_row point_dipole: the MinDif/MaxDif box in Afile is wider than the set of atoms
     # the file lists, so a spin can sit inside the box with no row of its own. Without
     # this the range-boundary run stops there and never reaches a coherence to compare.
     "hyperfine":{"dipole":"dft","contact":True,"file":f"{H}/Afile{suf}",
                  "missing_row":("stop" if suf else "point_dipole")},
     # The QD reader has no counterpart to missing_row, and Qfile's box is wide in
     # the same way, so the range-boundary case runs without the quadrupole. S4 covers
     # the QD reader over the vertex boundary; the QD range boundary stays uncovered.
     "quadrupole":({"structure":"plain","file":f"{H}/Qfile{suf}"} if suf else "off")}},
     "outfile":"./OUT"}
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    hBN VB: hyperfine dft+contact, quadrupole={'plain' if suf else 'off'} "
      f"rbath={d['rbath']} boundary={'vertex' if suf else 'MinDif/MaxDif'}")
PY
           ;;

      S5)  mkdir -p "$WORK/S5"
           python3 - "$WORK/S5" "$EX" <<'PY' || not_written
import json, re, sys
work, ex = sys.argv[1:3]
E=f"{ex}/tutorial/03_Diamond_multiNV_natab/NV2_13C2"
s=open(f"{E}/ccein_NV6_13C2_hahnecho.json").read()
s=re.sub(r'#[^\n]*','',s); s=re.sub(r'!\S*','',s)
d=json.loads(s)
# bathfile as a LIST. Given as a bare string it is read as zero bath files and the
# run segfaults with "Total 1-Cluster # : 0" -- v1.1.0 behaves the same way, so this
# is the example being stale, not anything this branch changed.
d["bathfile"]=[f"{E}/bath"]; d["bath_state_file"]=f"{E}/state"; d["gyrofile"]=f"{E}/gyro_13C"
d["nstep"]=20; d["outfile"]="./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    2 NV: nqubit={len(d['Qubit']['qubit'])} intmap={len(d['Qubit']['intmap'])} method={d['method']}")
PY
           ;;
      *) echo "  unknown system: $sys"; FAIL=$((FAIL+1)) ;;
    esac
    echo
done

echo "==========================================="
echo " legacy inputs : written for ${#WANT[@]} systems, FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
