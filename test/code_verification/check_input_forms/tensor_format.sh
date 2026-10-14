#!/usr/bin/env bash
#
# tensor_format (2026-10-13, (internal notes) "텐서 파일 형식 정리"): hyperfine tensor files of format v3
# (line 1 VASP / QUANTUMESPRESSO, line 2 "# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ...,
# S = ...") are read x 1; files without the provenance line keep the rule before it (VASP line -> x 1, otherwise
# x 0.5/S with S the qubit spin, DefectTotSpin until A47) and say so in the log. On the five h-BN cases of check_numerics (2 ranks, as recorded there):
# T1  the repository's v3 file (the tutorial's v1 numbers / 2, exactly) gives the bytes of the old untagged v1 file
#     (hBN_VB_flat_ATensor_HSE_CellPBE, x 0.5 by the old rule): the new file is the old one read correctly;
# T2  an old file is read as before: the VASP-tagged v1 file (the _tag file until 2026-10-13, made here as "VASP" + the
#     untagged file) gives the bytes of a v3 file holding the same numbers (v3 x 2) - the old rule is x 1 for a VASP
#     line - and, while reference/head-20261012-aa883d0d4 (recorded with that _tag file) exists, its bytes (a binary
#     built with MKL, A46: its _noDiv/dm within 1e-12 relative, as check_numerics judges such a binary); on a CI runner
#     (CI set) the reference half is skipped and says so (A56: the reference is one machine's digits);
# T3  the log of an old file has the warning, the log of a v3 file the format-v3 line and SpinFactor = 1.
# T4  values of 10 characters or more (1234.56789, -12345.6) in the A and Q files of SpinTensorReader's VASP and QE
#     writers, in their format now and in the one before 2026-10-13 (tab-glued fields), are read as written
#     (example/tools/SpinTensorReader/tests/check_wide_fields.py).
# T5  a Q file alone from a VASP OUTCAR whose magnetization is 0 or missing is written (the same Q file), while an A file
#     from it, or with a --spin it contradicts, is still refused (review 20261013-180049 T1;
#     example/tools/SpinTensorReader/tests/check_qfile_spin.py).
# The refused files (a tag that contradicts the source, S other than the qubit spin, a malformed line) are in
# ../check_refusals/.
#
# Usage:  bash test/code_verification/check_input_forms/tensor_format.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/tensor_format"
rm -rf "$W"; mkdir -p "$W"
NUM="$REPO/test/code_verification/check_numerics"
V3="$NUM/inputs/hbn/hBN_VB_flat_ATensor_HSE_CellPBE_v3"
V1="$EX/tutorial/07_hBN_VB_natab/bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE"
OLDREF="$NUM/reference/head-20261012-aa883d0d4"
{ echo VASP; cat "$V1"; } > "$W/v1_tagged"
# v3 with every hyperfine number x 2 (exact decimal doubling of the 6-decimal values): the numbers of v1_tagged
python3 - "$V3" "$W/v3_double" <<'PY'
import sys
from decimal import Decimal
out = []
for line in open(sys.argv[1]).read().split('\n'):
    t = line.split('\t')
    try:
        v = [Decimal(x.strip()) for x in t[1:]] if len(t) == 14 else None
    except Exception:
        v = None
    if v is None:
        out.append(line); continue
    out.append('\t'.join([t[0]] + t[1:4] + ['{:>11}'.format(str(y * 2)) for y in v[3:]]))   # Decimal keeps -0
open(sys.argv[2], 'w').write('\n'.join(out))
PY

run_case(){    # <case> <name> <tensor file>
    local d="$NUM/cases/$1" o="$W/$1/$2"; mkdir -p "$o"
    python3 -c 'import json,sys; j=json.load(open(sys.argv[1])); j["hamiltonian"]["common"]["hyperfine"]["file"]=sys.argv[3]; json.dump(j,open(sys.argv[2],"w"),indent=1)' \
        "$d/ccein.json" "$o.json" "$3"
    ( cd "$d" && "$MPIRUN" -n 2 "$BIN_DEV" -f "$o.json" -o "$o/out" > "$o.log" 2>&1 < /dev/null ) || { echo "    $1/$2 exited $?"; return 1; }
}
same_bytes(){ [ -n "$(ls "$1")" ] && [ "$(ls "$1")" = "$(ls "$2")" ] && diff -rq "$1" "$2" >/dev/null; }

for d in "$NUM"/cases/hbn_*/; do
    case=$(basename "$d"); ok=1
    for f in v3:"$V3" v1:"$V1" v1_tagged:"$W/v1_tagged" v3_double:"$W/v3_double"; do run_case "$case" "${f%%:*}" "${f#*:}" || ok=0; done
    same_bytes "$W/$case/v3" "$W/$case/v1" || { echo "    $case: the v3 file and the untagged v1 file give different outputs"; ok=0; }
    if [ $ok -eq 1 ]; then echo "  PASS  T1 $case: the v3 file gives the bytes of the untagged v1 file (x 0.5 by the old rule)"; PASS=$((PASS+1))
    else echo "  FAIL  T1 $case (see $W/$case)"; FAIL=$((FAIL+1)); fi
    ok=1
    same_bytes "$W/$case/v1_tagged" "$W/$case/v3_double" || { echo "    $case: the VASP-tagged v1 file is not read x 1"; ok=0; }
    if [ -d "$OLDREF/$case" ] && [ -n "${CI:-}" ]; then
        # A56: the reference holds one machine's digits (WSL); on a CI runner (CI set, as GitHub Actions sets it, like
        # evolution_update.sh E1) only the self-contained half is judged (CI run 38487968260: three of five cases differ)
        echo "  SKIP  T2 $case: the bytes of $(basename "$OLDREF") on a CI runner (one machine's digits)"
        what="(the bytes of $(basename "$OLDREF") skipped on a CI runner)"
    elif [ -d "$OLDREF/$case" ] && grep -q "MKL for large clusters : on" "$W/$case/v1_tagged.log"; then
        # A46: a binary built with MKL is held to the reference as check_numerics holds it (its README): the _noDiv
        # and density-matrix files within 1e-12 of the reference relative to it (absolute where it is 0), _wiDiv not judged
        python3 - "$OLDREF/$case" "$W/$case/v1_tagged" <<'PY' || ok=0
import os, sys
ref, run = sys.argv[1:3]
bad = 0
for f in sorted(os.listdir(ref)):
    if f.endswith('_wiDiv'): continue
    a = [complex(t) for l in open(os.path.join(run, f)) for t in l.split()]
    b = [complex(t) for l in open(os.path.join(ref, f)) for t in l.split()]
    w = max((abs(x - y) / (abs(y) or 1.0) for x, y in zip(a, b)), default=0.0) if len(a) == len(b) else float('inf')
    if not w <= 1e-12: print(f'    {os.path.basename(ref)}: {f} differs from the reference by {w:.1e} (relative, MKL build)'); bad = 1
sys.exit(bad)
PY
        what="and $(basename "$OLDREF") within 1e-12 relative on _noDiv/dm (MKL build)"
    elif [ -d "$OLDREF/$case" ]; then
        for f in "$OLDREF/$case"/*; do cmp -s "$f" "$W/$case/v1_tagged/$(basename "$f")" || { echo "    $case: $(basename "$f") differs from $(basename "$OLDREF")"; ok=0; }; done
        what="and the bytes of $(basename "$OLDREF")"
    else what="($(basename "$OLDREF") is gone: only the self-contained half)"; fi
    if [ $ok -eq 1 ]; then echo "  PASS  T2 $case: the VASP-tagged v1 file is read as before (x 1), $what"; PASS=$((PASS+1))
    else echo "  FAIL  T2 $case (see $W/$case)"; FAIL=$((FAIL+1)); fi
    ok=1
    for f in v1 v1_tagged; do grep -q "WARNING: this tensor file has no provenance line" "$W/$case/$f.log" || { echo "    $case/$f: no warning"; ok=0; }; done
    grep -q "Tensor file format v3 (provenance line): the hyperfine of the spin-1 defect, from VASP" "$W/$case/v3.log" &&
        grep -q -- "-->SpinFactor = 1 (read as written)" "$W/$case/v3.log" || { echo "    $case/v3: the log does not say format v3, x 1"; ok=0; }
    grep -q "WARNING: this tensor file has no provenance line" "$W/$case/v3.log" && { echo "    $case/v3: warned"; ok=0; }
    if [ $ok -eq 1 ]; then echo "  PASS  T3 $case: old files warn, the v3 file is logged as read x 1"; PASS=$((PASS+1))
    else echo "  FAIL  T3 $case"; FAIL=$((FAIL+1)); fi
done
out=$(python3 "$EX/tools/SpinTensorReader/tests/check_wide_fields.py" "$BIN_DEV" "$MPIRUN" "$W/wide_fields" 2>&1); rc=$?
echo "$out" | sed 's/^  PASS  /  PASS  T4 /; s/^  FAIL  /  FAIL  T4 /'
if [ $rc -eq 0 ]; then PASS=$((PASS+$(grep -c "^  PASS" <<<"$out"))); else FAIL=$((FAIL+1)); fi
out=$(python3 "$EX/tools/SpinTensorReader/tests/check_qfile_spin.py" "$W/qfile_spin" 2>&1); rc=$?
echo "$out" | sed 's/^  PASS  /  PASS  T5 /; s/^  FAIL  /  FAIL  T5 /'
if [ $rc -eq 0 ]; then PASS=$((PASS+$(grep -c "^  PASS" <<<"$out"))); else FAIL=$((FAIL+1)); fi
echo "tensor_format: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
