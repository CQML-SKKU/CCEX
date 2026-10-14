#!/usr/bin/env bash
#
# evolution_update (docs/user/InputFile_Backend.md, 2026-10-11; the default per path since A29, 2026-10-12):
# E1  an explicit "exact" gives the outputs the check_numerics references were recorded with (exact): its 25 cases
#     against its newest head-* reference by its compare.py, same tolerance (1e-9 on _noDiv and dm); skipped on a CI
#     runner (the reference is one machine's digits);
# E2  "recursive" (A16 b, the phases of the time grid carried over from the previous step) runs on the paths it covers
#     and stays within 1e-8 of exact in the raw _noDiv. NV + P1 10 ppm (inputs/Diamond_NV_P1_10ppm_{cce,gcce}, nstep 40:
#     the recurrence runs between steps 0, 16 and 32): cCCE2 and gCCE2 (evolution=vector) with one sampled state; since
#     A18 (2026-10-11) also the ensemble cCCE2 (nstate 0: Hahn, CPMG-8, and a Manual sequence of unequal intervals) and
#     gCCE2 evolution=matrix (the ensemble Hahn echo, and CPMG-8 with one state);
# E3  without the key those paths run recursive: the bytes of an explicit "recursive", and the log says recursive, chosen
#     by default (the value every time evolution dispatches on: simulator_cce.cpp, simulator_gcce.cpp, cce_batch_window.cpp,
#     small_cluster_batch.cpp read cnf->evolution_update and nothing else). Until A56 (2026-10-14) the path was judged
#     by the bytes of "exact" differing too, which on another machine they need not (CI run 38487968260: gcce_vector);
#     now that is only noted. Single-state and ensemble cCCE2, gCCE2 vector and matrix, on check_numerics cases where the
#     two differ in the written digits on WSL (the NV + P1 10 ppm inputs of E2 do not);
# E4  without the key the paths recursive does not cover run exact, silently (not refused): gCCE propagator expm gives
#     the bytes of an explicit "exact";
# E5  backend auto in this CPU-only binary computes every cluster on the CPU, so it takes the CPU default (A36,
#     2026-10-12; exact until then): without the key the bytes of backend cpu (recursive), and a given "recursive" runs
#     (refused until A36) with the same bytes. Single-state and ensemble cCCE2, gCCE2 vector. The per-path default of the
#     GPU build (bin/main_gpu.out) is checked in test/gpu_verification/check_gpu/check_options.py.
# The refusals ("recursive" given elsewhere, other values) are in ../check_refusals/.
#
# Usage:  bash test/code_verification/check_input_forms/evolution_update.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/evolution_update"
rm -rf "$W"; mkdir -p "$W"
NUM="$REPO/test/code_verification/check_numerics"

# run_case <case> <name> <json changes>: a check_numerics case with the changes, 2 ranks as its run.sh,
# into $W/<name>/out*; relative paths in the case's input are read from the case's folder
run_case(){
    local d="$NUM/cases/$1" o="$W/$2"; mkdir -p "$o"
    python3 -c 'import json,sys; j=json.load(open(sys.argv[1])); j.update(json.loads(sys.argv[3])); json.dump(j,open(sys.argv[2],"w"),indent=1)' \
        "$d/ccein.json" "$o.json" "$3"
    ( cd "$d" && "$MPIRUN" -n 2 "$BIN_DEV" -f "$o.json" -o "$o/out" > "$o.log" 2>&1 < /dev/null ) || { echo "    $2 exited $?"; return 1; }
}
# same_bytes <dirA> <dirB>: the same output files, byte for byte
same_bytes(){ [ -n "$(ls "$1")" ] && [ "$(ls "$1")" = "$(ls "$2")" ] && diff -rq "$1" "$2" >/dev/null; }

# E1: an explicit exact on the check_numerics cases against the reference recorded with exact. The reference is one
# machine's digits (like check_numerics itself, which the CI does not run), so on a CI runner (CI set, as GitHub
# Actions sets it) E1 is skipped and says so; another compiler moves the last digits past 1e-9 (CI run 38243036431)
if [ -n "${CI:-}" ]; then
    echo "  SKIP  E1 on a CI runner: the check_numerics reference holds one machine's digits (run it where it was recorded)"
else
mkdir -p "$W/numerics"
for d in "$NUM"/cases/*/; do run_case "$(basename "$d")" "numerics/$(basename "$d")" '{"evolution_update":"exact"}'; done
REF=""; for r in "$NUM"/reference/head-*/; do [ -d "$r" ] && REF="$r"; done
if python3 "$NUM/compare.py" "$W/numerics" "$REF" --expected "$NUM/expected.json" > "$W/numerics.txt"; then
    echo "  PASS  E1 evolution_update exact gives the check_numerics reference $(basename "$REF") ($(tail -1 "$W/numerics.txt"))"
    PASS=$((PASS+1))
else echo "  FAIL  E1 (see $W/numerics.txt)"; FAIL=$((FAIL+1)); fi
fi

# E3: without the key, recursive where it exists. The gCCE rows use hbn_gcce2_dm: since the h-BN hyperfine was halved
# (tensor format v3, 2026-10-13) exact and recursive write the same bytes on hbn_gcce2_oh, which then cannot tell them apart;
# the one-state cCCE row uses p1_cce2e2n for the same reason since A42 (the dim-2/4 Jacobi eigensystems, 2026-10-13: on
# p1_cce2_single, natab_cce2_single and vp_cce2 exact and recursive now write the same bytes)
# name  case  changes
SUPPORTED="cce_single p1_cce2e2n {}
cce_ensemble hbn_cce2_hf3_qd0 {}
gcce_vector hbn_gcce2_dm {}
gcce_matrix hbn_gcce2_dm {\"evolution\":\"matrix\"}"
while read -r name case change; do
    ok=1
    run_case "$case" "e3_$name/default" "$change" || ok=0
    for v in exact recursive; do
        run_case "$case" "e3_$name/$v" "$(python3 -c 'import json,sys; print(json.dumps(dict(json.loads(sys.argv[1]), evolution_update=sys.argv[2])))' "$change" $v)" || ok=0
    done
    d="$W/e3_$name"
    same_bytes "$d/default" "$d/recursive" || { echo "    $name: the default run is not an explicit recursive's bytes"; ok=0; }
    grep -Eq '^ +evolution_update +: +recursive' "$d/default.log" &&
        grep -Eq '^ +chosen by +: +default \(supported path\)' "$d/default.log" || { echo "    $name: the log does not say recursive by default"; ok=0; }
    grep -Eq '^ +evolution_update +: +exact' "$d/exact.log" || { echo "    $name: the log of the explicit exact does not say exact"; ok=0; }
    same_bytes "$d/default" "$d/exact" &&
        echo "    $name: exact and recursive write the same bytes on this machine (the path is judged by the log alone)"
    if [ $ok -eq 1 ]; then echo "  PASS  E3 $name ($case): without the key the run is recursive"; PASS=$((PASS+1))
    else echo "  FAIL  E3 $name (see $d)"; FAIL=$((FAIL+1)); fi
done <<< "$SUPPORTED"

# E4: without the key, exact where recursive does not exist (A29: silently, not refused)
# name  case  changes  what the log names
UNSUPPORTED="expm hbn_gcce2_oh {\"propagator\":\"expm\"} method=gCCE with propagator=expm"
while read -r name case change why; do
    ok=1
    run_case "$case" "e4_$name/default" "$change" || ok=0
    run_case "$case" "e4_$name/exact" "$(python3 -c 'import json,sys; print(json.dumps(dict(json.loads(sys.argv[1]), evolution_update="exact")))' "$change")" || ok=0
    d="$W/e4_$name"
    same_bytes "$d/default" "$d/exact" || { echo "    $name: the default run is not an explicit exact's bytes"; ok=0; }
    grep -Eq '^ +evolution_update +: +exact' "$d/default.log" && grep -Fq "default ($why)" "$d/default.log" ||
        { echo "    $name: the log does not say exact by default ($why)"; ok=0; }
    if [ $ok -eq 1 ]; then echo "  PASS  E4 $name ($case): without the key, $why runs exact"; PASS=$((PASS+1))
    else echo "  FAIL  E4 $name (see $d)"; FAIL=$((FAIL+1)); fi
done <<< "$UNSUPPORTED"

# E5: backend auto in the CPU build is the CPU default (A36): against the backend cpu runs of E3 (the same cases)
AUTO="cce_single p1_cce2e2n
cce_ensemble hbn_cce2_hf3_qd0
gcce_vector hbn_gcce2_dm"
while read -r name case; do
    ok=1
    run_case "$case" "e5_$name/default" '{"backend":"auto"}' || ok=0
    run_case "$case" "e5_$name/recursive" '{"backend":"auto","evolution_update":"recursive"}' || ok=0
    d="$W/e5_$name"
    same_bytes "$W/e3_$name/default" "$d/default" || { echo "    $name: backend auto is not backend cpu's bytes"; ok=0; }
    same_bytes "$d/default" "$d/recursive" || { echo "    $name: a given recursive with backend auto is not the default's bytes"; ok=0; }
    grep -Eq '^ +evolution_update +: +recursive' "$d/default.log" &&
        grep -Eq '^ +chosen by +: +default \(supported path\)' "$d/default.log" || { echo "    $name: the log does not say recursive by default"; ok=0; }
    if [ $ok -eq 1 ]; then echo "  PASS  E5 $name ($case): backend auto in the CPU build runs recursive, as backend cpu"; PASS=$((PASS+1))
    else echo "  FAIL  E5 $name (see $d)"; FAIL=$((FAIL+1)); fi
done <<< "$AUTO"

# E2: recursive against exact on the NV + P1 10 ppm variants
# name  input folder  input file  changes
VARIANTS="cce cce ccein_nstate1.json {}
gcce gcce ccein_nstate1.json {}
cce_ensemble cce ccein_nstate0.json {}
cce_ensemble_cpmg8 cce ccein_nstate0.json {\"pulse\":{\"count\":8}}
cce_ensemble_manual cce ccein_nstate0.json {\"pulse\":{\"count\":2,\"name\":\"Manual\",\"sequence\":[[0.2,\"X\",180],[0.5,\"X\",180],[0.3,\"I\",0]]}}
gcce_matrix gcce ccein_nstate0.json {}
gcce_matrix_cpmg8 gcce ccein_nstate1.json {\"evolution\":\"matrix\",\"pulse\":{\"count\":8}}"
NAMES=""
while read -r name sys file change; do
    NAMES="$NAMES $name"; mkdir -p "$W/$name/out"; cp -r "$REG/Diamond_NV_P1_10ppm_$sys/inputfiles" "$W/$name/"
    python3 - "$W/$name" "$file" "$change" <<'PY2'
import json, re, sys
d, file, change = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
j = json.loads(re.sub(r"(?m)^\s*#.*$", "", open(f"{d}/inputfiles/{file}").read()))
j.update(change)
for value in ("exact", "recursive"):
    json.dump(dict(j, outfile=f"./out/{value}", evolution_update=value), open(f"{d}/{value}.json", "w"), indent=1)
PY2
done <<< "$VARIANTS"
for sys in $NAMES; do
    cd "$W/$sys"; ok=1
    for v in exact recursive; do run_ccex "$BIN_DEV" $v.json $v.log || { echo "    $sys $v exited $?"; ok=0; }; done
    if [ $ok -eq 1 ] && grep -Eq '^ +evolution_update +: +recursive' recursive.log && python3 - <<'PY'
import cmath, glob, math, sys
def read(p): return [complex(l.split()[1].replace("+-", "-")) for l in open(p) if l.strip()]
def same(x, y):    # nan in both is the same; nan on one side only, or nan vs inf, is not (|nan - y| would be skipped by max)
    return all(a == b or (math.isnan(a) and math.isnan(b)) for a, b in ((x.real, y.real), (x.imag, y.imag)))
worst, n, bad = 0.0, 0, 0
for f in sorted(glob.glob("out/exact*_noDiv")):
    a, b = read(f), read(f.replace("exact", "recursive"))
    assert len(a) == len(b) == 40, (f, len(a), len(b))
    for x, y in zip(a, b):
        if cmath.isfinite(x) and cmath.isfinite(y): worst = max(worst, abs(x - y))
        elif not same(x, y): bad += 1
    n += 1
print(f"    recursive - exact, raw _noDiv: max |difference| {worst:.3e} over {n} file(s), {bad} non-finite mismatch(es)")
sys.exit(0 if n and worst <= 1e-8 and bad == 0 else 1)
PY
    then echo "  PASS  E2 $sys: evolution_update recursive runs and agrees with exact to 1e-8"; PASS=$((PASS+1))
    else echo "  FAIL  E2 $sys (see $W/$sys/recursive.log)"; FAIL=$((FAIL+1)); fi
done
echo "evolution_update: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
