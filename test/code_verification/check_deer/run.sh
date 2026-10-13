#!/usr/bin/env bash
# DEER bath pulses: the sequence merge (unit), then the closed form and the evolution paths (Task A3).
# Usage: bash test/code_verification/check_deer/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
PASS=0; FAIL=0
UNIT_ENV_FILE="$(mktemp)"
bash "$REPO/test/unit_build_env.sh" > "$UNIT_ENV_FILE" || { rm -f "$UNIT_ENV_FILE"; exit 1; }
unit_value(){ sed -n "s/^$1='\(.*\)'\$/\1/p" "$UNIT_ENV_FILE"; }
CXX="$(unit_value CCEX_UNIT_CXX)"; CXXFLAGS="$(unit_value CCEX_UNIT_CXXFLAGS) $(unit_value CCEX_UNIT_INCLUDES)"
UNIT_LIBS="$(unit_value CCEX_UNIT_LIBS)"; UNIT_OBJS="$(unit_value CCEX_UNIT_OBJS)"
for kv in $(unit_value CCEX_UNIT_ENV); do export "$kv"; done; rm -f "$UNIT_ENV_FILE"
"$CXX" "$HERE/test_merge_unit.cpp" $UNIT_OBJS -o "$W/test_merge_unit.out" $CXXFLAGS -w \
    -I"$REPO/external/eigen" -I"$REPO/external/uthash/include/" -I"$REPO/include/" $UNIT_LIBS \
    && "$W/test_merge_unit.out" && PASS=$((PASS+1)) || { echo "  FAIL  merge unit"; FAIL=$((FAIL+1)); }
# explicit "pulse.bath": {"count": 0} is the same run as no bath block, byte for byte (outputs and nothing else)
CASE="$REPO/test/code_verification/check_numerics/cases/p1_manual_hahn"
python3 - "$CASE/ccein.json" "$W/bn0.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); d["pulse"]["bath"] = {"count": 0}; json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
# nstate 1 with bath_state_file "../../inputs/p1_1ppm/state_": -s as well, and the outputs are <o>_state1_*
ST=../../inputs/p1_1ppm/state_
( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f ccein.json -s "$ST" -o "$W/bn_absent" > "$W/proc_bn_absent" 2>&1 < /dev/null )
( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$W/bn0.json" -s "$ST" -o "$W/bn_zero" > "$W/proc_bn_zero" 2>&1 < /dev/null )
if [ -s "$W/bn_absent_state1_noDiv" ] && cmp -s "$W/bn_absent_state1_noDiv" "$W/bn_zero_state1_noDiv" \
   && cmp -s "$W/bn_absent_state1_wiDiv" "$W/bn_zero_state1_wiDiv"; then
    echo "  PASS  pulse.bath.count 0 == absent (bytes)"; PASS=$((PASS+1))
else echo "  FAIL  pulse.bath.count 0 differs from absent (or no output)"; FAIL=$((FAIL+1)); fi
# --- Task A3 appends the run checks here ---
# six branch cases (adding_MECCE example/z.Tests/DEER): gCCE order 1 ensemble, closed form at order 1 (check_deer.py)
python3 - "$HERE/inputs" "$REPO/test/code_verification/bath_data/single_data/DiaP1_gyro" "$W" <<'PY'
import json, sys, copy
inp, gyro, w = sys.argv[1:4]

def merge(d, extra):   # nested blocks (pulse, Qubit, hamiltonian) merge key by key
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            merge(d[k], v)
        else:
            d[k] = copy.deepcopy(v)
    return d

# the branch's ccein.js in the A47 keys: avaaxfile -> Defect.axis_file, Defect list -> Defect.species
# (dfname -> name, naddspin -> n_nuclear_spins, navaax -> n_configurations; apprx was removed, "it changed nothing"),
# overhaus -> hamiltonian.common.qubit_bath_mean_field, npulse -> pulse.count, addsubclus -> add_subclusters;
# nk [0, 0], hf_readmode 0 and savemode normal are the defaults and are left out
base = {
 "method": "gCCE", "quantity": "coherence", "seed": 3528, "gyrofile": gyro,
 "Qubit": {"qubit": [{"name": "qubit1", "spin": 1.0, "gyro": -17608.597050, "xyz": [0.0, 0.0, 0.0],
            "alphams": -1, "betams": 0, "detuning": [0.0]}],
           "intmap": [{"between": ["qubit1", "qubit1"],
                       "tensor": [[-960000.0, 0, 0], [0, -960000.0, 0], [0, 0, 1920000.0]]}]},
 "hamiltonian": {"common": {"qubit_bath_mean_field": True}},
 "bathfile": [f"{inp}/bath_DEERtest_12spin.txt"],
 "order": 1, "bfield": 180, "rbath": 200, "rdip": 200, "deltat": 2.5e-04, "nstep": 41,
 "add_subclusters": False, "pulse": {"count": 1}, "nstate": 0,
 "Defect": {"axis_file": f"{inp}/avaax_DEERtest",
            "species": [{"name": "P1", "n_nuclear_spins": 0, "n_configurations": 8, "detuning": [
     [1, "e", 79.975539], [2, "e", 59.579552], [3, "e", 59.579552], [4, "e", 59.579552],
     [5, "e", -79.975539], [6, "e", -59.579552], [7, "e", -59.579552], [8, "e", -59.579552]]}]}}

def pump(tones, seq):
    return {"pulse": {"bath": {"count": 1, "defect": "P1", "spin": 0.5, "alphams": 0.5, "betams": -0.5,
                               "energy_shift": tones, "sequence": seq}}}
HAHN = [["0.50", "X", 180.0], ["0.50", "I", 0.0]]
cases = {
 "nopulse":    {},
 "onres":      pump([79.975539], HAHN),
 "offres":     pump([84.98], HAHN),
 "multi_hit":  pump([79.975539, 200], HAHN),   # the branch value: 79.98 is 28 rad*kHz off, outside the window
 "multi_miss": pump([200, 300], HAHN),
 "nearpulse":  pump([79.975539], [["0.505", "X", 180.0], ["0.495", "I", 0.0]]),
}
# the evolution paths on a pure state (nstate 1): matrix, vector, exact and recursive must agree with bath pulses.
# evolution_update is written out: its default is recursive on these paths since 15a136555 (general.cpp:539-548)
st = {"nstate": 1, "bath_state_file": f"{inp}/state_", "seed": 3528}
paths = {"path_matrix":           {"evolution": "matrix", "evolution_update": "exact"},
         "path_matrix_recursive": {"evolution": "matrix", "evolution_update": "recursive"},
         "path_vector":           {"evolution": "vector", "evolution_update": "exact"},
         "path_recursive":        {"evolution": "vector", "evolution_update": "recursive"},
         "path_nopump":           {"evolution": "vector", "evolution_update": "exact"}}
for path, extra in paths.items():
    cases[path] = merge(merge({}, cases["onres"] if path != "path_nopump" else {}), {**st, **extra})
cases["bad_bspin"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"spin": 1.0, "alphams": 1.0, "betams": 0.0}}})
# controller rulings A3 defect-name and ms-values: a defect name no spin has, and levels that are not two levels of S
cases["bad_defect"]   = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"defect": "P2"}}})
cases["bad_ms_equal"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"betams": 0.5}}})
cases["bad_ms_range"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"alphams": 1.5}}})
for name, extra in cases.items():
    d = merge(copy.deepcopy(base), extra); d["outfile"] = f"{w}/{name}"
    json.dump(d, open(f"{w}/{name}.json", "w"), indent=1)
PY
for c in nopulse onres offres multi_hit multi_miss nearpulse; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
for c in path_matrix path_matrix_recursive path_vector path_recursive path_nopump; do   # bath_state_file key AND -s flag (CLAUDE.md)
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" -s "$HERE/inputs/state_" -N 1 || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
if python3 "$HERE/check_deer.py" "$W" "$HERE/inputs"; then PASS=$((PASS+1)); else FAIL=$((FAIL+1)); fi
# a pump whose tones reach no spin of pulse.bath.defect runs (a control), with one warning (ruling A3 defect-name)
for c in offres multi_miss; do
    if grep -qF "Warning: pulse.bath drives no spin" "$W/proc_$c"; then echo "  PASS  $c : off-resonant pump warned"; PASS=$((PASS+1))
    else echo "  FAIL  $c : no off-resonant warning"; FAIL=$((FAIL+1)); fi
done
# refused <case> <message fragment> <label>: the run stops with that message, not later with a kron size error
refused(){
    run_ccex "$BIN_DEV" "$W/$1.json" "$W/proc_$1"; local rc=$?
    if [ $rc -ne 0 ] && grep -qF "$2" "$W/proc_$1"; then echo "  PASS  $3"; PASS=$((PASS+1))
    else echo "  FAIL  $3 (exit $rc)"; grep -m1 -i error "$W/proc_$1"; FAIL=$((FAIL+1)); fi
}
# pulse.bath.spin that is not the driven spin's spin
refused bad_bspin    "pulse.bath.spin 1 differs from the spin 0.5" "pulse.bath.spin mismatch refused"
refused bad_defect   'pulse.bath.defect "P2" names no bath spin' "pulse.bath.defect naming no spin refused"
refused bad_ms_equal "pulse.bath.alphams 0.5 and betams 0.5 must be two different levels" "pulse.bath alphams == betams refused"
refused bad_ms_range "pulse.bath.alphams 1.5 and betams -0.5 must be two different levels" "pulse.bath |alphams| > spin refused"
echo " DEER checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
