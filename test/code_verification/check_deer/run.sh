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
echo " DEER checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
