#!/usr/bin/env bash
#
# Every tutorial input (example/tutorial/**/ccein*.json, .jsonc) still runs with this build: each is run from its own
# folder, as its run.sh runs it, with nstep cut to 3 (the clusters and every reader are the same; only the time axis is
# short), and has to exit 0, write its outputs and leave no "Warning: input key" line (a key no part of the run reads or
# uses). An automatic-convergence input runs with auto_convergence.enabled false (the ordinary calculation of its keys,
# with rbath 30, rdip 8 and configuration 1, which the
# runner would choose);
# an adaptive input runs its refinement with max_rounds 1 (without -o, which the refinement refuses: its folder is
# cluster.adaptive.workdir) and may end with exit 2 (stopped without passing, the last run the result: round_limit,
# docs/user/OutputFiles.md); a cluster_only input writes its plan (.clusters.json) and no curve. An input whose bath or
# state files are made by a script of the tutorial (make_natab_baths.sh, the StateFiles tool in its run.sh) is
# reported as skipped when they are not there. In main.out (no CUDA) the note that execution.gpu_options is not used
# is expected, and not counted, as is the
# note that the auto_convergence block is not read when it is off. The curves are not judged: the
# tutorials' expected.png are compared by eye, and check_numerics holds 27 calculations cut from them to recorded
# numbers.
#
# Usage:  bash test/code_verification/check_tutorials/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
TUT="$REPO/example/tutorial"; SET="$REPO/example/tools/InputKeys/set_input_keys.py"
PASS=0; FAIL=0; SKIP=0
while IFS= read -r f; do
    dir=$(dirname "$f"); id=$(echo "${f#$TUT/}" | tr '/' '_'); id=${id%.json*}; out="$W/$id"; mkdir -p "$out"
    keys=(nstep=3)
    if grep -q '"auto_convergence"' "$f"; then
        keys+=(auto_convergence.enabled=false)
        # the runner finds rbath and rdip and takes the configurations k of the bath files, so the input leaves them
        # out; the ordinary run needs them (one configuration)
        keys+=(rbath=30 rdip=8 'configuration.k_range=[1, 1]')
    fi
    oflag=(-o "$out"); ok_rc=" 0 "
    if grep -q '"adaptive"' "$f"; then
        keys+=(cluster.adaptive.max_rounds=1 cluster.adaptive.workdir="$out/adaptive"); oflag=(); ok_rc=" 0 2 "
    fi
    if ! python3 "$SET" "$f" "$out/ccein.json" "${keys[@]}" > "$out/set.log" 2>&1; then
        echo "  FAIL  $id : set_input_keys.py ($(tail -1 "$out/set.log"))"; FAIL=$((FAIL+1)); continue
    fi
    # bath and state files a tutorial script writes (01, 09): skipped until made
    missing=$(cd "$dir" && python3 - "$out/ccein.json" <<'PY'
import json, glob, sys, re
c = json.load(open(sys.argv[1])).get('configuration', {})
k = c.get('k_range', [1, 1])[0] if isinstance(c.get('k_range'), list) else 1
s = c.get('s_range', [1, 1])[0] if isinstance(c.get('s_range'), list) else 1
for b in c.get('bathfile', []) + [x for x in c.get('bath_state_file', []) if x]:
    p = b.replace('{k}', str(k)).replace('{s}', str(s))
    if '{s}' not in b and b in c.get('bath_state_file', []): p += str(s)   # a name without {s} has the state appended
    if not glob.glob(p): print(p); break
PY
)
    if [ -n "$missing" ]; then echo "  SKIP  $id : no $missing (made by the tutorial's script)"; SKIP=$((SKIP+1)); continue; fi
    t0=$SECONDS
    ( cd "$dir" && timeout 600 "$MPIRUN" -n 2 "$BIN_DEV" -f "$out/ccein.json" "${oflag[@]}" > "$out/run.log" 2>&1 < /dev/null ); rc=$?
    warn=$(grep "Warning: input key" "$out/run.log" | grep -v 'gpu_options" is not used by this run: this main.out is built without CUDA' \
          | grep -vc 'input key "auto_convergence" is not read by this run')
    nout=$(find "$out" \( -name '*_wiDiv' -o -name '*_noise_correlation' -o -name '*.clusters.json' \) | wc -l)
    if [[ "$ok_rc" == *" $rc "* ]] && [ "$warn" -eq 0 ] && [ "$nout" -gt 0 ]; then
        echo "  PASS  $id ($((SECONDS - t0)) s)"; PASS=$((PASS+1))
    else
        echo "  FAIL  $id : exit $rc, $warn not-used warning(s), $nout output(s) (see $out/run.log)"; FAIL=$((FAIL+1))
    fi
done < <(find "$TUT" -name 'ccein*.json*' -not -path '*/results/*' | sort)
echo "check_tutorials: $PASS passed, $FAIL failed, $SKIP skipped"
[ $FAIL -eq 0 ]
