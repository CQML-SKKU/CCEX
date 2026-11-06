#!/bin/bash
# run.sh [NP]   Run every case with bin/main.out and compare it against every folder under reference/.
#               One table per reference; the exit status is 0 only if the newest head-* reference
#               (this line of development's own) shows no unexplained difference.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
# macOS has no timeout(1) (GNU coreutils): the same limit by perl's alarm, which every Mac has; exported for the runs
# started by xargs (check_rdip_ladder)
command -v timeout >/dev/null 2>&1 || { timeout(){ perl -e 'alarm shift; exec @ARGV or die "timeout: $ARGV[0]: $!\n"' "$@"; }; export -f timeout; }
# One run at a time per folder: a second run would delete the first one's work files (D29).
# A script this run starts (check_refusals runs legacy.sh) inherits the lock and does not take it again.
if [ "${CCEX_RUN_LOCK:-}" != "$HERE" ]; then
    exec 9>"$HERE/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $HERE (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$HERE"; export CCEX_RUN_LOCK
fi
BIN="${CCEX_BIN:-$REPO/bin/main.out}"; NP="${1:-2}"; WORK="$HERE/_work"
[ -x "$BIN" ] || { echo "no binary at $BIN (bash do_compile.sh wsl-gnu)"; exit 2; }
rm -rf "$WORK"; mkdir -p "$WORK"
# -o is the output folder since A63, and a case's input has no outfile: its copy named out.json gives the outputs their
# names (out_noDiv, ...), from the case's folder, where its relative paths are read
for d in "$HERE"/cases/*/; do n=$(basename "$d"); mkdir -p "$WORK/$n"; s=$(date +%s.%N)
  cp "$d/ccein.json" "$WORK/$n/out.json"
  ( cd "$d" && timeout 900 "${CCEX_MPI:-mpirun}" -n "$NP" "$BIN" -f "$WORK/$n/out.json" -o "$WORK/$n" > "$WORK/$n/run.log" 2>&1 ); e=$?
  printf '%-24s exit=%-3s %6.1fs\n' "$n" "$e" "$(echo "$(date +%s.%N) - $s" | bc)"
done | tee "$WORK/timing.txt"
# tags (v*) are history; the newest head-* folder is the reference this line of development is held to
status=0; deciding=""
for ref in "$HERE"/reference/head-*/; do [ -d "$ref" ] && deciding="$ref"; done
[ -n "$deciding" ] || for ref in "$HERE"/reference/*/; do deciding="$ref"; done
for ref in "$HERE"/reference/*/; do
  history=--history; [ "$ref" = "$deciding" ] && history=""    # a case newer than a reference of the past is "new"
  echo; python3 "$HERE/compare.py" "$WORK" "$ref" --expected "$HERE/expected.json" $history; rc=$?
  [ "$ref" = "$deciding" ] && status=$rc
done
# a case that did not finish (exit or timeout) fails the run whatever the tables say
if grep -qv ' exit=0 ' "$WORK/timing.txt"; then
  echo; echo "cases that did not finish:"; grep -v ' exit=0 ' "$WORK/timing.txt"; status=1
fi
echo; echo "decided by reference/$(basename "$deciding"): exit $status"; exit $status
