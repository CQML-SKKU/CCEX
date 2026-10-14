#!/bin/bash
# make_reference.sh TAG        Build the tagged version in a scratch worktree and run every case with it
#                              into reference/TAG/ (the outputs and each run.log).
# make_reference.sh current    Run bin/main.out into reference/head-<date>-<commit>/: the reference the
#                              commits after this one are held to (run.sh decides by the newest head-*).
#                              Delete the older head-* folder when the new one replaces it.
set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"; NP="${NP:-2}"
what="$1"
if [ "$what" = current ]; then
  # date and time in the name: two references of one day sort by time, not by commit hash
  name="head-$(date +%Y%m%d%H%M)-$(cd "$REPO" && git rev-parse --short HEAD)"; BIN="$REPO/bin/main.out"
else
  name="$what"; WT="${CCEX_WORKTREE_DIR:-/tmp}/ccex-$name"
  if [ ! -x "$WT/bin/main.out" ]; then
    rm -rf "$WT"; (cd "$REPO" && git worktree add -q "$WT" "$name")
    (cd "$WT" && bash do_compile.sh wsl-gnu > build.log 2>&1) || { echo "build of $name failed: see $WT/build.log"; exit 1; }
  fi
  BIN="$WT/bin/main.out"
fi
OUT="$HERE/reference/$name"; LOG="$HERE/_work/logs-$name"; rm -rf "$OUT" "$LOG"; mkdir -p "$OUT" "$LOG"
for d in "$HERE"/cases/*/; do n=$(basename "$d"); mkdir -p "$OUT/$n"
  ( cd "$d" && timeout 900 "${CCEX_MPI:-mpirun}" -n "$NP" "$BIN" -f ccein.json -o "$OUT/$n/out" > "$LOG/$n.log" 2>&1 ) || echo "$n: exit $? under $name (see $LOG/$n.log)"
done
echo "reference/$name: $(ls "$OUT" | wc -l) cases"
if [ "$what" != current ] && [ -z "${CCEX_KEEP_WORKTREE:-}" ]; then (cd "$REPO" && git worktree remove --force "$WT"); fi
