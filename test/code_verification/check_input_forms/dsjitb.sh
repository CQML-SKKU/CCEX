#!/usr/bin/env bash
#
# cluster.mode dsjitb at order 2 against the standard clusters, on the natab cCCE2 system of
# inputs/Diamond_NV_natab_cce (docs/methods/cce.md#cce_dsj):
#   D1  dsjitb order 2 == standard CCE2                        every value within 1e-12
#   D2  dsjitb lists every bath spin in its single clusters     the log's "Clusters per order" as standard
#   D3  the same on a bath sorted so that spin 0 has an inter-bath pair
# At order 2 the groups are pairs and single spins and every other connected pair is an inter-bath pair, so dsjitb is
# CCE2 written another way. Until 2026-11-05 the inter-bath pairs took their power from the spin before (a row index
# one off), and spin 0's pairs shortened the list of single spins.
#
# Usage:  bash test/code_verification/check_input_forms/dsjitb.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/dsjitb"
rm -rf "$W"; mkdir -p "$W"
cp -r "$REG/Diamond_NV_natab_cce/inputfiles" "$W/"
cd "$W"
python3 - <<'PY'
import json, re
src = open("inputfiles/ccein_nstate0.json").read()
j = json.loads(re.sub(r"(?m)^\s*#.*$", "", src))
# D3: the bath with its spins in order of distance from spin 1 of the file's list, so that the first spin of the
# selection sits among close neighbours (an inter-bath pair at order 2)
lines = open("inputfiles/bath_Diamond_NV_natab_1").read().splitlines()
head, rows = lines[0], [l for l in lines[1:] if l.strip()]
q = [float(x) for x in open("inputfiles/bath_Diamond_NV_natab_defect").read().split()[:3]]
pos = lambda l: [float(x) for x in l.split()[:3]]
d = lambda a, b: sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5
inside = [l for l in rows if d(pos(l), q) <= j["rbath"]]
dense = max(inside, key=lambda l: sum(1 for m in inside if 0 < d(pos(l), pos(m)) < j["rdip"]))
rest = sorted(rows, key=lambda l: d(pos(l), pos(dense)))
open("inputfiles/bath_sorted", "w").write("\n".join([head] + rest) + "\n")
def write(name, mode, bath=None):
    k = json.loads(json.dumps(j))
    if mode != "standard":
        k["cluster"] = {"mode": mode}
    if bath:
        k["configuration"]["bathfile"] = [bath]
    k["output"] = {"file": f"./out/{name}"}
    open(f"{name}.json", "w").write(json.dumps(k, indent=1))
write("standard", "standard"); write("dsjitb", "dsjitb")
write("standard_sorted", "standard", "./inputfiles/bath_sorted"); write("dsjitb_sorted", "dsjitb", "./inputfiles/bath_sorted")
PY
mkdir -p out
for n in standard dsjitb standard_sorted dsjitb_sorted; do
    run_ccex "$BIN_DEV" $n.json $n.log || { echo "  FAIL  $n did not run (see $W/$n.log)"; FAIL=$((FAIL+1)); }
done
close(){ # close <label> <prefixA> <prefixB>: every value of _wiDiv and _noDiv within 1e-12
    if python3 - "$2" "$3" <<'PY'
import sys
def load(p):
    return [[complex(t) for t in l.split()] for l in open(p) if l.strip() and not l.startswith("#")]
ok = True
for f in ("_noDiv", "_wiDiv"):
    a, b = load(sys.argv[1] + f), load(sys.argv[2] + f)
    if len(a) != len(b) or max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb)) > 1e-12:
        ok = False
sys.exit(0 if ok else 1)
PY
    then echo "  PASS  $1"; PASS=$((PASS+1)); else echo "  FAIL  $1"; FAIL=$((FAIL+1)); fi
}
counts(){ grep -o "Clusters per order : [^(]*" "$1"; }
close "D1 dsjitb order 2 == standard CCE2" out/standard out/dsjitb
if [ -n "$(counts standard.log)" ] && [ "$(counts standard.log)" = "$(counts dsjitb.log)" ]; then
    echo "  PASS  D2 dsjitb has the clusters of standard ($(counts dsjitb.log))"; PASS=$((PASS+1))
else
    echo "  FAIL  D2 the clusters differ: standard '$(counts standard.log)', dsjitb '$(counts dsjitb.log)'"; FAIL=$((FAIL+1))
fi
close "D3 the same on a bath whose first spin has close neighbours" out/standard_sorted out/dsjitb_sorted
echo "dsjitb: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
