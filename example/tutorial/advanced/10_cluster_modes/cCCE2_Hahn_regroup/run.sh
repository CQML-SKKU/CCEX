#!/bin/sh
# advanced/10_cluster_modes/cCCE2_Hahn_regroup: bundled cCCE2 of the NV in six P1 (five states) refined by regroup until every state is physical
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
# exit status 0: physical (every |L| <= 1.01); 2: stopped first (no culprit, no candidate, budget or change limit)
mpirun -n $NP $CCEX -f ccein.jsonc > results/ccein.log || echo "regroup stopped with status $?"
# the bundles before and after, the unphysical states (U) and the culprit pairs of each round, and the end
python3 -c 'import json; s = json.load(open("results/adaptive/summary.json")); g = s["regroup"]
for r in g["rounds"]:
    v = r.get("verification") or {}
    print("round %d: bundles %s, %d unphysical states (worst: state %d, time point %d, |L| - 1.01 = %.3g)" % (r["round"], r["units"], r["U"], r["s"], r["t"], r["excess"]))
    print("  culprit pairs (g > 0): %s" % "; ".join("%s-%s g=%.3g" % (a, b, f["g"]) for (a, b), f in zip(r["culprits"], r["factors"])))
    print("  accepted: %s -> rerun: %s unphysical states, excess %.3g" % (r.get("accepted"), v.get("U_actual"), v.get("excess_actual", float("nan"))))
print("status:", s["status"], "| changes:", g["changes"], "| final bundles:", s["final_bundles"])'
FINAL=$(python3 -c 'import json; print(json.load(open("results/adaptive/summary.json"))["final_directory"])')
python3 $PLOT results/coherence.png --unit us "state 4, first bundles=results/adaptive/0001-run/coherence_state4_wiDiv" \
  "state 4, final bundles=$FINAL/coherence_state4_wiDiv"
