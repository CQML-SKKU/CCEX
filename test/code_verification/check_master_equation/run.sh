#!/usr/bin/env bash
# ME-CCE / ME-gCCE: the superoperators (unit), then rate-0 limits, the branch's outputs and closed forms (B3, B4).
# Usage: bash test/code_verification/check_master_equation/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
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
"$CXX" "$HERE/test_superoperator_unit.cpp" $UNIT_OBJS -o "$W/test_superoperator_unit.out" $CXXFLAGS -w \
    -I"$REPO/external/eigen" -I"$REPO/external/uthash/include/" -I"$REPO/include/" $UNIT_LIBS \
    && "$W/test_superoperator_unit.out" && PASS=$((PASS+1)) || { echo "  FAIL  superoperator unit"; FAIL=$((FAIL+1)); }
# --- Task B3/B4 append the run checks here ---
# the branch's ME system (adding_MECCE test/verification): S=1/2 electron qubit, e-spin bath 1e16 cm^-3, B 300 G, CCE2, Hahn
mk(){ # mk <name> <python dict update>
python3 - "$HERE/inputs" "$W/$1.json" "$W/$1" <<PY
import json, sys
inp, out, of = sys.argv[1:4]
# A47 keys; savemode (coherence) and max_clusters_per_order (all 0; its length must be order+1) are defaults, left out
d = {"method": "cCCE", "quantity": "coherence",
     "Qubit": {"file": f"{inp}/qubit_espin", "spin": 0.5, "alphams": 0.5, "betams": -0.5},
     "gyrofile": f"{inp}/gyro_espin", "bathfile": [f"{inp}/bath_espin"],
     "order": 2, "bfield": 300, "rbath": 500, "rdip": 300, "deltat": 0.04, "nstep": 50, "nstate": 0,
     "add_subclusters": True, "pulse": {"count": 1}}
$2
d["outfile"] = of
json.dump(d, open(out, "w"), indent=1)
PY
}
MEJ(){ echo "d.setdefault('hamiltonian', {})['$1'] = {'jump_operators': $2}"; }   # MEJ <section> <python list>
RATE=0.00015915494309189535    # MHz = 1 rad/ms (the branch's value)
JUMPS="[{'bath_name':'e','operator':'+','rate':$RATE},{'bath_name':'e','operator':'-','rate':$RATE}]"
mk cce          ""
mk cce_ramsey   "d['pulse']={'count':0}"
mk me_r0        "d['method']='meCCE'"
mk me_r0_ramsey "d['method']='meCCE'; d['pulse']={'count':0}"
mk me_bath      "d['method']='meCCE'; $(MEJ meCCE "$JUMPS")"
mk me_bath_d01  "d['method']='meCCE'; $(MEJ meCCE "$JUMPS"); d['deltat']=0.1; d['nstep']=20"     # PyCCE grid (linspace(0,2,21))
# PRL2024 Eq. 6: one bath spin, order 1, FID, infinite temperature, L+ and L- at the same rate
RATE10=0.0015915494309189535   # MHz = 10 rad/ms
JUMPS10="[{'bath_name':'e','operator':'+','rate':$RATE10},{'bath_name':'e','operator':'-','rate':$RATE10}]"
ONE_Z="d['method']='meCCE'; d['bathfile']=[inp + '/bath_1e_z']; d['order']=1; d['rbath']=1000; d['pulse']={'count':0}; d['deltat']=0.05; d['nstep']=41"
mk me1_eq6_g1   "$ONE_Z; $(MEJ meCCE "$JUMPS")"
mk me1_eq6_g10  "$ONE_Z; $(MEJ meCCE "$JUMPS10")"
for c in cce cce_ramsey me_r0 me_r0_ramsey me_bath me_bath_d01 me1_eq6_g1 me1_eq6_g10; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
# E1: meCCE reads hamiltonian.cCCE as cCCE does. hfmedi moves the e-spin system above by 1e-10 only (A^2/omega_e), so
# the knight term carries the check, on check_numerics' natab_knight_shift (the term moves that case by 0.97, noDiv)
KCASE="$REPO/test/code_verification/check_numerics/cases/natab_knight_shift"
python3 -c 'import json, sys; d = json.load(open(sys.argv[1])); d["method"] = "meCCE"; json.dump(d, open(sys.argv[2], "w"), indent=1)' \
    "$KCASE/ccein.json" "$W/me_r0_knight.json"
for c in cce_knight me_r0_knight; do
    j=ccein.json; [ "$c" = cce_knight ] || j="$W/$c.json"
    ( cd "$KCASE" && run_ccex "$BIN_DEV" "$j" "$W/proc_$c" -o "$W/$c" ) || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
chk(){ # chk <label> <compare.py args>
    if python3 "$HERE/compare.py" "${@:2}"; then echo "  PASS  $1"; PASS=$((PASS+1)); else echo "  FAIL  $1"; FAIL=$((FAIL+1)); fi
}
REF="$HERE/reference/branch-2c63eac52"
PYREF="$HERE/reference/pycce"
chk "meCCE rate 0 == CCE (Hahn, noDiv)"          "$W/me_r0_noDiv" "$W/cce_noDiv" --tol 1e-6
chk "meCCE rate 0 == CCE (Hahn, wiDiv)"          "$W/me_r0_wiDiv" "$W/cce_wiDiv" --tol 1e-6
chk "meCCE rate 0 == CCE (Ramsey, complex)"      "$W/me_r0_ramsey_noDiv" "$W/cce_ramsey_noDiv" --tol 1e-6
chk "meCCE rate 0 == cCCE with hamiltonian.cCCE.knight true (natab, Hahn, noDiv)" "$W/me_r0_knight_noDiv" "$W/cce_knight_noDiv" --tol 1e-6
chk "meCCE bath jumps == branch 2c63eac52"       "$W/me_bath_noDiv" "$REF/ccex_mecce_noDiv" --tol 1e-6
chk "meCCE bath jumps == PyCCE (S=1/2)"         "$W/me_bath_d01_noDiv" "$PYREF/pycce_mecce.dat" --subset-b --tol 1e-6
python3 - "$W" <<'PY'
# PRL2024 Eq. 6 on the CCEX grid (deltat 0.05, nstep 41)
import sys, numpy as np
w = sys.argv[1]
HBAR, G = 1.054571729, -17608.597050      # include/common/utilities.h:72 H_BAR, :85 GAMMA_ELECTRON (= gyro_espin)
a = 2 * HBAR * G * G / 500.0 ** 3          # |A_zz| = 2 hbar g^2 / r^3 for the spin at (0,0,500 A) [rad/ms] = 5.23
for tag, g in (("g1", 1.0), ("g10", 10.0)):    # rad/ms: 2g < a (oscillating), 2g > a (overdamped)
    wt = np.sqrt(complex(4 * g * g - a * a))
    with open(f"{w}/eq6_{tag}", "w") as f:
        for k in range(41):
            t = k * 0.05
            L = np.exp(-g * t) * (np.cosh(wt * t / 2) + 2 * g / wt * np.sinh(wt * t / 2))
            f.write(f"{t:.10f} ({L.real:.15e}{L.imag:+.15e}j)\n")
PY
chk "meCCE 1 spin FID == PRL2024 Eq.6 (2g < a)"  "$W/me1_eq6_g1_noDiv"  "$W/eq6_g1"  --tol 1e-8
chk "meCCE 1 spin FID == PRL2024 Eq.6 (2g > a)"  "$W/me1_eq6_g10_noDiv" "$W/eq6_g10" --tol 1e-8
# B4, megCCE: the branch's 2_verify_megcce grid (deltat 0.1, nstep 20, the other keys as the template); the 0th cluster
# alone (order 0) with a qubit jump, S=1/2 and B || z: z and - both damp rho_ab at g/2, so L = exp(-g t/2)
mk gcce          "d['method']='gCCE'; d['deltat']=0.1; d['nstep']=20"
mk meg_r0        "d['method']='megCCE'; d['deltat']=0.1; d['nstep']=20"
mk meg_bath      "d['method']='megCCE'; d['deltat']=0.1; d['nstep']=20; $(MEJ megCCE "$JUMPS")"
mk meg0_qz       "d['method']='megCCE'; d['order']=0; $(MEJ megCCE "[{'bath_name':'qubit','operator':'z','rate':$RATE}]")"
mk meg0_qminus   "d['method']='megCCE'; d['order']=0; $(MEJ megCCE "[{'bath_name':'qubit','operator':'-','rate':$RATE}]")"
for c in gcce meg_r0 meg_bath meg0_qz meg0_qminus; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
chk "megCCE rate 0 == gCCE (noDiv)"              "$W/meg_r0_noDiv" "$W/gcce_noDiv" --tol 1e-6
chk "gCCE == conj(PyCCE) (S=1/2)"                "$W/gcce_noDiv" "$PYREF/pycce_gcce.dat" --conj-b --subset-b --tol 1e-6
chk "megCCE bath jumps == conj(PyCCE) (S=1/2)"   "$W/meg_bath_noDiv" "$PYREF/pycce_megcce.dat" --conj-b --subset-b --tol 1e-6
ONE="$W/one"; python3 -c "[print(f'{k*0.04:.10f} (1+0j)') for k in range(50)]" > "$ONE"   # template grid (deltat 0.04, nstep 50)
chk "megCCE order 0, qubit z == exp(-g t/2)"     "$W/meg0_qz_noDiv" "$ONE" --tol 1e-9 --scale-exp $RATE
chk "megCCE order 0, qubit - == exp(-g t/2)"     "$W/meg0_qminus_noDiv" "$ONE" --tol 1e-9 --scale-exp $RATE
echo " master-equation checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
