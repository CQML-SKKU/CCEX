#!/usr/bin/env bash
#
# method pseudospin (src/calculation/pseudospin.cpp) where it has to be exact: two bath spins on the qubit axis, the field
# along it, point dipoles. There the terms the pseudospin Hamiltonian drops (the pseudo-secular hyperfine, the non-secular
# dipolar) are zero, so
#  - pair "block" equals cCCE order 2, spins 1/2 and 3/2, free induction, Hahn echo and CPMG-2;
#  - pair "channel" equals it for spins 1/2 (one channel), and for spins 3/2 (no quadrupole: the channels of a chain
#    overlap, the model is off by design) equals L_1 L_2 P with P the closed form of its own channel file (Hahn echo).
# Two spins in the qubit's plane (z = the qubit's), the field along z: the pair coupling has its flip-flips
# |m1,m2> <-> |m1+1,m2+1> and nothing that changes one spin, the hyperfine no pseudo-secular part, so block equals cCCE order
# 2 (spins 1/2 and 3/2) and channel too for spins 1/2, whose flip-flop and flip-flip share no state.
# Two spins 1 with a quadrupole: an axial one along z gives the flip-flops alone and block equals cCCE order 2; an
# asymmetric one (eta, mixing m and m +- 2) adds the double-quantum channels |m1,m2> <-> |m1-2,m2+2>.
# Then the refusals (order 3, one sampled state, ranks_per_gpu).
#
# Usage:  bash test/code_verification/check_pseudospin/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
printf '  10.0 10.0 10.0\n' > "$W/qubit"
printf '2 0 0 0\n 10.0 10.0 13.0 X\n 10.0 10.0 15.5 X\n' > "$W/bath"
printf '2 0 0 0\n 13.0 10.0 10.0 X\n 10.0 12.5 10.0 X\n' > "$W/bath_plane"
qfile(){ # qfile <name> <EFG xx yy zz> (Hartree/Bohr^2, principal axes x y z): every spin takes the species row (no position row)
    local sep; sep=$(printf '%.0s-' $(seq 170))
    { echo "$sep"; printf '  AtomName\t r_1i_X\t r_1i_Y\t r_1i_Z\t EFGxx\t EFGxy\t EFGxz\t EFGyx\t EFGyy\t EFGyz\t EFGzx\t EFGzy\t EFGzz\t[Hatree/Bohr_radius^2]\n'
      echo "$sep"; printf 'X  eQ___:\t  2.044000\t[Q/milibarn]*(10^-1)\n'
      printf 'X  etc____:\t---\t---\t---\t%s\t0.0\t0.0\t0.0\t%s\t0.0\t0.0\t0.0\t%s\n' "$2" "$3" "$4"
      printf 'MinDif[A]:\t -0.1\t -0.1\t -0.1\nMaxDif[A]:\t 0.1\t 0.1\t 0.1\n'; echo "$sep"; } > "$W/$1"
}
qfile Qaxial -0.2 -0.2 0.4
qfile Qeta -0.1 -0.3 0.4
PASS=0; FAIL=0
input(){ # input <name> <method> <spin> <npulse> <order> <extra keys>
    printf 'X %s 6.728284\n' "$3" > "$W/gyro_$1"
    cat > "$W/$1.json" <<JSON
{"method": "$2", "quantity": "coherence", "Qubit": {"qubit": [{"file": "qubit"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
 "configuration": {"gyrofile": "gyro_$1", "bathfile": ["bath"], "seed": 1},
 "order": $5, "bfield": 50, "rbath": 20, "rdip": 10, "deltat": 0.05, "nstep": 200, $6
 "pulse": {"count": $4}, "output": {"file": "out_$1"}}
JSON
}
run(){ ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" ); }
for I in 0.5 1.5; do for np in 0 1 2; do
    input c_${I}_$np cCCE $I $np 2 ""
    input b_${I}_$np pseudospin $I $np 2 '"hamiltonian": {"pseudospin": {"pair": "block"}},'
    input h_${I}_$np pseudospin $I $np 2 '"hamiltonian": {"pseudospin": {"pair": "channel"}},'
    for k in c b h; do run ${k}_${I}_$np || { echo "  FAIL  ${k}_${I}_$np : CCEX exit"; FAIL=$((FAIL+1)); }; done
done; done
input o1 cCCE 1.5 1 1 ""; run o1 || { echo "  FAIL  o1 : CCEX exit"; FAIL=$((FAIL+1)); }
for I in 0.5 1.5; do for np in 0 1 2; do
    input pc_${I}_$np cCCE $I $np 2 ""
    input pb_${I}_$np pseudospin $I $np 2 '"hamiltonian": {"pseudospin": {"pair": "block"}},'
    input ph_${I}_$np pseudospin $I $np 2 '"hamiltonian": {"pseudospin": {"pair": "channel"}},'
    for k in pc pb ph; do sed -i 's/"bathfile": \["bath"\]/"bathfile": ["bath_plane"]/' "$W/${k}_${I}_$np.json"
                          run ${k}_${I}_$np || { echo "  FAIL  ${k}_${I}_$np : CCEX exit"; FAIL=$((FAIL+1)); }; done
done; done
for q in Qaxial Qeta; do
    Q='"common": {"quadrupole": {"structure": "plain", "file": "'$q'"}}'
    input qc_$q cCCE 1 1 2 "\"hamiltonian\": {$Q},"
    input qb_$q pseudospin 1 1 2 "\"hamiltonian\": {$Q, \"pseudospin\": {\"pair\": \"block\"}},"
    for k in qc qb; do run ${k}_$q || { echo "  FAIL  ${k}_$q : CCEX exit"; FAIL=$((FAIL+1)); }; done
done
python3 "$HERE/check.py" "$W"; rc=$?
if [ -f "$W/check.log" ]; then
    PASS=$((PASS+$(grep -c PASS "$W/check.log"))); FAIL=$((FAIL+$(grep -c FAIL "$W/check.log")))
fi
# check.py did not finish (no numpy, an output missing): a failure, not a pass
if [ $rc -ne 0 ] && ! grep -q FAIL "$W/check.log" 2>/dev/null; then echo "  FAIL  check.py exit $rc"; FAIL=$((FAIL+1)); fi
# refusals
refuse(){ # refuse <name> <expected fragment>
    run "$1"; local rc=$?
    if [ $rc -ne 0 ] && grep -qF -- "$2" "$W/proc_$1"; then echo "  PASS  refused $1"; PASS=$((PASS+1))
    else echo "  FAIL  $1 : not refused with \"$2\" (exit $rc)"; FAIL=$((FAIL+1)); fi
}
input r_order pseudospin 0.5 1 3 ""; refuse r_order "method pseudospin is a pair model"
input r_nstate pseudospin 0.5 1 2 ""; sed -i 's/"seed": 1}/"seed": 1, "s_range": [1, 1]}/' "$W/r_nstate.json"; refuse r_nstate "method pseudospin takes the ensemble only"
input r_gpu pseudospin 0.5 1 2 '"execution": {"backend": "auto", "gpu_options": {"ranks_per_gpu": 4}},'
refuse r_gpu "method pseudospin runs on the CPU only, leave out execution.gpu_options.ranks_per_gpu"
echo "check_pseudospin: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
