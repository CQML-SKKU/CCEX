#!/bin/bash
# Every check under test/code_verification/, at once: each check writes its own _work folder, except check_refusals,
# which writes check_input_forms/_work too, so those two run one after the other. The logs are printed afterwards in
# the order below, each framed as before. Exit status is the first failure's, in that order.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; status=0
CHECKS=(check_analytic check_deer check_master_equation check_numerics check_rotation check_refusals check_input_forms
        check_mean_field check_rdip_ladder check_noise check_pseudospin check_tutorials)
# check_input_forms/lib.sh's lock (D29) keeps a second run of this checkout out of check_input_forms/_work. The checks
# that source lib.sh take it too, so it is taken here once and handed down (CCEX_RUN_LOCK), as check_input_forms/run.sh
# hands it to its stages; within this run, the order of check_refusals and check_input_forms keeps them apart.
LOCKDIR="$(cd "$HERE/../.." && pwd)/test/code_verification/check_input_forms"
if [ "${CCEX_RUN_LOCK:-}" != "$LOCKDIR" ]; then
    exec 9>"$LOCKDIR/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $LOCKDIR; wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$LOCKDIR"; export CCEX_RUN_LOCK
fi
# check_rotation, check_deer and check_master_equation bring the build up to date before they compile their unit tests
# (test/unit_build_env.sh); a relink while the other checks run bin/main.out makes them fail ("Permission denied"), so
# it is done here once, first
bash "$HERE/../unit_build_env.sh" > /dev/null || { echo "ERROR: the build could not be brought up to date" >&2; exit 1; }
LOGS="$(mktemp -d)"; trap 'rm -rf "$LOGS"' EXIT
one(){ bash "$HERE/$1/run.sh" > "$LOGS/$1.log" 2>&1; echo $? > "$LOGS/$1.rc"; }
t0=$SECONDS
for c in "${CHECKS[@]}"; do
    case $c in
        check_refusals) ( one check_refusals; one check_input_forms ) & ;;
        check_input_forms) ;;
        *) one "$c" & ;;
    esac
done
wait
for c in "${CHECKS[@]}"; do
    echo; echo "################ $c/run.sh"; cat "$LOGS/$c.log"; rc=$(cat "$LOGS/$c.rc" 2>/dev/null || echo 1)
    echo "################ $c/run.sh: exit $rc"
    [ "$rc" -ne 0 ] && [ $status -eq 0 ] && status=$rc
done
echo; echo "all checks: $((SECONDS - t0)) s"
exit $status
