#!/usr/bin/env bash
#
# Legacy <-> new input form equivalence suite, under one binary.
#
#     bash test/code_verification/check_input_forms/run.sh          run everything, in order
#     bash test/code_verification/check_input_forms/run.sh axis2    just one stage
#
# Needs bin/main.out built (bash do_compile.sh wsl-gnu), or CCEX_BIN; no old binary (the
# comparison with v1.1.0 is ../check_numerics/). Everything written lands in _work/,
# which is gitignored.
#
# The stages are ordered by what depends on what, and the suite stops at the first
# failure: axis2 compares the legacy inputs legacy writes, and errcorr reads axis2's outputs.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The lock of lib.sh (D29), held here across the three stages so that no other run gets in between
# them (review 8fecd5 R2); the stages inherit it. The folder is named as lib.sh names it.
LOCKDIR="$(cd "$HERE/../../.." && pwd)/test/code_verification/check_input_forms"
if [ "${CCEX_RUN_LOCK:-}" != "$LOCKDIR" ]; then
    exec 9>"$LOCKDIR/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $LOCKDIR (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$LOCKDIR"; export CCEX_RUN_LOCK
fi
STAGES=(legacy axis2 errcorr gyro comments nodiv allfull evolution_update column_evolution gpu_options propagator tensor_format cli configuration autoconv units pcce_rest dsjitb)
WANT=("${@:-${STAGES[@]}}")
rc_all=0

for stage in "${WANT[@]}"; do
    echo
    echo "###########################################################"
    echo "# $stage"
    echo "###########################################################"
    case "$stage" in
        legacy)   bash "$HERE/legacy.sh" ;;
        axis2)    bash "$HERE/axis2.sh" ;;
        errcorr)  bash "$HERE/errcorr.sh" ;;
        gyro)     bash "$HERE/gyro.sh" ;;
        comments) bash "$HERE/comments.sh" ;;
        nodiv)    bash "$HERE/nodiv.sh" ;;
        allfull)  bash "$HERE/allfull.sh" ;;
        evolution_update) bash "$HERE/evolution_update.sh" ;;
        column_evolution) bash "$HERE/column_evolution.sh" ;;
        gpu_options) bash "$HERE/gpu_options.sh" ;;
        propagator) bash "$HERE/propagator.sh" ;;
        tensor_format) bash "$HERE/tensor_format.sh" ;;
        cli) bash "$HERE/cli.sh" ;;
        configuration) bash "$HERE/configuration.sh" ;;
        autoconv) bash "$HERE/autoconv.sh" ;;
        units) bash "$HERE/units.sh" ;;
        pcce_rest) bash "$HERE/pcce_rest.sh" ;;
        dsjitb) bash "$HERE/dsjitb.sh" ;;
        *) echo "unknown stage: $stage"; rc_all=1; continue ;;
    esac
    rc=$?
    [ $rc -ne 0 ] && { echo "STAGE FAILED: $stage (exit=$rc)"; rc_all=$rc; break; }
done

echo
echo "==========================================================="
[ $rc_all -eq 0 ] && echo " all stages passed" || echo " FAILED (exit=$rc_all)"
echo "==========================================================="
exit $rc_all
