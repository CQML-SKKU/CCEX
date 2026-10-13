#!/bin/bash
# Every check under test/code_verification/, in order. Exit status is the first failure's.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; status=0
for s in check_analytic/run.sh check_deer/run.sh check_numerics/run.sh check_rotation/run.sh check_refusals/run.sh check_input_forms/run.sh check_mean_field/run.sh check_rdip_ladder/run.sh; do
  echo; echo "################ $s"; bash "$HERE/$s"; rc=$?; echo "################ $s: exit $rc"
  [ $rc -ne 0 ] && [ $status -eq 0 ] && status=$rc
done
exit $status
