#!/usr/bin/env bash
#
# method noise (src/calculation/noise.cpp) against references outside CCEX: the OU closed forms and scipy quadrature
# for the model spectra, and semi-cce's mode_sum for a P1 bath (check.py has the cases and the tolerance).
#
# Usage:  bash test/code_verification/check_noise/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# Needs:  numpy and scipy on python3; the P1 half also semi-cce (SEMI_CCE, default ~/bin/code_develop/semi-cce),
#         without it that half is skipped and said so
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
SEMI="${SEMI_CCE:-$HOME/bin/code_develop/semi-cce}"
status=0
python3 "$HERE/check.py" model "$W/model" "$BIN_DEV" "$MPIRUN" || status=1
if [ -f "$SEMI/quantum-bath-noise/run_quantum.py" ]; then
    python3 "$HERE/check.py" semi "$W/semi" "$BIN_DEV" "$MPIRUN" "$SEMI" || status=1
else
    echo "  SKIP  P1 against semi-cce: no $SEMI/quantum-bath-noise (set SEMI_CCE)"
fi
[ $status -eq 0 ] && echo "check_noise: PASS" || echo "check_noise: FAIL"
exit $status
