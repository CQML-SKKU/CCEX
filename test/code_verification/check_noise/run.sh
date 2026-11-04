#!/usr/bin/env bash
#
# method noise (src/calculation/noise.cpp) against references outside CCEX: the OU closed forms and scipy quadrature
# for the model spectra, and an independent semi-classical CCE code (the reference) for P1 and 13C baths (check.py has the cases and the tolerance).
#
# Usage:  bash test/code_verification/check_noise/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# Needs:  numpy and scipy on python3; PyCCE for that comparison (skipped without it); the reference half: NOISE_REFERENCE_DIR, the folder of
#         the reference code (2851449 or later); without it that half is skipped and said so
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
REF="${NOISE_REFERENCE_DIR:-}"
status=0
python3 "$HERE/check.py" model "$W/model" "$BIN_DEV" "$MPIRUN" || status=1
python3 "$HERE/check.py" compute "$W/compute" "$BIN_DEV" "$MPIRUN" || status=1
if python3 -c "import os, pycce; os._exit(0)" > /dev/null 2>&1; then
    python3 "$HERE/check.py" pycce "$W/pycce" "$BIN_DEV" "$MPIRUN" || status=1
else
    echo "  SKIP  correlation against PyCCE: no pycce on python3"
fi
GPU_BIN="${CCEX_GPU_BIN:-$REPO/bin/main_gpu.out}"
if [ -x "$GPU_BIN" ] && { [ -e /dev/nvidia0 ] || [ -e /dev/dxg ]; }; then
    python3 "$HERE/check.py" gpu "$W/gpu" "$BIN_DEV" "$MPIRUN" "$GPU_BIN" || status=1
else
    echo "  SKIP  the GPU against the CPU: no $GPU_BIN or no card (set CCEX_GPU_BIN)"
fi
if [ -n "$REF" ] && [ -f "$REF/quantum-bath-noise/run_quantum.py" ]; then
    python3 "$HERE/check.py" reference "$W/reference" "$BIN_DEV" "$MPIRUN" "$REF" || status=1
else
    echo "  SKIP  the baths against the reference code: set NOISE_REFERENCE_DIR"
fi
[ $status -eq 0 ] && echo "check_noise: PASS" || echo "check_noise: FAIL"
exit $status
