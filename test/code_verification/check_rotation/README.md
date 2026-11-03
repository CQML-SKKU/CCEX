# check_rotation/

Checks `frame.rotation`, which turns the system so that the computational z axis is the qubit axis. Positions and
tensors must move with one matrix R about the first qubit, the rest must not move, and directional input without a
`frame` must be refused. The baths are toy systems written for this suite.

```
bash test/code_verification/check_rotation/run.sh     # bin/main.out (do_compile.sh wsl-gnu) and obj/*.o, or CCEX_BUILD_DIR=<CMake build>
bash test/code_verification/check_rotation/dft_legacy_check.sh <reference-binary> [<test-binary>]
```

`run.sh` wipes and rewrites `_work/`, which git ignores. The exit status is 0 only if every check passed.

## What it checks

The scripts rebuild R from the convention (ez = q, ex = q × b, ey = ez × ex, as rows). The log prints positions to
3 decimals, which sets the geometry tolerances. Some calls in `run.sh` pass their own `--tol`; below are the defaults.

| Case | Compared with | Tolerance |
|---|---|---|
| `test_rotation_unit.cpp`: matrix, fixed point, invariants | exact geometry | unit test |
| ¹³C bath; a spin on [111] lands at r0 + (0, 0, √3) | the unrotated log (`check_geometry.py`) | 5e-3 |
| hyperfine, quadrupole and intmap tensors | R T Rᵀ, or unchanged where `frame` says so (`check_tensors.py`) | 0.05 |
| P1 defect bath, identity rotation | the unrotated run | same bytes |
| P1 defect bath, [111] rotation | coherence must differ; `paxes` and `detuning` must not | exact |
| NV defect bath written five ways; on-site Defect units | each other; the legacy input | numerical agreement |
| `bfield` `[0,0,±500]` | Rᵀ [0,0,Bz] = Bz · q̂ (`check_srcfield.py`) | 1e-6 |
| multi-qubit positions and intmap provenance | R about the first qubit (`check_qubits.py`) | 3e-3 |
| A/B: the system written by hand in the qubit frame, no rotation | the original with rotation (`check_ab.py`) | coherence 1e-8 |
| gCCE turned about z with the state turned (`check_zrot_gcce.py`) | the unturned run | 5e-8 in \|L\|; with fixed mₛ it must differ by > 1e-5 |
| `"enabled": true`, no block | the default, `[disabled]` | same bytes |
| refusals: transverse field, no `frame`, removed keys, multi-qubit with a tensor file | non-zero exit, the message, no DFT reader line | — |

## Files

| File | Contents |
|---|---|
| `run.sh` | builds the unit test, writes and runs every case |
| `test_rotation_unit.cpp` | the geometry unit test |
| `check_geometry.py`, `check_tensors.py`, `check_qubits.py`, `check_srcfield.py` | compare an unrotated and a rotated log |
| `check_ab.py`, `make_manual_frame.py` | the A/B test |
| `check_zrot_gcce.py` | the gCCE turn about the qubit axis |
| `dft_legacy_check.sh` | the DFT tensor readers against a reference binary |
| `configuration_block.py` | puts the generators' flat keys into blocks |
| `inputfiles/` | toy baths, axis, gyro and qubit files |
