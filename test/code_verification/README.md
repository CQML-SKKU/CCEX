# test/code_verification/

Does the code still do what it did, and what it says. Nine checks and the data
they read; `run.sh` runs all nine in order and then exits with the status of the first that failed.

| folder | question | run |
|---|---|---|
| `check_analytic/` | does CCEX give the closed-form coherence of one 13C beside the NV (FID and Hahn echo, order 1, ensemble)? does `delta_A` bring the cCCE to the exact one-spin curve, and does `bath_bath_mean_field` "off" give what it must? compared with a closed form, not with CCEX (`check_deer/` is the other such check) | `bash check_analytic/run.sh` |
| `check_deer/` | do bath pump pulses (DEER) act as the closed form says, once each, on the gCCE CPU paths (propagator eigen/expm; sparse refused) and in megCCE, with `cluster_mode` standard and pcce, whatever the order of the spins in the bath file? | `bash check_deer/run.sh` |
| `check_master_equation/` | does meCCE (bath jumps) give the cCCE at rate 0 (Hahn, Ramsey, CPMG with 2 and 3 pulses, a pure bath state, the knight term), the adding_MECCE branch's and PyCCE's coherence with bath jumps (1e-6, S=1/2 electron qubit, electron-spin bath) and PRL 132, 250401 Eq. 6 for one bath spin (1e-8)? does megCCE give the gCCE at rate 0 (one qubit and two), conj(PyCCE) with bath jumps and with a qubit jump beside them (as gCCE does) and exp(-g t/2) for a qubit jump in the 0th cluster? with cluster_mode pcce, is gCCE with partition_size 1 the gCCE, and meCCE / megCCE at rate 0 the cCCE / gCCE on the same partition? the Liouville-space superoperators are unit-tested first | `bash check_master_equation/run.sh` |
| `check_numerics/` | are the numbers the ones earlier builds gave? 25 calculations against this branch's recorded run (largest difference of `_noDiv`/dm ≤ 1e-9; `_wiDiv` shown, not judged) and the tags v1.1.0 and v1.2.2 (with the known differences named) | `bash check_numerics/run.sh` |
| `check_rotation/` | does the coordinate-frame rotation do what it says? 149 checks: a unit test and runs against a hand-rotated input | `bash check_rotation/run.sh` |
| `check_refusals/` | are the 252 inputs that must be refused refused, most before anything is read, do eight that must run run, and is `--help` `-h`? | `bash check_refusals/run.sh` |
| `check_input_forms/` | do the new input forms say the same physics as the legacy forms? 32 comparisons under one binary, the analyzer reading the outputs, then the gyro, comments, nodiv, evolution_update, gpu_options, tensor_format and autoconv stages | `bash check_input_forms/run.sh` |
| `check_mean_field/` | does `bath_bath_mean_field` "xyz" (`vector_mean_field` until 2026-10-13) give the default z mean field's curve byte for byte where the answer is known (one P1 with its 14N and no other bath spin)? | `bash check_mean_field/run.sh` |
| `check_rdip_ladder/` | does each `rdip_ladder` rung equal a separate run at that rdip (≤ 1e-9), and are the inputs it does not cover refused? about 7 min | `bash check_rdip_ladder/run.sh` |
| `check_noise/` | does method `noise` give the OU closed forms, the quadrature of a 1/f spectrum and semi-cce's `mode_sum`, `cumulant` and `psd` for P1 and ¹³C baths (FID, Hahn, CPMG; ≤ 2e-10, the printed digits), and do the three `compute` values agree with each other and `correlation` with PyCCE? about 8 min; numpy, scipy, PyCCE (skipped without it), and semi-cce 2851449 or later for the half against it | `bash check_noise/run.sh` |
| `bath_data/` | not a check: the frozen baths, tensors and state files the checks read | — |

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`;
`check_numerics` needs a python3 with numpy, the analyzer stage one with numpy, scipy and
matplotlib (the analyzer's command line imports them; `test/requirements.txt`). What a
check writes goes to its `_work/`, which is ignored and wiped at each run.
Each folder's README says what it checks and how it decides; `check_mean_field/`, `check_rdip_ladder/` and
`check_noise/` have none and are described in `../README.md` (`check_noise/check.py` in its header).
