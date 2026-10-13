# test/code_verification/

Does the code still do what it did, and what it says. Nine checks and the data
they read; `run.sh` runs all nine in order and then exits with the status of the first that failed.

| folder | question | run |
|---|---|---|
| `check_analytic/` | does CCEX give the closed-form coherence of one 13C beside the NV (FID and Hahn echo, order 1, ensemble)? do `delta_A` false and `bath_bath_mean_field` "off" give what they must? compared with a closed form, not with CCEX (`check_deer/` is the other such check) | `bash check_analytic/run.sh` |
| `check_deer/` | do bath pump pulses (DEER) act as the closed form says, once each, on the gCCE CPU paths (propagator eigen/expm; sparse refused)? | `bash check_deer/run.sh` |
| `check_master_equation/` | does meCCE (bath jumps) give the cCCE at rate 0 (Hahn, Ramsey, the knight term), the adding_MECCE branch's and PyCCE's coherence with bath jumps (1e-6, S=1/2 electron qubit, electron-spin bath) and PRL 132, 250401 Eq. 6 for one bath spin (1e-8)? the Liouville-space superoperators are unit-tested first | `bash check_master_equation/run.sh` |
| `check_numerics/` | are the numbers the ones earlier builds gave? 25 calculations against this branch's recorded run (largest difference of `_noDiv`/dm ≤ 1e-9; `_wiDiv` shown, not judged) and the tags v1.1.0 and v1.2.2 (with the known differences named) | `bash check_numerics/run.sh` |
| `check_rotation/` | does the coordinate-frame rotation do what it says? 149 checks: a unit test and runs against a hand-rotated input | `bash check_rotation/run.sh` |
| `check_refusals/` | are the 172 inputs that must be refused refused, most before anything is read, do three that must run run, and is `--help` `-h`? | `bash check_refusals/run.sh` |
| `check_input_forms/` | do the new input forms say the same physics as the legacy forms? 32 comparisons under one binary, the analyzer reading the outputs, then the gyro, comments, nodiv, evolution_update, gpu_options, tensor_format and autoconv stages | `bash check_input_forms/run.sh` |
| `check_mean_field/` | does `bath_bath_mean_field` "xyz" (`vector_mean_field` until 2026-10-13) give the default z mean field's curve byte for byte where the answer is known (one P1 with its 14N and no other bath spin)? | `bash check_mean_field/run.sh` |
| `check_rdip_ladder/` | does each `rdip_ladder` rung equal a separate run at that rdip (≤ 1e-9), and are the inputs it does not cover refused? about 7 min | `bash check_rdip_ladder/run.sh` |
| `bath_data/` | not a check: the frozen baths, tensors and state files the checks read | — |

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`;
`check_numerics` needs a python3 with numpy, the analyzer stage one with numpy, scipy and
matplotlib (the analyzer's command line imports them; `test/requirements.txt`). What a
check writes goes to its `_work/`, which is ignored and wiped at each run.
Each folder's README says what it checks and how it decides; `check_mean_field/`,
`check_rdip_ladder/` and `check_master_equation/` have none and are described in `../README.md`.
