# test/code_verification/

Does the code still do what it did, and what it says? Twelve checks and the data they read. `run.sh` runs all twelve
at once (`check_refusals`, then `check_input_forms`, since both write `check_input_forms/_work`), prints their logs one
after another and exits with the status of the first that failed.

| Folder | Question | Run |
|---|---|---|
| `check_analytic/` | Does CCEX give the closed-form coherence of one ¹³C beside the NV (FID and Hahn echo, order 1, ensemble)? Does `delta_A` bring the cCCE to the exact one-spin curve, and does `bath_bath_mean_field` "off" give what it must? It compares with a closed form, not with CCEX (`check_deer/` is the other such check). | `bash check_analytic/run.sh` |
| `check_deer/` | Do bath pump pulses (DEER) act as the closed form says, once each, on the gCCE CPU paths (propagator eigen/expm, sparse refused) and in megCCE, with `cluster_mode` standard and pcce, whatever the order of the spins in the bath file? | `bash check_deer/run.sh` |
| `check_master_equation/` | Does meCCE (bath jumps) give the cCCE at rate 0 (Hahn, Ramsey, CPMG with 2 and 3 pulses, a pure bath state, the knight term), the adding_MECCE branch's and PyCCE's coherence with bath jumps (branch 1e-6, PyCCE 1e-7, S=1/2 electron qubit, electron-spin bath) and PRL 132, 250401 Eq. 6 for one bath spin (1e-8)? Does megCCE give the gCCE at rate 0 (one qubit and two), conj(PyCCE) with bath jumps and with a qubit jump beside them (as gCCE does) and exp(-g t/2) for a qubit jump in the 0th cluster? With cluster_mode pcce, is gCCE with partition_size 1 the gCCE, and meCCE / megCCE at rate 0 the cCCE / gCCE on the same partition? The Liouville-space superoperators are unit-tested first. | `bash check_master_equation/run.sh` |
| `check_numerics/` | Are the numbers the ones earlier builds gave? 27 calculations against this branch's recorded run (largest difference of `_noDiv`/dm ≤ 1e-9, `_wiDiv` shown, not judged) and the tags v1.1.0, v1.2.2 (with the known differences named) and v2.0.1 (1e-9). A case added after a reference was recorded is "new" there. | `bash check_numerics/run.sh` |
| `check_rotation/` | Does the coordinate-frame rotation do what it says? 143 checks (2026-11-03): a unit test and runs against a hand-rotated input. | `bash check_rotation/run.sh` |
| `check_refusals/` | Are the 396 inputs that must be refused (2026-11-03) refused, most before anything is read, do the eight inputs that must run finish, and is `--help` `-h`? | `bash check_refusals/run.sh` |
| `check_input_forms/` | Do the new input forms say the same physics as the legacy forms? The `axis2` stage's 22 comparisons under one binary (byte for byte, or numeric where a unit conversion sits between the forms, 2026-11-03), the analyzer reading the outputs, then the other stages listed in `check_input_forms/README.md`. | `bash check_input_forms/run.sh` |
| `check_mean_field/` | Does `bath_bath_mean_field` "xyz" (`vector_mean_field` until 2026-10-13) give the default z mean field's curve byte for byte where the answer is known (one P1 with its 14N and no other bath spin)? | `bash check_mean_field/run.sh` |
| `check_rdip_ladder/` | Does each `rdip_ladder` rung equal a separate run at that rdip (≤ 1e-9), and are the inputs it does not cover refused? About 7 min. | `bash check_rdip_ladder/run.sh` |
| `check_noise/` | Does method `noise` give the OU closed forms, the quadrature of a 1/f spectrum and an independent semi-classical CCE code (mode sum, cumulant, spectrum, with `NOISE_REFERENCE_DIR`) for P1 and ¹³C baths (FID, Hahn, CPMG, tolerance ≤ 2e-10, the printed digits), and do the three `compute` values agree with each other and `correlation` with PyCCE? About 8 min. Needs numpy, scipy, PyCCE (skipped without it), and the reference code (2851449 or later) for the half against it. | `bash check_noise/run.sh` |
| `check_pseudospin/` | Is method pseudospin exact where it has to be, and are order 1 and 3, an unknown `channels` value and sampled states refused? It is checked on two spins on the qubit axis, within 1e-9. Times cCCE order 1, it is cCCE order 2 for spins 1/2 and, with `channels` grouped, 3/2. It also meets Seo's eq. 23, and with `pairwise` its channel file's closed form for 3/2. | `bash check_pseudospin/run.sh` |
| `check_tutorials/` | Does every tutorial input (`example/tutorial/`) still run with this build, from its own folder, with `nstep` cut to 3: exit 0, its outputs written, no key reported as not used? Inputs whose bath or state files a tutorial script makes are skipped until they exist. The curves are not judged. | `bash check_tutorials/run.sh` |
| `bath_data/` | Not a check. The frozen baths, tensors and state files the checks read. | — |

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`. `check_numerics` needs a python3 with
numpy, the analyzer stage one with numpy, scipy and matplotlib (the analyzer's command line imports them, see
`test/requirements.txt`). What a check writes goes to its `_work/`, which is ignored and wiped at each run.
Each folder's README says what it checks and how it decides. `check_mean_field/`, `check_rdip_ladder/` and
`check_noise/` have none and are described in `../README.md` (`check_noise/check.py` in its header).
`check_pseudospin/` and `check_tutorials/` say it in their `run.sh`.
