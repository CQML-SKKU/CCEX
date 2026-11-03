# test/

The checks, in three folders by what they verify. Each check folder has its own README with the details.

## Which check for which change

| Change | Run |
|---|---|
| anything under `src/` or `include/` (except `clustering/algo_adaptive/`) | `bash test/code_verification/run.sh` and `make numeric-test cpu-test` |
| `src/clustering/algo_adaptive/` | `sh test/cluster_verification/check_adaptive/verify.sh unit`, `regress`, `inprocess` (and `sanitize` for memory-sensitive changes) |
| `src/clustering/bundle/` | the above, plus `make cpu-test` (its unit tests are there) |
| `src/calculation/gpu/` or the CUDA path | `make gpu-test` inside a GPU job, and the jobs of `gpu_verification/server/` |
| `example/tools/AutoConvergence/` | `python -m pytest -q -n 8 tests` in that folder |

## Requirements

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`. The Python parts need numpy, scipy and
matplotlib (`pip install -r test/requirements.txt`). Environment variables:

- `CCEX_BIN`: the binary to test, for example `CCEX_BIN=$PWD/build/bin/main.out` for a CMake build.
- `CCEX_BUILD_DIR`: the CMake build whose objects the unit tests link (`check_rotation`, `check_deer`,
  `check_master_equation`, the `check_pcce` drivers). Without it they use the `do_compile.sh` build.
  `unit_build_env.sh` brings that build up to date first.
- `CCEX_MPI`: the `mpirun` to use. Without it, the Intel MPI of the loaded module, else the one on `PATH`.

What a check writes goes to its `_work/` (or `cpp/out/`) folder, which git ignores and the check wipes. The reference
outputs under `code_verification/check_numerics/reference/` and the adaptive golden answers are committed on purpose
and change only in the commit that is meant to change them.

A build with MKL changes the last digits of large clusters. The checks read the `MKL for large clusters : on` line of
the run log and then compare relatively: `check_numerics` within 1e-12, `check_adaptive`'s `regress` within 1e-5 for
floating-point values (every decision stays exact).

`.github/workflows/ci.yml` runs the fast, machine-independent checks on every push: the unit tests, `check_analytic`,
`check_rotation`, `check_mean_field`, `check_deer`, `check_master_equation`, `check_noise` (without semi-cce and
PyCCE), `check_pseudospin`, `check_tutorials`, `check_refusals`, `check_input_forms` and `check_rdip_ladder` (quick
sizes). `check_numerics`, the long adaptive suites and the GPU checks are run by hand.

## code_verification/: does the code still give what it gave, and what it says

`bash test/code_verification/run.sh` runs all twelve at once (about 8 minutes on 28 cores) and prints their logs in
the order below. `check_refusals` and `check_input_forms` run one after the other, since both write
`check_input_forms/_work`.

| Check | Question |
|---|---|
| `check_analytic/` | Does CCEX match the closed-form coherence of one ¹³C beside the NV (FID and Hahn)? Also `delta_A` against gCCE order 1, and `bath_bath_mean_field` "off". |
| `check_deer/` | Do DEER bath pulses act as the closed form says, on every gCCE and megCCE path, whatever the order of the bath file? |
| `check_master_equation/` | Do meCCE and megCCE reduce to cCCE and gCCE at rate 0, and match PyCCE's curves and the closed forms with jumps? |
| `check_numerics/` | Are the numbers the ones recorded before? 27 calculations against this branch's newest reference (≤ 1e-9) and the releases v1.1.0, v1.2.2 and v2.0.1, with the known differences named in `expected.json`. |
| `check_rotation/` | Does the coordinate-frame rotation do what it says? A unit test and runs against hand-rotated inputs. |
| `check_refusals/` | Are the inputs that must be refused refused, most of them before anything is read? |
| `check_input_forms/` | Do the input forms of the current and older layouts give the same physics, and do the template (`-e`), the command line, units, comments and the other input stages work? |
| `check_mean_field/` | Does the vector mean field give the z mean field's curve where the two must agree? |
| `check_rdip_ladder/` | Does each `rdip_ladder` rung equal a separate run at that rdip? |
| `check_noise/` | Does method noise match the OU closed forms, quadrature, semi-cce and PyCCE? |
| `check_pseudospin/` | Is method pseudospin exact where it has to be (two spins on the qubit axis, Seo's eq. 23)? |
| `check_tutorials/` | Does every tutorial input still run (with `nstep` cut to 3)? |
| `bath_data/` | not a check: the frozen baths, tensors and state files the checks read |

## cluster_verification/: the clustering methods

| Folder | Contents | Run |
|---|---|---|
| `check_adaptive/` | the adaptive refinement (`src/clustering/algo_adaptive/`): unit tests, 33 recorded runs with their golden answers, fuzz and fault-injection scripts | `sh check_adaptive/verify.sh unit\|regress\|inprocess\|sanitize\|fuzz` |
| `check_bundle/` | the records of the bundle clustering experiments; the bundle code is checked by the unit tests in `make cpu-test` and by `check_adaptive` | — |
| `check_pcce/` | pCCE: the partition search, its cache, reproducibility across ranks | `python3 check_pcce/check_reproducibility.py --cpu bin/main.out --output DIR` and the others |

## gpu_verification/: the CUDA backend

`check_gpu/` has the unit tests (`make numeric-test` and `make cpu-test` run the ones that need no GPU, `make
gpu-test` the rest inside a GPU job) and the `check_*.py` scripts that run the CPU and CUDA builds on the same inputs.
`gpu_benchmarks/` times the CPU against the GPU. `server/` holds the Slurm jobs for idm4-gpu; its README records their
last runs and the commits they ran on.
