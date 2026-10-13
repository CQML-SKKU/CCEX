# test/

Three folders, named by what they verify; each has a README. Every suite needs
`bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`; the Python ones need what
`test/requirements.txt` lists - numpy, scipy (the analyzer imports it) and matplotlib -
on `python3` (on WSL: `PATH=~/anaconda3/bin:$PATH`; elsewhere `pip install -r test/requirements.txt`).
A CMake build keeps its binary in its own folder: `CCEX_BIN=<build>/bin/main.out` tells
the suites which binary to run, and `CCEX_BUILD_DIR=<build>` tells `check_rotation` and
the `check_pcce` unit-test drivers whose objects to link (the `do_compile.sh` build otherwise). Both ask
`unit_build_env.sh` here, which brings that build up to date and prints its compiler, flags and objects as the
build writes them: the `unit-env` target of the generated Makefile, or the `unit-build.env` of a CMake build of any
generator (until 2026-10-02 the drivers read `flags.make` and so needed the Unix Makefiles generator, D26).
`mpirun` is `CCEX_MPI` when set, else the Intel MPI of the loaded module, else the one
on `PATH`. `.github/workflows/ci.yml` runs the fast, machine-independent part on every
push and pull request (with `check_analytic`, the closed form); `check_numerics` (exact against
references recorded on one machine), the long adaptive suites and the GPU checks stay manual.
What a suite writes goes to a `_work/` or `cpp/out/` folder, which is gitignored;
the regression references under `code_verification/check_numerics/reference/` are committed on purpose.
They were recorded with builds without MKL. A binary built with MKL (A46; the build turns it on when it finds MKL,
and every run log's header says `MKL for large clusters : on` or `off`) changes the last digits of the large
clusters, so the checks read that line: `check_numerics` (and `check_input_forms`' tensor_format T2) then judge
within 1e-12 relative to the reference, and `check_adaptive`'s regress keeps every decision exact but lets a float
of the golden answers move by 1e-5 relative (their READMEs say why).

## code_verification/ — does the code still do what it did, and what it says

`bash test/code_verification/run.sh` runs the nine in order.

| Folder | Answers | Run |
|---|---|---|
| `check_analytic/` | CCEX gives the closed-form coherence of one 13C beside the NV (free induction and Hahn echo, order 1, ensemble), within the 10 decimals it writes; with the closed-form part of `check_deer` the only comparison here not of CCEX with CCEX; machine-independent, so the CI runs it. A self-test requires a 1 % wrong gyro to fail; `hfmedi` runs at zero field check degenerate qubit levels (finite where uncoupled - by S or by a zero hyperfine -, refused where coupled, the same decision after an energy shift of the ZFS); `hamiltonian.cCCE.delta_A` false and `hamiltonian.common.bath_bath_mean_field` "off" (A47) against what they must give | `bash check_analytic/run.sh` |
| `check_deer/` | DEER bath pulses (`pulse.bath`): the qubit/bath sequence merge (unit test), `pulse.bath.count` 0 byte-identical to no block, and the six cases of the adding_MECCE branch (12 P1, gCCE order 1, ensemble) against the closed form (2e-4) and its identities (off-resonant pump = no pump, a second far tone changes nothing, a pump 0.005 beside the qubit pi applied once); on one pure state the evolution vector = matrix, recursive = exact (1e-9) with the pump acting; refusals of a `pulse.bath` that does not fit the driven spins | `bash check_deer/run.sh` |
| `check_master_equation/` | ME-CCE (`method` meCCE, `hamiltonian.meCCE.jump_operators` on bath spins): the Liouville-space superoperators and jump operators (unit test); at rate 0 the cCCE curve (Hahn `_noDiv` and `_wiDiv`, Ramsey with its imaginary part, and the knight term on `check_numerics/cases/natab_knight_shift`; 1e-6); with S+ and S- jumps on an electron-spin bath (S=1/2 electron qubit, CCE2, Hahn) the adding_MECCE branch's output and PyCCE's (1e-6, `reference/`, md5 against the branch in `MD5SUMS`); one bath spin on the qubit's z axis, order 1, FID against PRL 132, 250401 Eq. 6 for 2g < a and 2g > a (1e-8). ME-gCCE (`method` megCCE): at rate 0 the gCCE curve, with bath jumps (and with a qubit `z` jump beside them) conj(PyCCE) and the gCCE itself against conj(PyCCE) (1e-6; PyCCE measures the pulse-turned qubit states), and the 0th cluster alone (order 0) with a qubit `z` or `-` jump against exp(-g t/2) (1e-9). Judged by `compare.py` (exit 2 when a file is missing or the grids differ) | `bash check_master_equation/run.sh` |
| `check_numerics/` | this build gives the numbers earlier builds gave: 25 calculations cut from the tutorials (CCE, gCCE, pCCE, dsj; ensemble and single approaches; Hahn, Ramsey, CPMG, manual sequences; ¹³C, P1, Vp, h-BN baths; coherence and density matrices), compared against outputs recorded on the same machine from this branch (`reference/head-*`, largest difference of `_noDiv`/dm ≤ 1e-9; `_wiDiv` is not judged) and from the tagged releases v1.1.0 and v1.2.2 (1e-4, the known differences named in `expected.json`). `make_reference.sh` records a reference, `make_cases.py` wrote the cases. Its `test_fixed_size.cpp` (in `make cpu-test`) holds the fixed-size single-state evolution of dimension 2 and 4 to the dynamic one's bits with `evolution_update: exact`, and its `test_cce_batch.cpp` (same target) holds the batched small-cluster evolution (A30) to the scalar path's bits on dumped clusters | `bash check_numerics/run.sh` |
| `check_rotation/` | the coordinate-frame rotation does what it says: the matrix and its fixed point (a unit test), then real runs against a hand-rotated input, the same transform for every bath file, the DFT tensor files moving with the positions; 149 checks | `bash check_rotation/run.sh` |
| `check_refusals/` | 172 input combinations that must be refused, most before anything is read (a tensor file without its frame, a tilted field with the rotation on, `reference_qubit` missing, `zfs` in D/E form without an axis, both rotation keys at once, ...), and three that must run, on the real inputs | `bash check_refusals/run.sh` |
| `check_mean_field/` | `hamiltonian.common.bath_bath_mean_field` "xyz" (`vector_mean_field` until 2026-10-13, D45) against the default z mean field where the answer is known: one P1 with its <sup>14</sup>N sub-spin and no other bath spin inside rbath has no mean field, so both modes must give the same curve byte for byte (it failed before review 158861 R1) | `bash check_mean_field/run.sh` |
| `check_rdip_ladder/` | `rdip_ladder` (auto-truncation-plan.md, improvements.md A23): one order-2 run writes the curve of every smaller rdip; each rung must equal a separate run at that rdip and the main output the run without the ladder (≤ 1e-9; cCCE P1 1 ppm with 2 states and <sup>13</sup>C ensemble; gCCE ensemble, vector and matrix on <sup>13</sup>C, h-BN and P1 10 ppm; pCCE2 centers and spins on <sup>13</sup>C, h-BN and P1 1 ppm), and inputs it does not cover (order 3 of cCCE and pCCE, a rung at or above rdip, unsorted rungs, savemode cluster_coherence, gCCE dm) are refused; about 7 min | `bash check_rdip_ladder/run.sh` |
| `check_input_forms/` | the new input forms this branch added say the same physics as the legacy forms (32 comparisons under one binary: 22 byte-identical, 8 numeric where a unit conversion sits between them, 2 that must differ; then the stages `gyro`, `comments`, `nodiv`, `allfull`, `evolution_update`, `gpu_options`, `tensor_format` and `autoconv`), and the analyzer reads the outputs. Its `inputs/` are the four legacy inputs the old regression suite used; `legacy.sh` writes every system's legacy form first | `bash check_input_forms/run.sh` |
| `bath_data/` | not a check: the frozen baths, DFT tensors and state files the checks above and `cluster_verification/check_pcce/check_baths.py` read (its README says who reads what) | — |

## gpu_verification/ — the CUDA backend

`check_gpu/` holds the checks: 15 `test_*.cpp` unit tests (`make numeric-test`
runs three on a CPU - the third is `test_small_cluster_cuda` without the card, `test_small_cluster_host`, `make cpu-test` four more on a CPU and on 1 and 4 ranks, `make gpu-test`
the numeric two and the CUDA ones inside a GPU job, among them `test_small_cluster_cuda`, the packed kernel against
CPU-dumped clusters (A26), and `test_warp_cluster_cuda`, the warp kernel likewise for every compiled dimension (A40); CMake has the same targets; `check_gpu/README.md`
lists them), and ten `check_*.py` that run CPU and CUDA builds on the same
inputs and compare (pulses, gCCE4 accuracy and batching, pCCE, seeds and sampled
states, options, host leaks, the packed batch on a host stand-in for the card). Most need a CUDA build; five run on CPU
builds (`check_gpu/README.md` says which and with what; `check_packed_emulation.py` builds its own binary from the CPU
build's objects).
`gpu_benchmarks/` times CPU against CUDA, `server/` holds the three standing Slurm
jobs for idm4-gpu and two measurement jobs (`min_dim`, `gpu_sharing`); their last runs and the commits they ran on are in `server/README.md` (it is the record; a
run on an older commit does not check a later change). The dated `validation-*.json` under `check_gpu/` are the recorded
results (internal notes) cites.

## cluster_verification/ — the clustering methods

`docs/methods/clustering-terminology.md` names the two methods this branch added.

| Folder | Method | What is here | Run |
|---|---|---|---|
| `check_adaptive/` | **adaptive refinement** (적응형): the partition is rebuilt from measured error by the merge policy in `src/clustering/algo_adaptive/` (the only one since 2026-10-10). The driver is `ccex-adaptive`; the input block is `"adaptive"` (`"Adaptive"` until 2026-10-13, `"Refine"` until 2026-09-23) | `cpp/` unit tests, `regression/` 33 recorded runs with their golden answers, `reference/` JSON and number-text cases, fuzz and fault-injection scripts | `sh check_adaptive/verify.sh unit\|regress\|inprocess\|sanitize\|fuzz` |
| `check_bundle/` | **bundle clustering** (가변형): sites tied into bundles once, by coupling strength (`cluster_mode=bundle`, `src/clustering/bundle/`) | the records of the 2026-09-12..16 experiments and five historical `check_*.py` that no longer run with the current build (its README); what a change here has to pass is the three bundle unit tests of `check_pcce/` and the recorded runs of `check_adaptive/` | `make cpu-test` (CMake: `--target cpu-test`) runs the three unit tests since 2026-10-02 (D31), then `sh check_adaptive/verify.sh regress\|inprocess` |
| `check_pcce/` | **pCCE**: partition cache, cycle skipping and search, reproducibility across ranks, the PL6 product-convergence diagnostics | eight `check_*.py`, seven `test_*.cpp` (four built by `make cpu-test`, three by the `check_*.py` through `unit_build.py`), three benchmarks, and (still) the sweep tools and records of the PL6 work | `python3 check_pcce/check_reproducibility.py --cpu bin/main.out --output DIR` and the others |

`check_adaptive/verify.sh` is what a change to `src/clustering/algo_adaptive/` has to pass;
`code_verification/run.sh` is what any change to the core has to pass.
