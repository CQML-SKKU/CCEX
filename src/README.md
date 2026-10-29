# src/ — what is where

The core is five directories, one per stage of a run, and `include/` has the same
shape: a header lives in the directory of the same name as its source. An optional
feature gets its own directory beside them. `main.cpp` is the entry point and calls
the groups below in this order: **read input → build the physical objects → make
clusters → calculate → write**.

## Core (`src/<group>/*.cpp`, `include/<group>/*.h`)

| Group | Files | Does |
|---|---|---|
| Input (`input/`) | `option.cpp`, `reader.cpp`, `json.cpp`, `general.cpp`, `isotope_table.cpp`, `example_input.cpp`, `mpi_seed.h` | parse the input JSON, read the bath / gyro / tensor files, hold the settings (`Config`); the built-in isotope table for species without a gyrofile entry (generated); the example input `main.out -e` writes |
| Physical objects (`physical/`) | `qubit.cpp`, `bath.cpp`, `defect.cpp`, `pulse.cpp`, `pulse_bath.cpp` | qubit, bath spins, defects, pulse sequence; the bath pulses of DEER (`pulse.bath`) merged into it |
| Clustering (`clustering/`) | `cluster.cpp`, `cluster_hash.cpp`, `cluster_pairs.cpp`, `cluster_pcce.cpp`, `cluster_dsjitb.cpp`, `pcce_partition_cache.cpp`, `pcce_cycle.h` | standard CCE clusters (order 2 without the hash: `cluster_pairs.cpp`), pCCE, the dsj/itb variants, and the optional on-disk cache of a pCCE partition; bundle clustering and the adaptive refinement (`cluster_mode` = `bundle`, the `adaptive` block) are in the two subfolders `clustering/bundle/` and `clustering/algo_adaptive/` below |
| Calculation (`calculation/`) | `simulator.cpp`, `simulator_cce.cpp`, `simulator_cce_recursive.cpp`, `simulator_cce_columns.cpp`, `cce_batch_window.cpp`, `cce_batch_kernel.cpp`, `small_eigensolver.cpp`, `cce_setup_small.cpp`, `simulator_gcce.cpp`, `simulator_gcce_vector.cpp`, `simulator_gcce_recursive.cpp`, `simulator_gcce_bath.cpp`, `simulator_mecce.cpp`, `simulator_megcce.cpp`, `pseudospin.cpp`, `superoperator.cpp`, `simulator_hamiltonian.cpp`, `small_hamiltonian.cpp`, `hamiltonian.cpp`, `memory_estimate.cpp`, `scaled_product.h` | build Hamiltonians, propagate, form the coherence; small cCCE clusters (bath dimension 2, 4) evolved several at once on the CPU (`cce_batch_kernel.cpp` is compiled with `-ffp-contract=off` in both builds, A30); their eigensystems by a fixed-sweep Jacobi, several at once (`small_eigensolver.cpp`, also `-ffp-contract=off`; `CCEX_SMALL_EIG=0` for Eigen, A42); large one-state clusters with the time points as matrix columns (`simulator_cce_columns.cpp`, A62, input block `column_evolution`, on by default from bath dimension 256); the bath pulses on the gCCE paths (`simulator_gcce_bath.cpp`); the master-equation methods meCCE and megCCE (`simulator_mecce.cpp`, `simulator_megcce.cpp`) on the Liouville-space superoperators and jump operators of `superoperator.cpp`; method pseudospin, the pair clusters in the pseudospin model (`pseudospin.cpp`); the memory check before the calculation (`memory_estimate.cpp`, A67: the largest cluster's peak per rank against what a rank of the node may use) |
| Shared tools (`common/`) | `utilities.cpp`, `memory.cpp`, `output.cpp`, `blas_option.cpp` | maths and string helpers, allocation, result files; the MKL calls (A46; empty in a build without MKL) |

## Optional features (a directory here, and one of the same name in `include/`)

`"cluster_mode": "bundle"` (`"adaptive"` until 2026-10-13) names the whole feature; `clustering/bundle/` and
`clustering/algo_adaptive/` are its two halves (both under `clustering/` since 2026-10-08).

| Directory | Does | Turned on by |
|---|---|---|
| `clustering/bundle/` | the calculation half of adaptive mode (was `cluster_bundle/` until 2026-09-23; it is a clustering method and lives with the others): spins are tied into bundles, clusters are made of bundles, and the expansion is formed over them. Bundles come from the input, or from a fixed coupling rule applied once (`bundle_partition`, `bundle_connection_check`). `bundle_work` is the file protocol the other half talks through | `"cluster_mode": "bundle"` and the `bundle` block in the input |
| `clustering/algo_adaptive/` | the algorithm half (was `algo_adaptive/` until 2026-10-08): run a calculation, measure its error, rebuild the bundles, repeat. Compiled separately, with its own flags, into the same binary | a `cluster.adaptive` block in the input (`docs/user/InputFile_Cluster.md`), needs `cluster.mode` `bundle` or `pcce` |
| `calculation/gpu/` | choice of compute backend, per-cluster timing, and the CUDA kernels (`gpu_backend.cu`, built into `main_gpu.out` by `do_compile.sh gpu` or by CMake where it finds nvcc); the packed kernel for small cCCE clusters (`small_cluster_gpu.cu`, the per-cluster arithmetic in `include/calculation/gpu/small_cluster_math.h`, the node feeder and windows in `small_cluster_batch.cpp`, A26); the warp kernel for clusters up to dimension 32 (`warp_cluster_gpu.cu`, `include/calculation/gpu/warp_cluster.h`, the same feeder, A40; on by default in the GPU build) | `"backend"` and the `"gpu_options"` block (`packed_max_dim` for the packed kernel, `warp_max_dim` for the warp kernel) in the input |
| `auto_convergence/` | the in-run half of the automatic convergence protocol ((internal notes); the driver between runs is the Python tool `example/tools/AutoConvergence/`): `rdip_ladder` assembles, from one order-2 run of cCCE, gCCE or pCCE, the curve of every smaller rdip | `"rdip_ladder"` in the input |

## Inside `clustering/algo_adaptive/`

One policy, merge, over four shared layers. Error-gated, budgeted, auto-partition and cost-preserving
(one file each, built on merge's `Driver`) were removed on 2026-10-10; they are in the history:
`git show 9251d93db:src/clustering/algo_adaptive/policy_<name>.cpp`.

| File | Holds |
|---|---|
| `adaptive_main.cpp` | the entry, `adaptive_main()` (both called `refine_main` until 2026-09-23): read the options, refuse a removed policy, run merge. (Its six-line `main()` is only for the test driver `ccex-adaptive`; CCEX's own `main.cpp` compiles it out.) |
| `options.cpp`, `options_table.h` | the options. The table is the one place an option is declared |
| `policy_merge` | `"merge"`, the only policy, and its `Driver` |
| `policy_regroup` | method D, merge's `--regroup` option (2026-10-11, (internal notes)): where the expansion gives \|L\| > 1, regroup the bundles around the culprit pair; its judgement as pure functions, and `Driver::run_regroup` |
| `base_bundles` | what the policy uses about bundles and curves: the plan, the expansion and its cost |
| `base_run` | running CCEX once: launcher, input manifest, the probe cache |
| `base_json` | the JSON value, reader and writer, and the text of numbers |
| `base_numeric` | error types, checked 64-bit integers, exact summation, complex abs and division |

The policies were written in Python first (`adaptive_refine.py` and the other
`adaptive_*.py` a comment may name) and ported to C++ one by one. The Python was
removed once the C++ had its own recorded answers; it is in the history:
`git show 3ba4cd829:test/cluster_verification/check_adaptive/oracle/<name>`.

`test/cluster_verification/check_adaptive/verify.sh` is what a change to `clustering/algo_adaptive/` has to pass: `unit`,
`regress` (33 recorded runs against `test/cluster_verification/check_adaptive/regression/golden/`), `inprocess`
and `sanitize`. When a result is meant to change, review the difference and then
rewrite the golden files with `test/cluster_verification/check_adaptive/regress.py WORKDIR --update`.

## Adding a file

A new source file goes in the directory of its group (`src/<group>/`, its header in
`include/<group>/`); both builds pick up every `*.cpp` there. `src/main.cpp` is the
only file at the top: `do_compile.sh` would also take a new `src/*.cpp`, CMake would not,
so do not put one there. A new directory has to be named in both builds:
`CORE_DIRS` in the Makefile template of `do_compile.sh`, with a compile rule beside the
others (one per directory, for the `obj/` and the `obj_gpu/` objects), and
`CCEX_CORE_DIRS` in `CMakeLists.txt`. `src/clustering/algo_adaptive/` is the exception: it is
compiled on its own, with its own floating-point flags, in both builds, and its standalone driver and unit tests
(`test/cluster_verification/check_adaptive/cpp/Makefile`, built by `verify.sh unit`) list the sources they link
by hand in `SHARED`: a new file there goes into that list too. The API reference lists its
directories one by one as well (`INPUT` in `docs/doxygen/Doxyfile`, which does not recurse).
