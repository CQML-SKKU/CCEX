# src/

The code is split into five core directories, one per stage of a run. `include/` has the same layout, and each header
sits in the directory of the same name as its source. `main.cpp` is the entry point and runs the stages in order: read
the input, build the physical objects, make the clusters, calculate, write the results.

## Core directories

| Directory | Contents |
|---|---|
| `input/` | The input reader (`option.cpp`, `json.cpp`, `reader.cpp`) and the settings (`general.cpp`, the `Config` struct). The bath, gyro and tensor file readers. The built-in isotope table (`isotope_table.cpp`, generated). The template that `main.out -e` writes (`example_input.cpp`). |
| `physical/` | The qubits, bath spins, defects and pulse sequence, with the DEER bath pulses merged into it. |
| `clustering/` | Standard CCE clusters (`cluster.cpp`, pairs in `cluster_pairs.cpp`), pCCE partitions and their cache, and the dsj/itb variants. Bundles and the adaptive refinement are in the subdirectories below. |
| `calculation/` | The Hamiltonians, the time evolution and the coherence of every method: cCCE (`simulator_cce*.cpp`, with batched small clusters and column evolution for large ones), gCCE (`simulator_gcce*.cpp`), meCCE and megCCE (`simulator_mecce.cpp`, `simulator_megcce.cpp`, `superoperator.cpp`), noise (`noise.cpp`) and pseudospin (`pseudospin.cpp`). The memory estimate made before a run (`memory_estimate.cpp`). |
| `common/` | Maths and string helpers, memory, output files, and the MKL calls (`blas_option.cpp`, empty without MKL). |

`calculation/cce_batch_kernel.cpp` and `calculation/small_eigensolver.cpp` are compiled with `-ffp-contract=off` in
both builds, so that their results do not depend on where the compiler places fused multiply-adds.

## Optional features

| Directory | Contents | Used when |
|---|---|---|
| `clustering/bundle/` | Bundle clustering: spins tied into bundles, clusters made of bundles, and the expansion over them. | `cluster.mode` `bundle` |
| `clustering/algo_adaptive/` | The adaptive refinement: run a calculation, measure its error, rebuild the bundles, repeat. It is compiled with its own floating-point flags into the same binary. | a `cluster.adaptive` block (with `cluster.mode` `bundle`) |
| `calculation/gpu/` | The compute backend and the CUDA kernels, built into `main_gpu.out`. The files are `gpu_backend.cu`, the packed kernel for small cCCE clusters (`small_cluster_gpu.cu`), the warp kernel for clusters up to dimension 32 (`warp_cluster_gpu.cu`) and method noise's mode sums (`noise_gpu.cu`). | `execution.backend` and `execution.gpu_options` |
| `auto_convergence/` | `rdip_ladder`: one order-2 run that also gives the curve of every smaller rdip, for the automatic convergence (whose driver is `example/tools/AutoConvergence/`). | set by the automatic convergence |

### Inside clustering/algo_adaptive/

| File | Contents |
|---|---|
| `adaptive_main.cpp` | The entry point `adaptive_main()`, which reads the options and runs the merge policy. |
| `options.cpp`, `options_table.h` | The options. The table is the one place an option is declared. |
| `policy_merge` | The merge policy and its `Driver`. |
| `policy_regroup` | The `regroup` option. Where the expansion gives \|L\| > 1, the bundles around the responsible pair are regrouped. |
| `base_bundles`, `base_run`, `base_json`, `base_numeric` | What the policy shares: bundle plans and their cost, running CCEX once, JSON, exact and checked arithmetic. |

A change here has to pass `test/cluster_verification/check_adaptive/verify.sh` (`unit`, `regress`, `inprocess`,
`sanitize`). When a result is meant to change, review the difference and then update the recorded answers with
`regress.py WORKDIR --update` in the same commit.

## Adding a file

Put a new source file in the directory of its group (`src/<group>/`) and its header in `include/<group>/`. Both builds
pick up every `*.cpp` there. Keep `src/main.cpp` the only file at the top. `do_compile.sh` would also compile a new
`src/*.cpp`, but CMake would not.

A new directory must be named in both builds: in `CORE_DIRS` of the Makefile template in `do_compile.sh` (with a
compile rule for `obj/` and `obj_gpu/` beside the others) and in `CCEX_CORE_DIRS` in `CMakeLists.txt`. It also has to
be added to `INPUT` in `docs/doxygen/Doxyfile`, which does not recurse.

`src/clustering/algo_adaptive/` is the exception. It is compiled on its own with its own floating-point flags in both
builds, and its standalone driver and unit tests (`test/cluster_verification/check_adaptive/cpp/Makefile`) list their
sources by hand in `SHARED`, so a new file there goes into that list too.
