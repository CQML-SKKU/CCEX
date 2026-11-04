# Getting started
<a id="GettingStarted"></a>

This page takes you through building CCEX, a first run and its log. The command-line options, exit codes, environment
variables and the order in which a run reads its input follow as reference.

<a id="Building"></a>
## Building

Build with one of two commands from the repository root:

```
bash do_compile.sh <preset>                     # wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac, or generic
cmake -B build && cmake --build build -j 8      # any machine: finds the MPI wrapper and, if there is one, nvcc
```

**What is built.**

- `do_compile.sh` writes a `Makefile` at the top, `obj/` and `bin/main.out`. The `Makefile` is generated: change
  `do_compile.sh`, never the `Makefile`.
- Where there is CUDA, it also writes `bin/main_gpu.out`, for `execution.backend` `cuda` or `auto`.
- CMake builds the same two binaries in `build/bin/` and keeps its objects in `build/`. `-DCCEX_BIN_DIR=$PWD/bin` puts
  the binaries in `bin/` instead.
- `bin/`, `obj/`, `build/` and the `Makefile` are build output and are not in git.

You need a C++17 compiler behind an MPI wrapper (`mpicxx`, `mpic++`, `mpiicpc`). The GPU binary also needs the CUDA
toolkit: nvcc, cuSOLVER, cuBLAS and cudart ([GPU execution](../methods/gpu.md)). Eigen, uthash and SQLite are bundled in
`external/`.

<a id="Prerequisites"></a>
**Prerequisites.**

| What | For | Notes |
|---|---|---|
| a C++17 compiler behind an MPI wrapper | `bin/main.out` | the recorded checks were run with g++ 13, icpx 2023.2 and icpx 2026.1, with Intel MPI and Open MPI 4.1.5 |
| `make` or CMake 3.18 or newer | the build | `do_compile.sh` writes a Makefile; `CMakeLists.txt` asks for CMake 3.18 |
| the CUDA toolkit (nvcc, cuSOLVER, cuBLAS) | `bin/main_gpu.out` only | compute capability 12.0 and above needs CUDA 13.2 or newer: 12.8 miscompiles a kernel for the RTX PRO 6000 (see `CHANGELOG.md`), and `do_compile.sh` and CMake refuse an older nvcc there. Below 12.0 any version works (`CUDA_ARCH`, CMake `-DCMAKE_CUDA_ARCHITECTURES`; `make gpu-test` checks it). The `gpu` preset uses 13.2 |
| Intel MKL | optional | used when the build finds it (below) |
| Python 3 with numpy and matplotlib | the tutorials' `run.sh` (`example/tutorial/plot_coherence.py`) | |
| scipy as well | the analyzer (`example/analysis/`) | |
| Python 3.10 or newer with numpy, scipy, matplotlib | the automatic convergence (`example/tools/AutoConvergence`) | `CCEX_PYTHON` names the interpreter |
| psutil (BathGenerator); pandas and psutil (SpinTensorReader) | the input tools ([Tools](Tools.md)) | |
| everything in `test/requirements.txt` | the checks under `test/` | |

Eigen, uthash and SQLite come with the code (`external/`). Nothing else has to be installed for `bin/main.out`.

<a id="Installing"></a>
**Installing on a new machine.**

- **Ubuntu (or WSL) without Intel oneAPI.** Install a compiler, Open MPI and CMake from the distribution. Then build
  with the `generic` preset (or CMake), which takes `mpicxx` from `PATH`:

  ```
  sudo apt install build-essential cmake openmpi-bin libopenmpi-dev python3-numpy python3-matplotlib python3-scipy
  bash do_compile.sh generic                       # or: cmake -B build && cmake --build build -j 8
  ```

  Without `$MKLROOT` the build prints `MKL off` and runs without MKL. With oneAPI (Intel MPI, `setvars.sh`) use
  `wsl-gnu` or `wsl-icpx` instead.
- **macOS.** Install an MPI wrapper from Homebrew (`brew install open-mpi`, or `mpich`). Then run
  `bash do_compile.sh mac`, or use CMake (`brew install cmake`).
  - The `mac` preset uses `mpicxx`, or `CXX`. On Apple Silicon it compiles with Apple clang at `-O3 -mcpu=native` and
    uses Accelerate for the large clusters.
  - The checks under `test/` run there too, with the system's bash 3.2; where `timeout` is missing they use perl.
  - `check_numerics` compares a build with references recorded on WSL (x86). A Mac differs from them in the last digits
    (Apple clang on arm64, where `long double` is `double`; below 1e-8 in noDiv). Record a reference of its own on the
    Mac (`make_reference.sh current`) and compare later Mac builds against it. Keep that reference out of commits,
    since the newest `head-*` folder decides for everyone.
  - Keep the checkout out of the folders iCloud Drive syncs (Desktop and Documents, when it keeps them). iCloud copies
    files the build rewrites to `obj/<name> 2.d`, which make reads as two targets, and `make unit-env` then fails (and
    so do the checks that compile unit tests with it). A CMake build folder elsewhere (`cmake -B ~/build/ccex`) also
    avoids it.
- **A cluster with environment modules.** Load a compiler and an MPI (`module load ...`). Then run
  `bash do_compile.sh generic`, with `CXX` set to the wrapper if it is not `mpicxx`
  (`CXX=mpiicpc bash do_compile.sh generic`). You can also add a preset for the machine to `do_compile.sh`, as the
  `idm4`, `gpu` and `nurion` presets do (each loads its modules itself).
- **A GPU.** `CUDA_HOME=/usr/local/cuda CUDA_ARCH=<compute capability without the dot> bash do_compile.sh generic`
  builds `bin/main_gpu.out` beside `bin/main.out` ([GPU execution](../methods/gpu.md)).

**Presets.** Each preset loads the modules of its machine itself. The header of `do_compile.sh` describes each one.

| preset | machine |
|---|---|
| `wsl-gnu` (old name `wsl`) | WSL/Ubuntu with Intel oneAPI: Intel MPI's `mpicxx` with g++, `I_MPI_FABRICS=shm` |
| `wsl-icpx` | the same with icpx behind Intel MPI (`apt install intel-oneapi-compiler-dpcpp-cpp`) |
| `idm4` | idm: module intel-23.2, `mpiicpc` with icpx |
| `gpu` | idm4-gpu: module `intel/oneapi-2023.2` (`mpiicpc` with icpx, MKL) + CUDA 13.2; also `bin/main_gpu.out`, compute capability 12.0 by default |
| `nurion` | KISTI Nurion: intel/19.0.5 + impi/19.0.5 |
| `mac` | macOS with an MPI wrapper (Homebrew open-mpi or mpich): Apple clang `-O3 -mcpu=native` (Apple Silicon), Apple's Accelerate for the large clusters |
| `generic` | any machine with an MPI C++ wrapper on `PATH`; `CXX` (default `mpicxx`), `CC`, `CXXFLAGS`, and for the GPU binary `CUDA_HOME` (or nvcc on `PATH`) and `CUDA_ARCH` from the environment |

`do_compile.sh` also reads these variables:

| variable | effect |
|---|---|
| `CCEX_JOBS` | parallel compile jobs, and nvcc `-split-compile` threads per CUDA file (default 5) |
| `CCEX_BLAS` | MKL or Accelerate on or off (below) |
| `CXXFLAGS_EXTRA`, `LIBRARY_EXTRA` | extra compiler flags and libraries |
| `CCEX_SANITIZE=1` | address and undefined-behaviour checks |

**CMake options.**

| option | effect |
|---|---|
| `-DCCEX_CUDA=OFF` | no CUDA binary |
| `-DCMAKE_CUDA_ARCHITECTURES=<cc>` | the compute capability |
| `-DMPI_CXX_COMPILER=<wrapper>` | the MPI wrapper |
| `-DCCEX_NATIVE=OFF` | no `-march=native` (`-mcpu=native` on Apple Silicon) |
| `-DCCEX_SANITIZE=ON` | address and undefined-behaviour checks |
| `-DCCEX_BLAS=OFF\|MKL\|ACCELERATE` | MKL or Accelerate (below) |
| `-DCCEX_BIN_DIR=<dir>` | where the binaries go |
| `-DCCEX_STRICT_WARNINGS=ON` | warnings on, with the format, sequence-point and misleading-indentation ones as errors; the CI builds with it |

**Intel MKL.** The build links MKL when it finds it: `$MKLROOT` set by setvars.sh or the intel module, with its
libraries there. The build prints `[build] MKL on (...)` or `MKL off (why)`, and every run log's header says which.

- MKL takes the large cCCE/gCCE clusters. zgemm is used in the recursive path of the cCCE ensemble Hahn echo from bath
  dimension 21 on. zheevd is used from dimension 63 (the bath dimension for cCCE, the whole qubit × bath dimension for
  gCCE).
- On idm it ran h-BN in 0.72–0.88 of the time and left P1 unchanged. The curves move in the last digits (up to 9e-8 at
  late times of h-BN gCCE).
- `CCEX_BLAS=off bash do_compile.sh <preset>` (CMake `-DCCEX_BLAS=OFF`) builds without it.
- `CCEX_BLAS=mkl` (`-DCCEX_BLAS=MKL`) stops when it is missing.

**On macOS, Apple's Accelerate** takes MKL's place, since there is no MKL for Apple Silicon. It comes with the system
and has the same zgemm and zheevd, so `do_compile.sh mac` and CMake link it by default (`[build] Accelerate on`). The
run log then says `MKL for large clusters : on (Apple Accelerate in its place; ...)`, and `Accelerate on` in its `Build`
line. `CCEX_BLAS=off` (`-DCCEX_BLAS=OFF`) leaves it out. `CCEX_BLAS=accelerate` (`-DCCEX_BLAS=ACCELERATE`) insists on
it.

**The binary of a CMake build** is in `build/bin/`, so point the tutorials and the checks at it:

- tutorials: `export CCEX_BIN=<repository>/build/bin/main.out`, and `CCEX_GPU_BIN` for `advanced/11_GPU`;
- checks: `CCEX_BIN=$PWD/build/bin/main.out CCEX_BUILD_DIR=$PWD/build`. The unit-test drivers take the compiler, flags
  and objects from the `unit-build.env` the build writes (any generator), through `test/unit_build_env.sh`.

<a id="CheckBuild"></a>
**Checking the build.**

- `make numeric-test cpu-test` (CMake: `cmake --build build --target numeric-test cpu-test`) runs the unit tests on the
  build tree.
- `bash test/code_verification/run.sh` compares results with recorded references. `test/README.md` lists every check
  and what it needs.
- The quickest check of a working install is the [first calculation](#FirstCalculation) below.

<hr>

<a id="FirstCalculation"></a>
## First calculation

This example computes the Hahn-echo coherence of an NV centre in natural-abundance ¹³C with cCCE to order 2. It takes
about a second on four ranks ([basic/01_first_run](../../example/tutorial/basic/01_first_run/README.md)).

```
bash do_compile.sh wsl-gnu                               # 1. build bin/main.out (or another preset, or CMake)
cd example/tutorial/basic/01_first_run
NP=4 bash run.sh                                         # 2. run: results/ccein_wiDiv, ccein_noDiv, ccein.log, coherence.png
```

1. `run.sh` runs `mpirun -n $NP bin/main.out -f ccein.jsonc -o results > results/ccein.log`. `CCEX_BIN` names another
   binary, e.g. `build/bin/main.out` of a CMake build. It then runs `example/tutorial/plot_coherence.py`, which needs
   numpy and matplotlib.
2. `results/` then holds:
   - `ccein_wiDiv`, the coherence ([Output files](OutputFiles.md));
   - `ccein_noDiv`, for the analyzer's error correction (it is not a coherence);
   - `ccein.log` and `coherence.png`.
3. In `results/ccein.log` ([The run log](#RunLog)) look at:
   - the `Build` line, which says which binary ran;
   - `Read BathFile[0]` and `Number of bath spins`;
   - `Clusters per order` and `Memory estimate`;
   - the last line, `JOB DONE`.

   A run that stops prints `Error: ...` on standard error instead. Here it goes to the terminal, since only standard
   output is redirected.
4. `run.sh` prints T2, the last time $|L| \ge 1/e$. Compare `results/coherence.png` with `expected.png` in the same
   folder: T2 = 1.068 ms.

**Your own input.**

1. `bin/main.out -e` writes `ccex_example.jsonc`. It holds the setting of this calculation, and every other key as a
   commented line with its allowed values and its default.
2. Set `configuration.bathfile` to a bath file. The template does not run without one.
   `example/tutorial/basic/01_first_run/bath/bath_Diamond_NV_natab_1` is one; [Tools](Tools.md) makes others.
3. Change the keys you need ([Input file](InputFile.md)).
4. Run `mpirun -n 4 bin/main.out -f ccex_example.jsonc`.

The [Examples](Examples.md) go on from here in learning order.

<hr>

<a id="RunLog"></a>
## The run log

Rank 0 writes the log to standard output. Refusals and errors (`Error: ...`) and the `[seed]` line go to standard
error. Without `-v` the log is short, in this order:

1. **The banner**: the start time, `MKL for large clusters : on (...)` or `off`, and the `Build` line, which says what
   binary ran:

   ```
   Build : CPU | g++ 13.3.0 | Intel(R) MPI Library 2021.18.0 for Linux* OS | MKL on | preset wsl-gnu | commit <hash>[-dirty]
   ```

   Its fields are `CPU` or `GPU (<card>)`, the compiler, the MPI library, MKL on or off, the build preset (`cmake` for
   a CMake build) and the commit it was built from. `-dirty` means the source differed from that commit; outside git
   the commit is `unknown`.
2. **`Warning: input key "..." is not read by this run and is ignored`**, once per key the run does not use (a key of
   another method, for example). A key that looks like a misspelling of a key the run reads stops the run instead
   (`Error: input key "..." is not read, did you mean "..."?`).
3. **The settings** (`Report all structures`): `< Structure Config >`, `< Structure Pulse >` and
   `< Structure Output >`.
   - Each setting is one `key : value` row under the input key's name (`cluster.mode`, `configuration.s_range`,
     `output.file`, ...). `< Structure Pulse >` also shows the pulse sequence as a table of start, end, axis and angle.
   - A `chosen by` row under a key says where a value the input did not give came from.
   - `configuration.seed` shows the seed when the run draws random numbers (sampled bath states, Defect axes and
     on-site states drawn without their files, the pCCE partition), and `(not used)` otherwise.
   - A seed the run generated itself is also printed as a `[seed]` line on standard error. Write that value into the
     input to replay the run.
4. **The files read**: `Read QubitFile`, `Read GyroFile` (or the built-in isotope table), `Read BathFile[i]` with
   `Number of bath spins`, and the hyperfine and quadrupole read modes.
5. **The clusters**:
   - `Clusters per order : 1: 967, 2: 3366 (total 4333, 980 bath spins, 4 ranks)` gives the clusters of each order
     over all ranks.
   - `Memory estimate : ... GB per rank (...), ... GB available per rank (...)` is the memory check
     ([the memory check](InputFile_Backend.md#memory_check)). A run that does not fit is refused right after this line.
6. **The calculation**:
   - `psia` / `psib` and which m<sub>s</sub> state they mostly are;
   - one progress line every 10 % of rank 0's clusters (`At rank 0, 109 of 1084 clusters (10%) are computed at 0.05 s`);
   - per state, `State done : ensemble, 4333 clusters, 0.16 s` (or `State done : 3 of 20, ...` for sampled states).
7. **`JOB DONE : Wall time = ... s`**, the last line of a run that finished.

With [configuration.k_range](InputFile_General.md#tag_configuration_k_range) each configuration starts with a
`Configuration k = <k> (configuration.k_range [first, last])` header. From the second configuration on, the log leaves
out the settings and reports that only repeat the first, up to the clusters. It keeps the files read and their spin
counts, the tagged lines and every warning, note and error.

**`-v`** prints everything:

- each block of the input as it is read;
- every structure: the qubits, the bath spins with their hyperfine and quadrupole tensors in rad/ms, the clusters of
  every order with their powers;
- the progress and wall time per cluster order;
- the steps of each state (product, gather, save) and their wall times;
- with `k_range`, one `[configuration]` line per rank with its memory after each configuration.

It is long: thousands of lines for a bath of a thousand spins.

<a id="RunLogTags"></a>
**Tagged lines.** A line that starts with a tag in brackets reports a decision or a measurement. These lines are kept
in the short log too.

| Tag | When | What it says |
|---|---|---|
| `[backend]` | `execution.backend` `auto` or `cuda` | which device each rank uses, or why the CPU is used (`CPU-only build`, a method or propagator the GPU does not take, CUDA initialisation failed) |
| `[gpu-packed]`, `[gpu-warp]` | GPU build, `auto` or `cuda` | whether the packed kernel (cCCE clusters of bath dimension 2 and 4) and the warp kernel (small cCCE and gCCE clusters, one per warp) are on, and why not |
| `[gpu-memory-fallback]` | `auto` | the GPU ran out of memory: this cluster and every later one of that dimension or more go to the CPU on that rank |
| `[gpu-batch]`, `[gpu-memory]` | GPU, `execution.profile` | the cluster batches per dimension and the GPU memory cache |
| `[cpu-batch]` | cCCE | the batched kernel of small clusters is off, and why (it checks itself at start against the one-by-one path) |
| `[dm]` | `quantity` dm | clusters put values into density-matrix elements the free evolution keeps at zero ([The density matrix](OutputFiles.md#out_dm)) |
| `[seed]` | a generated seed (standard error) | the seed to write into the input to replay the run or a configuration |
| `[partition]`, `[partition-search]`, `[partition-iteration]`, `[partition-detail]` | pCCE | the k-means partition search: mode, trials, the best trial, iterations |
| `[partition-cache]` | `cluster.pcce.cache` | `hit`, `miss` (with the reason) or `stored`, and the file |
| `[adaptive]`, `[adaptive-bundle]`, `[adaptive-link]`, `[adaptive-connections]`, `[adaptive-check]`, `[adaptive-expansion]` | `cluster.mode` bundle | how the sites were bundled, the connections kept and left out, the expansion |
| `[probe-budget]` | the probes of the adaptive refinement | the memory for the probes running at once, and the time a rank waited for it |
| `[cluster-contributions]` | `output.cluster_contributions` | the rows of the table and the time it took ([Cluster contributions](OutputFiles.md#out_contributions)) |
| `[cluster-plan]` | `output.cluster_only` | the sites, the clusters and the path of `.clusters.json` |
| `[rdip-ladder]` | `rdip_ladder` | the extra rdip values computed from the same run |
| `[profile]`, `[phase]`, `[small-eig]` | `execution.profile` | timings per method, backend and cluster dimension; the state loop split into phases; the small eigensystems solved by the Jacobi solver and those left to Eigen ([execution.profile](InputFile_Backend.md#tag_profile)) |
| `[sparse]` | `execution.propagator` sparse, `execution.profile` | the Krylov propagation per cluster: dimension, non-zeros, matrix-vector products |
| `[numerical-invalid]` | a cluster term is not finite, or a zero it would be divided by (standard error) | the first such value on a rank: the state, time step, cluster order, power and the raw value |
| `[auto]` | an `auto_convergence` block | the runner, the Python and the binary used |

<hr>

<a id="HowToRun"></a>
## How to run

Run the binary under MPI with the input file:

```
mpirun -n 4 bin/main.out -f ccein.json
```

- `-n 4`: the number of MPI processes.
- `bin/main.out`: the binary (`build/bin/main.out` from CMake, `bin/main_gpu.out` for the CUDA backend).
- `ccein.json`: the [input file](InputFile.md), JSON with `//`, `#` or `!` comments (`.json` or `.jsonc`).

**Command-line options**

| option | meaning |
|---|---|
| `-f FILE` | the [input file](InputFile.md); required for a calculation |
| `-o DIR` | the output folder, made when it is missing (below) |
| `-v` | verbose output (a longer log) |
| `-e` | write the input template `ccex_example.jsonc` and stop (below) |
| `-h`, `--help` | print a short usage and exit |

`-o DIR` puts the outputs in `DIR/<output.file>`. When the input has no [output.file](InputFile_Output.md#tag_outfile),
they are `DIR/<input file name without .json or .jsonc>`. An absolute `output.file` together with `-o` is refused.
Without `-o` the outputs go to [output.path](InputFile_Output.md#tag_output_path), else to the working directory. The
order is `-o DIR` > `output.path` > the working directory.

`-e` writes `ccex_example.jsonc` in the working directory (in `DIR` with `-o DIR`) and stops.

- Its uncommented part is the setting of `example/tutorial/basic/01_first_run`: NV-<sup>13</sup>C cCCE2 Hahn echo,
  ensemble, rbath 50 Å, 2 ms.
- It has a `bathfile` for you to fill in. The template does not run until it names a bath file; a ready one is
  `example/tutorial/basic/01_first_run/bath/bath_Diamond_NV_natab_1`.
- `gyrofile` is commented, so the run uses the built-in isotope table.
- Every other input key is a commented line with its allowed values and its default. The comments are `//`.
- An existing file is not overwritten. `-e` with `-f` is refused.

Every other setting of a run is a key of the input file. Older versions took the options below. They are now refused,
with the key that replaces them:

| removed option | input key |
|---|---|
| `-I FILE` | [configuration.bathfile](InputFile_General.md#tag_bathfile): `["FILE"]` |
| `-N N` | [configuration.s_range](InputFile_General.md#tag_s_range) `[1, N]` |
| `-m METHOD` | [method](InputFile_General.md#tag_method) (with [cluster.mode](InputFile_Cluster.md#tag_cluster_mode)) |
| `-q QUANTITY` | [quantity](InputFile_General.md#tag_quantity) |
| `-B BZ` | [bfield](InputFile_General.md#tag_bfield) (a number is Bz) |
| `-s FILE` | [configuration.bath_state_file](InputFile_General.md#tag_bath_state_file): `["FILE"]` |
| `-S FILE` | [configuration.defect_state_file](InputFile_General.md#tag_defect_state_file): `["FILE"]` |
| `-a FILE` | [configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file): `["FILE"]` |
| `-P` | [output.cluster_only](InputFile_Output.md#tag_cluster_only): `true` |
| `-o FILE` (the output name) | [output.file](InputFile_Output.md#tag_outfile); `-o` names the output folder |
| `--auto_conv` | an [auto_convergence](InputFile_AutoConvergence.md) block (`main.out -f FILE` starts it) |
| `--auto_status` | `PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status FILE` (the tool's own command) |
| `--adaptive ...` | a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block |

Two blocks of the input file replace the ordinary calculation with a run mode of their own:

- An `auto_convergence` block runs the automatic convergence
  ([Automatic convergence tags](InputFile_AutoConvergence.md)). It takes no option but `-f FILE`.
- An `adaptive` block runs the adaptive refinement ([cluster.adaptive](InputFile_Cluster.md#tag_adaptive)). `-o` does
  not apply to it; its steps go to the block's `output`.

<a id="ExitCodes"></a>
**Exit codes.** `main.out` exits with these codes:

| exit code | ordinary run |
|---|---|
| 0 | the calculation finished |
| 1 | the input is refused |
| nonzero, from `mpirun` | an error during the calculation |

- A refusal is decided on every rank, so all ranks stop together, before anything is calculated. Rank 0 prints one
  `Error: ...` line on standard error.
- An error during the calculation (a file that cannot be written, for example) aborts every rank through MPI. `mpirun`
  then ends with a nonzero code of its own.
- The adaptive refinement exits with 0 when its final status is a converged or physical one, and with 2 when it stopped
  at another status (its last step is then the result). It exits with 1 on an error
  ([Adaptive refinement](OutputFiles.md#out_adaptive)).
- The automatic convergence exits with 0 converged, 3 unresolved, 2 refused, 4 stopped, 1 error
  ([its results](InputFile_AutoConvergence.md#tag_ac_results)).

<a id="EnvironmentVariables"></a>
**Environment variables**

An ordinary run needs none of them. They switch an optimisation off for a comparison, point the automatic convergence
at its tool, or serve the checks. Set them in front of the command
(`CCEX_SMALL_EIG=0 mpirun -n 32 main.out -f ccein.json`) or export them in a job script.

| variable | when it is read | effect |
|---|---|---|
| `CCEX_CPU_BATCH=0` | every run | turns off the CPU batch of small cCCE clusters (bath dimension 2 and 4): each cluster is computed alone. The batch gives the same bits as that path. Where it would not (another compiler or its options), it checks itself at start and stays off by itself, with a line in the log |
| `CCEX_SMALL_EIG=0` | every run | turns off the fixed-sweep Jacobi eigensolver of the same small clusters: Eigen solves them as well (the reference for a comparison). Matrices outside a safe range or not converged already go to Eigen one by one |
| `CUDA_VISIBLE_DEVICES` | GPU build (`main_gpu.out`) | which card the run may use. Empty, or not matching any card, hides them all, and `execution.backend` `auto` then stays on the CPU |
| `CCEX_PYTHON` | a run with an `auto_convergence` block | the Python that runs the automatic-convergence tool (default `python3`; it needs numpy, scipy, matplotlib) |
| `CCEX_AUTOCONV_DIR` | a run with an `auto_convergence` block | the tool's folder when it is not found next to the binary (`example/tools/AutoConvergence`). A folder others can write to is refused |
| `CCEX_CHECK_DISORDER=1` | every run with `s_range` (sampled states) | for the checks only: also builds the static mean field of the bath on one rank the old way, and stops on any bit that differs from the split one |
| `CCEX_SMALL_CLUSTER_DUMP=FILE`, `CCEX_SMALL_CLUSTER_DUMP_N=N` | cCCE runs | for the checks only: writes the first N (400) small clusters with the CPU result to FILE (the CPU batch is off meanwhile) |
| `CCEX_WARP_CLUSTER_DUMP=FILE`, `CCEX_WARP_CLUSTER_DUMP_N=N` | gCCE runs | for the checks only: the same for the clusters of the GPU warp kernel |
| `CCEX_TEST_GPU_OOM_BYTES=N` | GPU build | for the checks only: every GPU allocation from N bytes fails as if the card were full, to exercise the fall-back to the CPU |

The variables `do_compile.sh` reads at build time are in [Building](#Building).

<hr>
<a id="Workflow"></a>
## Workflow

A run goes through these steps:

1. The command line: the input file (`-f`), the output folder (`-o`) and the log length (`-v`).
2. Every key of the input file, grouped as in [Input file](InputFile.md). A likely typo of a key, or a combination
   that is not supported, stops the run here, before any file is read. Any other key the run does not read is
   reported in the log and ignored.
3. The files the input names: the qubit file, the gyromagnetic ratios (or the built-in isotope table), the bath files,
   hyperfine and quadrupole tensor files, and the bath states and defect axes when they are given.
4. For each bath configuration: the clusters of the bath (`cluster.mode`, `order`), their propagation (CPU or
   `execution.backend` `cuda`) and the product of the cluster contributions.
5. The output files (`output.path`, `output.file`): the coherence or the density matrix as a function of time
   ([Output files](OutputFiles.md)).

<img src="CCEX_workflow.png" alt="CCEX workflow" width="600">

<hr>
<a id="Troubleshooting"></a>
## Troubleshooting

- **Another `mpirun` comes first on `PATH`** (Anaconda's, for example, when its `bin/` is on `PATH`). An `mpirun` of
  another MPI than the one the binary was built with starts the processes without joining them. Every process then
  runs the whole calculation as a run of its own (`Clusters per order` ends with `1 rank`). Check `which mpirun`
  against the MPI in the `Build` line, and put the right one first (or call it by its full path).
- **`Error: not enough memory: ...`** before the first cluster. The largest cluster does not fit what one rank may use.
  The message says what to change: fewer ranks per node, a smaller order or smaller clusters. See
  [the memory check](InputFile_Backend.md#memory_check).
- **No GPU.** With `execution.backend` `auto`, every cluster stays on the CPU when the binary is CPU-only, a card cannot
  be opened, or `CUDA_VISIBLE_DEVICES` hides all cards. The `[backend]` lines say why. With `cuda` the run is refused
  instead. The GPU binary is `bin/main_gpu.out` ([GPU execution](../methods/gpu.md)).
- **`Error: input key "..." is not read, did you mean "..."?`** is a misspelled key; the run stops before reading any
  file. A `Warning: input key ... is not read by this run` is only a key this run does not use.
- **A key of an older version** is refused with the name that replaced it ([old names](InputFile.md#OldNames)).
  `python3 example/tools/InputKeys/convert_input_keys.py OLD.json` rewrites such an input.
- **A curve above 1** ($|L| > 1$ in `_wiDiv`) comes from a division by a cluster term near zero. It is not physics. See
  [The coherence files](OutputFiles.md#out_coherence) and the analyzer's error correction ([Analysis](Analysis.md)).
