# Building and running
<a id="BuildingAndRunning"></a>

This page lists the options and rules of building and running CCEX, to look up when needed. To build and run CCEX
for the first time, follow [Getting started](GettingStarted.md).

1. [Building](#Building): the presets, CMake, MKL, CUDA and macOS.
2. [Command line](#HowToRun): the options, where the outputs go, and the other run modes.
3. [Exit codes](#ExitCodes): what `main.out` returns, in each run mode.
4. [The run log](#RunLog): the short log part by part, the long log, and the tagged lines.
5. [Environment variables](#EnvironmentVariables): the variables a run reads.
6. [Troubleshooting](#Troubleshooting): common problems.

<a id="Building"></a>
## 1. Building

### Presets

`bash do_compile.sh <preset>` sets the compiler and options for one kind of machine. It writes the top-level
`Makefile`, so change `do_compile.sh` rather than the `Makefile`.

| Preset | Machine |
|---|---|
| `generic` | any machine with an MPI C++ wrapper on `PATH`. It reads `CXX` (default `mpicxx`), `CC` and `CXXFLAGS`, and for the GPU binary `CUDA_HOME` (or nvcc on `PATH`) and `CUDA_ARCH` |
| `wsl-gnu` | WSL or Ubuntu with Intel oneAPI: Intel MPI's `mpicxx` with g++ |
| `wsl-icpx` | the same with icpx behind Intel MPI |
| `mac` | macOS with an MPI wrapper from Homebrew: Apple clang with `-O3` and `-mcpu=native` (`-march=native` on an Intel Mac), and Apple's Accelerate |

The header of `do_compile.sh` lists further presets for the clusters of the developers' group; each loads its own
modules. A preset for another machine can be added in the same way.

`do_compile.sh` also reads these variables:

| Variable | Effect |
|---|---|
| `CCEX_JOBS` | the number of parallel compile jobs (default 5) |
| `CCEX_BLAS` | `auto` (default), `mkl`, `accelerate` or `off` ([below](#BuildBlas)) |
| `CXXFLAGS_EXTRA`, `LIBRARY_EXTRA` | extra compiler flags and libraries |

### CMake

`cmake -B build && cmake --build build` puts the binaries in `build/bin/`. Point the tutorials at them with
`export CCEX_BIN=$PWD/build/bin/main.out` (and `CCEX_GPU_BIN` for Advanced 11).

| Option | Effect |
|---|---|
| `-DCCEX_CUDA=OFF` | no GPU binary |
| `-DCMAKE_CUDA_ARCHITECTURES=<cc>` | the compute capability |
| `-DMPI_CXX_COMPILER=<wrapper>` | the MPI wrapper |
| `-DCCEX_NATIVE=OFF` | no `-march=native` (`-mcpu=native` on Apple Silicon) |
| `-DCCEX_BLAS=AUTO\|MKL\|ACCELERATE\|OFF` | as `CCEX_BLAS` ([below](#BuildBlas)) |
| `-DCCEX_BIN_DIR=<dir>` | where the binaries go, for example `$PWD/bin` |

<a id="BuildBlas"></a>
### MKL and Accelerate

The build links Intel MKL when `$MKLROOT` is set (by `setvars.sh` or an intel module), and uses it for the large
clusters. On macOS, Apple's Accelerate takes its place; the `mac` preset and CMake link it by default. The results
change only in the last digits. The build says which one it linked, and so does the header of every run log.
`CCEX_BLAS=off` (CMake `-DCCEX_BLAS=OFF`) builds without either.

### CUDA

The GPU binary `bin/main_gpu.out` needs nvcc, cuSOLVER, cuBLAS and cudart. Set the compute capability with `CUDA_ARCH`
(CMake: `-DCMAKE_CUDA_ARCHITECTURES`). Compute capability 12.0 or above needs CUDA 13.2 or newer:

- `do_compile.sh` checks `CUDA_ARCH` (80 by default for `generic`). With an older nvcc it stops, or, when `generic`
  found nvcc only on `PATH`, warns and builds `bin/main.out` alone.
- CMake accepts an older nvcc only when `CMAKE_CUDA_ARCHITECTURES` names capabilities below 12.0 alone. Otherwise
  `CCEX_CUDA` `AUTO` (the default) skips the GPU binary with a warning, and `ON` stops.

### macOS

Do not build in a folder that iCloud Drive syncs, such as Desktop or Documents: iCloud makes copies of build files,
and the build fails. A CMake build folder elsewhere (`cmake -B ~/build/ccex`) also avoids this.

<a id="HowToRun"></a>
## 2. Command line

A run is `mpirun -n N bin/main.out -f INPUT`, with `build/bin/main.out` for a CMake build and `bin/main_gpu.out` for
the GPU. The input is JSON with `//`, `#` or `!` comments ([Input file](InputFile.md)).

| Option | Meaning |
|---|---|
| `-f FILE` | the [input file](InputFile.md); required for a calculation |
| `-o DIR` | the output folder, made when it is missing (below) |
| `-v` | verbose output (a longer log) |
| `-e` | write the input template `ccex_example.jsonc` and stop (below) |
| `-h`, `--help` | print a short usage and exit |

**Where the outputs go.** The folder is `-o DIR` if given, else [output.path](InputFile_Output.md#tag_output_path),
else the working directory. The file name is [output.file](InputFile_Output.md#tag_outfile), or the input file's name
without `.json` or `.jsonc`. An absolute `output.file` together with `-o` is refused.

**The template.** `-e` writes `ccex_example.jsonc` in the working directory, or in `DIR` with `-o DIR`, and stops
([Getting started](GettingStarted.md#OwnInput) shows how to use it). It does not overwrite an existing file, and `-e`
together with `-f` is refused.

**Removed options.** Every other setting is a key of the input file. The options of older versions (`-I`, `-N`, `-m`,
`-q`, `-B` and others) are refused with the name of the key that replaces them.

**Other run modes.** Two blocks of the input file replace the ordinary calculation:

- An `auto_convergence` block runs the automatic convergence
  ([Automatic convergence tags](InputFile_AutoConvergence.md)). It takes no option but `-f FILE`.
- A `cluster.adaptive` block runs the adaptive refinement ([cluster.adaptive](InputFile_Cluster.md#tag_adaptive)).
  `-o` and `output.path` are refused with it; its steps go to `cluster.adaptive.workdir`.

<a id="ExitCodes"></a>
## 3. Exit codes

| Code | Ordinary run | Adaptive refinement | Automatic convergence |
|---|---|---|---|
| 0 | finished | finished with a converged or physical status | converged |
| 1 | input refused | error | error |
| 2 | | stopped at another status; its last step is the result | refused |
| 3 | | | unresolved |
| 4 | | | stopped |
| other, from `mpirun` | error during the calculation | | |

A refusal is decided on every rank before anything is calculated, so all ranks stop together, and rank 0 prints one
`Error: ...` line on standard error. With `execution.backend` `cuda`, a method the GPU does not take, a CPU-only binary
or a card that cannot be opened stops the run differently: every rank prints `Error: GPU backend: ...` and aborts, so
the code comes from `mpirun`. An error during the calculation, such as a file that cannot be written, aborts
every rank through MPI, and `mpirun` returns a nonzero code of its own. The statuses are explained in
[Adaptive refinement](OutputFiles.md#out_adaptive) and in [Automatic convergence](InputFile_AutoConvergence.md#tag_ac_results).

<a id="RunLog"></a>
## 4. The run log

Rank 0 writes the log to standard output, apart from a few tagged lines that every rank writes (`[backend] rank=...`).
Errors (`Error: ...`), the `[seed]` line and `[numerical-invalid]` go to standard error.

### The short log

Without `-v`, the log has these parts in order:

1. **The banner**: the start time, `MKL for large clusters : on (...)` or `off`, and the `Build` line, which says what
   binary ran:

   ```
   Build : CPU | g++ 13.3.0 | Intel(R) MPI Library 2021.18.0 for Linux* OS | MKL on | preset wsl-gnu | commit <hash>[-dirty]
   ```

   It names the device, the compiler, the MPI library, MKL, the preset and the commit the binary was built from.
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
4. **The files read**: `Read GyroFile` (or the built-in isotope table), `Read BathFile[i]` with
   `Number of bath spins`, and the hyperfine and quadrupole read modes.
5. **The clusters**:
   - `Clusters per order : 1: 965, 2: 3349 (total 4314, 978 bath spins, 4 ranks)` gives the clusters of each order
     over all ranks.
   - `Memory estimate : ... GB per rank (...), ... GB available per rank (...)` is the memory check
     ([the memory check](InputFile_Backend.md#memory_check)). A run that does not fit is refused right after this line.
6. **The calculation**:
   - `psia` / `psib` and which m<sub>s</sub> state they mostly are;
   - one progress line every 10 % of rank 0's clusters (`At rank 0, 108 of 1080 clusters (10%) are computed at 0.05 s`);
   - per state, `State done : ensemble, 4314 clusters, 0.17 s` (or `State done : 3 of 20, ...` for sampled states).
7. **`JOB DONE : Wall time = ... s`**, which ends a run that finished.

With [configuration.k_range](InputFile_General.md#tag_configuration_k_range), each configuration starts with a
`Configuration k = <k> (configuration.k_range [first, last])` line. From the second configuration on, the settings and
the reports that repeat the first configuration are left out. The files read, their spin counts, the tagged lines and
every warning, note and error stay.

### The long log

`-v` adds:

- each block of the input as it is read;
- every structure: the qubits, the bath spins with their hyperfine and quadrupole tensors in rad/ms, the clusters of
  every order with their powers;
- the progress and wall time per cluster order;
- the steps of each state (product, gather, save) and their wall times;
- with `k_range`, one `[configuration]` line per rank with its memory after each configuration.

For a bath of a thousand spins this is thousands of lines.

<a id="RunLogTags"></a>
### Tagged lines

A line that starts with a tag in brackets reports a decision or a measurement. Both logs keep them.

| Tag | When | What it says |
|---|---|---|
| `[backend]` | `execution.backend` `auto` or `cuda` | which device each rank uses, or why the CPU is used |
| `[gpu-memory-fallback]` | `auto` | the GPU ran out of memory; this cluster and the larger ones go to the CPU on that rank |
| `[dm]` | `quantity` dm | clusters put values into density-matrix elements the free evolution keeps at zero ([The density matrix](OutputFiles.md#out_dm)) |
| `[seed]` | a generated seed (standard error) | the seed to write into the input to replay the run |
| `[numerical-invalid]` | a cluster term is not finite (standard error) | the first such value on a rank: the state, time step, cluster order and value |
| `[partition-cache]` | `cluster.pcce.cache` | `hit`, `miss` (with the reason) or `stored`, and the file |
| `[cluster-plan]` | `output.cluster_only` | the sites, the clusters and the path of `.clusters.json` |
| `[cluster-contributions]` | `output.cluster_contributions` | how many rows the table has and how many were added for it, the time and the memory per rank; the rows are in the `.tsv` file ([Cluster contributions](OutputFiles.md#out_contributions)) |
| `[auto]` | an `auto_convergence` block | the runner, the Python and the binary used |

Other tags report the details of pCCE partitions, bundles and the adaptive refinement, the GPU kernels, and the
timings of [execution.profile](InputFile_Backend.md#tag_profile); they matter mostly for performance studies.

<a id="EnvironmentVariables"></a>
## 5. Environment variables

An ordinary run needs none of them. Set them in front of the command
(`CUDA_VISIBLE_DEVICES=0 mpirun -n 4 bin/main_gpu.out -f ccein.json`) or export them in a job script.

| Variable | Read by | Effect |
|---|---|---|
| `CUDA_VISIBLE_DEVICES` | GPU build (`main_gpu.out`) | the cards the run may use. Empty, or matching no card, hides them all, and `execution.backend` `auto` then stays on the CPU |
| `CCEX_PYTHON` | a run with an `auto_convergence` block | the Python that runs the automatic-convergence tool (default `python3`; it needs numpy, scipy, matplotlib) |
| `CCEX_AUTOCONV_DIR` | a run with an `auto_convergence` block | the tool's folder, used in place of the `example/tools/AutoConvergence` next to the binary. A folder others can write to is refused |

The memory check also reads the memory Slurm gives the job (`SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU`,
`SLURM_CPUS_ON_NODE`), and a `MKL_CBWR` set by the user is kept. Variables that turn optimisations off for a comparison, or serve the checks, are listed in
[test/README.md](../../test/README.md). The build variables are in [section 1](#Building).

<a id="Troubleshooting"></a>
## 6. Troubleshooting

- **Every process runs the whole calculation alone** (`Clusters per order` ends with `1 rank`). Another MPI's
  `mpirun` comes first on `PATH`, for example Anaconda's. It starts the processes without joining them. Compare
  `which mpirun` with the MPI in the `Build` line of the log, and put the right one first or call it by its full path.
- **`Error: not enough memory: ...`** before the first cluster. The largest cluster does not fit in the memory of one
  rank. The message says what to change: fewer ranks per node, a lower order or smaller clusters
  ([the memory check](InputFile_Backend.md#memory_check)).
- **No GPU.** With `execution.backend` `auto`, every cluster stays on the CPU when the binary is CPU-only, a card cannot
  be opened, or `CUDA_VISIBLE_DEVICES` hides all cards. The `[backend]` lines say why. With `cuda` the run stops
  instead, with an `Error: GPU backend: ...` line from every rank. The GPU binary is `bin/main_gpu.out` ([GPU execution](../methods/gpu.md)).
- **`Error: input key "..." is not read, did you mean "..."?`** A misspelled key; the run stops before the calculation starts. `Warning: input key ... is not read by this run` only means that this run does not use the key.
- **A key of an older version** is refused with the name that replaced it ([old names](InputFile.md#OldNames)).
  `python3 example/tools/InputKeys/convert_input_keys.py OLD.json` rewrites such an input.
- **A curve above 1** ($|L| > 1$ in `_wiDiv`). A cluster term near zero was divided by; the value is an artefact.
  See [The coherence files](OutputFiles.md#out_coherence) and the analyzer's error correction ([Analysis](Analysis.md)).
