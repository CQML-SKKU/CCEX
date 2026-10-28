# CCEX

CCEX computes the decoherence of spin qubits in a bath of nuclear and electron spins by the cluster-correlation
expansion (CCE). From an input file that names the qubit(s), the bath and the pulse sequence, it builds the spin
Hamiltonian, cuts the bath into clusters, solves each cluster exactly and combines them into the coherence L(t) of the
qubit (or, with gCCE, its density matrix). It is a C++17/MPI code with an optional CUDA backend, by the CQML group at
Sungkyunkwan University (<https://github.com/CQML-SKKU/CCEX>); this is version 2.0.1, released 2026-10-22.

## Quick start

```
bash do_compile.sh wsl-gnu                 # or another preset (see Build), or: cmake -B build && cmake --build build -j 8
cd example/tutorial/00_quickstart && NP=4 bash run.sh             # results/ccein_wiDiv, _noDiv (about 1 s), coherence.png, prints T2
                                                                  # compare results/coherence.png with expected.png (T2 1.068 ms)
```

`bin/main.out -e` writes an input template, `ccex_example.jsonc`: the settings of 00_quickstart with every other key
commented, its default beside it; it runs once its `bathfile` names a bath file.

After a CMake build the binary is `build/bin/main.out`; the tutorials take it as `CCEX_BIN=<repository>/build/bin/main.out`.

## What runs where

`method` × `cluster.mode` × `execution.backend`. "CPU, GPU" means `cuda` and `auto` may put clusters on the card;
"CPU" means `cuda` is refused and `auto` runs on the CPU; "refused" means the run stops with a message.

| method | `standard` | `pcce` | `bundle` | `dsj`, `itb`, `dsjitb` |
|---|---|---|---|---|
| `cCCE` (default) | CPU, GPU | CPU, GPU | CPU, GPU | CPU |
| `gCCE` | CPU, GPU¹ | CPU | CPU, GPU¹ | refused |
| `meCCE` (cCCE + Lindblad jumps on bath spins) | CPU | CPU | refused | refused |
| `megCCE` (gCCE + jumps on bath spins and the qubit) | CPU | CPU | refused | refused |

¹ with `execution.propagator` `eigen` (the default); `expm` runs on the CPU. Also CPU only: DEER bath pulses
(`pulse.bath`, gCCE and megCCE with `standard` or `pcce`) and an ensemble (no `s_range`) `Manual` sequence. `quantity`
`dm` is for gCCE; meCCE and megCCE compute `coherence` only. The GPU needs the CUDA binary `bin/main_gpu.out`. Details:
[method](docs/user/InputFile_General.md#tag_method), [cluster.mode](docs/user/InputFile_Cluster.md#tag_cluster_mode),
[backend](docs/user/InputFile_Backend.md#tag_backend), [GPU](docs/methods/gpu.md),
[master equation](docs/methods/master-equation.md).

## The input at a glance

One JSON file (`.json` or `.jsonc`) with `//`, `#` and `!` comments; keys are read without regard to case, a missing key
takes its default, a misspelt or misplaced key is reported or refused. `bin/main.out -e` writes an input template with
every key, its allowed values and its default. The full reference is [docs/user/InputFile.md](docs/user/InputFile.md).

Top level: `method`, `quantity` (`coherence` or `dm`), `order`, `bfield`, `rbath`, `rdip`, `deltat`, `nstep`
([General](docs/user/InputFile_General.md)). A number is in the key's unit (`deltat` ms, `bfield` G, `rbath` and
`rdip` Å, ...); those keys also take `{"value": 5, "unit": "us"}` ([units](docs/user/InputFile.md#Units)).

| Block | What it holds | Page |
|---|---|---|
| `configuration` | bath and gyro files, sampled states and defect axes, `s_range`, `seed`, several configurations (`k_range`) | [General](docs/user/InputFile_General.md#tag_configuration) |
| `Qubit` | the qubit(s) (`qubit`: positions or a file, spins), their states (`alpha_state`, `beta_state`), interactions and zero-field splittings (`intmap`) | [Qubit](docs/user/InputFile_Qubit.md) |
| `Defect` | the species of a paramagnetic-defect bath (P1, ...) | [Defect](docs/user/InputFile_Defect.md) |
| `hamiltonian` | `common` (mean fields, hyperfine and quadrupole tensors) and one section per method: `cCCE` (`knight`, `hfmedi`, `delta_A`), `meCCE` / `megCCE` (`jump_operators`) | [General](docs/user/InputFile_General.md#tag_hamiltonian) |
| `pulse` | the sequence: Ramsey, Hahn, CPMG, WAHUHA or written out; DEER bath pulses | [Pulse](docs/user/InputFile_Pulse.md) |
| `cluster` | `mode`, the cluster selection, and the `pcce`, `bundle`, `adaptive` blocks | [Cluster](docs/user/InputFile_Cluster.md) |
| `execution` | `backend`, GPU options, propagator and evolution choices, profiling | [Backend](docs/user/InputFile_Backend.md) |
| `output` | the output folder and name, what is saved | [Output](docs/user/InputFile_Output.md) |
| `frame` | `rotation` (`qubit_axis`) from the frame of the bath file to the qubit frame; each tensor says its own `frame` (`"file"` or `"qubit"`) | [General](docs/user/InputFile_General.md#tag_frame) |
| `auto_convergence` | chooses rbath, rdip, the time window and the number of configurations, then runs | [AutoConvergence](docs/user/InputFile_AutoConvergence.md) |

Inputs written for version 1.2.2 or earlier use old key names; they are refused with the new name, and
`python3 example/tools/InputKeys/convert_input_keys.py OLD.json` rewrites them ([old and new names](docs/user/InputFile.md#OldNames)).

## Command line

| Option | |
|---|---|
| `-f FILE` | the input file; every setting of a run is a key in it |
| `-o DIR` | the output folder, made when missing |
| `-v` | a longer log |
| `-e` | write the input template `ccex_example.jsonc` (in `DIR` with `-o DIR`) and stop |
| `-h`, `--help` | a short usage |

Output folder: `-o DIR`, else `output.path`, else the working directory. The output name is `output.file`, else the
input file's name without `.json` or `.jsonc`. The options of earlier versions (`-I`, `-N`, `-m`, `-q`, `-B`, ...) are refused with
the key that replaces them ([Getting started](docs/user/GettingStarted.md#HowToRun)).

## Output files

| File | Holds |
|---|---|
| `<name>_wiDiv` | the CCE result: a time column (ms), then the coherence (or, with `dm`, the density matrix) as complex numbers |
| `<name>_noDiv` | the same product without the division by the sub-clusters; the analyzer uses it to repair diverged points (`output.savemode.nodiv` false leaves it out) |
| `<name>_state<s>_wiDiv`, `_noDiv` | one per sampled state s of `s_range` [first, last]; the ensemble (no `s_range`) has no `_state` |
| `<name>_k<k>_...` | one set per configuration with `configuration.k_range` (no average) |

`output.savemode` adds per-cluster files ([Output](docs/user/InputFile_Output.md#tag_savemode)). Averaging over
configurations and states, error correction and T2 fits are done by `example/analysis/` ([Analysis](docs/user/Analysis.md)).
The run log (stdout) starts with a `Build :` line: CPU or GPU, compiler, MPI, MKL on/off, preset and commit.
Every output file and its columns: [Output files](docs/user/OutputFiles.md); what the run log says:
[the run log](docs/user/GettingStarted.md#RunLog).

## Build

| Preset (`bash do_compile.sh <preset>`) | For |
|---|---|
| `wsl-gnu` | WSL/Ubuntu with Intel oneAPI: g++ behind Intel MPI |
| `wsl-icpx` | the same with icpx |
| `idm4` | the idm cluster: Intel 2023.2 module, mpiicpc with icpx |
| `gpu` | idm4-gpu: module `intel/oneapi-2023.2` (mpiicpc with icpx, MKL) + CUDA 13.2; also builds `bin/main_gpu.out` |
| `nurion` | KISTI Nurion: intel/19.0.5 + impi/19.0.5 |
| `mac` | macOS with an MPI wrapper (Homebrew open-mpi or mpich): Apple clang `-O3 -mcpu=native` (Apple Silicon), Apple's Accelerate for the large clusters |
| `generic` | any machine with an MPI C++ wrapper; `CXX`, `CXXFLAGS`, `CUDA_HOME`, `CUDA_ARCH` from the environment |

- `do_compile.sh` writes the top-level `Makefile`, `obj/` and `bin/main.out` (and `bin/main_gpu.out` with CUDA).
- CMake: `cmake -B build && cmake --build build -j 8` keeps everything in `build/` (`build/bin/main.out`).
- Needs: a C++17 compiler behind an MPI wrapper; Eigen, uthash and SQLite are bundled in `external/`.
- MKL: used for the large clusters when the build finds it (on macOS Apple's Accelerate in its place); `CCEX_BLAS=off`
  (CMake `-DCCEX_BLAS=OFF`) leaves it out.
- CUDA: requires CUDA 13.2 or newer (12.8 miscompiles a kernel for the RTX PRO 6000, see `CHANGELOG.md`;
  `do_compile.sh` and CMake refuse an older nvcc). The `gpu` preset, `generic` with `CUDA_HOME`, or CMake when it finds nvcc (`-DCCEX_CUDA=OFF`
  to skip).

Every option, environment variable and the MKL details: [Getting started](docs/user/GettingStarted.md#Building).

## Folders

| Folder | What | Read |
|---|---|---|
| `src/`, `include/` | the code: input, physical objects, clustering, calculation, common; `src/main.cpp` is the entry point | `src/README.md` |
| `example/` | `tutorial/` (twelve worked examples, `00_quickstart` … `11_noise`, each with `run.sh` and `expected.png`), `analysis/` (the coherence analyzer), `tools/` (automatic convergence, bath generator, input-key converter, ...) | `example/README.md` |
| `test/` | `code_verification/`, `cluster_verification/`, `gpu_verification/` | `test/README.md` |
| `docs/` | `user/` (the manual), `methods/` (the added methods in detail), `doxygen/` (the API site: `cd docs/doxygen && doxygen Doxyfile`) | `docs/README.md` |
| `external/` | Eigen, uthash and SQLite, so every machine builds without packages | `external/README.md` |

## Checks

All need `bin/main.out` and `mpirun`; the Python ones numpy, scipy and matplotlib (`test/requirements.txt`). For a
CMake build set `CCEX_BIN=$PWD/build/bin/main.out CCEX_BUILD_DIR=$PWD/build`. Details: [test/README.md](test/README.md).

| What | Command | Needs |
|---|---|---|
| the core: closed forms, numerics against recorded references, frame rotation, refusals, input forms, ... | `bash test/code_verification/run.sh` | `bin/main.out`, `mpirun`, Python |
| unit tests | `make numeric-test cpu-test` (CMake: the same targets) | the build tree |
| bundle clustering and the adaptive refinement | `sh test/cluster_verification/check_adaptive/verify.sh unit\|regress\|inprocess` | the build tree, Python |
| the automatic convergence runner | `cd example/tools/AutoConvergence && python -m pytest -q -n 8 tests` | pytest, pytest-xdist |
| the CUDA backend | `make gpu-test` and the jobs of `test/gpu_verification/server/` | a CUDA build, a GPU |

`.github/workflows/ci.yml` runs the fast, machine-independent part on every push and pull request.

## For contributors

Where to look: [src/README.md](src/README.md) (what each source file does, how to add one),
[test/README.md](test/README.md) (every check, what it answers and needs), [docs/README.md](docs/README.md) (the
documentation map) and [docs/user/](docs/user/) (the manual, one page per input group).

What keeps the results trustworthy:

- A change that is not meant to change the physics leaves the outputs byte-identical; run the checks above for the
  part you touched (`test/README.md` says which).
- Reference outputs (`test/code_verification/check_numerics/reference/`, the adaptive golden answers, `bath_data/`)
  change only on purpose, in the same commit as the change, with the reason; never re-record them to make a check pass.
- When an input key is added, renamed or moved, update together: the converter
  `example/tools/InputKeys/convert_input_keys.py`, the `-e` example `src/input/example_input.cpp` (held to the readers'
  key names by `test/code_verification/check_input_forms/cli.sh`), the page in `docs/user/`, and `CHANGELOG.md`.
- A new refusal gets a case in `test/code_verification/check_refusals/`.
- A new source file must be seen by both builds: `do_compile.sh` (never edit the generated `Makefile`) and
  `CMakeLists.txt`. `src/clustering/algo_adaptive/` keeps its own floating-point flags in both (`src/README.md`).

## Changes and citing

- What changed in 2.0.0 since v1.2.2, how to update old inputs, the measured speed-ups and the results that move:
  [CHANGELOG.md](CHANGELOG.md).
- The authors and the version to cite: [CITATION.cff](CITATION.cff).
