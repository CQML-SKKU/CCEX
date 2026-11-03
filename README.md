# CCEX

CCEX computes the decoherence of spin qubits in a bath of nuclear and electron spins with the cluster-correlation
expansion (CCE). The input file names the qubits, the bath and the pulse sequence. CCEX builds the spin Hamiltonian,
splits the bath into clusters, solves each cluster exactly and multiplies the cluster terms into the qubit coherence
L(t), or into the density matrix with gCCE. It is written in C++17 with MPI and an optional CUDA backend by the CQML
group at Sungkyunkwan University. This is version 2.1.0, released 2026-11-03.

Documentation: <https://cqml-skku.github.io/CCEX/>. The same pages are in [docs/user/](docs/user/).

## Features

- Methods: cCCE and gCCE, their master-equation versions meCCE and megCCE, the Gaussian phase-noise model `noise`,
  and the pseudospin model of nuclear pairs `pseudospin`.
- Clustering by order and distance, by partitions (pCCE) or by bundles, with an optional adaptive refinement.
- Nuclear spin baths with point-dipole or DFT hyperfine and quadrupole tensors, and electron-spin defects such as P1.
  The bath can be in its ensemble state or in sampled states, and one input can run several bath configurations.
- Ramsey, Hahn echo, CPMG, sequences written out pulse by pulse, and DEER pulses on the bath.
- MPI on CPUs, a CUDA backend for large clusters, and an automatic convergence of rbath, rdip, the time window and the
  number of configurations.

The table [What runs where](docs/user/Overview.md#WhatRunsWhere) shows which methods, cluster modes and backends work
together.

## Quick start

```
bash do_compile.sh wsl-gnu                     # build bin/main.out (presets below)
cd example/tutorial/00_quickstart
NP=4 bash run.sh                               # about 1 s: results/ccein_wiDiv and coherence.png (T2 1.068 ms)
```

To run your own input, use `mpirun -n 4 bin/main.out -f ccein.jsonc -o results`. `bin/main.out -e` writes
`ccex_example.jsonc`, a template that lists every key with its meaning and default.

## Build

`bash do_compile.sh <preset>` writes `bin/main.out`, and also `bin/main_gpu.out` when CUDA is enabled. The presets are
`wsl-gnu`, `wsl-icpx`, `idm4`, `gpu` (idm4-gpu with CUDA), `nurion`, `mac` and `generic` (any MPI C++ wrapper). CMake
also works: `cmake -B build && cmake --build build -j 8` puts everything in `build/`.

You need a C++17 compiler behind an MPI wrapper. Eigen, uthash and SQLite are included in `external/`, and MKL is used
when the build finds it. [Getting started](docs/user/GettingStarted.md#Building) covers the build options, the CUDA
versions and the environment variables.

## Documentation

- [Tutorials](example/tutorial/README.md): thirteen examples in learning order, from the quickstart to GPU runs, the
  automatic convergence and the noise and pseudospin models. Each has commented inputs, a `run.sh` and an
  `expected.png`.
- [Input file](docs/user/InputFile.md): every key by block, with its unit and default.
- [Output files](docs/user/OutputFiles.md) and the [analyzer](docs/user/Analysis.md).
- [Methods](docs/methods/Methods.md): cCCE and gCCE, pCCE, bundles, master-equation CCE, DEER and GPU execution.

## Folders

| Folder | Contents |
|---|---|
| `src/`, `include/` | the code ([src/README.md](src/README.md)) |
| `example/` | the tutorials, the analyzer and the tools (bath generator, tensor reader, automatic convergence, input keys) |
| `test/` | the checks ([test/README.md](test/README.md)) |
| `docs/` | the manual, the method notes and the website ([docs/README.md](docs/README.md)) |
| `external/` | Eigen, uthash and SQLite |

## Checks and contributing

`bash test/code_verification/run.sh` runs the core checks in about 8 minutes. [test/README.md](test/README.md) lists
the checks each kind of change needs. `.github/workflows/ci.yml` runs the fast checks on every push.

When you add a source file, add it to both builds ([src/README.md](src/README.md)). Change a reference output only on
purpose, in the same commit as the change that moves it, and give the reason. When you add or rename an input key,
update the converter, the `-e` template, its page in `docs/user/` and `CHANGELOG.md` in the same change.

## Changes and citing

[CHANGELOG.md](CHANGELOG.md) lists every change, and the [Changes](docs/user/Changes.md) page summarizes each version.
[CITATION.cff](CITATION.cff) gives the authors and the version to cite.
