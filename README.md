# CCEX

CCEX computes the coherence of spin qubits in a bath of nuclear and electron spins
by the cluster-correlation expansion (CCE). It is a C++ and MPI code of the CQML
group at Sungkyunkwan University: <https://github.com/CQML-SKKU/CCEX>. This is
version 2.0.0 (provisional; see `CITATION.cff`).

## What it computes

The input names a qubit (or several), the bath around it and a pulse sequence.
From the positions, gyromagnetic ratios and interaction tensors CCEX builds the
spin Hamiltonian, divides the bath into clusters, solves the dynamics of every
cluster exactly and combines them into the coherence `L(t)` of the qubit or, with
gCCE, its density matrix.

- Methods (`method`): conventional CCE (`cCCE`) and generalized CCE (`gCCE`); the names are read
  without regard to case.
- Clustering (`cluster_mode`): `standard` (the default), clusters of spins by order and distance;
  `pcce`, partition CCE (the `pcce` block); `dsj`, `itb`, `dsjitb`, the disjoint-cluster and inter-bath
  variants (these four with cCCE only); or `bundle`: strongly coupled spins tied into bundles and
  clusters made of bundles (the `bundle` block), with an adaptive refinement that rebuilds the bundles
  where the error is large (the `adaptive` block).
- Input keys were renamed and grouped on 2026-10-13; an older input is refused with the new name, and
  `example/tools/InputKeys/convert_input_keys.py` rewrites it ([the table](docs/user/InputFile.md#OldNames)).
- Sequences: Ramsey, Hahn echo, CPMG, WAHUHA, or any sequence written out pulse
  by pulse (the `pulse` block, `sequence`).
- Baths: nuclear spins with point-dipole or DFT tensors (hyperfine, quadrupole),
  electron-spin defects with sampled bath states; one configuration or an
  ensemble.
- A CUDA backend for the cluster propagation (`backend`, the `gpu_options`
  block), and a coordinate frame rotation that lets the bath file stay in the
  crystal frame (`coordinate_frame_rotation`).
- Automatic convergence: an `auto_convergence` block chooses rbath, rdip, the time
  window and the number of bath configurations (`example/tools/AutoConvergence/`).
- `rdip_ladder` writes the curves of smaller rdip values from one order-2 run;
  `evolution_update` (`recursive` by default where it exists, or `exact`) sets how
  the cluster phases are updated in time; `gyrofile` is optional (a built-in
  isotope table).

## Build

```
cmake -B build && cmake --build build -j 8      # any machine: finds the MPI wrapper and, if there is one, nvcc
bash do_compile.sh <preset>                     # the group's machines: wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac (or generic)
```

`do_compile.sh` gives `bin/main.out`, and, where there is CUDA, `bin/main_gpu.out` for
`backend` cuda or auto; CMake gives the same two in `build/bin/` (`-DCCEX_BIN_DIR=$PWD/bin`
puts them in `bin/`). Needed: a C++17 compiler behind an MPI wrapper
(`mpicxx`, `mpic++`, `mpiicpc`) and, for the GPU binary, the CUDA toolkit;
Eigen, uthash and SQLite are bundled in `external/`. Intel MKL is linked when the build finds it (`$MKLROOT`
from setvars.sh or the intel module, with its libraries there; the build prints `[build] MKL on (...)` or
`MKL off (why)`, and every run log's header says which): it takes the large cCCE/gCCE clusters (zgemm in the
recursive path of the cCCE ensemble Hahn echo from bath dimension 21; zheevd from dimension 63, the bath dimension for
cCCE and the whole qubit × bath dimension for gCCE). `CCEX_BLAS=off bash do_compile.sh <preset>` (CMake `-DCCEX_BLAS=OFF`) builds without
it, `CCEX_BLAS=mkl` (`-DCCEX_BLAS=MKL`) stops when it is missing. On idm it ran h-BN 0.72–0.88 of the time and left
P1 unchanged; the curves move in the last digits (up to 9e-8 at late times of h-BN gCCE). `CMakeLists.txt` takes `-DCCEX_CUDA=OFF`,
`-DCMAKE_CUDA_ARCHITECTURES=<cc>`, `-DMPI_CXX_COMPILER=<wrapper>`,
`-DCCEX_NATIVE=OFF`, `-DCCEX_SANITIZE=ON` and `-DCCEX_STRICT_WARNINGS=ON` (warnings on, format, sequence-point
and misleading-indentation ones as errors; the CI builds with it); its objects stay in `build/`.
`do_compile.sh` writes a `Makefile` at the top and loads the modules of the
group's machines itself; `generic` takes `CXX`, `CC`, `CXXFLAGS`, `CUDA_HOME`
and `CUDA_ARCH` from the environment. `bin/`, `obj/`, `build/` and the
`Makefile` are build output and not in git.

## Run

```
mpirun -n 4 bin/main.out -f ccein.json           # do_compile.sh build
mpirun -n 4 build/bin/main.out -f ccein.json     # CMake build
```

The input file is JSON with `#` and `!` comments. Its keys are the reference in
[docs/user/InputFile.md](docs/user/InputFile.md); the command-line options that
override them are in [docs/user/GettingStarted.md](docs/user/GettingStarted.md).
A run writes the coherence (or density matrix) per bath configuration into text
files, which the analyzer in `example/analysis/` error-corrects, averages over
the ensemble and plots. The quickest way in is a tutorial:

```
cd example/tutorial/02_Diamond_NV_natab/ensemble_CCE
NP=4 NCONF=1 bash run.sh                      # writes results/
```

The tutorials run `bin/main.out`; after a CMake build give them the binary with
`export CCEX_BIN=<repository>/build/bin/main.out` (and `CCEX_GPU_BIN` for
`09_GPU_backend`).

## Layout

| Folder | What | Read |
|---|---|---|
| `src/`, `include/` | the code: five core groups (input, physical, clustering, calculation, common) plus `clustering/bundle/`, `clustering/algo_adaptive/`, `calculation/gpu/` and `auto_convergence/` (the rdip ladder); `src/main.cpp` is the entry point | `src/README.md` |
| `example/` | `tutorial/` (twelve worked examples, one system or feature each, with `run.sh`), `analysis/` (the coherence analyzer and the convergence test), `tools/` (automatic convergence, bath generator, isotope table, job runner, DFT tensor reader) | `example/README.md` |
| `test/` | `code_verification/` (a closed form, numerics against recorded references, the frame rotation, refusals, input forms, the mean field, the rdip ladder), `cluster_verification/` (adaptive, bundle, pCCE), `gpu_verification/` (the CUDA backend, benchmarks, server jobs) | `test/README.md` |
| `docs/` | `user/` (the manual pages), `methods/` (the added methods in detail), `doxygen/` (the site) | `docs/README.md` |
| `external/` | Eigen 3.4.1, uthash and SQLite, kept in the tree so every server builds without packages | `external/README.md` |

## Documentation

- The manual: [Overview](docs/user/Overview.md), [Getting started](docs/user/GettingStarted.md),
  [Input file](docs/user/InputFile.md) (every key, by group) and [Examples](docs/user/Examples.md).
- The same pages with the API reference as a site: `cd docs/doxygen && doxygen Doxyfile`,
  then open `html/index.html` (Doxygen 1.9 and graphviz). GitHub Actions builds
  it from `.github/workflows/deploy_docs.yml`.
- The methods this branch added, in detail: `docs/methods/` (clustering
  terminology, bundle clustering, the `adaptive` block, the GPU backend, the
  frame rotation).

## Checks

```
bash test/code_verification/run.sh          # needs bin/main.out, mpirun, and test/requirements.txt (numpy, scipy, matplotlib)
```

runs the seven code checks in order; the cluster and GPU checks and what each
needs are in `test/README.md`. To check a CMake build, give the checks its binary
and its build folder: `CCEX_BIN=$PWD/build/bin/main.out CCEX_BUILD_DIR=$PWD/build` (the unit-test drivers take
the compiler, flags and objects from the `unit-build.env` the build writes, any generator; `test/unit_build_env.sh`).
`.github/workflows/ci.yml` builds with CMake and runs the fast checks on every push
and pull request.

## Changes

What changed since v1.2.2, with the measured speed-ups and the results that move, is in
[CHANGELOG.md](CHANGELOG.md).

## Citing

The authors and the version are in [CITATION.cff](CITATION.cff).
