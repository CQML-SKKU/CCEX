# Overview

CCEX computes the coherence of spin qubits in a bath of nuclear and electron
spins by the cluster-correlation expansion (CCE). It is a C++ and MPI code of the
CQML group at Sungkyunkwan University: <https://github.com/CQML-SKKU/CCEX>.

## What it computes

The input names a qubit (or several), the bath around it and a pulse sequence.
From the positions, gyromagnetic ratios and interaction tensors CCEX builds the
spin Hamiltonian, divides the bath into clusters, solves the dynamics of every
cluster exactly and combines them into the quantity asked for as a function of
time:

- the coherence function `L(t)` of the qubit under the pulse sequence, or with
  gCCE its density matrix (`quantity`);
- by conventional CCE (cCCE) or generalized CCE (gCCE) (`method`), with partition CCE (pCCE), bundle
  clustering and the disjoint-cluster and inter-bath variants dsj, itb and dsjitb (`cluster_mode`);
- for a Hahn echo, a Ramsey experiment, CPMG (`pulse.count`) or any sequence written
  out pulse by pulse (`sequence`);
- for one bath configuration or an ensemble of them, and for electron-spin baths
  with sampled bath states (`nstate`).

Baths can be nuclear spins with point-dipole or DFT tensors (hyperfine and
quadrupole), electron-spin defects such as P1 centres with their own on-site
spins and axes, or both. Beyond the standard clustering by order and distance
CCEX has bundle clustering (`cluster_mode` bundle, `bundle.max_dim`) and an adaptive
refinement of the clusters (the `adaptive` block), a CUDA backend for the
cluster propagation (`backend`), and a coordinate frame rotation that lets the
bath file stay in the crystal frame (`coordinate_frame_rotation`).

## Building

```
cmake -B build && cmake --build build -j 8      # any machine: finds the MPI wrapper and, if there is one, nvcc
bash do_compile.sh <preset>                     # the group's machines: wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac (or generic)
```

`do_compile.sh` writes `bin/main.out` and, where there is CUDA, `bin/main_gpu.out` for
`backend` cuda or auto; the CMake line writes the same two files under `build/bin/`
(pass `-DCCEX_BIN_DIR=$PWD/bin` to put them in `bin/` instead). The build needs a C++17 compiler behind an MPI wrapper
and, for the GPU binary, the CUDA toolkit. Eigen, uthash and SQLite are
bundled in `external/`. Intel MKL is linked when the build finds it (`$MKLROOT` from setvars.sh or the intel
module; the build prints `[build] MKL on/off`, the run log's header too): it takes the large cCCE/gCCE clusters
(zgemm in the recursive path of the cCCE ensemble Hahn echo from bath dimension 21; zheevd from dimension 63, the bath
dimension for cCCE and the whole qubit × bath dimension for gCCE), which helps runs whose time is in such clusters (h-BN: 0.72–0.88 of
the time on idm; P1 unchanged); the curves then differ in the last digits (up to 9e-8 at late times of h-BN gCCE).
`CCEX_BLAS=off bash do_compile.sh <preset>` (CMake `-DCCEX_BLAS=OFF`) builds without it. `CMakeLists.txt` takes
`-DCCEX_CUDA=OFF`, `-DCMAKE_CUDA_ARCHITECTURES=<cc>` and
`-DMPI_CXX_COMPILER=<wrapper>`; `do_compile.sh` loads the modules of the
group's machines itself, and its `generic` preset takes `CXX`, `CUDA_HOME` and
`CUDA_ARCH` from the environment.

## Running

```
mpirun -n 4 bin/main.out -f ccein.json
```

The input file is JSON with `#` and `!` comments; its keys are the reference in
[Input file](InputFile.md), and the command-line options that override them are
in [Getting started](GettingStarted.md). A run writes the coherence (or density
matrix) per bath configuration into text files, which the analyzer in
`example/analysis/` error-corrects, averages over the ensemble and plots.

## Where to go next

- [Getting started](GettingStarted.md): the command-line options, the input
  files (bath, gyromagnetic ratios, qubit, states, axes) and the order in which
  a run reads them.
- [Input file](InputFile.md): every key, by group.
- [Examples](Examples.md): twelve worked examples in `example/tutorial/`, from a
  hand-written bath to the GPU backend, the frame rotation, the automatic
  convergence and the rdip ladder.
- [Analysis](Analysis.md) and [Tools](Tools.md): what runs on the output, and what
  makes the inputs (`example/analysis/`, `example/tools/`).
- Data structures and Source files: the code itself, for developers.
- Methods: [clustering terminology](../methods/clustering-terminology.md),
  [adaptive cluster mode](../methods/adaptive-cluster-mode.md), [the `adaptive` block](../methods/adaptive-input-block.md),
  [GPU execution](../methods/gpu.md) and [the coordinate frame rotation](../methods/coordinate-frame-rotation.md),
  in detail (`docs/methods/` in the repository; partly in Korean).

## Citing

The authors and the version are in `CITATION.cff` at the repository root.
