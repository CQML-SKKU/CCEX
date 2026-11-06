# Getting started
<a id="GettingStarted"></a>

This page takes you from a fresh checkout to a first result. The build options, the command-line options and the run
log are described in full in [Building and running](BuildingAndRunning.md).

1. [What you need](#Prerequisites): a compiler with MPI, and Python for the plots.
2. [Building](#Installing): the commands for Ubuntu or WSL, macOS, a cluster, a GPU, and CMake.
3. [Checking the build](#CheckBuild): the unit tests and the reference comparison.
4. [A first calculation](#FirstCalculation): the Hahn echo of an NV centre, in a few seconds.
5. [Your own input](#OwnInput): a template with every key, and what to do when a run fails.

<a id="Prerequisites"></a>
## 1. What you need

| Needed | For |
|---|---|
| a C++17 compiler behind an MPI wrapper (`mpicxx`, `mpic++` or `mpiicpc`) | `bin/main.out` |
| `make` or CMake 3.18 or newer | the build |
| the CUDA toolkit (nvcc, cuSOLVER, cuBLAS) | only the GPU binary `bin/main_gpu.out` ([GPU execution](../methods/gpu.md)) |
| Intel MKL | optional, used when the build finds it |
| Python 3 with numpy and matplotlib | the figures of the tutorials |
| scipy as well | the analyzer ([Analysis](Analysis.md)) |
| Python 3.10 or newer with numpy, scipy and matplotlib | the automatic convergence |
| mpi4py and psutil | making bath files ([Tools](Tools.md)) |

Eigen, uthash and SQLite come with the code in `external/`. Nothing else has to be installed for `bin/main.out`.

<a id="Installing"></a>
## 2. Building

Run the build from the repository root. It writes `bin/main.out`, and also `bin/main_gpu.out` where CUDA is found.
Pick the case that fits your machine.

### Ubuntu or WSL

Install a compiler, Open MPI and CMake from the distribution, then use the `generic` preset, which takes `mpicxx` from
`PATH`:

```sh
sudo apt install build-essential cmake openmpi-bin libopenmpi-dev \
                 python3-numpy python3-matplotlib python3-scipy
bash do_compile.sh generic
```

With Intel oneAPI installed (`setvars.sh`), use `bash do_compile.sh wsl-gnu` or `wsl-icpx` instead.

### macOS

Install an MPI wrapper from Homebrew, then use the `mac` preset:

```sh
brew install open-mpi
bash do_compile.sh mac
```

Do not build in a folder that iCloud Drive syncs, such as Desktop or Documents, because the build fails there.

### A cluster with environment modules

Load a compiler and an MPI, then use `generic`. Set `CXX` when the wrapper is not called `mpicxx`:

```sh
module load <compiler> <mpi>
CXX=mpiicpc bash do_compile.sh generic
```

### With a GPU

Give the CUDA folder and the compute capability without the dot (86 for an RTX 30 series card):

```sh
CUDA_HOME=/usr/local/cuda CUDA_ARCH=86 bash do_compile.sh generic
```

### With CMake

CMake works on any machine, in place of `do_compile.sh`:

```sh
cmake -B build && cmake --build build -j 8
```

The binaries are then in `build/bin/`. Point the tutorials at them with `export CCEX_BIN=$PWD/build/bin/main.out`.

The other presets, the build variables, the CMake options and MKL are in
[Building and running](BuildingAndRunning.md#Building).

<a id="CheckBuild"></a>
## 3. Checking the build

The quickest check is the first calculation in the next section. For a full check, run the unit tests and compare
the results with the recorded references:

```sh
make numeric-test cpu-test                  # CMake: cmake --build build --target numeric-test cpu-test
bash test/code_verification/run.sh         # about 8 minutes on 28 cores
```

`test/README.md` lists every check and what it needs.

<a id="FirstCalculation"></a>
## 4. A first calculation

The first tutorial computes the Hahn echo of an NV centre in diamond with natural-abundance ¹³C. It takes a few
seconds on four ranks:

```sh
cd example/tutorial/basic/01_first_run
NP=4 bash run.sh
```

It prints T2 and writes `results/coherence.png`. The figure should match `expected.png` in the same folder. The log
is `results/ccein.log`. It ends with a `JOB DONE` line when the run finished. A run that stops prints `Error: ...` on the
terminal.

[Basic 01: NV in ¹³C, from a bath to T2](../../example/tutorial/basic/01_first_run/README.md) explains this run step by
step. The other tutorials are listed in [Examples](Examples.md).

<a id="OwnInput"></a>
## 5. Your own input

To write an input template, run

```sh
bin/main.out -e
```

This writes `ccex_example.jsonc`, the input of the first calculation, with every other key as a comment that gives its
meaning and default. To use it:

1. Set `configuration.bathfile` to a bath file. The template does not run without one.
   `example/tutorial/basic/01_first_run/bath/bath_Diamond_NV_natab_1` is one, and [Tools](Tools.md) makes others.
2. Change the keys you need ([Input file](InputFile.md)).
3. Run it: `mpirun -n 4 bin/main.out -f ccex_example.jsonc`.

When a run fails or gives a strange result, see [Troubleshooting](BuildingAndRunning.md#Troubleshooting).
