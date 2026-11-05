# Advanced 04: The GPU backend

Folder: `example/tutorial/advanced/04_GPU/`

CCEX can run the clusters on an NVIDIA GPU. The result is the same as on the CPU. The GPU saves time when the bath has
large clusters, of Hilbert dimension 256 and more in the cases below.

The GPU runs need `bin/main_gpu.out`, which `bash do_compile.sh gpu` makes together with `bin/main.out`
([Building with CUDA](../../../../docs/user/BuildingAndRunning.md#BuildCuda)). `CCEX_GPU_BIN` names another GPU binary.
Without one, each `run.sh` runs `auto` on the CPU and skips `cuda`. Every key not written in an input has its default
(`main.out -e`).

Run each case on a GPU node with `cd <folder> && NP=32 bash run.sh`, then compare `results/coherence.png` with
`expected.png`. The expected figures were made with 4 ranks on one GPU. With many ranks on one GPU, run the job under
CUDA MPS ([GPU execution](../../../../docs/methods/gpu.md#gpu_ranks)).

This page has three cases:

1. [Running on the GPU](#a04_auto) (`cCCE5_Hahn_autoOptions/`). One key in the input, and CCEX chooses the rest.
2. [The backends and the GPU options](#a04_cce5) (`cCCE5_Hahn_backends/`). `cpu`, `auto` and `cuda` on the same run,
   and the keys of `gpu_options`.
3. [When the GPU saves time](#a04_gcce4) (`gCCE4_Hahn_backends/`). The times of the three backends with gCCE.

<a id="a04_auto"></a>
## 1. Running on the GPU

Folder: `cCCE5_Hahn_autoOptions/`

<a id="a04_auto_input"></a>
### 1.1 The input

The system is the cCCE5 Hahn echo of the NV in five P1, with the `Defect` block of the P1 (¹⁵N, four Jahn-Teller
axes). `../bath/` holds the five P1 (`bath_p1_5`), the last four of them (`bath_p1_4`, section 3), the NV position
(`qubit_p1`) and the P1 gyromagnetic ratio (`gyro_P1`). One cluster has dimension 1024, and the run takes about a
minute of CPU.

To run on the GPU, the input needs one key, `execution.backend`:

```jsonc
"execution" : {"backend": "auto", "profile": true}
```

`auto` sends the large clusters to the GPU and keeps the rest on the CPU. `profile` is optional. It prints the timings
per backend and cluster dimension at the end of the run (`[profile]` lines).

<a id="a04_auto_run"></a>
### 1.2 Running it

Run the input with the GPU binary:

```sh
mpirun -n 32 <repository>/bin/main_gpu.out -f ccein.jsonc -o results
```

`NP=32 bash run.sh` does this, and runs the same input on `cpu` to compare.

For cCCE, CCEX chooses at the start the GPU options the input leaves out, and prints them in a `[gpu-auto]` line
([Execution tags](../../../../docs/user/InputFile_Backend.md#tag_gpu_options)). Writing a key fixes it.

<a id="a04_auto_result"></a>
### 1.3 The result

![cCCE5_Hahn_autoOptions: expected result](cCCE5_Hahn_autoOptions/expected.png)

`auto` is within 1.0e-10 of `cpu`. With 32 ranks on an RTX PRO 6000 the log line is
`[gpu-auto] ranks_per_gpu=4 min_dim=256 evolution_batch_size=32 cluster_batch_mb=16208`. With 4 ranks it has no
`ranks_per_gpu`, which is chosen only with more than four ranks a GPU. Section 2 explains these keys.

The bath has only five P1, so |L| does not decay to zero. The cases of this page are made to compare the backends on
a large cluster, and the P1 bath of [Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md) gives
a full decay.

Next, `"gpu_options": {"ranks_per_gpu": 0}` lets every rank open the GPU, to compare.

Input: [`ccein.jsonc`](cCCE5_Hahn_autoOptions/ccein.jsonc)

<a id="a04_cce5"></a>
## 2. The backends and the GPU options

Folder: `cCCE5_Hahn_backends/`

<a id="a04_cce5_input"></a>
### 2.1 The input

The system of section 1, with `backend` cpu and two GPU options written:

```jsonc
"execution" : {"backend": "cpu", "gpu_options": {"min_dim": 256, "evolution_batch_steps": 20}, "profile": true}
```

The input has `backend` cpu, the reference. `run.sh` runs it, then copies of it with `backend` auto and cuda.

`execution.backend` chooses where the clusters run:

- `cpu` (the default) runs every cluster on the CPU.
- `auto` runs on the GPU every cluster whose Hilbert dimension is at least `gpu_options.min_dim`, and the rest on the
  CPU. A binary built without CUDA runs everything on the CPU.
- `cuda` runs every cluster on the GPU, and is refused when CUDA is not available.

The clusters of `cluster.mode` pcce and bundle go to the GPU as those of standard cCCE do
([GPU execution](../../../../docs/methods/gpu.md#gpu_what)).

The keys of `execution.gpu_options` used most are below, and all of them are in
[Execution tags](../../../../docs/user/InputFile_Backend.md#tag_gpu_options).

| Key | Default | What it sets |
|---|---|---|
| `min_dim` | `768` (`256` with `ranks_per_gpu`) | With `auto`, the smallest cluster dimension sent to the GPU. |
| `ranks_per_gpu` | `0` (`4` with `auto`, cCCE standard or pcce, and more than four ranks a GPU) | How many ranks per GPU send clusters to it. With `0` every rank does. |
| `evolution_batch_size` | `8` (`32` for cCCE) | The most clusters of one dimension evolved together. |
| `evolution_batch_steps` | `64` | The most time points evolved together. |
| `cluster_batch_mb` | `64` (for cCCE, two thirds of the card shared by its ranks) | The GPU memory in MiB one batch of a rank may hold. |

<a id="a04_cce5_run"></a>
### 2.2 Running it

```sh
NP=32 bash run.sh
```

`run.sh` runs the input on `cpu` with `bin/main.out`, then with `backend` auto and cuda with `bin/main_gpu.out`. The
`cpu` run alone takes 24 s on 4 ranks of a desktop CPU.

<a id="a04_cce5_result"></a>
### 2.3 The result

![cCCE5_Hahn_backends: expected result](cCCE5_Hahn_backends/expected.png)

`auto` and `cuda` are within 1.0e-10 of `cpu`.

Next, change `gpu_options.min_dim`, which sets the clusters that go to the GPU under `auto`.

Input: [`ccein.jsonc`](cCCE5_Hahn_backends/ccein.jsonc)

<a id="a04_gcce4"></a>
## 3. When the GPU saves time

Folder: `gCCE4_Hahn_backends/`

<a id="a04_gcce4_input"></a>
### 3.1 The input

The gCCE4 Hahn echo of the NV, a spin 1 with its zero-field splitting, in four P1. gCCE evolves the density matrix of
the qubit and the bath on the GPU, another path than that of cCCE
([GPU execution](../../../../docs/methods/gpu.md#gpu_what)). The cluster of all four P1 has dimension 768 with the
qubit:

```jsonc
"method"        : "gCCE",
"configuration" : {"gyrofile": "../bath/gyro_P1", "bathfile": ["../bath/bath_p1_4"], "seed": 1},
"order"         : 4,
"execution"     : {"backend": "cpu", "gpu_options": {"min_dim": 256, "evolution_batch_steps": 64}, "profile": true}
```

The input has `backend` cpu, the reference. `run.sh` runs it, then copies of it with `backend` auto and cuda.

<a id="a04_gcce4_run"></a>
### 3.2 Running it

```sh
NP=32 bash run.sh
```

The `cpu` run alone takes 85 s on 4 ranks of a desktop CPU. On the RTX PRO 6000 with 32 ranks the three backends
took:

| `cpu` | `auto` | `cuda` |
|---|---|---|
| 33.8 s | 2.1 s | 4.7 s |

<a id="a04_gcce4_result"></a>
### 3.3 The result

![gCCE4_Hahn_backends: expected result](gCCE4_Hahn_backends/expected.png)

Both GPU curves are within 1.0e-10 of `cpu`. `auto` is the fastest, 16 times faster than `cpu`. `cuda` sends every
cluster to the GPU, the qubit alone included, and is slower than `auto` (4.7 s against 2.1 s).

Input: [`ccein.jsonc`](gCCE4_Hahn_backends/ccein.jsonc)
