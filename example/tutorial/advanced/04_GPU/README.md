# Advanced 04: The GPU backend

Folder: `example/tutorial/advanced/04_GPU/`

The CUDA backend runs the clusters on the GPU. `execution.backend` chooses where:

- `cpu` (the default) runs every cluster on the CPU.
- `auto` runs on the GPU every cluster whose Hilbert dimension is at least `execution.gpu_options.min_dim`, the rest on
  the CPU, and everything on the CPU in a binary built without CUDA.
- `cuda` runs every cluster on the GPU, and is refused when CUDA is not available.

In the GPU build, `auto` and `cuda` with cCCE also send the clusters of bath dimension 2 and 4 to a packed kernel
(`gpu_options.packed_max_dim` 4). [GPU execution](../../../../docs/methods/gpu.md) has the details, and
`test/gpu_verification/check_gpu/` the checks.

Each folder has one input with `backend` cpu. `run.sh` runs it with `bin/main.out`, and then with `backend` auto and
cuda (`set_input_keys.py`) with `bin/main_gpu.out` when that exists. `bash do_compile.sh gpu` makes both binaries, and
`CCEX_GPU_BIN` names another GPU binary, for example one for a laptop card built with
`CUDA_HOME=/usr/local/cuda-13.2 CUDA_ARCH=86 bash do_compile.sh generic`. Without it `auto` runs on the CPU and `cuda`
is skipped. With many ranks on one GPU, start CUDA MPS in the job (`test/gpu_verification/server/mps.sh`). Every key
not written in an input has its default (`main.out -e`).

`bath/` is the P1 bath of `test/gpu_verification/check_gpu/examples/`. It holds five P1 (`bath_p1_5`, an electron and
a ¹⁵N each), the last four of them (`bath_p1_4`), the NV position (`qubit_p1`) and the P1 gyromagnetic ratio
(`gyro_P1`).

The expected figures were made on an RTX PRO 6000 with 32 ranks and MPS, with the binaries of `do_compile.sh gpu`.

1. [cCCE5_Hahn_backends](#a04_cce5): one large cCCE cluster on each backend.
2. [cCCE5_Hahn_autoOptions](#a04_auto): the `gpu_options` chosen at the start.
3. [gCCE4_Hahn_backends](#a04_gcce4): gCCE on each backend.
4. [pCCE2_Hahn_backends](#a04_pcce2): pCCE on each backend.

<a id="a04_cce5"></a>
## 1. cCCE5_Hahn_backends: one large cCCE cluster on each backend

- **What it shows:** the cCCE5 Hahn echo of the NV in the five P1. One cluster has dimension 1024, the case the GPU is
  for, and takes about a minute of CPU for its 100 time steps. The same run is made on each backend, with the
  `[profile]` timings and the curves.
- **Keys that matter:** `execution.backend`, `execution.gpu_options.min_dim` 256, `execution.profile` true, `order` 5,
  and the `Defect` block (¹⁵N, four Jahn-Teller axes).
- **How to run:** on a GPU node, `NP=32 bash run.sh`. The CPU part alone takes 24 s on 4 ranks of a desktop CPU.
- **Expected result:** `expected.png`. `auto` and `cuda` are within 1.0e-10 of `cpu`. The three runs and the plot took
  200 s.
- **What to change next:** `gpu_options.min_dim`, which sets the clusters that go to the GPU.

![cCCE5_Hahn_backends: expected result](cCCE5_Hahn_backends/expected.png)

Input: [`ccein.jsonc`](cCCE5_Hahn_backends/ccein.jsonc)

<a id="a04_auto"></a>
## 2. cCCE5_Hahn_autoOptions: the gpu_options chosen at the start

- **What it shows:** the input of `cCCE5_Hahn_backends` with `backend` `auto` and no `gpu_options`. For cCCE with
  `auto` or `cuda`, `bin/main_gpu.out` sets the keys an input leaves out at the start and reports them in a
  `[gpu-auto]` line, which `run.sh` prints ([Execution tags](../../../../docs/user/InputFile_Backend.md#tag_backend)):
  - `ranks_per_gpu` 4, with `auto` and more than 4 ranks per GPU on a node.
  - `min_dim` 256, with `ranks_per_gpu`.
  - `evolution_batch_size` 32.
  - `cluster_batch_mb`, two thirds of the smallest card over the ranks that open it.

  On an RTX PRO 6000 with 32 ranks these values were as fast as the best ones chosen by hand on a larger bath.
  The curve is compared with the same input on `cpu`.
- **Keys that matter:** `execution.backend` `auto` and `execution.gpu_options` left out. Writing a key fixes it.
- **How to run:** on a GPU node, `NP=32 bash run.sh`. Without `bin/main_gpu.out` it runs `auto` on the CPU and says
  so.
- **Expected result:** `expected.png`, with the log line
  `[gpu-auto] ranks_per_gpu=4 min_dim=256 evolution_batch_size=32 cluster_batch_mb=16208`. `auto` is within 1.0e-10
  of `cpu`, and the two runs took 82 s.
- **What to change next:** `"gpu_options": {"ranks_per_gpu": 0}`, which lets every rank open the GPU, to compare.

![cCCE5_Hahn_autoOptions: expected result](cCCE5_Hahn_autoOptions/expected.png)

Input: [`ccein.jsonc`](cCCE5_Hahn_autoOptions/ccein.jsonc)

<a id="a04_gcce4"></a>
## 3. gCCE4_Hahn_backends: gCCE on each backend

- **What it shows:** the gCCE4 Hahn echo of the NV (a spin 1 with its zero-field splitting) in four P1. The cluster of
  all four has dimension 768 with the qubit.

  | Machine | `cpu` | `auto` | `cuda` |
  |---|---|---|---|
  | RTX PRO 6000, 32 ranks, MPS (median of three runs) | 33.8 s | 2.1 s | 4.7 s (5.8 s without MPS) |
  | RTX 3050, 4 ranks | 40.9 s | 16.8 s | 16.8 s |

  Both GPU curves are within 1.0e-10 of `cpu`.
- **Keys that matter:** `method` gCCE, `execution.backend`, `gpu_options.min_dim` 256 and `evolution_batch_steps` 64.
- **How to run:** on a GPU node, `NP=32 bash run.sh`. The CPU part alone takes 85 s on 4 ranks of a desktop CPU.
- **Expected result:** `expected.png`. `auto` and `cuda` are within 1.0e-10 of `cpu`, and the three runs took 45 s.
- **What to change next:** compare `auto` with `cuda`. `cuda` sends every cluster to the GPU, the qubit alone
  included, and was slower than `auto` here (4.7 s against 2.1 s).

![gCCE4_Hahn_backends: expected result](gCCE4_Hahn_backends/expected.png)

Input: [`ccein.jsonc`](gCCE4_Hahn_backends/ccein.jsonc)

<a id="a04_pcce2"></a>
## 4. pCCE2_Hahn_backends: pCCE on each backend

- **What it shows:** pCCE2 (partitions of two P1) on the five P1. The pair of partitions is one cluster of dimension
  256, so it goes to the GPU under `auto` too (`min_dim` 256). On the RTX PRO 6000 with 32 ranks and MPS it took
  0.22 s on `cpu`, 0.29 s on `auto` and 2.8 s on `cuda` (3.8 s without MPS), with no difference in the curves.
- **Keys that matter:** `cluster.mode` pcce, `cluster.pcce.partition_size` 2 and `execution.backend`.
- **How to run:** on a GPU node, `NP=32 bash run.sh`. The CPU part alone takes 1 s on 4 ranks of a desktop CPU.
- **Expected result:** `expected.png`. `auto` and `cuda` give the `cpu` curve exactly, and the three runs took 17 s.
- **What to change next:** `partition_size` 3, for larger clusters and more of them on the GPU.

![pCCE2_Hahn_backends: expected result](pCCE2_Hahn_backends/expected.png)

Input: [`ccein.jsonc`](pCCE2_Hahn_backends/ccein.jsonc)
