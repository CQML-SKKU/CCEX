# 08_GPU

The CUDA backend. `execution.backend` chooses it: `cpu` (the default), `auto` (CUDA for every cluster whose Hilbert
dimension is at least `execution.gpu_options.min_dim`, the CPU for the rest, and the CPU everywhere in a binary built
without CUDA) or `cuda` (every cluster; refused when CUDA is not available). In the GPU build, `auto` and `cuda` with
cCCE also send clusters of bath dimension 2 and 4 to a packed kernel (`gpu_options.packed_max_dim` 4).
[GPU execution](../../../docs/methods/gpu.md) has the details and `test/gpu_verification/check_gpu/` the checks. Each folder has one input with
`backend` cpu; `run.sh` runs it with `bin/main.out`, then with `backend` auto and cuda (`set_input_keys.py`) with
`bin/main_gpu.out` when that exists (`bash do_compile.sh gpu` on the GPU server makes both; `CCEX_GPU_BIN` names another GPU
binary, e.g. a laptop card: `CUDA_HOME=/usr/local/cuda-13.2 CUDA_ARCH=86 bash do_compile.sh generic`). Without it
`auto` runs on the CPU and `cuda` is skipped. With many ranks on one GPU, start CUDA MPS in the job
(`test/gpu_verification/server/mps.sh`). Every key not written in an input has its default (`main.out -e`).

`bath/` is the P1 bath of `test/gpu_verification/check_gpu/examples/`: five P1 (`bath_p1_5`, electron + ¹⁵N each), the
last four of them (`bath_p1_4`), the NV position (`qubit_p1`) and the P1 gyromagnetic ratio (`gyro_P1`).

## cCCE5_Hahn_backends

- **What it shows:** cCCE5 Hahn echo of the NV in the five P1: one cluster of dimension 1024, the case the GPU is for
  (about a minute of CPU for its 100 time steps); the same run per backend, the `[profile]` timings and the curves.
- **Keys that matter:** `execution.backend`, `execution.gpu_options.min_dim` 256, `execution.profile` true, `order` 5,
  the `Defect` block (¹⁵N, four Jahn-Teller axes).
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 24 s on 4 ranks, WSL).
- **Expected result:** no `expected.png` yet (the run needs a GPU node); the backends agree to about 1e-10.
- **What to change next:** `gpu_options.min_dim` (which clusters go to the GPU).

Input: [`ccein.jsonc`](cCCE5_Hahn_backends/ccein.jsonc)

## cCCE5_Hahn_autoOptions

- **What it shows:** the input of `cCCE5_Hahn_backends` with `backend` `auto` and no `gpu_options`: for cCCE with
  `auto` or `cuda`, `bin/main_gpu.out` sets the keys an input leaves out at the start and says so in a `[gpu-auto]`
  line, which `run.sh` prints: `ranks_per_gpu` 4 when a node has more than 4 ranks a GPU (`auto`, `cluster.mode`
  standard), `min_dim` 256 with it, `evolution_batch_size` 32 and `cluster_batch_mb` two thirds of the smallest card
  over the ranks that open it ([Execution tags](../../../docs/user/InputFile_Backend.md#tag_backend)). On an RTX PRO
  6000 with 32 ranks these values were as fast as the best ones chosen by hand on a larger bath (487 s against
  490 s, `CHANGELOG.md` 2.1.0). The curve is compared with the same input on `cpu`.
- **Keys that matter:** `execution.backend` `auto`, `execution.gpu_options` left out (write a key to fix it).
- **How to run:** on a GPU node, `NP=32 bash run.sh`; without `bin/main_gpu.out` it runs `auto` on the CPU and
  says so.
- **Expected result:** no `expected.png` yet (the run needs a GPU node); `auto` within about 1e-10 of `cpu`.
- **What to change next:** `"gpu_options": {"ranks_per_gpu": 0}` (every rank may open the GPU) to compare.

Input: [`ccein.jsonc`](cCCE5_Hahn_autoOptions/ccein.jsonc)

## gCCE4_Hahn_backends

- **What it shows:** gCCE4 Hahn echo of the NV (spin 1 with its zero-field splitting) in four P1: the cluster of all
  four has dimension 768 with the qubit. Measured on a GPU server (32 ranks + RTX PRO 6000, CUDA 13.2, MPS, median of
  three runs): cpu 33.8 s, auto 2.1 s, cuda 4.7 s (5.8 s without MPS), both GPU curves within 1.0e-10 of cpu; WSL,
  4 ranks + RTX 3050: 40.9 / 16.8 / 16.8 s.
- **Keys that matter:** `method` gCCE, `execution.backend`, `gpu_options.min_dim` 256, `evolution_batch_steps` 64.
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 85 s on 4 ranks, WSL).
- **Expected result:** no `expected.png` yet (the run needs a GPU node).
- **What to change next:** `cuda` sends every cluster, the qubit alone included, to the GPU and was slower than `auto`
  here (4.7 s against 2.1 s).

Input: [`ccein.jsonc`](gCCE4_Hahn_backends/ccein.jsonc)

## pCCE2_Hahn_backends

- **What it shows:** pCCE2 (partitions of two P1) on the five P1: the pair of partitions is one cluster of dimension
  256, on the GPU under `auto` too (`min_dim` 256). Measured on the same server (32 ranks, MPS): cpu 0.22 s,
  auto 0.29 s, cuda 2.8 s (3.8 s without MPS), difference 0.
- **Keys that matter:** `cluster.mode` pcce, `cluster.pcce.partition_size` 2, `execution.backend`.
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 1 s on 4 ranks, WSL).
- **Expected result:** no `expected.png` yet (the run needs a GPU node).
- **What to change next:** `partition_size` 3 (larger clusters, more of them on the GPU).

Input: [`ccein.jsonc`](pCCE2_Hahn_backends/ccein.jsonc)
