# GPU execution

CCEX has an optional NVIDIA CUDA backend. It runs in the GPU binary `bin/main_gpu.out`, and the CPU binary
`bin/main.out` runs every cluster on the CPU. The keys are in [Execution tags](../user/InputFile_Backend.md), and
[Advanced 11: The GPU backend](../../example/tutorial/advanced/11_GPU/README.md) runs it.

1. [What runs on the GPU](#gpu_what): the methods and what stays on the CPU.
2. [Building](#gpu_build): the GPU binary.
3. [Choosing the backend](#gpu_backend): `cpu`, `auto` and `cuda`, and what happens when the GPU fails.
4. [The options chosen at the start](#gpu_auto): the automatic `gpu_options` of cCCE.
5. [Kernels for small clusters](#gpu_small): the packed and the warp kernel.
6. [Several ranks and GPUs](#gpu_ranks): the device of each rank, MPS and `ranks_per_gpu`.
7. [Memory](#gpu_memory): what a cluster batch holds.
8. [How the GPU computes](#gpu_how): the cCCE and gCCE paths and their accuracy.

<a id="gpu_what"></a>
## 1. What runs on the GPU

| Method | Ensemble (no `configuration.s_range`) | Sampled states (`configuration.s_range`) |
|---|---|---|
| cCCE, `cluster.mode` standard, pcce or bundle | The normalized bath trace of each cluster. | The two conditional state vectors. |
| gCCE with `propagator` eigen, `cluster.mode` standard or bundle | The density-matrix evolution. | The vector or density-matrix evolution. |
| noise with `backend` `cuda` | The sums over the bath's modes. | The same. |

Everything else runs on the CPU: gCCE with `propagator` expm, meCCE and megCCE, pseudospin, DEER bath pulses, and the
model spectra of noise. The CPU also finds the clusters and the pCCE partitions, builds every Hamiltonian, and forms
the product over the clusters. The GPU diagonalizes the cluster Hamiltonians and evolves them in time.

<a id="gpu_build"></a>
## 2. Building

`bash do_compile.sh gpu` builds both binaries with one `make`, from the same sources into separate object folders.
The GPU half needs a CUDA toolkit with nvcc, cuSOLVER, cuBLAS and cudart. `CUDA_HOME` gives the toolkit and `CUDA_ARCH`
the compute capability without the decimal point (for example `86` for an RTX 3050):

```bash
CUDA_HOME=/usr/local/cuda CUDA_ARCH=86 bash do_compile.sh gpu
```

After changing only `CUDA_HOME` or `CUDA_ARCH`, run `make clean` first. Keep the toolkit's shared libraries available
when the binary runs. The build options, CMake and the CUDA versions are in
[Building and running](../user/BuildingAndRunning.md#BuildCuda).

<a id="gpu_backend"></a>
## 3. Choosing the backend

| `execution.backend` | Where the clusters are computed |
|---|---|
| `cpu` (default) | All on the CPU. |
| `auto` | The clusters of dimension `gpu_options.min_dim` and above on the GPU, the rest on the CPU. A run the GPU cannot take stays on the CPU. |
| `cuda` | All on the GPU. A run the GPU cannot take is refused when the input is read. |

The dimension counts every spin of the cluster, the on-site spins of a defect and, with gCCE, the qubit, so it is
larger than the order suggests. `min_dim` is a starting point and does not guarantee a speedup. The `[backend]` lines
of the log say what was chosen, and `execution.profile` prints the times per backend and dimension to compare.

`auto` also runs in `bin/main.out` and when CUDA cannot be initialized, and then every cluster is computed on the CPU.
Once the GPU has started, an error of the solver or of a kernel stops the run. A cluster is never computed again on
the CPU after a partial GPU result, with one exception. When the GPU runs out of memory under `auto`, the cluster is
freed from the card and computed on the CPU from the start, and that rank sends no later cluster of that dimension or
larger to the GPU. The log has one `[gpu-memory-fallback]` line per rank and dimension. With `cuda`, running out of
memory stops the run.

<a id="gpu_auto"></a>
## 4. The options chosen at the start

For cCCE (`cluster.mode` standard or pcce) with `backend` `auto` or `cuda`, the `gpu_options` that the input leaves out
are chosen at the start from the number of GPUs, their memory and the ranks on the node. The log prints them in one
`[gpu-auto]` line.

| Option | Chosen value |
|---|---|
| `ranks_per_gpu` | 4, with `backend` `auto` and more than four ranks per GPU on the node. |
| `min_dim` | 256 when `ranks_per_gpu` is set. |
| `evolution_batch_size` | 32. |
| `cluster_batch_mb` | Two thirds of the smallest card's memory, shared by the ranks that use it, at least 64. |

A value written in the input is kept. gCCE uses the fixed defaults of [gpu_options](../user/InputFile_Backend.md#tag_gpu_options).

<a id="gpu_small"></a>
## 5. Kernels for small clusters

A cluster of the large-cluster path has a fixed cost of its own solver calls, copies and synchronizations. A P1 bath
has millions of clusters of bath dimension 2 and 4 per bath state, for which this cost is far larger than the CPU's
time. Two kernels take such clusters many at a time. They share one CUDA context per node. Every rank builds its
clusters' Hamiltonians on the CPU and puts them into node shared memory, and node rank 0 sends them to the first GPU of
the node and writes the coherences back. Until the GPU is ready, the ranks compute their clusters on the CPU.

<a id="the-packed-kernel-for-small-ccce-clusters-gpu_optionspacked_max_dim"></a>
**The packed kernel** (`gpu_options.packed_max_dim`, 0, 2 or 4) runs one cCCE cluster per thread, for clusters of bath
dimension up to the value. Its default is 4 in `bin/main_gpu.out` with `backend` `auto` or `cuda` and cCCE standard or
bundle with `quantity` coherence, and 0 otherwise. It takes at most 64 pulses with at most four distinct interval
lengths, and both values of `evolution_update`.

<a id="the-warp-kernel-for-clusters-up-to-dimension-32-gpu_optionswarp_max_dim"></a>
**The warp kernel** (`gpu_options.warp_max_dim`, 0 to 32) runs one cluster per warp, for the gCCE clusters of an NV-like
qubit (dimension 3 times bath dimension 2, 3, 4, 6, 7, 8 or 9) and the cCCE clusters of bath dimension 8 and 16. Its
default is 27 in `bin/main_gpu.out` with `backend` `auto` or `cuda` and cCCE (`quantity` coherence) or gCCE, and 0
otherwise. It switches itself off, with one `[gpu-warp] off: ...` line, for `gcce_refinement`, `pulse.per_qubit`, a
pulse detuning, `savemode` `cluster_density_matrix`, a qubit of another dimension or more than one qubit, and more than
four distinct pulse intervals.

Both kernels send a cluster back to the CPU when its eigensystem does not converge or is not finite, or when its
coherence is not finite. A GPU that cannot be opened turns them off under `auto`, and stops the run with `cuda` or with
an explicit `packed_max_dim` under `cpu`. With either kernel on, the large-cluster batching of `evolution_batch_size`
is off.

<a id="gpu_ranks"></a>
## 6. Several ranks and GPUs

```bash
mpirun -n 4 ./bin/main_gpu.out -f input.json
```

Each rank takes the GPU `node-local rank % visible GPUs`. When the scheduler gives each rank one GPU, every rank uses
its device 0. The `[backend]` log records the mapping. The GPU returns scalar coherences (cCCE) or small reduced
density matrices (gCCE) to the host before the products and the MPI reductions, so CUDA-aware MPI is not needed.

When many ranks share one GPU, run the job under CUDA MPS, with `nvidia-cuda-mps-control -d` at its start and
`echo quit | nvidia-cuda-mps-control` at its end.

`gpu_options.ranks_per_gpu` lets only that many ranks per GPU of a node use the card. Written in the input it needs
cCCE standard, and chosen at the start it also applies to pcce. The other ranks
stay on the CPU, and the clusters are split between the two kinds of rank by their measured cost. The `[gpu-split]` log
line gives how many ranks use the GPU. With 0 (the default unless chosen at the start) every rank may use it.

<a id="gpu_memory"></a>
## 7. Memory

The GPU memory of one time batch of `evolution_batch_steps` points is about:

- gCCE with `evolution` matrix: two complex-double buffers of 16 D² × `evolution_batch_steps` bytes;
- cCCE ensemble: three such buffers;
- a pure state: buffers of D × `evolution_batch_steps` elements.

The eigensystems, the solver workspace and the cache of `cache_mb` come on top, per rank. `cluster_batch_mb` limits
what one batch of clusters may hold, but it is an estimate and not a cap on the card's memory. When memory is tight,
lower `evolution_batch_steps`, `cluster_batch_mb` or `cache_mb`, or use fewer ranks per GPU. gCCE `savemode`
`cluster_density_matrix` needs extra density-matrix storage on the card and on the host.

<a id="gpu_how"></a>
## 8. How the GPU computes

**cCCE.** The CPU builds the two conditional Hamiltonians of each cluster. cuSOLVER diagonalizes them, for a batch of
clusters of one dimension in one call (one cluster alone in its own call). The GPU forms the overlap of the two
eigenbases and evolves many time points as the columns of one matrix with cuBLAS. A reduction returns the overlap of
the conditional states (sampled state) or its normalized trace (ensemble). The sampled-state path evolves two state
vectors instead of the full propagators of the CPU path.

**gCCE.** The CPU builds each Hamiltonian and pulse operator. cuSOLVER diagonalizes the Hamiltonian once per cluster,
and the eigensystem, the states and the pulse data stay on the GPU through the time grid. cuBLAS evolves a batch of
time points, and kernels apply the phases, the qubit pulses and the partial trace over the bath. With
`gcce_refinement` on, the eigensystem is refined on the GPU as on the CPU.

All the arithmetic is complex double, without reduced precision or fast-math. The coherence products keep a separate
exponent as on the CPU. The GPU checks hold the GPU results within 10⁻⁸ + 2×10⁻⁶·|L| of the CPU's. Where the product
diverges (|L| > 1), the two can differ more, since the rounding of the eigensolver is amplified there.

`test/gpu_verification/` holds the GPU checks and the server jobs that run them.
