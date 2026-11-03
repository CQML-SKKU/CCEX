# GPU execution

CCEX has an optional NVIDIA CUDA backend for conventional CCE, pCCE and gCCE, in both ensemble and single-sample
calculations:

- cCCE: `method=cCCE`;
- pCCE: `method=cCCE` with `cluster.mode=pcce`;
- gCCE: `method=gCCE` with `propagator=eigen`.

The CPU backend is the default. gCCE with `propagator=expm` and all other methods run on the CPU only.

The keys on this page (`backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `gcce_refinement`,
`profile`) belong to the `execution` block, so `backend=auto` means `"execution": {"backend": "auto"}`. The other
short names are:

| Short name | Full key |
| --- | --- |
| `pcce.*` | `cluster.pcce.*` |
| `savemode` | `output.savemode` |
| `probe_memory_gb` | `cluster.adaptive.probe_memory_gb` |

What the GPU computes for each method:

| Method | Ensemble (no `configuration.s_range`) | Single (`configuration.s_range`) |
| --- | --- | --- |
| Conventional CCE, also called cCCE | GPU normalized bath trace | GPU conditional state vectors |
| Partitioned CCE, `cluster.mode=pcce` | GPU normalized bath trace | GPU conditional state vectors |
| gCCE, `propagator=eigen` | GPU density-matrix evolution | GPU vector or density-matrix evolution |
| Bundle clustering, `cluster.mode=bundle` (cCCE, or gCCE with `propagator=eigen`) | as cCCE or gCCE above, per cluster of bundles | as cCCE or gCCE above |

For cCCE and pCCE, `propagator` and `evolution` do not select an algorithm. cCCE diagonalizes the two conditional
Hamiltonians, and its GPU single path uses state vectors. `evolution` is a gCCE option: for a cCCE input, leave it at
the default or set `matrix`, because the parser rejects an explicit `evolution=vector` outside gCCE.

pCCE uses the cCCE coherence kernels after the bath spins are partitioned. The partition search, the cluster
construction and the Hamiltonian construction stay on the CPU. CUDA does the diagonalization of the conditional
Hamiltonians and the time evolution. pCCE profiles are labelled `method=pcce`.

Coherence products keep a separate exponent (see Implementation below), which avoids an intermediate
overflow when the final product is representable. It says nothing about the convergence of the CCE.

The adaptive refinement ([cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive)) uses the same CPU and CUDA
engines to test and improve groups of variable size.

## Build

`do_compile.sh gpu` builds two executables with one `make`, from the same sources into separate object directories:

- `bin/main.out`: CPU only, what every other preset builds;
- `bin/main_gpu.out`: compiled with `-DCCEX_USE_CUDA` and linked with the CUDA kernels.

The build uses MPI and the bundled Eigen headers. It does not require Intel MKL. The CUDA half needs a toolkit with
nvcc, cuSOLVER, cuBLAS and cudart.

Tested toolchains:

| CUDA | GCC | GPU |
| --- | --- | --- |
| 12.8.1 | 13.3 | RTX 3050 |
| 13.2.2 | 11.5 | idm4-gpu |

On idm4-gpu (compute capability 12.0), CUDA 12.8 miscompiles the pulse kernel `applyPulse` (see `CHANGELOG.md`), so
the `gpu` preset uses CUDA 13.2.

```bash
bash do_compile.sh gpu              # idm4-gpu: loads gcc, intel/oneapi-2023.2 (mpiicpc with icpx, MKL) and CUDA 13.2 itself
CCEX_JOBS=4 bash do_compile.sh gpu  # fewer parallel compiles on a shared node
```

On a machine without those modules, put your MPI wrapper (`mpicxx`) on `PATH` and give the toolkit location and the
card:

```bash
CUDA_HOME=/usr/local/cuda CUDA_ARCH=86 bash do_compile.sh gpu
```

`CUDA_ARCH` is the CUDA compute capability without the decimal point:

- `120`: the RTX PRO 6000 Blackwell on idm4-gpu (the default);
- `86`: the RTX 3050.

The binary embeds both device code and PTX for that architecture. Choose an architecture that both your GPU and your
toolkit support. Switching preset cleans the objects before building. After changing only `CUDA_HOME` or
`CUDA_ARCH`, run `make clean` first.

Keep the toolkit's shared libraries available when running the binary. On idm4-gpu, load the same three modules in the
job script (`gcc/gcc-11.5.0 intel/oneapi-2023.2 cuda/cuda-13.2`). Without them the binary stops at start-up on the
system's older libstdc++.

gCCE eigen propagation can refine the eigenpairs, on CPU and CUDA, to resolve small bath energies next to large
central-spin energies. Writing the `gcce_refinement` block turns it on; it is off by default. On the NV–¹³C gCCE4 case
it was made for, it changed the error-corrected final coherence by at most 6.7e-7, while it took most of a cluster's
setup time. The raw `wiDiv` diverges at some points whatever the setting (`../user/InputFile_Backend.md`).

When refinement is on, it is strict: a cluster that does not converge within `gcce_refinement.max_sweeps` (default 32)
stops the run, with the dimension and the sweep limit in the message. It never accepts the unrefined eigensystem
silently. Compensated phase products and scaled coherence accumulation are active either way. These options apply
only to gCCE with `propagator=eigen`.

When `configuration.seed` is omitted, MPI rank 0 generates a seed, broadcasts it to every rank and prints it for
replay.

pCCE prints at most N trial details from rank 0 per partition search with `pcce.print_trials=N` (default 0, none).
The global partition summary still aggregates every rank.

`cluster.pcce.partition` selects the tie rule of the search:

- `deterministic`: larger silhouette throughout;
- `legacy`: the original asymmetric rule, larger silhouette within a rank and smaller silhouette between ranks at
  equal inertia.

The C RNG is initialized once in `main`, from the resolved seed. It drives:

- the sampled bath and subspin states;
- random defect axes;
- legacy pCCE searches.

Deterministic pCCE partition trials use their own seeded RNG. Fixed-bath ensemble calculations without random defect
axes or a legacy partition search do not use the C RNG. Partition cache contents do not depend on it either, because
deterministic partitions do not use it. Raw coherence caches must match the executing binary and the sampled states.
The regression checks different seeds, an independently generated C RNG state sequence, CPU/CUDA agreement and MPI
replay.

Cluster batch admission accounts for the difference between state vectors and density matrices. The budget is a
conservative estimate. It does not guarantee the solver workspace or the total device memory.

Device pointer tables and batch output buffers persist with the first work object, and unchanged tables are not
uploaded again. A change in membership, buffer address or density time count refreshes the corresponding tables. With
profiling, `[gpu-memory] pointer_uploads` counts the table uploads.

## Input and execution

Add this block to an existing CCE, pCCE or gCCE input:

```json
{
    "execution": {
        "backend": "auto",
        "gpu_options": {
            "min_dim": 768,
            "evolution_batch_steps": 64,
            "evolution_batch_size": 8,
            "cluster_batch_mb": 64,
            "cache_mb": 64
        },
        "propagator": "eigen",
        "profile": true
    }
}
```

This fragment is not a complete input. Complete examples for all four combinations, run from the repository root:

```bash
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_cce5_ensemble.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_cce5_single.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_gcce4.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_gcce4_single.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/pcce_ensemble.json
# one state: a copy of pcce_ensemble.json with "configuration.s_range": [1, 1] and "output.file": "pcce_gpu_single"
./bin/main_gpu.out -f pcce_single.json -o bin
```

The cCCE examples use five P1 sites (bath dimension 1024). The gCCE examples use four (qubit-plus-bath dimension 768).
Each writes its own output names under `bin/`. To compare the same inputs between CPU and GPU, use the regression
suite below.

The full list of these keys, with types and defaults, is in [the backend tags](../user/InputFile_Backend.md) and, for
pCCE, [the cluster tags](../user/InputFile_Cluster.md). Summary:

| Option | Behavior |
| --- | --- |
| `execution.backend=cpu` | Default; no GPU computation for the clusters the packed kernel does not take (with `packed_max_dim` > 0, given explicitly, node rank 0 still opens CUDA for the packed clusters) |
| `execution.backend=cuda` | Require CUDA for every CCE/pCCE/gCCE cluster the packed kernel does not take; fail if unsupported or unavailable |
| `execution.backend=auto` | Of the clusters the packed kernel does not take, use CUDA for supported clusters with dimension at least `gpu_options.min_dim`; otherwise use CPU |
| `execution.gpu_options.min_dim` | Positive integer, default 768 (the GPU wins from about 670-750 on, every path); applies only to `auto` |
| `execution.gpu_options.packed_max_dim` | 0 (off), 2 or 4; default 4 in `main_gpu.out` with `backend` `cuda`/`auto` and CCE, else 0: cCCE clusters of bath dimension up to it go to [the packed kernel](#the-packed-kernel-for-small-ccce-clusters-gpu_optionspacked_max_dim) |
| `execution.gpu_options.warp_max_dim` | 0 (off) to 32; default 27 in `main_gpu.out` with `backend` `cuda`/`auto` and CCE or gCCE, else 0: clusters of a compiled total dimension up to it go to [the warp kernel](#the-warp-kernel-for-clusters-up-to-dimension-32-gpu_optionswarp_max_dim) |
| `execution.gpu_options.evolution_batch_steps` | Positive integer, default 64; maximum time points processed together |
| `execution.gpu_options.evolution_batch_size` | Positive integer, default 8; maximum consecutive compatible clusters evolved together; 1 uses the sequential cluster path |
| `execution.gpu_options.cluster_batch_mb` | Positive integer, default 64 MiB per rank; estimated active evolution-buffer budget used when grouping clusters |
| `execution.gpu_options.cache_mb` | Nonnegative integer, default 64 MiB per rank; reusable free GPU allocations; 0 disables the cache |
| `execution.propagator` | gCCE: `eigen` (default; one eigendecomposition per cluster, required for CUDA) or `expm` (a matrix exponential per step, CPU only) |
| `execution.evolution` | gCCE: `vector` (pure initial state, `s_range`, `propagator=eigen`) or `matrix` (density matrix); default `vector` where it can run, else `matrix` |
| `execution.gcce_refinement` | Block, off by default; `enabled` (true when the block is given) turns on the compensated refinement of the eigensystem, CPU and CUDA (below) |
| `execution.gcce_refinement.max_sweeps` | Positive integer, default 32; a cluster not converged within it stops the run |
| `cluster.adaptive.probe_memory_gb` | Nonnegative, default 0 (automatic: 80 % of the job's node memory and, when clusters may go to the GPU, of the smallest GPU); memory the adaptive probes running at once on a node may hold |
| `cluster.pcce.print_trials` | pCCE: rank 0 prints the details of at most this many partition trials (default 0) |
| `cluster.pcce.connections` | pCCE: `centers` (default; distance between partition centroids) or `spins` (any spin pair within `rdip`) |
| `cluster.pcce.partition` | `deterministic` (default) or `legacy`; pCCE partition RNG and trial distribution policy |
| `cluster.pcce.cache.mode` | `off` (default), `readwrite` (reuse or compute/store), `read` (require a matching stored partition); deterministic pCCE only |
| `cluster.pcce.cache.file` | Explicit partition cache path, relative to the working directory or absolute; parent directory must exist |
| `execution.profile` | Boolean, default false; emit timings grouped by rank, method, backend and Hilbert dimension |

A fresh pCCE search evaluates distances only to centers whose groups still have room. It keeps the spin order, the
ascending tie-breaking on centers, the trial RNG and the convergence criteria. This is a CPU optimization: it applies to
all three backends, with the cache off or on, and removes the N-by-K distance table from the search. Partitioning
stays on the CPU.

A fresh search also detects exactly repeated center/assignment states and skips whole periods. It then executes the
remaining iterations, so that it ends at the same `pcce.max_iterations` endpoint. It does not declare a cycle
converged and does not change the trial selection. This applies to deterministic and legacy pCCE on every backend.

The root `[partition-search]` log counts converged, cycling and limit-only trials, then the logical, executed and
skipped iterations summed over ranks. `pcce.print_trials=true` adds each trial's period and counts.

The partition cache algorithm identity is `pcce-deterministic-20260912-v3`. With `readwrite`, a cache file of the old
v1/v2 identities is recomputed once; `read` rejects it until it is regenerated. The on-disk cache format is version 1.

CUDA initialization happens in two places:

- The packed kernel's context (below) is opened by node rank 0 at the start whenever `packed_max_dim` > 0. In the GPU
  build that is the default for CCE with `cuda`/`auto`.
- The per-rank initialization of the large-cluster path. `auto` defers it until a rank meets its first eligible
  cluster, so a rank with only smaller clusters never initializes it. Ranks with no eligible cluster do not wait for
  the ones that initialize; only the node-rank metadata discovery is collective. `backend=cuda` initializes on every
  rank at start-up.

`auto` also works in a CPU-only executable and falls back to the CPU when CUDA cannot be initialized. Once GPU
computation has started, a solver or execution failure aborts the MPI job with an error. CCEX does not rerun a
partially completed cluster silently.

Running out of GPU memory is the exception under `auto`. (Without it, 20 ranks holding dimension-2197 ensemble clusters
on one 96 GB card at once aborted the run and the MPS server with it.) On an out-of-memory error:

- what the cluster had on the card is freed;
- the cluster is computed on the CPU from the start;
- that rank sends no later cluster of that dimension or larger to the GPU.

The log has one `[gpu-memory-fallback] rank=R dim=D` line per rank and dimension, and a count per rank at the end. The
curves for those clusters are the CPU path's. With `cuda`, an out-of-memory error still stops the run.

`CCEX_TEST_GPU_OOM_BYTES=<n>` is a test hook: it fails every device allocation of n bytes or more, as a full card
would. `check_options.py` uses it for cCCE and gCCE, one cluster at a time and in cluster batches.

An explicit `cuda` request with gCCE `expm` or with an unsupported method is rejected once the input is read.

`cluster.mode=bundle` adds automatic grouping of variable size to `method=cCCE` and `method=gCCE`, with the selected
CPU, CUDA or auto engine. The dimension cap is `cluster.bundle.max_dim` (required). Bundle orders of 3 and more are
supported. The refinement on top of it is in [cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive).

The dimension threshold is a starting point. It does not guarantee a speedup. Check the profile for your spin species,
defect spins, state mode, time grid and GPU. The matrix dimension includes all attached defect spins and, for gCCE, the
qubit, so it is larger than the CCE order suggests.

GPU memory use, per time batch:

- gCCE matrix evolution: two main complex-double buffers of `16 * D² * batch_steps` bytes each;
- cCCE ensemble: three such buffers;
- pure-state evolution: buffers of `D * batch_steps` elements each.

Eigensystems, solver workspace and other buffers come on top. gCCE `savemode` `cluster_density_matrix` needs extra
density-matrix storage on the device and output storage on the host. If memory is tight, reduce the time batch size.

For pCCE, `order` counts groups and `pcce.partition_size` is the number of bath sites per group. The largest cluster
can hold `order * partition_size` sites plus the attached defect spins. This is a different set of clusters from all
cCCE clusters of that spin order.

The partitioner drops the final `nspin % partition_size` sites, identically on CPU and CUDA. Inputs must have:

- a positive integer `pcce.partition_size` no greater than the loaded bath size;
- positive `pcce.trials` and `pcce.max_iterations`.

With `cluster.pcce.partition` = `deterministic`:

- `pcce.trials` is the exact global trial count;
- trial `t` has its own RNG stream derived from `(seed,t)` and runs on rank `t % ranks`;
- the search does not consume the bath-state C RNG;
- equal candidates are ordered by inertia, then silhouette, then global trial index;
- with an explicit seed and the same executable and input, changing the MPI size keeps the candidate set and the
  selected partition.

Floating-point products can still differ in their last bits, because the MPI reduction order changes. Bitwise identity
across compilers or architectures is not promised.

The deterministic default changes the random search and the sampled bath states relative to old pCCE runs. To
reproduce an existing result, set `cluster.pcce.partition` to `legacy` and keep the original seed and MPI size. Each
rank then uses `max(1, pcce.trials / ranks)` trials and the shared C RNG, as before. Legacy mode keeps its dependence
on the MPI size and its possibly repeated search streams.

<a id="the-packed-kernel-for-small-ccce-clusters-gpu_optionspacked_max_dim"></a>
## The packed kernel for small cCCE clusters (gpu_options.packed_max_dim)

The paths above give each cluster its own eigensolver calls, copies and synchronizations. A P1 bath has millions of
clusters of bath dimension 2 and 4 per bath state. For them this fixed cost is about 400 times the CPU's time per
cluster, and batching the time evolution (`evolution_batch_size`) does not remove it.
`"execution": {"gpu_options": {"packed_max_dim": 4}}` (or 2) moves those clusters to the packed kernel.

One thread per cluster (`src/calculation/gpu/small_cluster_gpu.cu`, arithmetic in
`include/calculation/gpu/small_cluster_math.h`):

- cyclic complex Jacobi for both conditional Hamiltonians (lower triangle, as Eigen reads it), then
  $M = V_\uparrow^\dagger V_\downarrow$;
- for every time point, the chain of `evolveSingleState`;
- the ensemble as the average of that chain over the basis states;
- the Hahn ensemble by the CPU's identity $\sum_{ab} |X_{ab}|^2\, \overline{D_{\uparrow,a}}\, D_{\uparrow,b} / D$ with $X = M D_\downarrow M^\dagger$.

Supported:

- pulse plans with at most 64 pulses (`maxIntervals` = `pulse.count` + 1 = 65) and at most four distinct interval
  lengths;
- phases `exact` or `recursive`, as the run's `evolution_update`.

These limits are applied after the general input checks. Before them, an ensemble `Manual` sequence is refused with
`backend=cuda` and switched to `cpu` with `auto`.

The Jacobi rotations run on the Hamiltonian scaled by a power of two (exact) to a largest component in [0.5, 1), with
`hypot` for the magnitudes, so that no finite Hamiltonian overflows or underflows in them. A cluster goes back to the
CPU path when:

- its Jacobi does not converge (off-diagonal above 1e-12 of the diagonal);
- its eigenvalues, eigenvectors or coherences are not finite;
- its Hamiltonian is not finite and Hermitian.

One CUDA context per node:

- every rank builds its clusters' Hamiltonians with the CPU path's own function and packs them into its part of a node
  shared-memory window (MPI-3, page-locked once);
- node rank 0 runs a feeder thread that takes the submitted windows of all ranks (2048 clusters each, two per rank)
  into combined batches in two lanes, and writes the coherences back;
- the ranks assemble the product as before.

With one context per rank, 32 ranks starting CUDA at once took 3.5 s, and tearing it down about 2 s. The context and
the shared memory are therefore started in `computeBackendInit`, in a thread, while the bath is read and the clusters
are made. Until they are up, the ranks compute their clusters on the CPU.

When no rank has a cluster of a dimension the packed or the warp kernel takes, that start is released before the
calculation and no feeder runs. (Each rank looks at its clusters until it finds one; one Allreduce.) A node with
several GPUs uses its first one.

A window holds fewer clusters when their coherences would pass 6 MiB per slot (19 at 20,000 time points). When one
cluster's coherences alone pass it, or no rank has a cluster, the batch is off and nothing is allocated past the
start. The batch runs threads beside MPI, so `main.out` asks MPI for `MPI_THREAD_FUNNELED` (only the main thread calls
MPI); the batch is refused if the MPI library provides less.

A GPU that cannot be opened:

- with `backend` = `auto`, the batch is turned off and every cluster is computed on the CPU (the `auto` rule above);
  the log says `[gpu-packed] off: the GPU could not be opened`;
- with `cuda`, or an explicit `packed_max_dim` with `cpu`, the run stops. This holds also when the ranks had finished
  their clusters before the opening failed (the feeder's state is read again after it is joined), so the outcome does
  not depend on the timing.

Errors after the GPU is up stop the run under every `backend`.

Every loop over the matrix indices has a compile-time trip count and is unrolled. Without that, nvcc 12.8's NVVM
optimizer (cicc -O3) miscompiled the dynamically indexed local arrays: the device read a diagonal from before the
first rotation, and dimension-4 eigenvalues came out wrong by up to 3.4e6.

The standing check `test_small_cluster_cuda` (`make gpu-test`) runs the kernel on CPU-dumped clusters: 11 files
covering single state and ensemble, Hahn, CPMG, Ramsey, `exact` and `recursive`, with 128 production P1 clusters. The
tolerance is 1e-8 + 2e-6 |CPU| per point.

<a id="the-warp-kernel-for-clusters-up-to-dimension-32-gpu_optionswarp_max_dim"></a>
## The warp kernel for clusters up to dimension 32 (gpu_options.warp_max_dim)

The warp kernel is on by default in the GPU build with backend `auto`/`cuda` (default 27, every compiled dimension).
It is 0 in the CPU build and with backend `cpu`.

gCCE clusters with an NV-like qubit (dimension 3) have dimension 3 × bath dimension. That gives 6 and 12 for a P1 bath
(dimension 12 is 83–87 % of the P1 gCCE2 time), and 9, 12, 21, 27 in h-BN. These are too large for the packed kernel's
registers and far below `min_dim`. `"execution": {"gpu_options": {"warp_max_dim": 27}}` (the default) sends the
clusters of a compiled total dimension up to it to `src/calculation/gpu/warp_cluster_gpu.cu` (plan and sizes in
`include/calculation/gpu/warp_cluster.h`).

One cluster per warp:

- the matrices are in shared memory;
- the eigensystem comes from a parallel cyclic Jacobi with round-robin ordering: each round rotates disjoint index pairs
  at once, and the 32 lanes share the column and row updates. Rotation, threshold, 20-sweep limit, convergence test,
  power-of-two scaling and finiteness checks are those of the packed kernel. A cluster whose eigensystem or coherence
  is not finite goes back to the CPU;
- lane l computes the time points l, l + 32, ..., each the whole interval chain of `propagateStateGcce` (gCCE) or
  `evolveSingleState` (cCCE), with the eigenvectors read from shared memory.

What it computes:

- gCCE ensemble (`evolution=matrix`, $\rho_0 = \rho_q \otimes 1/d_\text{bath}$): the same chain averaged over the
  $d_\text{bath}$ starts $\psi_q \otimes e_m$;
- a single state with `evolution=matrix`: the pure start $\psi_q \otimes \psi_\text{bath}$;
- quantity coherence or dm (the qdim × qdim reduced density matrix);
- cCCE ensemble Hahn by the CPU's identity, as in the packed kernel;
- phases exact (with the tail of `refinedPhases` for gCCE).

Compiled dimensions: gCCE qubit dimension 3 × bath 2, 3, 4, 6, 7, 8, 9; cCCE bath 8 and 16.

The warp kernel uses the packed kernel's node feeder. Its clusters go into the same windows of node shared memory, in
an input area of their own (6 MiB per window; a window closes early when it is full), with their outputs after the
packed ones. The feeder launches one kernel per dimension. The gCCE plan (pulses, ψa, ψb, ψq) is taken in the first
state (`beginState`) and must not change afterwards; it cannot, since those are set once per run.

What it does not take runs as with `warp_max_dim` 0, with one `[gpu-warp] off: ...` line (the default switches itself
off):

- `gcce_refinement`;
- `pulse.per_qubit`;
- a pulse detuning;
- savemode `cluster_density_matrix`;
- a qubit of another dimension, or more than one qubit;
- more than 4 distinct pulse intervals.

With the key on, the large-cluster batching of `evolution_batch_size` is off, as with the packed kernel: the windows
and that batching both take clusters ahead.

A thread-per-cluster variant was miscompiled at `-O3` on sm_120 (Jacobi unconverged for 97 % of the clusters; correct
on sm_86 and with `-Xcicc -O0`), the same NVVM-optimizer problem with dynamically indexed local arrays as above. The
warp kernel therefore keeps the matrices in shared memory and selects its phase slot without a run-time index into the
per-lane arrays.

The standing check `test_warp_cluster_cuda` (`make gpu-test`) runs the kernel through the batch API on CPU-dumped
clusters in `test/gpu_verification/check_gpu/data/warp_cluster/`. Its 15 files cover every compiled dimension, single
state and ensemble, coherence and dm, Hahn, CPMG, Ramsey, and 64 production P1 CR60 clusters. The tolerance is
1e-8 + 2e-6 |CPU| per value, and both lanes must give the same bits. A deliberately broken rotation (the conjugate in
one column update) fails it with 32 FAIL lines. The h-BN clusters cannot see that change, because their off-diagonals
are real.

What to expect (32 ranks + MPS, time relative to the GPU build without the warp kernel):

| Case | Time ratio |
| --- | --- |
| P1 CR0 gCCE2 | 0.69 |
| P1 CR60 σ4 | 0.60 |
| h-BN gCCE2 | 0.86 |

The curves agree within the tolerance, except where the product diverges (\|L\| > 1 and between). There CPU `exact`
and `recursive` agree within it too, so the difference is the eigensolver's rounding, amplified.

## MPI and multiple GPUs

### Reusing a deterministic pCCE partition

The partition search runs on the CPU even with `backend=cuda`. For repeated calculations with the same selected bath
and partition settings, turn on an explicit cache:

```json
{
  "method": "cCCE",
  "cluster": {
    "mode": "pcce",
    "pcce": {
      "partition": "deterministic",
      "cache": {"mode": "readwrite", "file": "/path/to/existing-directory/nv-sk3.partition"}
    }
  },
  "configuration": {"seed": 73}
}
```

This is a fragment; the usual bath, `pcce.partition_size`, order and dynamics options are still required.

The cache modes:

- `readwrite`: the first run does the normal search and stores the group centers and the spin-to-group assignment.
  Later matching runs skip the search. A missing, stale or corrupt file triggers a new search and an atomic
  replacement.
- `read`: requires an existing, valid, matching file and never changes it. A missing, stale or corrupt file fails on
  all ranks together.
- `off`: normal execution; the cache is not accessed, also in legacy mode.

Other I/O failures are reported as errors. A cache with legacy mode or with another method is rejected once the input
is read, because skipping the legacy search would change the random stream used for the sampled bath states.

The identity of a cache file is checked against:

- the selected, transformed, ordered sites (names, spins, binary64 coordinates, including an excluded remainder);
- the effective seed, partition_size, trials, max_iterations, kmeans_pp;
- the partition algorithm version, the host compiler and the selected arithmetic compatibility markers.

It is compared in full; a filename or hash is not accepted as proof of a match. Centers are stored without rounding.
A file checksum detects accidental corruption, and group indices, group sizes and the finiteness of the centers are
checked before use. The checksum is not authentication: these are user-controlled simulation data files.

These are not partition-search inputs and may change while the same partition is reused: CPU/CUDA/auto, MPI count,
order, rdip, pulses, state mode, time points and output options. Connectivity, cluster weights, Hamiltonians, sampled
states and time evolution are computed afresh for the current input.

The cache is invalidated by:

- radius or rotation changes that alter the selected coordinates;
- different compiler or arithmetic markers (a conservative choice).

A cache hit keeps the deterministic sampled-state RNG behavior. It does not restore a previous coherence result.
Changing the partition algorithm requires bumping its version in `src/clustering/pcce_partition_cache.cpp`.

Rank zero reads or writes the validated partition and broadcasts it, after all ranks have verified identical input
identities. A write goes to a unique temporary file beside the destination, is flushed, and is renamed atomically, so
concurrent writers cannot expose half-written data. Use distinct files for different baths, partition sizes or
settings to avoid repeated replacement.

Log lines:

- hit: `[partition-cache] status=hit` and `source=cache total_trials=0`;
- miss or store: lookup, search and write times.

First-run timing includes the full search and the cache write. Reuse timing includes validating and reading the file
plus the rest of the simulation. Compare CPU/CUDA/auto under the same cache mode, and report cold and reuse timings
separately.

### Running with MPI

```bash
mpirun -n 4 ./bin/main_gpu.out -f input.json
```

Allocate the GPUs through your site's scheduler before launching the job. Start with one MPI rank per GPU. The MPI
cluster partitioning is unchanged.

Each rank selects `node_local_rank % visible_GPU_count`. When the scheduler exposes one GPU per rank, each rank selects
its visible device 0. The `[backend]` log records the mapping. Several ranks can share a visible GPU, but that can
reduce performance; check the scheduler's visibility and binding.

When many ranks share one GPU, run the job under CUDA MPS: `nvidia-cuda-mps-control -d` in the job,
`echo quit | nvidia-cuda-mps-control` at its end (`test/gpu_verification/server/mps.sh` does both). On idm4-gpu,
32 ranks on one GPU took 0.49-0.62 of the time without MPS (gCCE4 and cCCE4, `auto` and `cuda`), with byte-identical
curves.

The GPU returns scalar cluster coherences (cCCE) or small reduced density matrices (gCCE) to the host before the
cluster products and MPI reductions. CUDA-aware MPI is not required.

Up to eight ranks sharing one RTX 3050 were tested. Physical multi-GPU and multi-node runs still need validation on the
target server.

## Implementation

### Conventional CCE

1. The CPU constructs `Halpha` and `Hbeta`.
2. cuSOLVER diagonalizes both once per cluster.
3. The GPU forms the basis overlap `M = V_alpha^H V_beta`. Since the eigenvectors are unitary, this equals the CPU's
   `V_alpha.inverse()*V_beta` without a general inverse.
4. Evolution alternates between the two eigenbases, with cuBLAS handling many time points as columns of one matrix.
5. A CUDA reduction returns the overlap of the conditional states (single) or its normalized Frobenius inner product
   (ensemble).

The cCCE semantics are those of the CPU:

- ensemble uses Ramsey/CPMG intervals from `pulse.count`, including the CPU's float half-interval rounding and the
  even-pulse complex-conjugation convention;
- single follows the configured interval sequence and its delay-cache indices;
- cCCE toggles the conditional Hamiltonians, as on the CPU. The GPU path adds no finite-angle or finite-duration pulse
  physics to cCCE; such qubit-pulse evolution belongs to gCCE.

The GPU single path evolves two state vectors instead of constructing the full propagators of the CPU cCCE path.
Speedups therefore include this algorithmic change as well as the GPU hardware. The CPU path is kept for direct
regression comparison.

### Generalized CCE

The CPU constructs each Hamiltonian and pulse operator. cuSOLVER diagonalizes the complex-double Hermitian Hamiltonian
once per cluster. Eigenvectors, eigenvalues, states and pulse data stay on the GPU through the time grid.

cuBLAS matrix products evolve a batch of time points. CUDA kernels apply the phases, the qubit pulses and the bath
partial traces. A pulse is applied with the qubit operator directly, without forming its tensor product with the bath
identity.

The host supplies the same repeated-addition time grid and delay-cache indices as the CPU path. Pure-state and
mixed-state evolution mean the same as on the CPU. Hamiltonians and evolution buffers are complex double; no
reduced-precision or fast-math mode is enabled.

With `gcce_refinement` on (the block written; off by default, see the backend tags), the CPU and CUDA eigen paths both
refine the complete eigensystem:

- compensated residuals are projected into the computed basis;
- small Jacobi rotations correct eigenvectors and energies (on the GPU, with parallel round-robin Jacobi sweeps);
- the small energy correction and the phase-product rounding remainder are kept separately.

This targets propagation errors amplified by negative cluster powers in the NV–13C order-4 regression. Refinement adds
setup work; include its cost when comparing with older binaries.

Scalar coherence products in cCCE, pCCE and gCCE keep a normalized complex mantissa and a separate binary exponent,
also during the MPI reduction and, for gCCE, the final normalization. The host mantissa is `long double`. A finite
final result is then not lost when an intermediate local product underflows or overflows. This applies to
`quantity=coherence`; gCCE `quantity=dm` and the DSJ/ITB variants keep their previous product arithmetic. It does not
repair a zero denominator or an unconverged CCE approximation.

Small-system tests compare against an independent general matrix exponential, including complex degenerate
eigenspaces, noncommuting pulses, full density matrices and trace conservation.

Consecutive clusters of the same dimension, state mode and pulse schedule can share batched cuBLAS evolution calls.
This applies to cCCE, pCCE and gCCE reduced outputs. Hamiltonians and eigensystems are still prepared per cluster, and
clustering and Hamiltonian construction stay on the CPU. Cluster order and product order are preserved.
`savemode=cluster_density_matrix` uses sequential gCCE cluster evolution.

Grouping uses a conservative buffer estimate and the two cluster-batch limits. The budget is not a hard cap on total
device memory: a single cluster, the cuSOLVER workspace, retained free buffers and CUDA runtime allocations also use
memory. Reduce the time or cluster batch sizes, or the cache size, if needed. Every rank has its own cache and active
buffers. The cache reuses allocations on the same CUDA stream, frees excess cached buffers, and is released at
finalization.

With profiling on:

- `[gpu-memory]` records allocations, reuses and cached bytes;
- `[gpu-batch]` records actual multi-cluster batches. A shared batch's evolution time is divided among its clusters for
  the dimension totals.

`profile` reports Hamiltonian construction (`build_s`), cluster preparation (`setup_s`) and time evolution
(`evolution_s`). On the GPU, preparation and evolution include the synchronization and transfers of those stages. The
stage totals leave out application start-up, clustering, MPI communication and output; judge end-to-end performance by
the benchmark's wall times.

The gCCE full-density output allocates its label array for the qubit-plus-bath labels, including the qubit entries.

Library references: [cuSOLVER Hermitian eigensolver](https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-syevd),
[cuBLAS](https://developer.nvidia.com/cublas),
[NVIDIA redistributable installation](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#tarball-and-zip-archive-deliverables).

## Validation and reproduction

On a GPU server, the CCE2 Slurm jobs (`test/gpu_verification/server/README.md`) build the code on the server and
compare all four combinations on the original inputs.

Which backend was faster on the tested server:

- small NV–¹³C and PL6 matrices: CPU;
- large MoS₂ pCCE matrices: CUDA/auto, with about four MPI ranks per GPU;
- P1 sK2: CUDA at two MPI ranks, but CPU at 16 ranks was faster than the best CUDA configuration.

The dimension-only `auto` threshold is a heuristic. Use bath-specific timings to choose the backend and the MPI count.

The Python tools use only the standard library, read the existing regression inputs, and require a new output
directory, so previous results are kept.

```bash
# Independent GPU numerical oracle (use your build variables as above).
make gpu-test        # after: bash do_compile.sh gpu

# gCCE4 accuracy, MPI consistency, cluster batching and deferred initialization.
python3 test/gpu_verification/check_gpu/check_acceleration.py --cpu bin/main.out --cuda bin/main_gpu.out \
    --ranks 2 --output /tmp/ccex-acceleration-new

# Deterministic pCCE across MPI 1, 2, 3 and 8; optional old binary checks legacy mode.
python3 test/cluster_verification/check_pcce/check_reproducibility.py --cpu bin/main.out --cuda bin/main_gpu.out \
    --output /tmp/ccex-reproducibility-new

# Small pCCE ensemble/single, pulse, partition and backend checks (Open MPI or Intel MPI).
python3 test/gpu_verification/check_gpu/check_pcce.py --cpu bin/main.out --cuda bin/main_gpu.out \
    --ranks 2 --output /tmp/ccex-pcce-check-new

# Eight existing cases: CCE/gCCE x ensemble/sampled x NV–13C/NV–P1.
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out \
    --reference bin/main.out --method all --repeats 3 \
    --output /tmp/ccex-gpu-regression

# Larger cCCE system; --nearest-sites also sets CCE order to five.
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out \
    --reference bin/main.out --method cce --case P1 --nearest-sites 5 \
    --steps 20 --batch 20 --repeats 3 --output /tmp/ccex-cce-gpu-large

# Larger, controlled system: nearest four P1 sites and their attached 15N spins.
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out \
    --reference bin/main.out --method gcce --case P1 --nearest-sites 4 --steps 20 \
    --output /tmp/ccex-gpu-large

# Backend errors and CPU fallback with GPUs hidden.
python3 test/gpu_verification/check_gpu/check_options.py --cpu bin/main.out \
    --cuda bin/main_gpu.out --output /tmp/ccex-gpu-options

# cCCE Ramsey, odd/even CPMG, manual single intervals and time-batch tails.
python3 test/gpu_verification/check_gpu/check_cce_pulses.py --cpu bin/main.out \
    --cuda bin/main_gpu.out --output /tmp/ccex-cce-gpu-pulses
```

Useful benchmark options:

- `--mpirun /path/to/mpirun --ranks 2 --backend auto --min-dim 16 --batch 7`: exercises MPI, mixed CPU/GPU selection
  and batch tails;
- `--method gcce --savemode cluster_density_matrix`: checks the exported full density matrices as well as the
  coherence;
- `--repeats 3`: alternates the execution order, reports median wall times and the largest observed errors, and keeps
  every run's timings and profile.

Exit code 77 from the standalone CUDA numerical test means the GPU could not be initialized. It is not a numerical
pass.
