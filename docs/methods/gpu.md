# GPU execution

CCEX has an optional NVIDIA CUDA backend for **conventional CCE, pCCE and gCCE** (`method=cCCE`; `method=cCCE` with
`cluster.mode=pcce`; `method=gCCE` with `propagator=eigen`), both ensemble and single-sample calculations.
The keys on this page (`backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `gcce_refinement`,
`profile`) are keys of the `execution` block: `backend=auto` below means
`"execution": {"backend": "auto"}`; `pcce.*` is `cluster.pcce.*`, `savemode` is `output.savemode`, `probe_memory_gb` is
`cluster.adaptive.probe_memory_gb`.
The original CPU backend remains the default. gCCE with `propagator=expm` and
other method names still run on the CPU.

| Method | Ensemble (no `configuration.s_range`) | Single (`configuration.s_range`) |
| --- | --- | --- |
| Conventional CCE, also called cCCE | GPU normalized bath trace | GPU conditional state vectors |
| Partitioned CCE, `cluster.mode=pcce` | GPU normalized bath trace | GPU conditional state vectors |
| gCCE, `propagator=eigen` | GPU density-matrix evolution | GPU vector or density-matrix evolution |
| Bundle clustering, `cluster.mode=bundle` (cCCE, or gCCE with `propagator=eigen`) | as cCCE or gCCE above, per cluster of bundles | as cCCE or gCCE above |

For cCCE and pCCE, `propagator` and `evolution` do not select an algorithm: the existing
cCCE method diagonalizes two conditional Hamiltonians, and its GPU single path
uses vectors. The `evolution` option belongs to gCCE; use the default or `matrix`
for a conventional CCE input, as the existing parser rejects an explicit
`evolution=vector` outside gCCE.

pCCE uses the same conventional coherence kernels after partitioning bath spins.
Its partition search, cluster construction and Hamiltonian construction stay on
the CPU. CUDA accelerates the conditional Hamiltonian diagonalization and time
evolution. pCCE profiles are labelled `method=pcce`.

Coherence products in cCCE, pCCE and gCCE retain a separate exponent through
local accumulation and MPI reduction. This prevents intermediate overflow when
the final product is representable; it does not establish CCE convergence.
The optional adaptive refinement ([cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive)) uses the
same CPU/CUDA engines to test and improve variable-size groups.

## Build

`do_compile.sh gpu` builds two executables with one `make`, from the same sources
into separate object directories: `bin/main.out` (CPU only, what every other
preset builds) and `bin/main_gpu.out` (compiled with `-DCCEX_USE_CUDA` and linked
with the CUDA kernels). MPI and the bundled Eigen headers are used; this build
does not require Intel MKL. The CUDA half needs a toolkit with nvcc, cuSOLVER,
cuBLAS and cudart. CUDA 12.8.1 with GCC 13.3 (RTX 3050) and CUDA 13.2.2 with GCC
11.5 (idm4-gpu) were tested. On idm4-gpu's compute capability 12.0, CUDA 12.8
miscompiles the pulse kernel `applyPulse` (A97, see `CHANGELOG.md`), so the
`gpu` preset uses 13.2.

```bash
bash do_compile.sh gpu              # idm4-gpu: loads gcc, intel/oneapi-2023.2 (mpiicpc with icpx, MKL) and CUDA 13.2 itself
CCEX_JOBS=4 bash do_compile.sh gpu  # fewer parallel compiles on a shared node
```

On a machine without those modules, put your MPI wrapper (`mpicxx`) on `PATH` and
say where the toolkit is and what the card is:

```bash
CUDA_HOME=/usr/local/cuda CUDA_ARCH=86 bash do_compile.sh gpu
```

`CUDA_ARCH` is the target CUDA compute capability without a decimal point (`120`
for the RTX PRO 6000 Blackwell on idm4-gpu, which is the default; `86` for the
RTX 3050 used in early validation). Both device code and PTX for that
architecture are embedded. Choose an architecture supported by both your GPU and
toolkit. Switching preset cleans the objects before building; after changing only
`CUDA_HOME` or `CUDA_ARCH`, run `make clean` first.

Keep the chosen toolkit's shared libraries available when running the binary. On
idm4-gpu that means loading the same three modules in the job script
(`gcc/gcc-11.5.0 intel/oneapi-2023.2 cuda/cuda-13.2`): without them the binary stops
at start-up on the system's older libstdc++.

gCCE eigen propagation can refine eigenpairs on both CPU and CUDA to resolve
small bath energies beside large central-spin energies: writing the `gcce_refinement` block turns it on.
It is off by default (on the NV–¹³C gCCE4 case it was made for
it changed the error-corrected final coherence by at most 6.7e-7 — the predeclared bound 1e-6 on every raw `wiDiv`
point was not met, because `wiDiv` diverges at some points whatever the setting (all-point differences up to 12.7, and
2.75e11 on P1 gCCE3; (internal notes)) — while it took most of a cluster's setup
time; `../user/InputFile_Backend.md`). When on, it is strict: nonconvergence within
`gcce_refinement.max_sweeps` (default 32) stops the run and reports the dimension
and sweep limit, and it never silently accepts the unrefined eigensystem.
Compensated phase products and scaled coherence accumulation are active either way. These options
apply only to gCCE with `propagator=eigen`.

An omitted `configuration.seed` is generated on MPI rank 0, broadcast to every rank, and
printed for replay. pCCE `pcce.print_trials=N` prints at most N
trial details from rank 0 per partition search (default 0, none).
The global partition summary still aggregates every rank. `cluster.pcce.partition` = `legacy`
preserves the original asymmetric tie rule: larger silhouette within a rank,
smaller silhouette between ranks at equal inertia. Deterministic mode uses
larger silhouette throughout.

The C RNG is initialized once in `main`, from the resolved seed. It drives the
sampled bath/subspin states, random defect axes and legacy pCCE searches;
deterministic pCCE partition trials use their own seeded RNG. The regression
checks different seeds and an independently generated C RNG state sequence as
well as CPU/CUDA and MPI replay. Fixed-bath ensemble calculations without random
defect axes or legacy partition search do not use the C RNG. Partition cache
contents do not depend on it either, because deterministic partitions do not use
the C RNG; raw coherence caches must match the executing binary and sampled states.

Cluster batch admission accounts for state vectors versus density matrices.
The budget remains a conservative estimate, not a guarantee for solver workspace
or total device memory. Device pointer tables and batch output buffers persist
with the first work object; unchanged tables are not re-uploaded. With profiling,
`[gpu-memory] pointer_uploads` counts table uploads. Membership, buffer-address,
or density time-count changes refresh the corresponding tables.

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

This fragment is not a complete simulation input. Complete examples for all four
combinations are available; run them from the repository root:

```bash
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_cce5_ensemble.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_cce5_single.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_gcce4.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/p1_gcce4_single.json
./bin/main_gpu.out -f test/gpu_verification/check_gpu/examples/pcce_ensemble.json
# one state: a copy of pcce_ensemble.json with "configuration.s_range": [1, 1] and "output.file": "pcce_gpu_single"
./bin/main_gpu.out -f pcce_single.json -o bin
```

The cCCE examples use five P1 sites (bath dimension 1024); the gCCE examples use
four (qubit-plus-bath dimension 768). Each is a controlled size experiment, with
distinct output names under `bin/`. Use the regression suite below to compare
the same inputs between CPU and GPU.

The complete list of these keys, with their types and defaults, is
[the backend tags](../user/InputFile_Backend.md) (and, for pCCE,
[the cluster tags](../user/InputFile_Cluster.md)); the table below is a summary.

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

pCCE also accelerates fresh searches by evaluating distances only to centers
whose groups still have room, retaining spin order, ascending center tie-breaking,
trial RNG and convergence criteria. This CPU optimization applies to all three
backends, with cache off or on. It removes the N-by-K distance table from the
search path; it does not move partitioning to CUDA.
Fresh searches also detect exact repeated center/assignment states and skip
whole periods, then execute the remaining iterations to preserve the original
`pcce.max_iterations` endpoint. They do not declare a cycle converged or change trial
selection. This applies to deterministic and legacy pCCE on every backend.
The root `[partition-search]` log counts converged, cycling and limit-only
trials, followed by logical, executed and skipped iterations summed over ranks.
`pcce.print_trials=true` adds each trial's period and counts.
The partition cache algorithm identity is now `pcce-deterministic-20260912-v3`:
`readwrite` recomputes old v1/v2 identities once; `read` rejects them until
regenerated. The on-disk cache format remains version 1.

The rest of this paragraph is about the per-rank CUDA initialization of the
large-cluster path; the packed kernel's context (below) is opened separately by
node rank 0 at the start whenever `packed_max_dim` > 0, which in the GPU build is
the default for CCE with `cuda`/`auto`. `auto` defers the per-rank initialization
until a rank encounters its first eligible cluster, so a run containing only
smaller clusters does not initialize it on that rank. Ranks
with no eligible clusters do not wait for those that initialize: only node-rank
metadata discovery is collective. `backend=cuda` still requires initialization
on every rank at startup.

`auto` also works in a CPU-only executable and falls back to CPU when CUDA
initialization is unavailable. Once GPU computation starts, solver or execution
failures abort the MPI job with an error; they do not silently rerun a partially
completed cluster. Running out of GPU memory is the exception under `auto` (without it, 20 ranks
holding dimension-2197 ensemble clusters on one 96 GB card at once aborted the run
and the MPS server with it): what the cluster had on
the card is freed, the cluster is computed on the CPU from the start, and that rank
sends no later cluster of that dimension or larger to the GPU. The log has one
`[gpu-memory-fallback] rank=R dim=D` line per rank and dimension, and a count per
rank at the end; the curves are the CPU path's for those clusters. With `cuda` an
out-of-memory still stops the run. `CCEX_TEST_GPU_OOM_BYTES=<n>` (a test hook) fails
every device allocation of n bytes or more, as a full card would; `check_options.py`
uses it for cCCE and gCCE, one by one and in cluster batches. An explicit `cuda` request with gCCE `expm` or an
unsupported method is rejected once the input is read.

`cluster.mode=bundle` adds variable-size automatic grouping to `method=cCCE`
and `method=gCCE` while retaining the selected CPU/CUDA/auto engine. The dimension
cap is `cluster.bundle.max_dim` (required); bundle orders of 3 and more are supported.
See [cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive) for the refinement on top of it.

The dimension threshold is a configurable starting point, **not a guarantee of
speedup**. Inspect the profile for your spin species, defect spins, state mode,
time grid and GPU. Matrix dimension includes all attached defect spins and, for
gCCE, the qubit. It is not just the CCE order. Reduce the time batch size if GPU
memory is tight. gCCE matrix evolution has two main complex-double buffers of
`16 * D² * batch_steps` bytes each; cCCE ensemble uses three. The corresponding
pure-state buffers hold `D * batch_steps` elements each. Eigensystems, solver
workspace and other buffers add to this. gCCE `savemode` `cluster_density_matrix` requires additional
density-matrix storage and host output storage.

For pCCE, `order` counts groups and `pcce.partition_size` is the number of bath sites per group.
The largest cluster can contain `order * partition_size` sites plus attached defect spins.
This grouping is different from enumerating every conventional CCE cluster of
that spin order. The current partitioner drops the final `nspin % partition_size` sites;
this existing behavior is identical for CPU and CUDA. Supply a positive integer
`pcce.partition_size` no greater than the loaded bath size, and positive `pcce.trials`/`pcce.max_iterations`.
With `cluster.pcce.partition` = `deterministic`, `pcce.trials` is the exact global trial count.
Trial `t` has its own RNG stream derived from `(seed,t)` and is assigned to rank
`t % ranks`. The search does not consume the bath-state C RNG. Equal candidates
are resolved by inertia, silhouette and then global trial index. With an explicit
seed and the same executable/input, changing MPI size preserves the candidate
set and selected partition. Floating-point products can still differ in their
last bits because MPI reduction order changes; this is not a promise of bitwise
identity across compilers or architectures.

This default changes the random search and sampled bath states relative to old
pCCE runs. Set `cluster.pcce.partition` to `legacy` when reproducing an existing result. In
that mode retain the original seed **and MPI size**; each rank uses the historical
`max(1, pcce.trials / ranks)` trials and shared C RNG behavior. Legacy mode preserves
its MPI dependence and potentially repeated search streams.

<a id="the-packed-kernel-for-small-ccce-clusters-gpu_optionspacked_max_dim"></a>
## The packed kernel for small cCCE clusters (gpu_options.packed_max_dim)

The paths above give each cluster its own eigensolver calls, copies and synchronizations. For the bath-dimension 2 and 4
clusters of a P1 bath (millions per bath state) that fixed cost is about 400 times the CPU's time per cluster
((internal notes)), and batching the time evolution
(`evolution_batch_size`) does not remove it. `"execution": {"gpu_options": {"packed_max_dim": 4}}` (or 2) takes those clusters off
that path:

- **one thread per cluster** (`src/calculation/gpu/small_cluster_gpu.cu`, arithmetic in
  `include/calculation/gpu/small_cluster_math.h`): cyclic complex Jacobi for both conditional Hamiltonians (lower
  triangle, as Eigen reads it), $M = V_\uparrow^\dagger V_\downarrow$, then for every time point the chain of `evolveSingleState` - the
  ensemble as the average of that chain over the basis states, the Hahn ensemble by the CPU's identity
  $\sum_{ab} |X_{ab}|^2\, \overline{D_{\uparrow,a}}\, D_{\uparrow,b} / D$ with $X = M D_\downarrow M^\dagger$. Pulse plans with at most 64 pulses (`maxIntervals` = `pulse.count` + 1 = 65)
  and at most four distinct interval lengths, chosen after the general input checks (an ensemble `Manual` sequence is
  refused with `backend=cuda` and switched to `cpu` with `auto` before this); phases
  `exact` or `recursive`, as the run's `evolution_update`. The Jacobi rotations run on the Hamiltonian scaled by a
  power of two (exact) to a largest component in [0.5, 1), with `hypot` for the magnitudes, so that no finite
  Hamiltonian overflows or underflows in them. A cluster whose Jacobi does not converge (off-diagonal above
  1e-12 of the diagonal), whose eigenvalues, eigenvectors or coherences are not finite, or whose Hamiltonian is not
  finite and Hermitian, goes back to the CPU path;
- **one CUDA context per node**: every rank builds its clusters' Hamiltonians with the CPU path's own function and packs
  them into its part of a node shared-memory window (MPI-3, page-locked once); node rank 0 runs a feeder thread that
  takes the submitted windows (2048 clusters each, two per rank) of all ranks into combined batches in two lanes and
  writes the coherences back. The ranks assemble the product as before. With one context per rank, 32 ranks starting
  CUDA at once took 3.5 s and tearing it down about 2 s; the context and the shared memory are therefore
  started in `computeBackendInit`, in a thread, while the bath is read and the clusters are made, and until they are
  up the ranks compute their clusters on the CPU. When no rank has a cluster of a dimension the packed or the warp
  kernel takes (each rank looks at its clusters until it finds one; one Allreduce), that start is released before the
  calculation and no feeder runs. A node with several GPUs uses its first one. A window holds fewer
  clusters when their coherences would pass 6 MiB per slot (19 at 20,000 time points); when one cluster's alone pass it,
  or no rank has a cluster, the batch is off and nothing is allocated past the start. The batch runs threads beside
  MPI, so `main.out` asks MPI for `MPI_THREAD_FUNNELED` (only the main thread calls MPI) and the batch is refused if
  the MPI provides less;
- **a GPU that cannot be opened**: with `backend` = `auto` the batch is turned off and every cluster is computed on the
  CPU (the `auto` contract above; the log says `[gpu-packed] off: the GPU could not be opened`); with `cuda`, or an
  explicit `packed_max_dim` with `cpu`, the run stops - also when the ranks had finished their clusters before the
  opening failed (the feeder's state is read again after it is joined), so the outcome does not depend on the timing.
  Errors after the GPU is up stop the run under every `backend`.

Every loop over the matrix indices has a compile-time trip count and is unrolled. The step-0 prototype without that
unrolling got dimension-4 eigenvalues wrong on the device (up to 3.4e6): after the first rotation the device read the
diagonal from before it. `-G` and `-Xcicc -O0` gave the host's values and `ptxas -O0` did not, so it is nvcc 12.8's
NVVM optimizer (cicc -O3) on the dynamically indexed local arrays
((internal notes)). The standing check `test_small_cluster_cuda`
(`make gpu-test`) runs the kernel on CPU-dumped clusters (11 files: single state and ensemble, Hahn, CPMG, Ramsey,
`exact` and `recursive`, 128 production P1 clusters) against 1e-8 + 2e-6 |CPU| per point.

<a id="the-warp-kernel-for-clusters-up-to-dimension-32-gpu_optionswarp_max_dim"></a>
## The warp kernel for clusters up to dimension 32 (gpu_options.warp_max_dim)

On by default in the GPU build with backend `auto`/`cuda` (default 27 = every compiled dimension; 0 in the CPU build
and with backend `cpu`). gCCE clusters with an NV-like qubit (dimension 3) are 3 × bath dimension: 6, 12 for a
P1 bath (dimension 12 is 83–87 % of P1 gCCE2 time), 9, 12, 21, 27 in h-BN — too large for the packed kernel's registers
and far below `min_dim`. `"execution": {"gpu_options": {"warp_max_dim": 27}}` (the default) sends the clusters of a compiled total
dimension up to it to
`src/calculation/gpu/warp_cluster_gpu.cu` (plan and sizes in `include/calculation/gpu/warp_cluster.h`):

- **one cluster per warp**: the matrices in shared memory; the eigensystem by a parallel cyclic Jacobi (round-robin
  ordering: each round rotates disjoint index pairs at once, the 32 lanes sharing the column and row updates; the
  rotation, threshold, 20-sweep limit, convergence test, power-of-two scaling and finiteness checks of the packed kernel;
  a cluster whose eigensystem or coherence is not finite goes back to the CPU); then lane l computes the time points
  l, l + 32, ..., each the whole interval chain of `propagateStateGcce` (gCCE) or of `evolveSingleState` (cCCE) with the
  eigenvectors read from shared memory. gCCE ensemble (`evolution=matrix`, $\rho_0 = \rho_q \otimes 1/d_\text{bath}$) is the same chain averaged
  over the $d_\text{bath}$ starts $\psi_q \otimes e_m$; a single state with `evolution=matrix` uses the pure start $\psi_q \otimes \psi_\text{bath}$. Quantity
  coherence or dm (the qdim × qdim reduced density matrix). cCCE ensemble Hahn by the CPU's identity, as the packed
  kernel. Phases exact (with the tail of `refinedPhases` for gCCE). Compiled: gCCE qubit dimension 3 × bath 2, 3, 4, 6, 7,
  8, 9; cCCE bath 8 and 16. The step-0 prototype ((internal notes)) of this scheme ran
  P1 dimension-12 clusters 1.6 times as fast as 32 CPU cores;
- **the packed kernel's node feeder**: the clusters go into the same windows of node shared memory, in an input area of
  their own (6 MiB per window; a window closes early when it is full), their outputs after the packed ones; the feeder
  launches one kernel per dimension. The gCCE plan (pulses, ψa, ψb, ψq) is taken in the first state (`beginState`) and
  must not change afterwards (it cannot: those are set once per run).

What it does not take runs as with `warp_max_dim` 0 (one `[gpu-warp] off: ...` line; the default switches itself off): `gcce_refinement`, `pulse.per_qubit`, a pulse
detuning, savemode `cluster_density_matrix`, a qubit of another dimension, more than one qubit, more than 4 distinct pulse intervals.
With the key on, the large-cluster batching of `evolution_batch_size` is off, as with the packed kernel (the windows and
that batching both take clusters ahead). The step-0 thread-per-cluster variant was wrong at `-O3` on sm_120 (Jacobi
unconverged for 97 % of the clusters; right on sm_86 and with `-Xcicc -O0`), the same NVVM-optimizer pattern as the packed
kernel's prototype (dynamically indexed local arrays). The warp kernel keeps the matrices in shared memory and selects its phase slot
without a run-time index into the per-lane arrays. The standing check `test_warp_cluster_cuda` (`make gpu-test`) runs it
through the batch API on CPU-dumped clusters (`test/gpu_verification/check_gpu/data/warp_cluster/`, 15 files: every
compiled dimension, single state and ensemble, coherence and dm, Hahn, CPMG, Ramsey, 64 production P1 CR60 clusters)
against 1e-8 + 2e-6 |CPU| per value and both lanes the same bits. Measured (32 ranks + MPS, against the GPU
build without it): P1 CR0 gCCE2 r 0.69, CR60 σ4 0.60, h-BN gCCE2 0.86; the curves within the tolerance except 14 of 2416
P1 CR60 values in the stretch where the product diverges (|L| > 1 and between), which the acceptance leaves out (there CPU exact and recursive agree within it, so the difference is the eigensolver's rounding amplified); a deliberately broken rotation (the conjugate in one
column update) fails it with 32 FAIL lines (the h-BN clusters, whose off-diagonals are real, cannot see that change).

## MPI and multiple GPUs

### Reusing a deterministic pCCE partition

Partition search is a CPU stage even with `backend=cuda`. For repeated calculations
with the same selected bath and partition settings, enable an explicit cache:

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

This is a configuration fragment; the usual bath, `pcce.partition_size`, order and dynamics options
are still required. In `readwrite`, the first run performs the unchanged search
and stores the group centers and spin-to-group assignment. Later matching runs
skip the search. `read` requires an existing, valid, matching file and never
changes it; missing, stale or corrupt files fail collectively. In `readwrite`,
those files trigger a new search and atomic replacement. Other I/O failures are
reported as errors. `off` preserves normal execution and does not access the
cache, including in legacy mode. Enabling a cache for legacy or another method
is rejected once the input is read, because skipping legacy search would
change the random stream used for sampled bath states.

Identity is checked against the actual selected, transformed, ordered sites
(names, spins, binary64 coordinates, including an excluded remainder), effective
seed, partition_size, trials, max_iterations, kmeans_pp, partition algorithm version, host compiler
and selected arithmetic compatibility markers. It is compared in full, rather
than accepting a filename or hash as proof of a match. Centers are stored without
rounding. A file checksum detects accidental corruption, and group indices,
group sizes and center finiteness are checked before use. The checksum is not
authentication: these are user-controlled simulation data files.

CPU/CUDA/auto, MPI count, order, rdip, pulses, state mode, time points and output
options are not partition-search inputs. They may change while reusing the same
partition. Subsequent connectivity, cluster weights, Hamiltonians, sampled
states and time evolution are calculated afresh for the current input. Radius
or rotation changes that alter the selected coordinates invalidate the cache.
Different compiler/arithmetic markers conservatively invalidate it too. Cache
hits preserve deterministic sampled-state RNG behavior. They do not restore a
previous coherence result. Changing the partition algorithm requires bumping its
version in `src/clustering/pcce_partition_cache.cpp`.

Rank zero reads/writes and broadcasts the validated partition. All ranks first
verify identical input identities. Writes use a unique temporary file beside the
destination, flush it and rename it atomically. Concurrent writers cannot expose
half-written data; using distinct files for different baths/partition sizes/settings avoids
repeated replacement. A hit reports `[partition-cache] status=hit` and
`source=cache total_trials=0`. A miss/store reports lookup, search and write times.

First-run timing includes the full search and cache write. Reuse timing includes
file validation/read plus the full remaining simulation. Compare CPU/CUDA/auto
under the same cache mode and report cold and reuse timings separately.

### Running with MPI

```bash
mpirun -n 4 ./bin/main_gpu.out -f input.json
```

Allocate the GPUs through your site's scheduler before launching the job. One MPI
rank per GPU is the initial recommended configuration. The existing MPI cluster
partitioning is retained. Each rank selects `node_local_rank % visible_GPU_count`;
when a scheduler exposes one GPU per rank, each rank selects its visible device 0.
The `[backend]` log records this mapping. Multiple ranks can share a visible GPU,
but that can reduce performance. Check the scheduler's visibility and binding.
When many ranks share one GPU, run the job under CUDA MPS (`nvidia-cuda-mps-control -d` in the job, `echo quit |
nvidia-cuda-mps-control` at its end; `test/gpu_verification/server/mps.sh` does both): on idm4-gpu, 32 ranks on one
GPU took 0.49-0.62 of the time without it (gCCE4 and cCCE4, `auto` and `cuda`), with byte-identical curves.

The GPU returns scalar cluster coherences (cCCE) or small reduced density matrices
(gCCE) to the host before the existing cluster products and MPI reductions.
CUDA-aware MPI is not required.
Up to eight ranks sharing one RTX 3050 were tested; physical multi-GPU and multi-node
execution still need validation on the target server.

## Implementation

### Conventional CCE

The CPU constructs `Halpha` and `Hbeta`. cuSOLVER diagonalizes both once per
cluster. The GPU forms the basis overlap `M = V_alpha^H V_beta`; since the
eigenvectors are unitary, this represents the CPU's `V_alpha.inverse()*V_beta`
without a general inverse. Evolution alternates between the two eigenbases,
with cuBLAS handling many time points as columns of one matrix. A CUDA reduction
returns the overlap of the conditional states (single) or its normalized
Frobenius inner product (ensemble).

This retains the existing cCCE semantics: ensemble uses Ramsey/CPMG intervals
from `pulse.count`, including the CPU's float half-interval rounding and even-pulse
complex-conjugation convention. Single follows the configured interval sequence
and its delay-cache indices. As on the CPU, cCCE toggles conditional Hamiltonians;
this change does not add general finite-angle or finite-duration pulse physics
to cCCE. Such qubit-pulse evolution belongs to gCCE.

The GPU single implementation evolves two state vectors instead of constructing
the full propagators used by the current CPU cCCE implementation. Reported
speedups therefore include this algorithmic improvement as well as GPU hardware.
The CPU calculation path has been retained for direct regression comparison.

### Generalized CCE

The CPU constructs each Hamiltonian and pulse operator. cuSOLVER diagonalizes the
complex-double Hermitian Hamiltonian once per cluster. Eigenvectors, eigenvalues,
states and pulse data remain on the GPU through the time grid. cuBLAS matrix
products evolve a batch of time points; CUDA kernels apply phases, qubit pulses
and bath partial traces. Pulse application uses the qubit operator directly,
without materializing its tensor product with the bath identity.

The host supplies the same repeated-addition time grid and delay-cache indices as
the CPU path. Pure-state and mixed-state evolution preserve their original
meanings. Hamiltonians and evolution buffers use complex double; no
reduced-precision or fast-math mode is enabled. With `gcce_refinement` on (the block written; off by
default, see the backend tags) both CPU and CUDA eigen paths
refine the complete eigensystem: compensated residuals are projected into the
computed basis and small Jacobi rotations correct eigenvectors and energies
(on the GPU with parallel round-robin Jacobi sweeps).
The small energy correction and phase-product rounding remainder are retained
separately. This addresses propagation errors amplified by negative cluster
powers in the NV–13C order-4 regression. Refinement adds setup work; its cost
must be included in comparisons with older binaries.

Scalar coherence products in cCCE, pCCE and gCCE keep a normalized complex
mantissa and separate binary exponent, including during MPI reduction and,
for gCCE, final normalization. The
host mantissa uses `long double`. This prevents a finite final result from being
lost when an intermediate local product underflows or overflows. It applies to
`quantity=coherence`; gCCE `quantity=dm` and DSJ/ITB variants retain their prior
product arithmetic. This does not repair a genuinely zero denominator or an
unconverged CCE approximation. Small-system tests compare against an independent
general matrix exponential, including complex degenerate eigenspaces,
noncommuting pulses, full density matrices and trace conservation.

Consecutive clusters of the same dimension, state mode and pulse schedule can
share batched cuBLAS evolution calls. This applies to cCCE, pCCE and gCCE reduced
outputs. Hamiltonians and eigensystems are still prepared per cluster; clustering
and Hamiltonian construction remain on the CPU. Cluster order and product order
are preserved. `savemode=cluster_density_matrix` uses sequential gCCE cluster evolution.

Grouping uses a conservative buffer estimate and the two cluster-batch limits.
The budget is **not a hard cap on total device memory**: a single cluster,
cuSOLVER workspace, retained free buffers and CUDA runtime allocations also use
memory. Reduce time/cluster batch sizes or cache size if needed. Every rank has
its own cache and active buffers. The cache reuses allocations on the same CUDA
stream, frees excess cached buffers, and releases the cache at finalization.

With profiling enabled, `[gpu-memory]` records allocations, reuses and cached
bytes; `[gpu-batch]` records actual multi-cluster batches. A shared batch's
evolution time is divided among its clusters for the dimension totals.

`profile` reports Hamiltonian construction (`build_s`), cluster preparation
(`setup_s`) and time evolution (`evolution_s`). GPU preparation and evolution
include synchronization and transfers in those stages. These stage totals do not
include the entire application startup, clustering, MPI communication or output;
use the benchmark's wall times to judge end-to-end performance.

During full-density output validation, the existing gCCE label-array allocation
was also corrected: the writer reads qubit-plus-bath labels, so the allocation
must include the qubit entries.

Library references: [cuSOLVER Hermitian eigensolver](https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-syevd),
[cuBLAS](https://developer.nvidia.com/cublas),
[NVIDIA redistributable installation](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#tarball-and-zip-archive-deliverables).

## Validation and reproduction

On a GPU server, CCE2 Slurm jobs (`test/gpu_verification/server/README.md`) build
the code on the server and compare all four combinations on the original inputs.
See the server CCE2–4 measurements ((internal notes)) and
the test summary ((internal notes), in Korean). CCE2 and CCE3 passed; CCE4
stopped on a numerical comparison failure after six passing cases; the full-bath CCE5 case was not run.
The accuracy fix, MPI reproducibility and cluster batching that followed were validated locally and on the server.
On the tested server, the small NV–¹³C and PL6 matrices favored CPU execution;
large MoS₂ pCCE matrices favored CUDA/auto with about four MPI ranks per GPU.
P1 sK2 favored CUDA at two MPI ranks, but CPU at 16 ranks was faster than the
best CUDA configuration. The dimension-only `auto` threshold is a heuristic;
use the recorded bath-specific timings to choose the backend and MPI count.
The pCCE validation report ((internal notes)) covers the subsequent pCCE
extension and the supplied PL6 archive comparison.

The Python tools use only the standard library, read the existing regression
inputs, and require a new output directory so previous results are preserved.

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

Add `--mpirun /path/to/mpirun --ranks 2 --backend auto --min-dim 16 --batch 7`
to exercise MPI, mixed CPU/GPU selection and batch tails. `--method gcce --savemode cluster_density_matrix`
checks exported full density matrices as well as coherence. With `--repeats 3`,
the benchmark alternates execution order, reports median wall times and the
largest observed errors, and retains every run's timings and profile. Exit code 77 from
the standalone CUDA numerical test means GPU initialization was unavailable;
it is not a numerical pass.

See cCCE measurements ((internal notes)) and gCCE measurements ((internal notes))
for accuracy, timing and the limits of the tested hardware and inputs.
