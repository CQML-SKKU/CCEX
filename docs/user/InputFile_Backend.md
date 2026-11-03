# Execution tags
Where and how the cluster dynamics is computed, and what the run reports about it: the CPU or CUDA backend, the gCCE
propagation and the timing. These keys are the `execution` block; at the top level they are refused with their place.
Two exceptions: `probe_memory_gb` is in `cluster.adaptive`, and `output.cluster_contributions` is with the
[Output tags](InputFile_Output.md#tag_cluster_contributions).

```
"execution": {"backend": "auto", "gpu_options": {"min_dim": 768}, "propagator": "eigen", "evolution": "vector",
              "evolution_update": "recursive", "gcce_refinement": {"max_sweeps": 32},
              "column_evolution": {"min_dim": 256, "columns": 64}, "profile": false}
```

- [execution.backend](#tag_backend)
- [execution.gpu_options](#tag_gpu_options) (min_dim, packed_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb, warp_max_dim, ranks_per_gpu)
- [execution.propagator](#tag_propagator)
- [execution.evolution](#tag_evolution)
- [execution.evolution_update](#tag_evolution_update)
- [execution.column_evolution](#tag_column_evolution) (on by default; `enabled`, `min_dim`, `columns`)
- [execution.gcce_refinement](#tag_gcce_refinement) (writing the block turns it on; `enabled`, [max_sweeps](#tag_gcce_refinement_max_sweeps))
- [execution.profile](#tag_profile)
- [cluster.adaptive.probe_memory_gb](#tag_probe_memory_gb)
- [output.cluster_contributions](InputFile_Output.md#tag_cluster_contributions) (in [Output tags](InputFile_Output.md))
- [The memory check](#memory_check) (no key)

<a id="memory_check"></a>
## The memory check before the calculation (no input key)
After the clusters are made and before the first one is calculated, every rank estimates its peak memory. The estimate
is what the rank holds then (VmRSS) plus what its largest cluster needs on the path the run takes for it. D is the
cluster's dimension (with the qubit for gCCE and megCCE), and one dense complex matrix is 16 D² bytes:

| Path | Matrices at the peak |
|---|---|
| cCCE, one state (`s_range`) | 15, and with `column_evolution` (D ≥ its `min_dim`) the column blocks: 2 per distinct interval and 4 more, D × `columns` each |
| cCCE ensemble | Ramsey 16, Hahn echo 18, CPMG 25, `Manual` 19 + 4 (n + 1) |
| gCCE | `evolution` `matrix` 18 + n + (n + 1), `vector` 11 + n; `savemode` `cluster_density_matrix` nstep + 2 more |
| meCCE, megCCE | 14 and 10 + n + (distinct intervals) Liouville-space matrices of 16 D⁴ bytes |
| `propagator` `sparse` | the Krylov basis and the Hamiltonians in compressed rows (no dense matrix) |
| GPU clusters | the host side: 3 (cCCE) or 5 (gCCE), with `cluster_density_matrix` nstep + the time batch more |

Notes on the counts:

- n = `pulse.count`.
- Clusters of bath dimension 2 and 4 count nothing.
- The libraries add 8 MB.
- `quantity` `dm` adds four nstep × qdim² buffers.

The counts round up the peak resident memory (VmHWM) measured in local runs of every method; the estimates came out
3-16 % above the measured peak. <br>
The estimate is compared with the memory one rank of the node may use:

- under Slurm, `SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE`, divided by the ranks on the node;
- otherwise `MemAvailable` of `/proc/meminfo` divided by the ranks on the node, plus what the rank already holds.

The run log prints one line,

```
      Memory estimate    : 0.316 GB per rank (rank 0: 0.0544 GB held + its largest cluster, dimension 1024), 20.4 GB available per rank (MemAvailable / 1 rank on the node)
```

A run whose estimate exceeds what a rank may use stops there, on every rank, before anything is calculated:
`Error: not enough memory: the largest cluster (dimension D) needs about X GB per rank, and Y GB is available per rank
(...). Run fewer ranks per node, a smaller order or smaller clusters (rdip, cluster.pcce or cluster.bundle max_dim)`.
Where they apply, the message adds `quantity coherence instead of dm` and
`output.savemode coherence instead of cluster_density_matrix`. <br>
Not checked: `output.cluster_only` (no dynamics) and the adaptive probe and factors steps, which have
[their own budget](#tag_probe_memory_gb). The GPU path checks the GPU's own memory and falls back to the CPU under
`backend` `auto`.

<a id="tag_backend"></a>
## execution.backend
<b>Description:</b><br>
Where the dynamics of each cluster is computed. <br>
`cuda` and `auto` need the CUDA build, `bin/main_gpu.out` from `bash do_compile.sh gpu`
([GPU execution](../methods/gpu.md)). In a CPU-only binary `auto` runs everything on the CPU and `cuda` is refused. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `cpu` | the CPU for every cluster the packed kernel does not take (all of them while `packed_max_dim` is 0, its default here) |
| `cuda` | CUDA for every cCCE (`cluster.mode` standard, pcce or bundle) and gCCE (`propagator` eigen) cluster the packed kernel does not take; the run stops if the device or the build is missing |
| `auto` | of the clusters the packed kernel does not take, CUDA for those whose Hilbert dimension is at least [gpu_options](#tag_gpu_options).`min_dim`, the CPU for the rest; the CPU also when the GPU cannot be opened, and when it runs out of memory (that cluster and, on that rank, every later one of its dimension or larger; a `[gpu-memory-fallback]` line in the log) |

Whatever the backend, the packed kernel ([gpu_options](#tag_gpu_options).`packed_max_dim`) takes its small cCCE
clusters first; the table applies to the clusters left. <br>
An ensemble (no `s_range`) `Manual` sequence is refused with `cuda` and runs as `cpu` with `auto`. This is decided
before the packed plan.

<b>Default:</b> `cpu`

<b>Example:</b> ``` "execution" : {"backend": "auto"} ```

<a id="tag_gpu_options"></a>
## execution.gpu_options
<b>Description:</b><br>
The GPU settings, in one block. Every key is optional; an unknown key is refused. The flat keys of older inputs are
refused with a message naming the key in this block ([old names](InputFile.md#OldNames)). <br>
Without the block every default below holds. The GPU build with `backend` `cuda` or `auto` then runs:

- the small cCCE clusters in the packed kernel (`packed_max_dim`, below);
- the gCCE and cCCE clusters of dimension up to 27 in the warp kernel (`warp_max_dim`, below).

<b>Automatic values.</b> These apply to cCCE (`method` cCCE, any `cluster.mode`) in `bin/main_gpu.out` with `backend`
`auto` or `cuda` and a GPU visible. The keys `min_dim`, `evolution_batch_size`, `cluster_batch_mb` and `ranks_per_gpu`
that the input leaves out are set at the start of the run:

- `ranks_per_gpu` 4 when a node has more than 4 ranks a GPU (`backend` `auto`, `cluster.mode` standard);
- `min_dim` 256 with `ranks_per_gpu`, given or automatic. The timing of `ranks_per_gpu` then keeps on the CPU every
  dimension the GPU is slower at;
- `evolution_batch_size` 32;
- `cluster_batch_mb` two thirds of the smallest card's memory over the ranks that open a card (at least 64).

A `[gpu-auto]` line in the log gives the values. A key the input gives is used as given. gCCE keeps the fixed defaults;
the defaults in the table are those of all other runs. <br>
On one cCCE input (RTX PRO 6000, 32 ranks, MPS, one run each) `ranks_per_gpu` 4 with batches of 32 in 16384 MiB
took 490 s, against 616 s on the CPU alone and 887 s for `backend` `auto` with the fixed defaults; `ranks_per_gpu` 8
and 16 were slower (526 and 608 s). On a GH200 node `ranks_per_gpu` 4 with batches of 32 in 8192 MiB took 136 s
against 170 s with the fixed batch defaults. <br>

<b>Data Type:</b> Object

<b>Keys:</b><br>
| Key | Type, default | Meaning |
|---|---|---|
| `min_dim` | positive integer, `768` (cCCE with `ranks_per_gpu`: `256`, automatic values above) | with `backend` = `auto`, the smallest cluster Hilbert dimension that goes to the GPU; smaller clusters stay on the CPU. Ignored by `cpu` and `cuda`. 768: on an RTX PRO 6000 with 32 ranks the GPU was faster from about 670-750 on in every path (cCCE and gCCE, single state and ensemble) and slower at 648 and below |
| `packed_max_dim` | `0`, `2` or `4`; `4` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE (`cluster.mode` standard or bundle), else `0` | the packed kernel (below): cCCE clusters of bath dimension up to this value on the GPU, many at once. `0` is off |
| `evolution_batch_size` | positive integer, `8` (cCCE on the GPU: `32`, automatic values above) | the largest number of clusters of one dimension taken together on the GPU (`backend` `cuda`/`auto`): cCCE diagonalizes and evolves them together, gCCE evolves them together; `1` takes one cluster at a time. With `cluster_batch_mb` above 64, a rank that uses the GPU goes through its clusters by size, so that clusters of one dimension follow each other |
| `evolution_batch_steps` | positive integer, `64` | the largest number of time points a GPU evolution processes together |
| `cluster_batch_mb` | positive integer, `64` (MiB per rank; cCCE on the GPU: two thirds of the card over its ranks, automatic values above) | the budget that limits how many clusters one GPU batch may hold: the evolution buffers and, for cCCE, the eigensolver's workspace, about 14 D² complex numbers a cluster. With 64 a cCCE batch holds 2 clusters at D = 256 and one from about D = 310 on |
| `cache_mb` | nonnegative integer, `64` (MiB per rank) | GPU memory that freed allocations are kept in for reuse; `0` disables the cache |
| `warp_max_dim` | `0` to `32`; `27` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE or gCCE, else `0` | the warp kernel (below): gCCE and cCCE clusters of a compiled total dimension up to this value on the GPU, one per warp. `0` is off |
| `ranks_per_gpu` | nonnegative integer, `0` (cCCE, `cluster.mode` standard, `backend` `auto`, more than 4 ranks a GPU on a node: `4`, automatic values above) | with `backend` = `auto` and `method` cCCE (`cluster.mode` standard), only this many ranks per GPU open it, and the clusters of `min_dim` and above go to them (below). `0`: every rank may open the GPU |

<b>The packed kernel (`packed_max_dim`).</b> The `cuda` and `auto` backends send each cluster to the GPU on its own. For
the dimension-2 and 4 clusters of a P1 bath that is about 400 times slower than the CPU. <br>
With `packed_max_dim` = 4 (or 2), every cCCE cluster of bath dimension 2 and 4 (or 2) is computed by one GPU thread:
both eigensystems by Jacobi rotations and the whole time chain, thousands of clusters per kernel. How the work is
shared:

- Each rank still builds its clusters' Hamiltonians and multiplies their coherences into the product.
- One process per node (node rank 0, in a thread of its own) hands the packed clusters of all the node's ranks to the
  GPU through node shared memory. The packed kernel thus opens one CUDA context per node.
- Node rank 0 opens it as soon as `packed_max_dim` > 0, also with `backend` = `cpu`. Ranks that send larger clusters to
  the GPU under `cuda`/`auto` open their own.
- Larger clusters follow `backend` (and `evolution_batch_size`).

A GPU that cannot be opened turns the packed kernel off under `backend` = `auto`, and every cluster runs on the CPU.
With `cuda`, or an explicit `packed_max_dim` with `cpu`, the run stops. <br>
The packed plan is chosen after the general input checks, for instance the refusal of an ensemble `Manual` sequence
with `cuda` above. <br>
It covers:

- single state and ensemble;
- every pulse sequence with at most 64 pulses and at most 4 distinct interval lengths (Hahn, CPMG, Ramsey, most
  `Manual` ones); other sequences stay on the CPU;
- `evolution_update` `exact` and `recursive`; the kernel follows the run's value.

It needs `bin/main_gpu.out` and `method` = `cCCE` with `cluster.mode` standard or bundle; otherwise it is refused. The
results agree with the CPU within the GPU tolerance $10^{-8} + 2\times10^{-6}\,|\mathrm{CPU}|$, not bit for bit,
since the eigensolver differs. The standing check is `test_small_cluster_cuda` in `make gpu-test`. <br>
<b>Default:</b> on (`4`) in `bin/main_gpu.out` with `backend` = `cuda` or `auto` and `method` cCCE (quantity
coherence). Off (`0`) with `backend` = `cpu`, with any other method, and always in the CPU build `bin/main.out`; there
nothing goes to the GPU and only an explicit nonzero value is refused. An explicit value works with any `backend`; with
`cpu` the packed clusters are then the only GPU work. The run log prints the value and why
(`chosen by : default (execution.backend=auto, method=CCE)`, `default (CPU build)`, ... or `given`). <br>
<b>When it pays.</b> Starting it costs a fixed part of a second: the CUDA context and page-locking the node shared
memory. This is done while the bath is read and the clusters are made, and released when no rank has a cluster it
takes; until it is up the ranks compute on the CPU. A run of a few seconds therefore gains nothing. On production P1
cCCE2 inputs (32 ranks + one RTX PRO 6000) the packed kernel brought the whole run to 0.59-0.62 of the CPU binary's
time for CR60 σ4, 0.84 for CR0 with 30 states and 0.95 for the CR0 ensemble, a 1.8 s run. Runs under 1 s finished
before the GPU was up and ran entirely on the CPU. <br>
<b>When they stay off.</b> This applies when neither `packed_max_dim` nor `warp_max_dim` is given. The packed and warp
kernels are then not started when the cluster evaluations they would take are fewer than 20000 per rank (640000 for 32
ranks). Evaluations are clusters × states: each cluster once per state, with the states of `s_range` or 1 for the
ensemble, counted over all ranks from order 1 up. The reason: on a full node starting the card costs about 0.5 s, and
the card saves about 4.4 rank-µs per packed evaluation and 21 per warp one. <br>
The log says so (`[gpu-packed] off: N cluster evaluations (clusters x states) for the packed kernel, fewer than the
threshold ...`). Every cluster is then computed on the CPU path, byte-identical to the run with both keys set to 0. <br>
The card is opened early, beside the bath reading, before the count is known, and is closed again at that point. The
early start is skipped only where the input alone rules the kernels out:

- the packed kernel's plan, the ensemble rule below, or the warp kernel's run conditions such as `gcce_refinement`,
  `propagator`, `pulse.per_qubit`, `output.savemode` `cluster_density_matrix`;
- cCCE whenever the input turns the packed kernel off (`packed_max_dim` 0, its plan, the ensemble rule). The warp
  kernel's clusters (bath dimension 8 and 16) are rare there, so the card starts only when some rank has one (the log
  says `[gpu-warp] no early start: ...`).

A given value of either key starts the card whatever the count. The warp kernel's part of the node shared memory is
allocated only when some rank has a cluster of a warp dimension; a cCCE run with clusters of bath dimension 2 and 4
keeps only the packed kernel's memory. <br>
<b>The packed kernel and the CPU batch.</b> The CPU binary computes the cCCE clusters of bath dimension 2 and 4 several
at a time in SIMD lanes, with the same bits as one at a time (the "CPU batch"). In `bin/main_gpu.out` the CPU batch
also takes every such cluster the GPU does not: the ones met while the card starts, those of a run whose packed kernel
is off, and the ones the card gives back. <br>
The CPU batch takes a run with cCCE, one qubit, coherence, and either one state with at most 16 intervals or the
ensemble Hahn echo; not `backend` `cuda`. For such a run, with `packed_max_dim` not given, the packed kernel is used
only where it beats the CPU batch:

- **one state**: when the clusters of bath dimension 2 and 4 are at least 40000 per rank per state (each state
  evaluates them all; summed over the ranks and divided by their number), else off;
- **the ensemble**: off.

On a GPU node with 32 ranks and MPS the card costs about 0.105 s per state and rank and saves 2.65 µs per evaluation,
a break-even of about 39.6 k; the ensemble gained nothing at 25-178 k. <br>
A run the CPU batch does not take keeps the threshold of the paragraph above. The warp kernel is unchanged, and a given
`packed_max_dim` is kept as given. The log says which case applies:

- `[gpu-packed] on: one state, at least 40000 evaluations ...`
- `[gpu-packed] off: one state, N evaluations ... fewer than 40000 - the CPU batch (...) computes them`
- `[gpu-packed] off: the ensemble - the CPU batch ...`
- `[gpu-packed] the CPU batch does not take this run (why)`

Where the packed kernel is off this way, `main_gpu.out` writes the bytes `main.out` writes. <br>
For one state in a standard cCCE run (not `cluster.mode` bundle, no `rdip_ladder`), the card is not opened before the
bath is read. Once the bath is in, the bath spins and the pairs closer than `rdip` of bath dimension 2 and 4 (the
clusters this threshold counts) are counted, and the card starts only if they reach the threshold. Below it the log
says `[gpu-packed] no early start: one state, N clusters ...`. <br>
<b>The CPU batch checks itself.</b> At the first run that could use it, each process computes a few made-up clusters
through the CPU batch and one at a time: bath dimension 2 and 4, one state and ensemble, exact and recursive, Ramsey,
Hahn echo and CPMG. This takes about 1 ms. If any bit differs, the batch is off in that run on every rank, and rank 0
prints `[cpu-batch] off: ... does not reproduce the scalar path in this build (...)`. Causes are a build with flags the
batch cannot see, such as `CXXFLAGS_EXTRA=-ffp-contract=off`, or another CPU's tuning. The results are then those of
the one-at-a-time path. `CCEX_CPU_BATCH=0` turns the batch off by hand, for comparisons. <br>

<b>The warp kernel (`warp_max_dim`).</b> For clusters too large for the packed kernel's one-thread-per-cluster
registers but far below `min_dim`. One cluster runs per warp of 32 GPU threads, its matrices in shared memory: the
eigensystem by a parallel Jacobi, then each thread a different time point. It uses the same node feeder and shared
memory as the packed kernel. <br>
Compiled dimensions; the clusters of a compiled dimension up to `warp_max_dim` go to it:

- gCCE with a qubit of dimension 3: total dimension 6, 9, 12, 18, 21, 24, 27 (bath dimension 2, 3, 4, 6, 7, 8, 9);
- cCCE of bath dimension 8 and 16 (2 and 4 stay with the packed kernel).

Supported runs:

- gCCE: quantity coherence or dm, single state (`evolution` `vector` or `matrix`) and ensemble (`matrix`), `propagator`
  eigen, one qubit;
- cCCE: as the packed kernel;
- pulse plans with at most 4 distinct interval lengths (Hahn, CPMG, Ramsey).

Phases are always `exact` there; the log's `chosen by` of `evolution_update` says "warp kernel included". The CPU's
default `recursive` differs from `exact` by about $10^{-11}$. <br>
<b>Default:</b> `27` (every compiled dimension) in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE
(quantity coherence) or gCCE. `0` with `backend` `cpu`, another method, and always in the CPU build. The log prints the
value and why (`chosen by : default (execution.backend=auto, method=gCCE)`, `default (CPU build)`, ... or `given`). <br>
In a run the kernel does not take, it switches itself off with one `[gpu-warp] off: ...` line and the run goes on
without it. Such runs have `gcce_refinement`, `pulse.per_qubit`, a pulse detuning, `output.savemode`
`cluster_density_matrix`, another qubit dimension, or more than 4 distinct intervals. <br>
A given nonzero value is refused in the CPU build, with `backend` `cpu`, with a method other than cCCE or gCCE, and
when no compiled dimension is up to it. <br>
<b>Results.</b> Within the GPU tolerance $10^{-8} + 2\times10^{-6}\,|\mathrm{CPU}|$ of the CPU, not bit-identical to
the same clusters on the CPU path. The standing check is `test_warp_cluster_cuda` in `make gpu-test` (CPU-dumped
clusters of every compiled dimension). <br>
Inputs whose qubit splitting is many orders above the bath couplings amplify the eigensolver's rounding through the
cluster powers. Two cases:

- an NV at zero-field splitting with a 13C bath, CCE3: the CPU itself differs from its refined eigensystem
  (`gcce_refinement`) by about $5\times10^{-7}$, and the warp kernel by as much;
- P1 CR60 σ4 gCCE2: 14 of 2416 values differ by $2$–$5.5\times10^{-6}$ relative, all in the stretch where the product
  diverges ($|L| > 1$ and the points between), which the acceptance leaves out.

On a GPU node (32 ranks + one RTX PRO 6000, MPS) the warp kernel brought the whole run to 0.69 (P1 CR0 gCCE2), 0.60
(CR60 σ4) and 0.86 (h-BN gCCE2) of its time without it. <br>

<b>GPU ranks (`ranks_per_gpu`).</b> Under `backend` = `auto` every rank that meets a cluster of `min_dim` or above
opens a CUDA context of its own. Without MPS the ranks on a GPU take turns on it. A node may also admit only a few
contexts per GPU (16 on a GH200 node with four GPUs), which held a 288-core node to 64 ranks. <br>
With `ranks_per_gpu` = k, on every node the first k × (visible GPUs) ranks are GPU ranks (node rank r on GPU r mod
GPUs, k per GPU). The other ranks never start CUDA, so the run keeps every core. <br>
How the clusters are split:

1. Before the clusters are split over the ranks, each dimension of `min_dim` and above is timed on a made-up cluster of
   that dimension: a random Hermitian pair and state through the run's own CPU and GPU paths, with the run's settings
   and pulses. The load is as in the run. On the CPU every rank times at once (on a GH200 node a D = 900 cluster took
   1.9 s on an idle node and 5.8 s with 288 ranks). On the GPU one dimension is timed at a time, every GPU rank on it
   at once from a barrier, k on a card. Without MPS they take turns and each waits about k times as long; with MPS
   they run together. The time is thus what such a cluster holds a GPU rank in the run, with or without MPS. The log
   prints `[gpu-split] dim=... cpu=... s gpu=... s gpu/cpu=...`.
2. Largest first, each such cluster goes to the least-loaded GPU rank, counted at its cost times GPU/CPU time, or to
   the least-loaded CPU rank at its cost, whichever would finish it first.
3. Every other cluster goes to the least-loaded rank of all, GPU ranks included.

A dimension that was as slow or slower on a GPU rank, or that the GPU could not take in the timing (out of its memory),
stays on the CPU on every rank, GPU ranks included. <br>
The timings differ from run to run, and so does the split; the product of the clusters is then formed in another
order. The outputs agree within rounding, not bit for bit. The timing adds about the time of the largest such cluster
on a loaded CPU to the run. <br>
The packed and warp kernels already use one context per node (node rank 0, a GPU rank) and are unchanged. Needs
`bin/main_gpu.out`, `backend` = `auto` and `method` cCCE with `cluster.mode` standard; refused otherwise (bundle and
pCCE split their work themselves). <br>

<b>Batches of large cCCE clusters.</b> A GPU batch of cCCE clusters goes through the eigensolver in one call
(`cusolverDnXsyevBatched`), so the transfer, launch and wait are paid once instead of per matrix. On a GH200 node with
4 processes on a card, 32 random Hermitian matrices per call ran at 8.7 times the rate of one at a time at D = 256,
falling to 2.3 at D = 1296. <br>
It needs many clusters of one dimension and the memory for them. On that node (120 GB a card, 4 ranks a card)
`evolution_batch_size` 32 with `cluster_batch_mb` 8192 took two cCCE cases from 477 s to 328 s and from 155 s to
137 s. <br>
A given `cluster_batch_mb` is the user's to weigh: ranks on a card × `cluster_batch_mb` should stay below the card's
memory. A batch that runs out of GPU memory leaves its clusters to the CPU. <br>

<b>Example:</b> ``` "execution" : {"gpu_options": {"min_dim": 256, "packed_max_dim": 4}} ```, ``` "execution" : {"gpu_options": {"warp_max_dim": 12}} ```, ``` "execution" : {"backend": "auto", "gpu_options": {"min_dim": 256, "ranks_per_gpu": 4}} ```, ``` "execution" : {"backend": "auto", "gpu_options": {"min_dim": 256, "ranks_per_gpu": 4, "evolution_batch_size": 32, "cluster_batch_mb": 8192}} ```

<a id="tag_propagator"></a>
## execution.propagator
<b>Description:</b><br>
How the free evolution of a cluster is propagated. `eigen` and `sparse` exist for cCCE (`method` cCCE, any
`cluster.mode`) and gCCE. `expm` is for gCCE only and refused with any other method. Must be a JSON string or an
object. <br>
`sparse` applies $e^{-iH\tau}$ to the state vector without diagonalizing H:

- The cluster Hamiltonian is assembled in compressed rows from its one- and two-spin terms; no dense matrix of the
  cluster is formed.
- The exponential is a Lanczos (Krylov) expansion of dimension 30 (the whole space up to dimension 64), with an
  a-posteriori error estimate per substep. The interval is split into as many substeps as the estimate requires.
- It runs for clusters of Hilbert dimension `min_dim` and above (gCCE: qubit times bath dimension). Smaller clusters
  keep the `eigen` path, so `"sparse"` with the default `min_dim` gives the bytes of `eigen` with `evolution_update` =
  `exact` wherever no cluster reaches 4096.
- It covers one sampled state (`s_range`) of cCCE and of gCCE with `evolution` = `vector` and one qubit, with `backend`
  = `cpu` and `evolution_update` = `exact`, without `hamiltonian.cCCE.knight` and `hfmedi`.
- Refused at the start: the ensemble (no `s_range`), gCCE `evolution` = `matrix`, `backend` = `cuda` or `auto`,
  `evolution_update` = `recursive`, `knight`, `hfmedi` and more than one qubit.
- Without an `evolution_update` key `sparse` runs `exact`; the default `recursive` does not apply to it, and an
  explicit `recursive` is refused.

**When `sparse` pays:** only when the dense eigendecomposition of a cluster does not fit in memory. Its memory is
O(dim) instead of the O(dim²) matrices of `eigen`: at dim 16384 (14 spin-½) `eigen` needs about 36 GB, `sparse` about
0.2 GB. <br>
**When it does not:** at high field, or with a wide spectrum and long times. The cost of an expansion grows with the
spectral width of the Hamiltonian times the interval. The second and later intervals of a sequence start from a new
state at every time point, so the work also grows with the square of the number of time points. On the P1 probe at
500 G (Hahn echo, 151 points to 15 µs, one core, extrapolated from 11 points) dim 4096 takes about 18 h with `sparse`
against 5 min with `eigen`, and dim 16384 about 5 days against 6 h. At 0 G the same dim-1024 cluster takes as long as
with `eigen`. <br>

<b>Data Type:</b> String, or Object `{"method": ..., "min_dim": ..., "tol": ...}`

<b>Options:</b><br>
| Option | Description |
|---|---|
| `eigen` | one eigendecomposition of the cluster Hamiltonian; every step is $e^{-iH\tau} = M\,\mathrm{diag}(e^{-i\lambda\tau})\,M^\dagger$ |
| `expm` | gCCE only: a matrix exponential per time step, the original CCEX path |
| `sparse` | Krylov action of $e^{-iH\tau}$ on the state, H sparse, for clusters of dimension ≥ `min_dim` |

| Object key | Description | Default |
|---|---|---|
| `method` | `eigen`, `expm` or `sparse` (required) | |
| `min_dim` | `sparse` only: the smallest cluster dimension that takes it; a positive integer | `4096` |
| `tol` | `sparse` only: the estimated error of one evolution interval, relative to the norm of the state; 0 < `tol` < 1 | `1e-12` |

<b>Default:</b> `eigen`

<b>Example:</b> ``` "execution" : {"propagator": "eigen"} ```, ``` "execution" : {"propagator": {"method": "sparse", "min_dim": 4096, "tol": 1e-12}} ```

<a id="tag_evolution"></a>
## execution.evolution
<b>Description:</b><br>
gCCE only. What is propagated. <br>
`vector` needs a pure initial state (`s_range`), `propagator` = `eigen` or `sparse` and `method` = `gCCE`. When the key
is absent the run picks `vector` where those hold and `matrix` otherwise. <br>
An explicit `vector` that cannot run stops with an error. With any other method an explicit `vector` is refused; an
explicit `matrix` is accepted and changes nothing, since the other methods do not read the key. <br>
The two agree to rounding but are not bit-identical, so a converged series should not switch between them. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `vector` | propagate the state vector; matrix-vector products, O(n²) per step |
| `matrix` | propagate the density matrix $\rho(t) = U\rho_0 U^\dagger$, O(n³) per step; works for every configuration |

<b>Default:</b> `vector` where it can run, else `matrix`

<b>Example:</b> ``` "execution" : {"evolution": "matrix"} ```

<a id="tag_evolution_update"></a>
## execution.evolution_update
<b>Description:</b><br>
How the phases $e^{-iE\tau}$ of the eigenvalues are formed at each time point of the uniform grid (`deltat`). Must be
a JSON string.

- `exact` computes them with sin and cos at every point.
- `recursive` does so at points 0, 16, 32, … In between it multiplies the previous point's phases by one step's factor
  $e^{-iEf\,\mathrm{deltat}}$ (f the interval's fraction of the sequence), so the rounding of at most 15 products
  accumulates before the next exact point.

`recursive` exists on:

- the CPU time evolution of cCCE (`method` cCCE, any `cluster.mode`; one sampled state, `s_range`, and the ensemble, no
  `s_range`);
- the packed kernel of `bin/main_gpu.out` (`gpu_options.packed_max_dim`, cCCE clusters of bath dimension 2 and 4);
- the CPU time evolution of gCCE with `propagator` = `eigen` (`evolution` = `vector` and `matrix`).

Every pulse sequence these paths run is covered. The ensemble cCCE carries over the phases of its τ = t/2 (Hahn), t
(Ramsey), t/(2 `pulse.count`) and 2τ (CPMG), or of each `Manual` interval; the others those of each interval of the
sequence. The GPU large-cluster path (`bin/main_gpu.out`: every other cluster with `backend` = `cuda`, the clusters of
`gpu_options.min_dim` and above with `auto`) computes every phase with sin and cos, as `exact`. <br>
Without the key each path takes `recursive` where it has it and `exact` elsewhere, without a message (per path in the
GPU build):

- `backend` = `cpu`, and the CPU build `bin/main.out` with `cpu` or `auto` (all on the CPU there): `recursive` (gCCE
  `expm`: `exact`);
- `bin/main_gpu.out` with `auto`: the clusters on the CPU (below `min_dim`) and the packed kernel run `recursive`, the
  clusters on the GPU `exact`;
- `bin/main_gpu.out` with `cuda`: the packed kernel runs `recursive` and every other cluster `exact` on the GPU. With
  the packed kernel off, as for gCCE and pCCE, the value is `exact`.

The run log prints the value and why, per path in the GPU build (`chosen by : default (supported path)`,
`default (clusters on the CPU and packed kernel: recursive, clusters on the GPU (execution.backend=auto): exact)`,
`default (execution.backend=cuda: every cluster on the GPU, exact there)`, … or `given`). <br>
A given `recursive` is refused, rather than run as `exact`, where a cluster may reach a path without it:

- gCCE with `propagator` = `expm` (no phases to carry over);
- `bin/main_gpu.out` with `backend` = `cuda` or `auto` and a method the GPU takes (cCCE with `cluster.mode` standard,
  pcce or bundle, gCCE with `eigen`). Leave the key out there to get `recursive` on the CPU clusters and the packed
  kernel.

`backend` = `cuda` in the CPU build is refused in any case. <br>
On the production P1 inputs (CR0 and CR60 σ4, cCCE2 and gCCE2, 30 states, 64 ranks) `recursive` takes about 0.73–0.76
of the time of `exact` for cCCE2 and 0.82–0.84 for gCCE2. The `wiDiv` curves moved by at most $2.2\times10^{-10}$ for
cCCE2 (the outputs have 10 decimals) and $1.3\times10^{-8}$ for gCCE2 at points with $|L| \le 1$. The results are not
bit-identical to `exact`, so a converged series should not switch between the two. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `exact` | sin and cos at every time point |
| `recursive` | sin and cos every 16th time point, the previous point's phases times one step's factor in between |

<b>Default:</b> `recursive` where it exists (CPU cCCE, the packed kernel, CPU gCCE with `propagator` = `eigen`), `exact` elsewhere (the GPU large-cluster path, gCCE `expm`; meCCE, megCCE: `exact`, and a given `recursive` is refused); per path in `bin/main_gpu.out`

<b>Example:</b> ``` "execution" : {"evolution_update": "exact"} ```

<a id="tag_column_evolution"></a>
## execution.column_evolution
<b>Description:</b><br>
The one-state cCCE time evolution of the large clusters with the time points as matrix columns. <br>
The usual path carries one row and one column vector per time point through the intervals of the sequence: two
matrix-vector products per interval and time point, each reading the cluster's D × D matrix again. With this block on,
`columns` time points go together: each interval is one matrix-matrix product (the library's zgemm in a build with
MKL), and the matrix is read once per block. <br>
The multiplications are the same but their summation order differs, so the outputs move in the last bits
($10^{-15}$ relative in the unit check). The phases follow the usual path's rule, `exact` and `recursive` alike
(`recursive`: sin and cos at steps 0, 16, 32, ..., the products carried across the blocks). The carried products may
differ in the last bit where the compiler fuses the complex multiplication otherwise. Every pulse sequence of the
one-state cCCE (Ramsey, Hahn echo, CPMG, `Manual`) is covered. <br>
Where it applies: on the CPU only, for one sampled state (`s_range`) of `method` cCCE (every `cluster.mode`), clusters of
bath dimension `min_dim` and above. Unchanged whatever the block says: bath dimension 2 and 4 (the fixed-size path and
the CPU batch), the ensemble (no `s_range`), `propagator` `sparse`, gCCE, and the clusters on the GPU. With the block
`enabled` false every output is the bytes of the run without it. <br>
Extra memory per rank: the phases of one block for each distinct interval and four D × `columns` matrices,
(2 × intervals + 4) × D × `columns` × 16 bytes (D 2048, CPMG-8, 64 columns: about 12 MB). <br>
**On by default, `min_dim` 256 and `columns` 64:** on a Xeon Gold 6530 node (64 ranks; pCCE clusters of dimension
1024, 512 and 256; Hahn and CPMG-8, nstep 200 and 1000; 12 inputs) it took with MKL 0.40 of the time (geometric mean;
up to 8.5 times faster for dimension 1024, 0.65–0.97 for dimension 256), with Eigen alone (`CCEX_BLAS` off) 0.61. The
curves agree with the run without it to $10^{-10}$, the last printed digit. Peak memory per rank moved by −8 to
+11 MB. `min_dim` 384 left the dimension-256 cases unchanged (1.01) and 768 also the dimension-512 ones (1.00);
`columns` 32 and 128 were within a few per cent of 64, 128 with up to +27 MB. <br>
How the block is read:

- without the block it is on with these values;
- a block without `enabled` is on with the values it gives (the others as above);
- `{"enabled": false}` turns it off (the usual one-vector path).

`enabled` must be a JSON boolean, `min_dim` and `columns` positive integers. The run log prints
`column_evolution.min_dim` and `.columns` when it is on. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "min_dim": Integer, "columns": Integer}`

<b>Keys:</b><br>
| Key | Type, default | Meaning |
|---|---|---|
| `enabled` | boolean, `true` (also when the block is absent) | on or off |
| `min_dim` | positive integer, `256` | the smallest bath dimension of a cluster that takes it |
| `columns` | positive integer, `64` | the time points per block |

<b>Default:</b> on, `{"enabled": true, "min_dim": 256, "columns": 64}`

<b>Example:</b> ``` "execution" : {"column_evolution": {"enabled": false}} ``` (off)

<a id="tag_probe_memory_gb"></a>
## cluster.adaptive.probe_memory_gb
<b>Description:</b><br>
The memory, in GB per node, that the probes running at once may hold, in the probe step of the adaptive clustering
(`cluster.mode` = `bundle` with the `cluster.adaptive` block). It is a key of the `cluster.adaptive` block; the
refinement driver copies it into the step inputs beside their `work_file`. <br>
The probes of a bath state go to the ranks in turn, and the ranks do not wait for each other between states, so the
large probes of several states can run together. A probe waits while the probes already running on its node would
exceed this budget; a probe larger than the budget runs alone. <br>
A probe is estimated on the path it takes (`estimate` in `src/clustering/bundle/bundle_work.cpp`), in dense dim × dim
complex matrices of 16 dim² bytes each (n = `pulse.count`):

- CPU, gCCE: 6 + n with `evolution` = `vector`, 7 + 2n with `matrix`;
- CPU, cCCE: 8 for a sampled state, 23 for the ensemble, 10 + 4(n + 1) for an ensemble `Manual` sequence;
- GPU: what the device allocates for the time batch and evolution;
- with `output.savemode` = `cluster_density_matrix`, also every time's full density matrix.

Probes under 64 MiB do not count. <br>
0 takes 80 % of the memory the job may use on the node (its RAM, or less if a cgroup limits the job), and, when the
run may put clusters on the GPU, 80 % of the smallest GPU. The run prints the budget (`[probe-budget]`) and how long
each rank waited. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> `0` (automatic)

<b>Example:</b> ``` "cluster" : {"adaptive": {..., "probe_memory_gb": 60}} ```

<a id="tag_gcce_refinement"></a>
## execution.gcce_refinement
<b>Description:</b><br>
gCCE with `propagator` = `eigen`: refine the eigenvalues and eigenvectors the solver returns by compensated correction
sweeps before they are used. <br>
A block `{"enabled": true|false, "max_sweeps": n}`; a plain boolean is refused. Writing the block turns the refinement
on: a block without `enabled` is on, as an `auto_convergence` block is. `{"enabled": false}` keeps a written block off.
`enabled` must be a JSON boolean. <br>
Off by default, because it takes most of the setup time of a gCCE cluster (dim 768: about 4 s of 4.8 s on the CPU)
and changes little. On the NV–¹³C gCCE4 case it was made for, and on P1 gCCE3, it changed the error-corrected final
coherence by at most $6.7\times10^{-7}$. The raw `wiDiv` differs by more at points where it diverges with either
setting (up to 12.7 on NV–¹³C and $2.75\times10^{11}$ on P1 gCCE3), so a $10^{-6}$ bound on every raw point was not
met. <br>
Turn it on where the bath energies are very small beside the qubit's and the curve is read to better than $10^{-6}$. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "max_sweeps": Integer}`

<b>Default:</b> off (absent)

<b>Example:</b> ``` "execution" : {"gcce_refinement": {}} ``` (on: the block is written)

<a id="tag_gcce_refinement_max_sweeps"></a>
## execution.gcce_refinement.max_sweeps
<b>Description:</b><br>
The largest number of refinement sweeps per cluster. A cluster that has not converged within it stops the run with an
error naming its dimension. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `32`

<b>Example:</b> ``` "execution" : {"gcce_refinement": {"max_sweeps": 32}} ```

<a id="tag_profile"></a>
## execution.profile
<b>Description:</b><br>
Print timings at the end of the run. Must be a JSON boolean.

- Every rank writes one `[profile]` line per method, backend and cluster dimension, with the number of clusters and the
  seconds spent building, setting up and evolving them.
- With CUDA also `[gpu-batch]` and `[gpu-memory]` lines.
- In a plain cCCE/gCCE run rank 0 also prints one `[phase]` line per rank: the wall seconds of the state loop split
  into state preparation, cluster creation, cluster calculation, product assembly, freeing, waiting at the barriers,
  the reduction and saving.

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "execution" : {"profile": true} ```
