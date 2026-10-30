# Execution tags
Where and how the cluster dynamics is computed, and what the run reports
about it: the CPU or CUDA backend, the gCCE propagation and the timing. These keys are the `execution` block (at the top level they are refused with
their place), except `probe_memory_gb` (in `cluster.adaptive`); `output.cluster_contributions` is with the
[Output tags](InputFile_Output.md#tag_cluster_contributions):

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
After the clusters are made and before the first one is calculated, every rank estimates its peak memory: what it holds
then (VmRSS) and what its largest cluster needs on the path the run takes for it, with D the cluster's dimension (with
the qubit for gCCE and megCCE) and one dense complex matrix 16 D² bytes:

| Path | Matrices at the peak |
|---|---|
| cCCE, one state (`s_range`) | 15, and with `column_evolution` (D ≥ its `min_dim`) the column blocks: 2 per distinct interval and 4 more, D × `columns` each |
| cCCE ensemble | Ramsey 16, Hahn echo 18, CPMG 25, `Manual` 19 + 4 (n + 1) |
| gCCE | `evolution` `matrix` 18 + n + (n + 1), `vector` 11 + n; `savemode` `cluster_density_matrix` nstep + 2 more |
| meCCE, megCCE | 14 and 10 + n + (distinct intervals) Liouville-space matrices of 16 D⁴ bytes |
| `propagator` `sparse` | the Krylov basis and the Hamiltonians in compressed rows (no dense matrix) |
| GPU clusters | the host side: 3 (cCCE) or 5 (gCCE), with `cluster_density_matrix` nstep + the time batch more |

(n = `pulse.count`; clusters of bath dimension 2 and 4 count nothing; 8 MB more for the libraries; `quantity` `dm` adds
four nstep × qdim² buffers.) The counts were measured against the peak resident memory (VmHWM) of local runs with the MKL
build (`do_compile.sh wsl-gnu`) and round up: the estimates came out 3-16 % above the measured peak (cCCE D 512-2048,
gCCE D 384-1536, pCCE with 17 intervals, meCCE and megCCE D 16-32). The estimate
is compared with the memory one rank of the node may use: `SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE`,
divided by the ranks on the node; without Slurm's, `MemAvailable` of `/proc/meminfo` divided by the ranks on the node, plus
what the rank already holds. The run log prints one line,

```
      Memory estimate    : 0.316 GB per rank (rank 0: 0.0544 GB held + its largest cluster, dimension 1024), 20.4 GB available per rank (MemAvailable / 1 rank on the node)
```

and a run whose estimate passes what a rank may use stops there, every rank, before anything is calculated:
`Error: not enough memory: the largest cluster (dimension D) needs about X GB per rank, and Y GB is available per rank
(...). Run fewer ranks per node, a smaller order or smaller clusters (rdip, cluster.pcce or cluster.bundle max_dim)`, with
`quantity coherence instead of dm` and `output.savemode coherence instead of cluster_density_matrix` where they apply.
Not checked: `output.cluster_only` (no dynamics) and the adaptive probe and factors steps (they have
[their own budget](#tag_probe_memory_gb)). The GPU's own memory is checked by the GPU path, which falls back to the CPU
under `backend` `auto`.

<a id="tag_backend"></a>
## execution.backend
<b>Description:</b><br>
Where the dynamics of each cluster is computed. <br>
`cuda` and `auto` need the CUDA build, `bin/main_gpu.out` from `bash do_compile.sh gpu`
([GPU execution](../methods/gpu.md)); in a CPU-only binary `auto` runs everything on the CPU and `cuda`
is refused. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `cpu` | the CPU for every cluster the packed kernel does not take (all of them while `packed_max_dim` is 0, its default here) |
| `cuda` | CUDA for every cCCE (`cluster.mode` standard, pcce or bundle) and gCCE (`propagator` eigen) cluster the packed kernel does not take; the run stops if the device or the build is missing |
| `auto` | of the clusters the packed kernel does not take, CUDA for those whose Hilbert dimension is at least [gpu_options](#tag_gpu_options).`min_dim`, the CPU for the rest; the CPU also when the GPU cannot be opened, and when it runs out of memory (that cluster and, on that rank, every later one of its dimension or larger; a `[gpu-memory-fallback]` line in the log) |

The packed kernel ([gpu_options](#tag_gpu_options).`packed_max_dim`) takes its small cCCE clusters first, whatever the
backend; the table is the choice for the clusters left. An ensemble (no `s_range`) `Manual` sequence is refused with
`cuda` and runs as `cpu` with `auto`; this is decided before the packed plan.

<b>Default:</b> `cpu`

<b>Example:</b> ``` "execution" : {"backend": "auto"} ```

<a id="tag_gpu_options"></a>
## execution.gpu_options
<b>Description:</b><br>
The GPU settings, in one block (the flat keys of older inputs are refused with a message naming the key in this
block; [old names](InputFile.md#OldNames)). Every key is optional; an unknown key is refused. Without the block every
default below holds: the GPU build with `backend` `cuda` or `auto` runs the small cCCE clusters in the packed kernel
(`packed_max_dim`, below) and the gCCE and cCCE clusters of dimension up to 27 in the warp kernel (`warp_max_dim`,
below). <br>

<b>Data Type:</b> Object

<b>Keys:</b><br>
| Key | Type, default | Meaning |
|---|---|---|
| `min_dim` | positive integer, `768` | with `backend` = `auto`, the smallest cluster Hilbert dimension that goes to the GPU; smaller clusters stay on the CPU. Ignored by `cpu` and `cuda`. 768: on an RTX PRO 6000 with 32 ranks the GPU was faster from about 670-750 on in every path (cCCE and gCCE, single state and ensemble) and slower at 648 and below |
| `packed_max_dim` | `0`, `2` or `4`; `4` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE (`cluster.mode` standard or bundle), else `0` | the packed kernel (below): cCCE clusters of bath dimension up to this value on the GPU, many at once. `0` is off |
| `evolution_batch_size` | positive integer, `8` | the largest number of consecutive clusters of one dimension evolved together on the GPU (`backend` `cuda`/`auto`); `1` evolves one cluster at a time |
| `evolution_batch_steps` | positive integer, `64` | the largest number of time points a GPU evolution processes together |
| `cluster_batch_mb` | positive integer, `64` (MiB per rank) | the evolution-buffer budget that limits how many clusters one GPU batch may hold |
| `cache_mb` | nonnegative integer, `64` (MiB per rank) | GPU memory that freed allocations are kept in for reuse; `0` disables the cache |
| `warp_max_dim` | `0` to `32`; `27` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE or gCCE, else `0` | the warp kernel (below): gCCE and cCCE clusters of a compiled total dimension up to this value on the GPU, one per warp. `0` is off |
| `ranks_per_gpu` | nonnegative integer, `0` | with `backend` = `auto` and `method` cCCE (`cluster.mode` standard), only this many ranks per GPU open it, and the clusters of `min_dim` and above go to them (below). `0`: every rank may open the GPU |

<b>The packed kernel (`packed_max_dim`).</b> The `cuda`/`auto` backends send a cluster to the GPU on its own, which
for the dimension-2 and 4 clusters of a P1 bath costs far more than the CPU (about 400 times slower).
With `packed_max_dim` = 4 (or 2), every cCCE cluster of bath dimension 2 and 4 (or 2) is instead computed by one GPU
thread — both eigensystems by Jacobi rotations and the whole time chain — thousands of clusters per kernel. Each rank
still builds its clusters' Hamiltonians and multiplies their coherences into the product as usual; one process per node
(node rank 0, in a thread of its own) hands the packed clusters of all the node's ranks to the GPU through node shared
memory, so the packed kernel opens one CUDA context per node (node rank 0 opens it as soon as `packed_max_dim` > 0,
also with `backend` = `cpu`; ranks that send larger clusters to the GPU under `cuda`/`auto` open their own). Larger
clusters follow `backend` (and `evolution_batch_size`). With `backend` = `auto` a GPU that cannot be opened
turns the packed kernel off and every cluster runs on the CPU; with `cuda` (or an explicit `packed_max_dim` with `cpu`)
the run stops. The packed plan is chosen after the general input checks (for instance the refusal
of an ensemble `Manual` sequence with `cuda` above). Single state and ensemble, every pulse sequence with at most 64
pulses and at most 4 distinct interval lengths (Hahn, CPMG, Ramsey, most `Manual` ones; others stay on the CPU), `evolution_update` `exact` and `recursive` (the kernel follows the run's value). Needs `bin/main_gpu.out` and
`method` = `cCCE` with `cluster.mode` standard or bundle (refused otherwise). The results agree with the CPU within the GPU tolerance 1e-8 + 2e-6 |CPU|, not bit
for bit (the eigensolver differs); the standing check is `test_small_cluster_cuda` in `make gpu-test`.
<b>Default:</b> on (`4`) in `bin/main_gpu.out` with `backend` = `cuda` or `auto` and `method` cCCE
(quantity coherence); off (`0`) with `backend` = `cpu`, with any other method, and always in the CPU build `bin/main.out`,
where nothing goes to the GPU and only an explicit nonzero value is refused. An explicit value works with any `backend`
(`cpu` included: then the packed clusters are the only GPU work). The run log prints the value and why
(`chosen by : default (execution.backend=auto, method=CCE)`, `default (CPU build)`, ... or `given`). <br>
<b>When it pays.</b> Starting it costs a fixed part of a second (the CUDA context and page-locking the node shared memory,
done while the bath is read and the clusters are made, and released when no rank has a cluster it takes; until it is up
the ranks compute on the CPU), so a run of a few
seconds gains nothing. Measured on a GPU node (32 ranks + one RTX PRO 6000, the production P1 cCCE2 inputs, whole
run against the CPU binary on the same node, median of AB BA AB): CR60 σ4 30 states 70.7 → 44.2 s (0.62), its ensemble
6.5 → 3.9 s (0.59), CR0 30 states 13.6 → 11.4 s (0.84), but the CR0 ensemble, a 1.8 s run, 1.80 → 1.71 s (0.95);
(internal notes). In that measurement the 12 small fixture runs (under 1 s) finished before the GPU was
up and ran entirely on the CPU. <br>
<b>When they stay off.</b> With neither `packed_max_dim` nor `warp_max_dim` given, the packed
and warp kernels are not started at all when the cluster evaluations they would take (clusters × states: each cluster
once per state, the states of `s_range` or 1 for the ensemble), counted over all ranks from order 1 up, are fewer than 20000 per
rank (640000 for 32 ranks): on a full node starting the card costs about 0.5 s, and the card saves about 4.4 rank-µs per
packed evaluation and 21 per warp one, so such runs gain little or lose (the reasoning and the numbers:
(internal notes)). The log says so (`[gpu-packed] off: N cluster evaluations
(clusters x states) for the packed kernel, fewer than the threshold ...`), and every cluster is computed on the CPU path,
byte-identical to the run with both set to 0. The card's early start (opened beside the bath reading, before the count is
known) is closed again at that point; it is not avoided, except where the input alone rules the kernels out
(the packed kernel's plan, the ensemble rule below, the warp kernel's run conditions such as
`gcce_refinement`, `propagator`, `pulse.per_qubit`, `output.savemode` `cluster_density_matrix`), and for cCCE whenever the packed
kernel is off by the input (`packed_max_dim` 0, its plan, the ensemble rule): the warp kernel's clusters (bath
dimension 8 and 16) are rare there, so the card starts only when some rank has one (the log says
`[gpu-warp] no early start: ...`). A given value of either key starts the card whatever the count. The warp kernel's part of the
node shared memory is allocated only when some rank has a cluster of a warp dimension (a cCCE run with clusters of bath
dimension 2 and 4 only keeps the packed kernel's memory: P1 CR0 cCCE2 was 4.7 % slower with it). <br>
<b>The packed kernel and the CPU batch.</b> The CPU binary computes the cCCE clusters of
bath dimension 2 and 4 several at a time in SIMD lanes (the same bits as one at a time; the "CPU batch"). In
`bin/main_gpu.out` the CPU batch also takes every such cluster the GPU does not: the ones met while the card starts,
those of a run whose packed kernel is off, and the ones the card gives back. With `packed_max_dim` not given and a run the
CPU batch takes (cCCE, one qubit, coherence, one state with at most 16 intervals or the ensemble Hahn echo; not `backend`
`cuda`), the packed kernel is used only where it beats the CPU batch: **one state** - when the clusters of bath dimension
2 and 4 are at least 40000 per rank per state (each state evaluates them all; summed over the ranks and divided by their
number), else off; **the ensemble** - off. Measured on a GPU node with 32 ranks and MPS
((internal notes)): per state and rank the card loses about 0.105 s and saves 2.65 µs per
evaluation (break-even about 39.6 k; 36 k: +0.06 s, 54 k: -0.81 s); the ensemble gained nothing at 25-178 k. A run the
CPU batch does not take keeps the threshold of the paragraph above; the warp kernel is unchanged, and a given
`packed_max_dim` is kept as given. The log says which: `[gpu-packed] on: one state, at least 40000 evaluations ...`,
`[gpu-packed] off: one state, N evaluations ... fewer than 40000 - the CPU batch (...) computes them`,
`[gpu-packed] off: the ensemble - the CPU batch ...`, or `[gpu-packed] the CPU batch does not take this run (why)`.
Where the packed kernel is off this way, `main_gpu.out` writes the bytes `main.out` writes. For one state in a standard
cCCE run (not `cluster.mode` bundle, no `rdip_ladder`), the card is not opened before the bath is read: once the bath
is in, the bath spins and the pairs closer than `rdip` of bath dimension 2 and 4 (the clusters this threshold counts) are
counted, and the card starts then only if they reach it (P1 CR0 30 states had paid about 0.2 s for a card it did not
use). Below it the log says `[gpu-packed] no early start: one state, N clusters ...`. <br>
<b>The CPU batch checks itself.</b> At the first run that could use it, each process computes a few made-up clusters
(bath dimension 2 and 4, one state and ensemble, exact and recursive, Ramsey, Hahn echo and CPMG) through the CPU batch and
one at a time; if any bit differs - a build with flags the batch cannot see, such as `CXXFLAGS_EXTRA=-ffp-contract=off`,
or another CPU's tuning - the batch is off in that run (every rank alike) and rank 0 prints `[cpu-batch] off: ... does not
reproduce the scalar path in this build (...)`. About 1 ms. Results are then those of the one-at-a-time path, as
everywhere. `CCEX_CPU_BATCH=0` turns the batch off by hand (A/B). <br>

<b>The warp kernel (`warp_max_dim`).</b> Clusters too large for the packed kernel's
one-thread-per-cluster registers but far below `min_dim`: one cluster per warp of 32 GPU threads, its matrices in shared
memory, the eigensystem by a parallel Jacobi, then each thread a different time point. Through the same node feeder and
shared memory as the packed kernel. Compiled for gCCE with a qubit of dimension 3 — total dimension 6, 9, 12, 18, 21, 24,
27 (bath dimension 2, 3, 4, 6, 7, 8, 9) — and for cCCE of bath dimension 8 and 16 (2 and 4 stay with the packed kernel);
the clusters of a compiled dimension up to `warp_max_dim` go to it. gCCE: quantity coherence or dm, single state
(`evolution` `vector` or `matrix`) and ensemble (`matrix`), `propagator` eigen, one qubit; cCCE: as the packed kernel. Pulse
plans with at most 4 distinct interval lengths (Hahn, CPMG, Ramsey). Phases are always `exact` there (the log's `chosen by`
of `evolution_update` says "warp kernel included"); the CPU's default `recursive` differs from `exact` by about 1e-11.
<b>Default:</b> `27` (every compiled dimension) in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE
(quantity coherence) or gCCE; `0` with `backend` `cpu`, another method, and always in the CPU build. The log prints the
value and why (`chosen by : default (execution.backend=auto, method=gCCE)`, `default (CPU build)`, ... or `given`). Where the run is
one the kernel does not take (`gcce_refinement`, `pulse.per_qubit`, a pulse detuning, `output.savemode` `cluster_density_matrix`, another qubit
dimension, more than 4 distinct intervals) it switches itself off, with one `[gpu-warp] off: ...` line, and the run goes
on without it. A given nonzero value is refused in the CPU build, with `backend` `cpu`, with a method other than cCCE or
gCCE, and when no compiled dimension is up to it. <br>
<b>Results.</b> Within the GPU tolerance 1e-8 + 2e-6 |CPU| of the CPU, not bit-identical to the same clusters on the CPU
path. The standing check is `test_warp_cluster_cuda` in `make gpu-test` (CPU-dumped clusters of
every compiled dimension). Inputs whose qubit splitting is many orders above the bath couplings amplify the eigensolver's
rounding through the cluster powers: an NV at zero-field splitting with a 13C bath, CCE3, where the CPU itself differs from
its refined eigensystem (`gcce_refinement`) by about 5e-7 and the warp kernel by as much; and in P1 CR60 σ4 gCCE2
14 of 2416 values differ by 2–5.5e-6 relative - all in the stretch where the product diverges (|L| > 1 and the points
between), which the acceptance leaves out. Measured on a GPU node (32 ranks + one RTX PRO 6000, MPS, whole run against
the warp kernel off, AB BA AB): P1 CR0 gCCE2 14.2 → 9.8 s (0.69), CR60 σ4 60.0 → 36.0 s (0.60),
h-BN gCCE2 10.9 → 9.3 s (0.86); (internal notes). <br>

<b>GPU ranks (`ranks_per_gpu`).</b> Under `backend` = `auto` every rank that meets a cluster of `min_dim` or above opens a
CUDA context of its own. Without MPS the ranks on a GPU take turns on it, and a node may admit only a few contexts per GPU
(16 on a GH200 node with four GPUs), which held a 288-core node to 64 ranks. With `ranks_per_gpu` = k, on every node
the first k × (visible GPUs) ranks are GPU ranks (node rank r on GPU r mod GPUs, k per GPU) and the others never start
CUDA, so the run keeps every core. Before the clusters are split over the ranks, each dimension of `min_dim` and above is
timed on a made-up cluster of that dimension (a random Hermitian pair and state through the run's own CPU and GPU
paths, the run's settings and pulses), as loaded as in the run: on the CPU by every rank at once (on a GH200 node a
D = 900 cluster took 1.9 s on an idle node and 5.8 s with 288 ranks), and on the GPU one dimension at a time, every
GPU rank on it at once from a barrier, k on a card - without MPS they take turns on it and each waits about k times as long, with MPS they run together - so the time
is what such a cluster holds a GPU rank in the run, with or without MPS (the log prints
`[gpu-split] dim=... cpu=... s gpu=... s gpu/cpu=...`). Then, largest first, such a cluster goes to the least-loaded GPU
rank, counted at its cost times GPU/CPU time, or to the least-loaded CPU rank at its cost, whichever would finish it
first; every other
cluster goes to the least-loaded rank of all, GPU ranks included. The timings differ from run to run, so the split does
and the product of the clusters is formed in another order: the outputs agree within rounding, not bit for bit. A
dimension that was as slow or slower on a GPU rank, or that the GPU could not take in the timing (out of its memory),
stays with every rank on the CPU, GPU ranks included. The timing adds about the time of the largest such cluster on a loaded
CPU to the run. The packed and warp kernels already use one
context per node (node rank 0, a GPU rank) and are unchanged. Needs `bin/main_gpu.out`, `backend` = `auto` and `method`
cCCE with `cluster.mode` standard (refused otherwise; bundle and pCCE split their work themselves). <br>

<b>Example:</b> ``` "execution" : {"gpu_options": {"min_dim": 256, "packed_max_dim": 4}} ```, ``` "execution" : {"gpu_options": {"warp_max_dim": 12}} ```, ``` "execution" : {"backend": "auto", "gpu_options": {"min_dim": 256, "ranks_per_gpu": 4}} ```

<a id="tag_propagator"></a>
## execution.propagator
<b>Description:</b><br>
How the free evolution of a cluster is propagated. `eigen` and `sparse` exist for cCCE (`method` cCCE, any `cluster.mode`) and gCCE, `expm` for gCCE only (refused with any other method). <br>
`sparse` applies exp(-i H τ) to the state vector without diagonalizing H: the cluster Hamiltonian is
assembled in compressed rows from its one- and two-spin terms (no dense matrix of the cluster is formed) and the
exponential is a Lanczos (Krylov) expansion of dimension 30 (the whole space up to dimension 64) with an a-posteriori
error estimate per substep; the interval is split into as many substeps as the estimate requires. It runs for clusters of
Hilbert dimension `min_dim` and above (gCCE: qubit times bath dimension); smaller clusters keep the `eigen` path, so
`"sparse"` with the default `min_dim` gives the bytes of `eigen` with `evolution_update` = `exact` wherever no cluster reaches 4096. It covers one sampled
state (`s_range`) of the cCCE and gCCE with `evolution` = `vector` and one qubit, with `backend` = `cpu` and
`evolution_update` = `exact`, without `hamiltonian.cCCE.knight` and `hfmedi`; the ensemble (no `s_range`), gCCE `evolution` = `matrix`,
`backend` = `cuda` or `auto`, `evolution_update` = `recursive`, `knight`, `hfmedi` and more than one qubit are refused at
the start. <br>
**When `sparse` pays:** only when the dense eigendecomposition of a cluster does not fit in memory. Memory is O(dim)
instead of the O(dim²) matrices of `eigen`: at dim 16384 (14 spin-½) `eigen` needs about 36 GB, `sparse` about 0.2 GB.
**When it does not:** at high field or with a wide spectrum and long times. The cost of an expansion grows with the
spectral width of the Hamiltonian times the interval, and the second and later intervals of a sequence start from a new
state at every time point, so the work also grows with the square of the number of time points. On the P1 probe at
500 G (Hahn echo, 151 points to 15 µs) dim 4096 takes about 18 h with `sparse` (11 points: 6 min, × (151/11)²) against 5 min with `eigen`, and dim 16384 about
5 days against 6 h (11 points: 40 min; one core); at 0 G the same
dim-1024 cluster takes as long as `eigen` ((internal notes)). Without an
`evolution_update` key `sparse` runs `exact` (the default `recursive` does not apply to it); an explicit `recursive` is
refused. Must be a JSON string or an object. <br>

<b>Data Type:</b> String, or Object `{"method": ..., "min_dim": ..., "tol": ...}`

<b>Options:</b><br>
| Option | Description |
|---|---|
| `eigen` | one eigendecomposition of the cluster Hamiltonian; every step is exp(-i H τ) = M diag(exp(-i λ τ)) M^† |
| `expm` | gCCE only: a matrix exponential per time step, the original CCEX path |
| `sparse` | Krylov action of exp(-i H τ) on the state, H sparse, for clusters of dimension ≥ `min_dim` |

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
gCCE only. What is propagated. `vector` needs a pure initial state (`s_range`),
`propagator` = `eigen` or `sparse` and `method` = `gCCE`. When the key is absent the run picks `vector`
where those hold and `matrix` otherwise; an explicit `vector` that cannot run stops with an
error, and with any other method an explicit `vector` is refused (an explicit `matrix` is accepted and
changes nothing: the other methods do not read the key). The two agree to rounding but
are not bit-identical, so a converged series should not switch between them. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `vector` | propagate the state vector; matrix-vector products, O(n²) per step |
| `matrix` | propagate the density matrix ρ(t) = U ρ₀ U^†, O(n³) per step; works for every configuration |

<b>Default:</b> `vector` where it can run, else `matrix`

<b>Example:</b> ``` "execution" : {"evolution": "matrix"} ```

<a id="tag_evolution_update"></a>
## execution.evolution_update
<b>Description:</b><br>
How the phases exp(-i E τ) of the eigenvalues are formed at each time point of the uniform grid (`deltat`). `exact`
computes them with sin and cos at every point. `recursive` does so at points 0, 16, 32, … and in between multiplies the
previous point's phases by one step's factor exp(-i E f `deltat`) (f the interval's fraction of the sequence), so the
rounding of at most 15 products accumulates before the next exact point. <br>
`recursive` exists on the CPU time evolution of cCCE (`method` cCCE, any `cluster.mode`; one
sampled state, `s_range`, and the ensemble, no `s_range`), in the packed kernel of
`bin/main_gpu.out` (`gpu_options.packed_max_dim`, cCCE clusters of bath dimension 2 and 4) and on the CPU time evolution
of gCCE with `propagator` = `eigen` (`evolution` = `vector` and `matrix`). Every pulse sequence these paths
run is covered: the ensemble cCCE carries over the phases of its τ = t/2 (Hahn), t (Ramsey), t/(2 `pulse.count`) and 2τ
(CPMG) or of each `Manual` interval, the others those of each interval of the sequence. The GPU large-cluster path
(`bin/main_gpu.out`: every other cluster with `backend` = `cuda`, the clusters of `gpu_options.min_dim` and above with
`auto`) computes every phase with sin and cos, as `exact`. <br>
Without the key (per path in the GPU build) each path takes `recursive` where it has
it and `exact` elsewhere, without a message: with `backend` = `cpu`, and in the CPU build `bin/main.out` with `cpu` or
`auto` (all on the CPU there), the value is `recursive` (gCCE `expm`: `exact`); in `bin/main_gpu.out` with `auto` the
clusters on the CPU (below `min_dim`) and the packed kernel run `recursive` and the clusters on the GPU `exact`; with
`cuda` the packed kernel runs `recursive` and every other cluster `exact` on the GPU (with the packed kernel off, as for
gCCE and pCCE, the value is `exact`). The run log prints the value and why, per path in the GPU build
(`chosen by : default (supported path)`, `default (clusters on the CPU and packed kernel: recursive, clusters on the GPU
(execution.backend=auto): exact)`, `default (execution.backend=cuda: every cluster on the GPU, exact there)`, … or `given`). <br>
A given `recursive` is refused rather than run as `exact` where a cluster may reach a path without it: gCCE with
`propagator` = `expm` (no phases to carry over), and in `bin/main_gpu.out` `backend` = `cuda` or `auto` with a method
the GPU takes (cCCE with `cluster.mode` standard, pcce or bundle, gCCE with `eigen`) — leave the key out there to get `recursive` on the CPU clusters and the packed
kernel. `backend` = `cuda` in the CPU build is refused in any case. <br>
Measured on the production P1 inputs (CR0 and CR60 σ4, cCCE2 and gCCE2, 30 states; 64 ranks): `recursive` takes
about 0.73–0.76 of the time of `exact` for cCCE2 and 0.82–0.84 for gCCE2. The `wiDiv`
curves moved by at most 2.2e-10 for cCCE2 (the outputs have 10 decimals) and 1.3e-8 for gCCE2 at points with |L| ≤ 1
((internal notes)). The results are not bit-identical to `exact`, so a converged
series should not switch between the two. Must be a JSON string. <br>

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
The one-state cCCE time evolution of the large clusters with the time points as matrix columns. The usual path carries one row and one column vector per
time point through the intervals of the sequence: two matrix-vector products per interval and time point, each reading
the cluster's D × D matrix again. With this block on, `columns` time points go together: each interval is one
matrix-matrix product (the library's zgemm in a build with MKL), the matrix is read once per block. The multiplications
are the same; the summation order is not, so the outputs move in the last bits (1e-15 relative in the unit check). The
phases follow the usual path's rule, `exact` and `recursive` alike (`recursive`: sin and cos at steps 0, 16, 32, ..., the
products carried across the blocks; the carried products may differ in the last bit where the compiler fuses the complex
multiplication otherwise), and every
pulse sequence of the one-state cCCE (Ramsey, Hahn echo, CPMG, `Manual`) is covered. <br>
Read on the CPU only, for one sampled state (`s_range`) of `method` cCCE (every `cluster.mode`), clusters of bath
dimension `min_dim` and above. Bath dimension 2 and 4 (the fixed-size path and the CPU batch), the ensemble
(no `s_range`), `propagator` `sparse`, gCCE and the clusters on the GPU are unchanged whatever the block says; with the block
`enabled` false every output is the bytes of the run without it. Extra memory per rank: the phases of one block
for each distinct interval and four D × `columns` matrices, (2 × intervals + 4) × D × `columns` × 16 bytes (D 2048,
CPMG-8, 64 columns: about 12 MB). <br>
**On by default, `min_dim` 256 and `columns` 64:** measured on a Xeon Gold 6530 node (
64 ranks, median of three passes, wall time with / without; pCCE sK5 with clusters of dimension 1024, sK3o3 with 512,
sK4 with 256; Hahn and CPMG-8, nstep 200 and 1000; 12 inputs): with MKL the geometric mean is 0.40 of the time (sK5
CPMG-8 nstep 1000: 114.3 s → 13.4 s, 8.5 times faster; sK4: 0.65–0.97), with Eigen alone (`CCEX_BLAS` off) 0.61. The
curves agree with the run without it to 1e-10 (the last printed digit; most are identical). Peak memory per rank moved by
−8 to +11 MB (MKL: 344–346 → 352–356 MB for sK5). `min_dim` 384 left sK4 unchanged (1.01) and 768 also sK3o3 (1.00);
`columns` 32 and 128 were within a few per cent of 64, 128 with up to +27 MB. Results and scripts:
(internal notes).
Without the block it is on with these values; a block without `enabled` is on with the values it gives (the others as above); `{"enabled": false}`
turns it off (the usual one-vector path). `enabled` must be a JSON boolean, `min_dim` and `columns` positive
integers. The run log prints `column_evolution.min_dim` and `.columns` when it is on. <br>

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
In the `cluster.adaptive` block (the refinement driver copies it into the step
inputs beside their `work_file`). The probe step of the adaptive clustering (`cluster.mode` = `bundle` with the
`cluster.adaptive` block): the memory, in GB per node, that the probes
running at once may hold. The probes of a bath state go to the ranks in turn and the ranks do not wait
for each other between states, so the large probes of several states can run together; a probe waits while the probes
already running on its node would pass this budget, and a probe larger than the budget runs alone. A probe is estimated
on the path it takes (`estimate` in `src/clustering/bundle/bundle_work.cpp`), in dense dim × dim complex matrices of 16 dim²
bytes each: on the CPU, gCCE 6 + n with `evolution` = `vector` and 7 + 2n with `matrix`; cCCE 8 for a sampled state, 23 for
the ensemble and 10 + 4(n + 1) for an ensemble `Manual` sequence (n = `pulse.count`); on the GPU what the device allocates for
the time batch and evolution; with `output.savemode` = `cluster_density_matrix` also every time's full
density matrix. Probes under 64 MiB do not count. 0 takes 80 % of the memory the job may use on the node (its RAM, or
less if a cgroup limits the job) and, when the run may put clusters on the GPU, of the smallest GPU. The run prints the budget (`[probe-budget]`) and how long
each rank waited. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> `0` (automatic)

<b>Example:</b> ``` "cluster" : {"adaptive": {..., "probe_memory_gb": 60}} ```

<a id="tag_gcce_refinement"></a>
## execution.gcce_refinement
<b>Description:</b><br>
gCCE with `propagator` = `eigen`: refine the eigenvalues and eigenvectors the solver returns by
compensated correction sweeps before they are used. A block `{"enabled": true|false, "max_sweeps": n}` (a plain boolean is refused).
Writing the block
turns the refinement on (a block without `enabled` is on, as an `auto_convergence` block is); `{"enabled": false}` keeps a
written block off. `enabled` must be a JSON boolean. <br>
Off by default: on the NV–¹³C gCCE4 case the refinement was made for, and on P1 gCCE3, it
changed the error-corrected final coherence by at most 6.7e-7
((internal notes); the raw `wiDiv` differs by more at points where it diverges with either
setting — up to 12.7 on NV–¹³C and 2.75e11 on P1 gCCE3 — so the predeclared 1e-6 bound on every raw point was not met), while it took most of the setup time of a gCCE
cluster (dim 768: about 4 s of 4.8 s on the CPU). Turn it on where the bath energies are very small beside the
qubit's and the curve is read to better than 1e-6. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "max_sweeps": Integer}`

<b>Default:</b> off (absent)

<b>Example:</b> ``` "execution" : {"gcce_refinement": {}} ``` (on: the block is written)

<a id="tag_gcce_refinement_max_sweeps"></a>
## execution.gcce_refinement.max_sweeps
<b>Description:</b><br>
The largest number of refinement sweeps per cluster. A cluster that
has not converged within it stops the run with an error naming its dimension. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `32`

<b>Example:</b> ``` "execution" : {"gcce_refinement": {"max_sweeps": 32}} ```

<a id="tag_profile"></a>
## execution.profile
<b>Description:</b><br>
Print timings at the end of the run. Every rank writes one `[profile]` line per method,
backend and cluster dimension with the number of clusters and the seconds spent building,
setting up and evolving them; with CUDA also `[gpu-batch]` and `[gpu-memory]` lines. In a plain cCCE/gCCE run
rank 0 also prints one `[phase]` line per rank: the wall seconds of the state loop split into state preparation,
cluster creation, cluster calculation, product assembly, freeing, waiting at the barriers, the reduction and saving.
Must be a JSON boolean. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "execution" : {"profile": true} ```
