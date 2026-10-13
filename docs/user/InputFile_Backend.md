# Backend tags
Where and how the cluster dynamics is computed, and what the run reports
about it: the CPU or CUDA backend, the gCCE propagation, the timing and the
per-cluster diagnostics.

- [backend](#tag_backend)
- [gpu_options](#tag_gpu_options) (min_dim, packed_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb, warp_max_dim)
- [probe_memory_gb](#tag_probe_memory_gb)
- [propagator](#tag_propagator)
- [evolution](#tag_evolution)
- [evolution_update](#tag_evolution_update)
- [gcce_refinement](#tag_gcce_refinement) (`enabled`, [max_sweeps](#tag_gcce_refinement_max_sweeps))
- [profile](#tag_profile)
- [cluster_contributions](#tag_cluster_contributions) (`enabled`, [max_mb](#tag_cluster_contribution_max_mb))

<a id="tag_backend"></a>
## backend
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
| `cuda` | CUDA for every CCE, pCCE and gCCE cluster the packed kernel does not take; the run stops if the device or the build is missing |
| `auto` | of the clusters the packed kernel does not take, CUDA for those whose Hilbert dimension is at least [gpu_options](#tag_gpu_options).`min_dim`, the CPU for the rest; the CPU also when the GPU cannot be opened, and when it runs out of memory (that cluster and, on that rank, every later one of its dimension or larger; a `[gpu-memory-fallback]` line in the log) |

The packed kernel ([gpu_options](#tag_gpu_options).`packed_max_dim`) takes its small cCCE clusters first, whatever the
backend; the table is the choice for the clusters left. An ensemble (`nstate` 0) `Manual` sequence is refused with
`cuda` and runs as `cpu` with `auto`; this is decided before the packed plan.

<b>Default:</b> `cpu`

<b>Example:</b> ``` "backend" : "auto" ```

<b>Related variable:</b> [Config](#Config)::backend

<a id="tag_gpu_options"></a>
## gpu_options
<b>Description:</b><br>
The GPU settings, in one block (since 2026-10-12, A26; the flat keys before it — `gpu_min_dim`, `gpu_batch_steps`,
`gpu_batch_clusters`, `gpu_cluster_batch_mb`, `gpu_cache_mb` — are refused with a message naming the new form, e.g.
`gpu_min_dim is now gpu_options.min_dim`). Every key is optional; an unknown key is refused. Without the block every
default below holds; the outputs are those of the earlier flat defaults except that the GPU build with `backend` `cuda`
or `auto` now runs the small cCCE clusters in the packed kernel (`packed_max_dim`, below) and, since 2026-10-13, the gCCE
and cCCE clusters of dimension up to 27 in the warp kernel (`warp_max_dim`, below). <br>

<b>Data Type:</b> Object

<b>Keys:</b><br>
| Key | Was | Type, default | Meaning |
|---|---|---|---|
| `min_dim` | `gpu_min_dim` | positive integer, `768` (`256` until 2026-10-13) | with `backend` = `auto`, the smallest cluster Hilbert dimension that goes to the GPU; smaller clusters stay on the CPU. Ignored by `cpu` and `cuda`. 768: on idm4-gpu (A37, job 699; RTX PRO 6000, 32 ranks) the GPU was faster from about 670-750 on in every path (cCCE and gCCE, single state and ensemble) and slower at 648 and below |
| `packed_max_dim` | (new) | `0`, `2` or `4`; `4` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE (`cluster_mode` standard or bundle), else `0` | the packed kernel (below): cCCE clusters of bath dimension up to this value on the GPU, many at once. `0` is off |
| `evolution_batch_size` | `gpu_batch_clusters` | positive integer, `8` | the largest number of consecutive clusters of one dimension evolved together on the GPU (`backend` `cuda`/`auto`); `1` evolves one cluster at a time |
| `evolution_batch_steps` | `gpu_batch_steps` | positive integer, `64` | the largest number of time points a GPU evolution processes together |
| `cluster_batch_mb` | `gpu_cluster_batch_mb` | positive integer, `64` (MiB per rank) | the evolution-buffer budget that limits how many clusters one GPU batch may hold |
| `cache_mb` | `gpu_cache_mb` | nonnegative integer, `64` (MiB per rank) | GPU memory that freed allocations are kept in for reuse; `0` disables the cache |
| `warp_max_dim` | (new, A40) | `0` to `32`; `27` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` cCCE or gCCE, else `0` | the warp kernel (below): gCCE and cCCE clusters of a compiled total dimension up to this value on the GPU, one per warp. `0` is off |

<b>The packed kernel (`packed_max_dim`).</b> The `cuda`/`auto` backends send a cluster to the GPU on its own, which
for the dimension-2 and 4 clusters of a P1 bath costs far more than the CPU (A12, A26 step 0: about 400 times slower).
With `packed_max_dim` = 4 (or 2), every cCCE cluster of bath dimension 2 and 4 (or 2) is instead computed by one GPU
thread — both eigensystems by Jacobi rotations and the whole time chain — thousands of clusters per kernel. Each rank
still builds its clusters' Hamiltonians and multiplies their coherences into the product as before; one process per node
(node rank 0, in a thread of its own) hands the packed clusters of all the node's ranks to the GPU through node shared
memory, so the packed kernel opens one CUDA context per node (node rank 0 opens it as soon as `packed_max_dim` > 0,
also with `backend` = `cpu`; ranks that send larger clusters to the GPU under `cuda`/`auto` open their own). Larger
clusters follow `backend` (and `evolution_batch_size`) as before. With `backend` = `auto` a GPU that cannot be opened
turns the packed kernel off and every cluster runs on the CPU; with `cuda` (or an explicit `packed_max_dim` with `cpu`)
the run stops. The packed plan is chosen after the general input checks (for instance the refusal
of an ensemble `Manual` sequence with `cuda` above). Single state and ensemble, every pulse sequence with at most 64
pulses and at most 4 distinct interval lengths (Hahn, CPMG, Ramsey, most `Manual` ones; others stay on the CPU), `evolution_update` `exact` and `recursive` (the kernel follows the run's value). Needs `bin/main_gpu.out` and
`method` = `cCCE` with `cluster_mode` standard or bundle (refused otherwise). The results agree with the CPU within the GPU tolerance 1e-8 + 2e-6 |CPU|, not bit
for bit (the eigensolver differs); the standing check is `test_small_cluster_cuda` in `make gpu-test`.
<b>Default (user, 2026-10-12):</b> on (`4`) in `bin/main_gpu.out` with `backend` = `cuda` or `auto` and `method` cCCE
(quantity coherence); off (`0`) with `backend` = `cpu`, with any other method, and always in the CPU build `bin/main.out`,
where nothing goes to the GPU and only an explicit nonzero value is refused. An explicit value works with any `backend`
(`cpu` included: then the packed clusters are the only GPU work). The run log prints the value and why
(`chosen by : default (backend=auto, method=CCE)`, `default (CPU build)`, ... or `given`). <br>
<b>When it pays.</b> Starting it costs a fixed part of a second (the CUDA context and page-locking the node shared memory,
done while the bath is read and the clusters are made, and released when no rank has a cluster it takes; until it is up
the ranks compute on the CPU), so a run of a few
seconds gains nothing. Measured on idm4-gpu (32 ranks + one RTX PRO 6000, job 691, the production P1 cCCE2 inputs, whole
run against the CPU binary on the same node, median of AB BA AB): CR60 σ4 30 states 70.7 → 44.2 s (0.62), its ensemble
6.5 → 3.9 s (0.59), CR0 30 states 13.6 → 11.4 s (0.84), but the CR0 ensemble, a 1.8 s run, 1.80 → 1.71 s (0.95);
(internal notes). In that job the 12 small fixture runs (under 1 s) finished before the GPU was
up and ran entirely on the CPU. <br>
<b>When they stay off (A40 follow-up, 2026-10-13).</b> With neither `packed_max_dim` nor `warp_max_dim` given, the packed
and warp kernels are not started at all when the cluster evaluations they would take (clusters × states: each cluster
once per state, `nstate` states or 1 for the ensemble), counted over all ranks from order 1 up, are fewer than 20000 per
rank (640000 for 32 ranks): on a full node starting the card costs about 0.5 s, and the card saves about 4.4 rank-µs per
packed evaluation and 21 per warp one, so such runs gain little or lose (the reasoning and the numbers:
(internal notes)). The log says so (`[gpu-packed] off: N cluster evaluations
(clusters x states) for the packed kernel, fewer than the threshold ...`), and every cluster is computed on the CPU path,
byte-identical to the run with both set to 0. The card's early start (opened beside the bath reading, before the count is
known) is closed again at that point; it is not avoided. A given value of either key starts the card whatever the count. The warp kernel's part of the
node shared memory is allocated only when some rank has a cluster of a warp dimension (a cCCE run with clusters of bath
dimension 2 and 4 only keeps the packed kernel's memory: P1 CR0 cCCE2 was 4.7 % slower with it, job 708). <br>

<b>The warp kernel (`warp_max_dim`; A40, user 2026-10-13).</b> Clusters too large for the packed kernel's
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
value and why (`chosen by : default (backend=auto, method=gCCE)`, `default (CPU build)`, ... or `given`). Where the run is
one the kernel does not take (`gcce_refinement`, `pulse.per_qubit`, a pulse detuning, savemode `cluster_density_matrix`, another qubit
dimension, more than 4 distinct intervals) it switches itself off, with one `[gpu-warp] off: ...` line, and the run goes
on as before. A given nonzero value is refused in the CPU build, with `backend` `cpu`, with a method other than CCE or
gCCE, and when no compiled dimension is up to it. <br>
<b>Results.</b> Within the GPU tolerance 1e-8 + 2e-6 |CPU| of the CPU, not bit-identical to the GPU build before it (those
clusters were on the CPU path). The standing check is `test_warp_cluster_cuda` in `make gpu-test` (CPU-dumped clusters of
every compiled dimension). Inputs whose qubit splitting is many orders above the bath couplings amplify the eigensolver's
rounding through the cluster powers: an NV at zero-field splitting with a 13C bath, CCE3, where the CPU itself differs from
its refined eigensystem (`gcce_refinement`) by about 5e-7 and the warp kernel by as much; and in P1 CR60 σ4 gCCE2 (job 703)
14 of 2416 values differ by 2–5.5e-6 relative - all in the stretch where the product diverges (|L| > 1 and the points
between), which the acceptance leaves out (user 2026-10-13). Measured on idm4-gpu (32 ranks + one RTX PRO 6000, MPS, job
703, whole run against the GPU build without it, AB BA AB): P1 CR0 gCCE2 14.2 → 9.8 s (0.69), CR60 σ4 60.0 → 36.0 s (0.60),
h-BN gCCE2 10.9 → 9.3 s (0.86); (internal notes). <br>

<b>Example:</b> ``` "gpu_options" : {"min_dim": 256, "packed_max_dim": 4} ```, ``` "gpu_options" : {"warp_max_dim": 12} ```

<b>Related variables:</b> [Config](#Config)::gpu_min_dim, gpu_packed_max_dim, gpu_batch_clusters, gpu_batch_steps,
gpu_cluster_batch_mb, gpu_cache_mb, gpu_warp_max_dim

<a id="tag_propagator"></a>
## propagator
<b>Description:</b><br>
How the free evolution of a cluster is propagated. `eigen` and `sparse` exist for cCCE (`method` cCCE, any `cluster_mode`) and gCCE, `expm` for gCCE only (refused with any other method). <br>
`sparse` (A2, 2026-10-12) applies exp(-i H τ) to the state vector without diagonalizing H: the cluster Hamiltonian is
assembled in compressed rows from its one- and two-spin terms (no dense matrix of the cluster is formed) and the
exponential is a Lanczos (Krylov) expansion of dimension 30 (the whole space up to dimension 64) with an a-posteriori
error estimate per substep; the interval is split into as many substeps as the estimate requires. It runs for clusters of
Hilbert dimension `min_dim` and above (gCCE: qubit times bath dimension); smaller clusters keep the `eigen` path, so
`"sparse"` with the default `min_dim` gives the bytes of `eigen` with `evolution_update` = `exact` wherever no cluster reaches 4096. It covers one sampled
state (`nstate` > 0) of the cCCE and gCCE with `evolution` = `vector` and one qubit, with `backend` = `cpu` and
`evolution_update` = `exact`, without `hamiltonian.cCCE.knight` and `hfmedi`; the ensemble (`nstate` = 0), gCCE `evolution` = `matrix`,
`backend` = `cuda` or `auto`, `evolution_update` = `recursive`, `knight`, `hfmedi` and more than one qubit are refused at
the start. <br>
**When `sparse` pays:** only when the dense eigendecomposition of a cluster does not fit in memory. Memory is O(dim)
instead of the O(dim²) matrices of `eigen`: at dim 16384 (14 spin-½) `eigen` needs about 36 GB, `sparse` about 0.2 GB (idm job 17303/17304).
**When it does not:** at high field or with a wide spectrum and long times. The cost of an expansion grows with the
spectral width of the Hamiltonian times the interval, and the second and later intervals of a sequence start from a new
state at every time point, so the work also grows with the square of the number of time points. On the P1 probe at
500 G (Hahn echo, 151 points to 15 µs) dim 4096 takes about 18 h with `sparse` (11 points: 6 min, × (151/11)²) against 5 min with `eigen`, and dim 16384 about
5 days against 6 h (11 points: 40 min; idm, one core, jobs 17303/17304); at 0 G the same
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

<b>Example:</b> ``` "propagator" : "eigen" ```, ``` "propagator" : {"method": "sparse", "min_dim": 4096, "tol": 1e-12} ```

<b>Related variable:</b> [Config](#Config)::propagator

<a id="tag_evolution"></a>
## evolution
<b>Description:</b><br>
gCCE only. What is propagated. `vector` needs a pure initial state (`nstate` > 0),
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

<b>Example:</b> ``` "evolution" : "matrix" ```

<b>Related variable:</b> [Config](#Config)::evolution

<a id="tag_evolution_update"></a>
## evolution_update
<b>Description:</b><br>
How the phases exp(-i E τ) of the eigenvalues are formed at each time point of the uniform grid (`deltat`). `exact`
computes them with sin and cos at every point. `recursive` does so at points 0, 16, 32, … and in between multiplies the
previous point's phases by one step's factor exp(-i E f `deltat`) (f the interval's fraction of the sequence), so the
rounding of at most 15 products accumulates before the next exact point (A16 (b), 2026-10-11). <br>
`recursive` exists on the CPU time evolution of cCCE (`method` cCCE, any `cluster_mode`; one
sampled state, `nstate` > 0, and the ensemble, `nstate` = 0, since 2026-10-11 (A18)), in the packed kernel of
`bin/main_gpu.out` (`gpu_options.packed_max_dim`, cCCE clusters of bath dimension 2 and 4) and on the CPU time evolution
of gCCE with `propagator` = `eigen` (`evolution` = `vector`, and `matrix` since A18). Every pulse sequence these paths
run is covered: the ensemble cCCE carries over the phases of its τ = t/2 (Hahn), t (Ramsey), t/(2 `pulse.count`) and 2τ
(CPMG) or of each `Manual` interval, the others those of each interval of the sequence. The GPU large-cluster path
(`bin/main_gpu.out`: every other cluster with `backend` = `cuda`, the clusters of `gpu_options.min_dim` and above with
`auto`) computes every phase with sin and cos, as `exact`. <br>
Without the key (since 2026-10-12, A29; per path in the GPU build since A36) each path takes `recursive` where it has
it and `exact` elsewhere, without a message: with `backend` = `cpu`, and in the CPU build `bin/main.out` with `cpu` or
`auto` (all on the CPU there), the value is `recursive` (gCCE `expm`: `exact`); in `bin/main_gpu.out` with `auto` the
clusters on the CPU (below `min_dim`) and the packed kernel run `recursive` and the clusters on the GPU `exact`; with
`cuda` the packed kernel runs `recursive` and every other cluster `exact` on the GPU (with the packed kernel off, as for
gCCE and pCCE, the value is `exact`). The run log prints the value and why, per path in the GPU build
(`chosen by : default (supported path)`, `default (clusters on the CPU and packed kernel: recursive; clusters on the GPU
(backend=auto): exact)`, `default (backend=cuda: every cluster on the GPU, exact there)`, … or `given`). Until A36
(2026-10-12) `backend` `cuda` or `auto` made the whole run `exact`; `"evolution_update": "exact"` gives those runs, and
the runs before A29, bit for bit. <br>
A given `recursive` is refused rather than run as `exact` where a cluster may reach a path without it: gCCE with
`propagator` = `expm` (no phases to carry over), and in `bin/main_gpu.out` `backend` = `cuda` or `auto` with a method
the GPU takes (CCE, pCCE, gCCE `eigen`) — leave the key out there to get `recursive` on the CPU clusters and the packed
kernel. `backend` = `cuda` in the CPU build is refused in any case. <br>
Measured on the production P1 inputs (CR0 and CR60 σ4, cCCE2 and gCCE2, 30 states; idm job 17217, 64 ranks): the run
took 0.63–0.76 of the time of the build before A16 (geometric mean 0.70; this includes the byte-identical A16 (a)
changes, and against those alone the same job gives about 0.73–0.76 for cCCE2 and 0.82–0.84 for gCCE2). The `wiDiv`
curves moved by at most 2.2e-10 for cCCE2 (the outputs have 10 decimals) and 1.3e-8 for gCCE2 at points with |L| ≤ 1
((internal notes)). The results are not bit-identical to `exact`, so a converged
series should not switch between the two. Must be a JSON string. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `exact` | sin and cos at every time point |
| `recursive` | sin and cos every 16th time point, the previous point's phases times one step's factor in between |

<b>Default:</b> `recursive` where it exists (CPU cCCE, the packed kernel, CPU gCCE with `propagator` = `eigen`), `exact` elsewhere (the GPU large-cluster path, gCCE `expm`; meCCE, megCCE: `exact`, and a given `recursive` is refused); per path in `bin/main_gpu.out` since A36 (the whole run `exact` with `backend` `cuda`/`auto` until then; `exact` everywhere until A29, 2026-10-12)

<b>Example:</b> ``` "evolution_update" : "exact" ```

<b>Related variable:</b> [Config](#Config)::evolution_update

<a id="tag_probe_memory_gb"></a>
## probe_memory_gb
<b>Description:</b><br>
The probe step of the adaptive clustering (`cluster_mode` = `bundle` with the `adaptive` block): the memory, in GB per node, that the probes
running at once may hold. Since 2026-10-02 the probes of a bath state go to the ranks in turn and the ranks do not wait
for each other between states, so the large probes of several states can run together; a probe waits while the probes
already running on its node would pass this budget, and a probe larger than the budget runs alone. A probe is estimated
on the path it takes (`estimate` in `src/clustering/bundle/bundle_work.cpp`), in dense dim × dim complex matrices of 16 dim²
bytes each: on the CPU, gCCE 6 + n with `evolution` = `vector` and 7 + 2n with `matrix`; cCCE 8 for a sampled state, 23 for
the ensemble and 10 + 4(n + 1) for an ensemble `Manual` sequence (n = `pulse.count`); on the GPU what the device allocates for
the time batch and evolution; with `savemode` = `cluster_density_matrix` also every time's full
density matrix. Probes under 64 MiB do not count. 0 takes 80 % of the memory the job may use on the node (its RAM, or
less if a cgroup limits the job) and, when the run may put clusters on the GPU, of the smallest GPU. The run prints the budget (`[probe-budget]`) and how long
each rank waited. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> `0` (automatic)

<b>Example:</b> ``` "probe_memory_gb" : 60 ```

<b>Related variable:</b> [Config](#Config)::probe_memory_gb

<a id="tag_gcce_refinement"></a>
## gcce_refinement
<b>Description:</b><br>
gCCE with `propagator` = `eigen`: refine the eigenvalues and eigenvectors the solver returns by
compensated correction sweeps before they are used. A block `{"enabled": true|false, "max_sweeps": n}` (until 2026-10-13
the boolean `gcce_refinement` and the key `gcce_refinement_max_sweeps`; the boolean is refused now). A block without
`enabled` is on, as an `auto_convergence` block is; `enabled` must be a JSON boolean. <br>
Off by default since 2026-10-02: on the NV–¹³C gCCE4 case the refinement was made for, and on P1 gCCE3, it
changed the error-corrected final coherence by at most 6.7e-7 (idm job 16869,
(internal notes); the raw `wiDiv` differs by more at points where it diverges with either
setting — up to 12.7 on NV–¹³C and 2.75e11 on P1 gCCE3 — so the predeclared 1e-6 bound on every raw point was not met), while it took most of the setup time of a gCCE
cluster (dim 768: about 4 s of 4.8 s on the CPU). Turn it on where the bath energies are very small beside the
qubit's and the curve is read to better than 1e-6. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "max_sweeps": Integer}`

<b>Default:</b> off (absent; on until 2026-10-02)

<b>Example:</b> ``` "gcce_refinement" : {"enabled": true} ```

<b>Related variable:</b> [Config](#Config)::gcce_refinement

<a id="tag_gcce_refinement_max_sweeps"></a>
## gcce_refinement.max_sweeps
<b>Description:</b><br>
The largest number of refinement sweeps per cluster (`gcce_refinement_max_sweeps` until 2026-10-13). A cluster that
has not converged within it stops the run with an error naming its dimension. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `32`

<b>Example:</b> ``` "gcce_refinement" : {"enabled": true, "max_sweeps": 32} ```

<b>Related variable:</b> [Config](#Config)::gcce_refinement_max_sweeps

<a id="tag_profile"></a>
## profile
<b>Description:</b><br>
Print timings at the end of the run. Every rank writes one `[profile]` line per method,
backend and cluster dimension with the number of clusters and the seconds spent building,
setting up and evolving them; with CUDA also `[gpu-batch]` and `[gpu-memory]` lines. In a plain cCCE/gCCE run
rank 0 also prints one `[phase]` line per rank: the wall seconds of the state loop split into state preparation,
cluster creation, cluster calculation, product assembly, freeing, waiting at the barriers, the reduction and saving.
Must be a JSON boolean. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "profile" : true ```

<b>Related variable:</b> [Config](#Config)::profile

<a id="tag_cluster_contributions"></a>
## cluster_contributions
<b>Description:</b><br>
Record the contribution of every cluster to the coherence over all time points
(cCCE and gCCE, `quantity` = `coherence`). For each bath state the run writes
`<outfile>.contributions.state-<s>.tsv` (one row per cluster: sites, bundles, production
power, dimension, the largest |log L̃|, undefined points, worst step) and once
`<outfile>.contributions.sites.tsv` (site, name, spin, gyromagnetic ratio, position).
A diagnostic; the raw table per rank is bounded by
[max_mb](#tag_cluster_contribution_max_mb). A block `{"enabled": true|false, "max_mb": n}` (until 2026-10-13 a boolean,
refused now, and the key `cluster_contribution_max_mb`); a block without `enabled` is on. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "max_mb": Integer}`

<b>Default:</b> off (absent)

<b>Example:</b> ``` "cluster_contributions" : {"enabled": true} ```

<b>Related variable:</b> [Config](#Config)::cluster_contributions

<a id="tag_cluster_contribution_max_mb"></a>
## cluster_contributions.max_mb
<b>Description:</b><br>
The memory per rank the contribution table of `cluster_contributions` may take (`cluster_contribution_max_mb` until
2026-10-13). The run
stops before the propagation if the table would exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Unit:</b> MiB per rank

<b>Default:</b> `256`

<b>Example:</b> ``` "cluster_contributions" : {"enabled": true, "max_mb": 512} ```

<b>Related variable:</b> [Config](#Config)::cluster_contribution_max_mb
