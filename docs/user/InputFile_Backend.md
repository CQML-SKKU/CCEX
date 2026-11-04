# Execution tags

Where and how the cluster dynamics is computed: the CPU or the GPU, how a cluster is evolved in time, and what the run
reports about it. These keys are the `execution` block; at the top level they are refused with their place. Each key
has a small table (its type, default and unit), what it does, and an example. Two keys of other blocks are here too:
`cluster.adaptive.probe_memory_gb`, and the memory check, which has no key.

```json
"execution": {"backend": "auto", "gpu_options": {"min_dim": 768}, "propagator": "eigen", "profile": false}
```

1. [Where it runs](#ex_where): `backend`, `gpu_options`.
2. [How a cluster is evolved](#ex_evolve): `propagator`, `evolution`, `evolution_update`, `column_evolution`,
   `gcce_refinement`.
3. [Memory and timing](#ex_report): the memory check, `cluster.adaptive.probe_memory_gb`, `profile`.

The GPU backend itself is described in [GPU execution](../methods/gpu.md).

<a id="ex_where"></a>
## 1. Where it runs

<a id="tag_backend"></a>
### execution.backend

| Type | Default |
|---|---|
| string | `"cpu"` |

Where the dynamics of each cluster is computed. `cuda` and `auto` need the GPU binary `bin/main_gpu.out`
([GPU execution](../methods/gpu.md)).

| Value | Where |
|---|---|
| `cpu` | every cluster on the CPU |
| `cuda` | every cluster the GPU takes on the GPU: cCCE (`cluster.mode` standard, pcce or bundle) and gCCE with `propagator` eigen (not with `cluster.mode` pcce) |
| `auto` | the clusters of Hilbert dimension [gpu_options](#tag_gpu_options) `min_dim` and above on the GPU, the rest on the CPU |

With `auto`, a run stays on the CPU when the binary is CPU-only, no card can be opened or the method is one the GPU does
not take; a cluster that runs out of GPU memory goes to the CPU, and so do the later ones of its dimension or larger on
that rank (a `[gpu-memory-fallback]` line). The `[backend]` lines of the log say what was chosen.

With `cuda`, an unsupported method, a CPU-only binary or a card that cannot be opened stops the run: every rank prints
`Error: GPU backend: ...` and aborts. An ensemble (no `s_range`) `Manual` sequence with cCCE is refused with `cuda` and
runs on the CPU with `auto`.

In `bin/main_gpu.out` the small cCCE and gCCE clusters also go to two GPU kernels of their own, the packed and the warp
kernel ([gpu_options](#tag_gpu_options)), whatever the dimension.

```json
"execution": {"backend": "auto"}
```

<a id="tag_gpu_options"></a>
### execution.gpu_options

| Type | Default |
|---|---|
| object | every key at its default (below) |

The GPU settings. Every key is optional; an unknown key is refused. They are read by `bin/main_gpu.out` with `backend`
`cuda` or `auto`.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `min_dim` | positive integer | `768` (`256` with `ranks_per_gpu`) | with `auto`, the smallest cluster Hilbert dimension that goes to the GPU |
| `packed_max_dim` | `0`, `2` or `4` | `4` for cCCE with `cluster.mode` standard or bundle, else `0` | the packed kernel: the cCCE clusters of bath dimension up to this value, thousands at once; `0` is off |
| `warp_max_dim` | `0` to `32` | `27` for cCCE and gCCE, else `0` | the warp kernel: the cCCE clusters of bath dimension 8 and 16 and the gCCE clusters of total dimension 6 to 27 (one NV qubit), one per GPU warp; `0` is off |
| `evolution_batch_size` | positive integer | `8` (`32` for cCCE) | the largest number of clusters of one dimension taken together on the GPU |
| `evolution_batch_steps` | positive integer | `64` | the largest number of time points one GPU evolution takes together |
| `cluster_batch_mb` | positive integer | `64` (cCCE: two thirds of the card over its ranks) | MiB per rank that one GPU batch of clusters may hold |
| `cache_mb` | nonnegative integer | `64` | MiB per rank of freed GPU memory kept for reuse; `0` turns the cache off |
| `ranks_per_gpu` | nonnegative integer | `0` (`4` for cCCE, `cluster.mode` standard, `auto`, more than 4 ranks a GPU) | with `auto` and cCCE standard, only this many ranks per GPU use it, and the clusters of `min_dim` and above go to them; `0`: every rank may use it |

For cCCE with a card visible, the values the input leaves out of `min_dim`, `evolution_batch_size`, `cluster_batch_mb`
and `ranks_per_gpu` are set at the start of the run, as in the parentheses above; a `[gpu-auto]` line gives them. A
value the input gives is used as given. Ranks on a card × `cluster_batch_mb` should stay below the card's memory; a
batch that runs out of GPU memory leaves its clusters to the CPU.

The packed and warp kernels:

- They start only when the run has enough of their clusters to pay for opening the card (a `[gpu-packed]` or
  `[gpu-warp]` line says why not). A given `packed_max_dim` or `warp_max_dim` starts them whatever the count.
- The packed kernel takes pulse sequences of at most 64 pulses and 4 distinct intervals; the warp kernel one qubit,
  `propagator` eigen and at most 4 distinct intervals. In other runs they switch themselves off with a log line.
- In a cCCE run that the CPU's own batch of small clusters takes, the packed kernel is used only where it is faster:
  for one state with many small clusters, never for the ensemble.
- A given nonzero value is refused in the CPU binary, with `backend` `cpu` for `warp_max_dim`, and with a method they
  do not take.

With `ranks_per_gpu`, the run times each cluster dimension of `min_dim` and above on the CPU and on the GPU at the
start (`[gpu-split]` lines), and gives each large cluster to a GPU rank or a CPU rank, whichever finishes it first. The
split differs between runs, so the outputs agree within rounding, not bit for bit. It needs `auto` and cCCE with
`cluster.mode` standard.

The GPU results agree with the CPU within 10⁻⁸ + 2×10⁻⁶·|L|, not bit for bit.

```json
"execution": {"backend": "auto", "gpu_options": {"min_dim": 256, "ranks_per_gpu": 4}}
"execution": {"backend": "cuda", "gpu_options": {"warp_max_dim": 0}}
```

<a id="ex_evolve"></a>
## 2. How a cluster is evolved

<a id="tag_propagator"></a>
### execution.propagator

| Type | Default |
|---|---|
| string, or `{"method": string, "min_dim": integer, "tol": number}` | `"eigen"` |

How the free evolution of a cluster is computed.

| Value | Method |
|---|---|
| `eigen` | one diagonalization of the cluster Hamiltonian, then $e^{-iH\tau} = M\,\mathrm{diag}(e^{-i\lambda\tau})\,M^\dagger$ for every step; cCCE and gCCE |
| `expm` | a matrix exponential per time step; gCCE only |
| `sparse` | the action of $e^{-iH\tau}$ on the state by a Krylov expansion, without forming the dense matrix; for clusters of dimension `min_dim` and above, the smaller ones on `eigen` |

`sparse` is for clusters whose dense matrices do not fit in memory: at dimension 16384 `eigen` needs about 36 GB,
`sparse` about 0.2 GB. It is much slower at high field and for long times, since its work grows with the spectral width
times the interval. It runs for one sampled state (`s_range`) of cCCE, or of gCCE with `evolution` vector and one
qubit, on the CPU with `evolution_update` exact. Refused with it: the ensemble, gCCE `evolution` matrix, `backend` cuda
or auto, `evolution_update` recursive, `hamiltonian.cCCE.knight` and `hfmedi`, and more than one qubit.

| Object key | Meaning | Default |
|---|---|---|
| `method` | `eigen`, `expm` or `sparse` (required) | |
| `min_dim` | `sparse`: the smallest cluster dimension that takes it | `4096` |
| `tol` | `sparse`: the estimated error of one interval, relative to the norm of the state (0 < tol < 1) | `1e-12` |

```json
"execution": {"propagator": "eigen"}
"execution": {"propagator": {"method": "sparse", "min_dim": 4096}}
```

<a id="tag_evolution"></a>
### execution.evolution

| Type | Default |
|---|---|
| string | `vector` where it can run, else `matrix` |

gCCE only: what is propagated. `vector` needs a pure initial state (`s_range`) and `propagator` eigen or sparse. An
explicit `vector` that cannot run stops the run; with another method an explicit `vector` is refused. The two agree to
rounding, not bit for bit.

| Value | Propagates |
|---|---|
| `vector` | the state vector |
| `matrix` | the density matrix $\rho(t) = U\rho_0 U^\dagger$; works in every case |

```json
"execution": {"evolution": "matrix"}
```

<a id="tag_evolution_update"></a>
### execution.evolution_update

| Type | Default |
|---|---|
| string | `recursive` where it exists, `exact` elsewhere |

How the phases $e^{-iE\tau}$ are formed at each time point.

| Value | Phases |
|---|---|
| `exact` | sin and cos at every time point |
| `recursive` | sin and cos at every 16th point, and in between the previous point's phases times one step's factor |

`recursive` exists on the CPU for cCCE and for gCCE with `propagator` eigen, and in the GPU's packed kernel. The GPU's
other paths, gCCE `expm`, meCCE and megCCE always use `exact`. Without the key each path takes `recursive` where it has
it; the log says which. A given `recursive` is refused where a cluster may reach a path without it: gCCE `expm`, and
`bin/main_gpu.out` with `backend` cuda or auto. `recursive` is faster (about 0.75 of the time for cCCE2) and moves the
curves by about 10⁻¹⁰; a converged series should not switch between the two.

```json
"execution": {"evolution_update": "exact"}
```

<a id="tag_column_evolution"></a>
### execution.column_evolution

| Type | Default |
|---|---|
| `{"enabled": boolean, "min_dim": integer, "columns": integer}` | on, `min_dim` 256, `columns` 64 |

The one-state cCCE time evolution of large clusters with `columns` time points at once, as matrix-matrix products. It
applies on the CPU to one sampled state (`s_range`) of cCCE, for clusters of bath dimension `min_dim` and above, and
can be several times faster there. The outputs move only in the last digits. It needs (2 × intervals + 4) ×
D × `columns` × 16 bytes more per rank (about 12 MB for D = 2048, CPMG-8). A block without `enabled` is on;
`{"enabled": false}` turns it off.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `enabled` | boolean | `true` | on or off |
| `min_dim` | positive integer | `256` | the smallest bath dimension of a cluster that takes it |
| `columns` | positive integer | `64` | the time points per block |

```json
"execution": {"column_evolution": {"enabled": false}}
```

<a id="tag_gcce_refinement"></a>
### execution.gcce_refinement

| Type | Default |
|---|---|
| `{"enabled": boolean, "max_sweeps": integer}` | off |

gCCE with `propagator` eigen: refines the eigenvalues and eigenvectors of each cluster by correction sweeps before they
are used. Writing the block turns it on; `{"enabled": false}` keeps a written block off, and a plain boolean is refused.
It takes most of the setup time of a gCCE cluster and changed the error-corrected coherence by less than 10⁻⁶ in the
cases tried. Turn it on when the bath energies are very small beside the qubit's and the curve is needed to better
than 10⁻⁶.

```json
"execution": {"gcce_refinement": {}}
```

<a id="tag_gcce_refinement_max_sweeps"></a>
### execution.gcce_refinement.max_sweeps

| Type | Default |
|---|---|
| positive integer | `32` |

The largest number of refinement sweeps per cluster. A cluster that has not converged by then stops the run with an
error naming its dimension.

```json
"execution": {"gcce_refinement": {"max_sweeps": 64}}
```

<a id="ex_report"></a>
## 3. Memory and timing

<a id="memory_check"></a>
### The memory check (no key)

After the clusters are made and before the first one is calculated, every rank estimates its peak memory: what it holds
then plus what its largest cluster needs on its path. The estimate is compared with the memory one rank may use: under
Slurm the job's memory (`SLURM_MEM_PER_NODE`, or `SLURM_MEM_PER_CPU` × `SLURM_CPUS_ON_NODE`), otherwise `MemAvailable`,
divided by the ranks on the node. The log prints one line:

```
      Memory estimate    : 0.316 GB per rank (rank 0: 0.0544 GB held + its largest cluster, dimension 1024), 20.4 GB available per rank (MemAvailable / 1 rank on the node)
```

A run whose estimate is larger stops there on every rank, before anything is calculated, with
`Error: not enough memory: ...`. The message says what to change: fewer ranks per node, a lower order or smaller
clusters (`rdip`, `cluster.pcce` or `cluster.bundle` `max_dim`), and where they apply `quantity` coherence instead of
dm or `output.savemode` coherence instead of `cluster_density_matrix`. A dense cluster matrix is 16 D² bytes, so the
memory grows fast with the order.

Not checked: `output.cluster_only` (no dynamics) and the probe steps of the adaptive refinement, which have their own
budget ([probe_memory_gb](#tag_probe_memory_gb)). The GPU checks its own memory.

<a id="tag_probe_memory_gb"></a>
### cluster.adaptive.probe_memory_gb

| Type | Default | Unit |
|---|---|---|
| nonnegative number | `0` (automatic) | GB per node |

A key of the `cluster.adaptive` block: the memory that the probes of the adaptive refinement running at once may hold
on a node. A probe waits while the probes already running on its node would exceed it; a probe larger than the budget
runs alone. `0` takes 80 % of the memory the job may use on the node, and, when clusters may go to the GPU, 80 % of the
smallest GPU. The log prints the budget and how long each rank waited (`[probe-budget]`).

```json
"cluster": {"adaptive": {"workdir": "adaptive", "probe_memory_gb": 60}}
```

<a id="tag_profile"></a>
### execution.profile

| Type | Default |
|---|---|
| boolean | `false` |

Prints timings at the end of the run: one `[profile]` line per rank, method, backend and cluster dimension, with the
number of clusters and the seconds spent building, setting up and evolving them; with the GPU also `[gpu-batch]` and
`[gpu-memory]` lines; and in a plain cCCE or gCCE run one `[phase]` line per rank, the state loop split into its phases.

```json
"execution": {"profile": true}
```
