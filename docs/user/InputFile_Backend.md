# Execution tags

The `execution` block says where each cluster is computed, on the CPU or the GPU, and how its evolution in time is
computed. Every key has a default, and most runs need only `backend`. The memory check, which has no key, is described
here too. How the GPU backend works is in [GPU execution](../methods/gpu.md).

1. [Where it runs](#ex_where): `backend` and `gpu_options`.
2. [How a cluster is evolved](#ex_evolve): `propagator`, `evolution` and three switches.
3. [Memory and timing](#ex_report): the memory check, `cluster.adaptive.probe_memory_gb` and `profile`.

```json
"execution": {"backend": "auto"}
```

<a id="ex_where"></a>
## 1. Where it runs

<a id="tag_backend"></a>
### execution.backend

| Type | Default |
|---|---|
| string | `"cpu"` |

| Value | Where the clusters are computed |
|---|---|
| `cpu` | All on the CPU. |
| `auto` | The clusters of dimension `min_dim` ([gpu_options](#tag_gpu_options)) and above on the GPU, the rest on the CPU. A run the GPU cannot take stays on the CPU. |
| `cuda` | All on the GPU. A run the GPU cannot take stops with `Error: GPU backend: ...`. |

The GPU takes cCCE with `cluster.mode` standard, pcce or bundle, and gCCE with standard or bundle and `propagator`
eigen. `auto` and `cuda` need the GPU binary `bin/main_gpu.out`. With `bin/main.out`, `auto` runs on the CPU and
`cuda` stops. The `[backend]` lines of the log say what was chosen.

<a id="tag_gpu_options"></a>
### execution.gpu_options

| Type | Default |
|---|---|
| object | every key at its default |

Fine settings of the GPU, read by `bin/main_gpu.out`. The defaults suit most runs. For cCCE, `min_dim`,
`evolution_batch_size`, `cluster_batch_mb` and `ranks_per_gpu` are chosen at the start of the run when they are left
out (the values in brackets).

| Key | Default | Meaning |
|---|---|---|
| `min_dim` | `768` (`256` with `ranks_per_gpu`) | With `auto`, the smallest cluster dimension sent to the GPU. |
| `ranks_per_gpu` | `0` (`4` with `auto`, cCCE standard or pcce, and more than four ranks per GPU) | How many ranks per GPU send clusters to it. With `0` every rank does. Written in the input it needs `backend` auto and cCCE standard. |
| `packed_max_dim` | `4` for cCCE standard or bundle, else `0` | Small cCCE clusters up to this bath dimension go to a kernel that packs many of them together. `0` turns it off. |
| `warp_max_dim` | `27` for cCCE and gCCE (not pcce), else `0` | Small clusters up to this dimension go to a kernel that runs one cluster per warp. `0` turns it off. |
| `evolution_batch_size` | `8` (`32` for cCCE) | The most clusters of one dimension evolved together. |
| `evolution_batch_steps` | `64` | The most time points evolved together. |
| `cluster_batch_mb` | `64` (for cCCE, two thirds of the card shared by its ranks) | The GPU memory in MiB one batch of a rank may hold. |
| `cache_mb` | `64` | The freed GPU memory in MiB a rank keeps for reuse. `0` turns it off. |

```json
"execution": {"backend": "auto", "gpu_options": {"min_dim": 256, "ranks_per_gpu": 4}}
```

The GPU checks hold the GPU results to within 10⁻⁸ + 2×10⁻⁶·|L| of the CPU's.

<a id="ex_evolve"></a>
## 2. How a cluster is evolved

<a id="tag_propagator"></a>
### execution.propagator

| Type | Default |
|---|---|
| string, or `{"method", "min_dim", "tol"}` | `"eigen"` |

| Value | How the free evolution is computed |
|---|---|
| `eigen` | The cluster Hamiltonian is diagonalized once and reused for every step. For cCCE and gCCE. |
| `expm` | A matrix exponential per time step. For gCCE only. |
| `sparse` | A Krylov expansion acts on the state without the dense matrix. It is used for clusters of dimension `min_dim` (default 4096) and above, to the error `tol` (default 1e-12) per interval. |

Use `sparse` only when a cluster's dense matrices do not fit in memory. One dense matrix of dimension D takes 16 D²
bytes, 4.3 GB at D = 16384, and `eigen` holds several. `sparse` is much slower at high field and long times. It runs on
the CPU for one sampled state (`s_range`) and one qubit, with cCCE (`quantity` coherence) or with gCCE and `evolution`
vector. It needs `evolution_update` exact and no `knight`, `hfmedi` or `delta_A`.

```json
"execution": {"propagator": {"method": "sparse", "min_dim": 4096}}
```

<a id="tag_evolution"></a>
### execution.evolution

| Type | Default |
|---|---|
| `"vector"` or `"matrix"` | `vector` where it can run, else `matrix` |

What gCCE propagates. `vector` is the state vector and needs a pure initial state (`s_range`) and `propagator` eigen
or sparse. `matrix` is the density matrix and works in every case. The two agree to rounding.

```json
"execution": {"evolution": "matrix"}
```

<a id="tag_evolution_update"></a>
<a id="tag_column_evolution"></a>
<a id="tag_gcce_refinement"></a>
<a id="tag_gcce_refinement_max_sweeps"></a>
### Three switches

| Key | Type | Default | Meaning |
|---|---|---|---|
| `evolution_update` | `"exact"` or `"recursive"` | `recursive` on the CPU with cCCE or gCCE and `eigen`, else `exact` | `recursive` computes the phases exactly at every 16th point and multiplies in between. It is faster, and the curves move by about 10⁻¹⁰. Keep one choice through a convergence series. |
| `column_evolution` | `{"enabled", "min_dim", "columns"}` | on, `min_dim` 256, `columns` 64 | For one sampled state of cCCE on the CPU, clusters of bath dimension `min_dim` and above evolve `columns` time points at once. It is several times faster. |
| `gcce_refinement` | `{"enabled", "max_sweeps"}` | off (writing the block turns it on) | Refines each cluster's eigensystem, for gCCE with `eigen`. It helps when the bath energies are very small beside the qubit's. A cluster not converged after `max_sweeps` (default 32) sweeps stops the run. |

```json
"execution": {"evolution_update": "exact", "column_evolution": {"enabled": false}, "gcce_refinement": {"max_sweeps": 64}}
```

<a id="ex_report"></a>
## 3. Memory and timing

<a id="memory_check"></a>
### The memory check (no key)

Before the first cluster, every rank estimates its peak memory and compares it with what one rank may use. That is the
job's memory under Slurm, otherwise the free memory, divided by the ranks on the node. The log prints
`Memory estimate : ...`. A run that does not fit stops before any calculation with `Error: not enough memory: ...`,
and the message says what to change, such as fewer ranks per node, a lower order or smaller clusters. A dense cluster
matrix takes 16 D² bytes, so the memory grows fast with the order.

<a id="tag_probe_memory_gb"></a>
<a id="tag_profile"></a>
### Other keys

| Key | Type | Default | Meaning |
|---|---|---|---|
| `execution.profile` | boolean | `false` | Prints timings at the end of the run, per rank, method, backend and cluster dimension (`[profile]` lines). |
| `cluster.adaptive.probe_memory_gb` | number, 0 or more (GB per node) | `0`, which means 80 % of the job's memory on the node (and of the smallest GPU when clusters may go there) | The memory the probes of the adaptive refinement may hold at once. A probe waits until it fits. |

```json
"execution": {"profile": true}
```
