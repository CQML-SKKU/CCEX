# Cluster tags
How the bath is cut into clusters: the `cluster` block - the [mode](#tag_cluster_mode), the
`selection` every mode reads, and the blocks of the `pcce` and `bundle` modes and of the adaptive refinement. `order`
and `method` stay at the top. Older key names and places are in the table of [InputFile.md](InputFile.md#OldNames).
How the methods work: [CCE and dsj/itb](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundle clustering and adaptive refinement](../methods/bundle-adaptive.md).

```
"cluster": {"mode": "pcce", "selection": {"add_subclusters": true, "max_per_order": [0, 0, 0]},
            "pcce": {"partition_size": 2}}
```

- [order](InputFile_General.md#tag_order), [method](InputFile_General.md#tag_method) (top level)
- [cluster.mode](#tag_cluster_mode)
- cluster.selection: [add_subclusters](#tag_add_subclusters), [max_per_order](#tag_max_clusters_per_order)
- [cluster.pcce](#tag_pcce):[partition_size](#tag_pcce_partition_size), [trials](#tag_pcce_trials),
  [max_iterations](#tag_pcce_max_iterations), [kmeans_pp](#tag_pcce_kmeans_pp), [print_trials](#tag_pcce_print_trials),
  [partition](#tag_pcce_partition), [connections](#tag_pcce_connections),
  [cache](#tag_pcce_cache) (`mode`, `file`)
- [cluster.bundle](#tag_bundle): [max_dim](#tag_bundle_max_dim), [merge_cutoff](#tag_bundle_merge_cutoff),
  [connections](#tag_bundle_connections), [connection_budget](#tag_bundle_connection_budget),
  [max_clusters](#tag_bundle_max_clusters), [max_work](#tag_bundle_max_work)
- [cluster.adaptive](#tag_adaptive): [the keys](#tag_adaptive_keys), [the results](#tag_adaptive_results),
  [probe_memory_gb](#tag_adaptive_probe_memory_gb), and its internal [work_file](#tag_adaptive_work_file)

<a id="tag_cluster_mode"></a>
## cluster.mode
<b>Description:</b><br>
How the bath is divided into clusters. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description | method |
|---|---|---|
| `standard` | clusters of up to `order` spins within `rdip` | cCCE, gCCE, meCCE, megCCE |
| `pcce` | partition CCE: the bath is cut into partitions of [pcce.partition_size](#tag_pcce_partition_size) spins, a cluster is up to `order` partitions | cCCE, gCCE, meCCE, megCCE |
| `bundle` | bundle clustering: strongly coupled spins are tied into bundles first, a cluster is up to `order` bundles | cCCE, gCCE |
| `dsj` | disjoint clusters: the bath split into groups of at most `order` spins (strongest pairs within `rdip` first), each group computed whole, no other cluster | cCCE |
| `itb` | inter-bath clusters: only the pair corrections L̃ of the connected pairs between different groups (the factor `dsjitb` adds to `dsj`) | cCCE |
| `dsjitb` | disjoint + inter-bath (order 2 gives the same result as standard CCE2; above order 2 only pairs connect the groups) | cCCE |

`pcce` runs with cCCE, gCCE, meCCE and megCCE (the gCCE and master-equation forms
on the CPU only, and [pcce.cache](#tag_pcce_cache) stays with cCCE). gCCE with `dsj`, `itb` or `dsjitb` is refused (not
implemented yet); meCCE and megCCE with `bundle`, `dsj`, `itb` or `dsjitb`
likewise. Old method names are refused with their current form ([old names](InputFile.md#OldNames)).
The [pcce](#tag_pcce) block belongs to `pcce`, the [bundle](#tag_bundle) block to `bundle`, the [adaptive](#tag_adaptive)
block to `pcce` and `bundle`; each is refused with another mode. `"adaptive"` as a value is refused with the name `bundle`.

With `bundle`, [order](InputFile_General.md#tag_order) counts bundles per cluster and
[bundle.max_dim](#tag_bundle_max_dim) (required) caps the Hilbert dimension of one bundle. With
[rdip](InputFile_General.md#tag_rdip) given, bundles are connected by distance; with `rdip` omitted,
connections are chosen by a coherence check and the weak ones are dropped. Adding the [adaptive](#tag_adaptive) block turns
the adaptive refinement on top of it. The words: a *site* is one bath spin (or one defect that cannot be split); a
*bundle* is a set of sites that always enters a cluster whole and is treated exactly inside; a *cluster* is the union of
up to `order` connected bundles, the unit whose Hamiltonian is built and propagated; a *partition* is the whole assignment
of every site to exactly one bundle. Bundle clustering builds one partition and computes it once; the adaptive refinement
measures the error and rebuilds the partition. With every bundle one site, `bundle` is standard CCE.
[dsj, itb and dsjitb](../methods/cce.md#cce_dsj), [pCCE](../methods/pcce.md) and
[bundle](../methods/bundle-adaptive.md) are described in Methods. <br>

<b>Default:</b> `standard`

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8}} ```

<a id="tag_add_subclusters"></a>
## cluster.selection.add_subclusters
<b>Description:</b><br>
Order ≥ 3. With `true`, every sub-cluster of every kept cluster of order ≥ 3 (each subset of 2 or more of its spins,
or partitions with `pcce`) is put into the cluster list before the powers of the expansion are computed, among them
sub-clusters whose spins are not all within `rdip` of each other (the two ends of a chain) and sub-clusters that
[max_per_order](#tag_max_clusters_per_order) cut at a lower order. With `false` they are not listed up front; a
sub-cluster that the powers need and that is missing is added while the powers are computed, so both give the same
result; `false` can be a little faster (at order 2 it builds the same list by a shorter path). Not read with
`cluster.mode` bundle (its expansion is always closed under subsets). <br>
<b>Data Type:</b> Boolean (`true` or `false`)

<b>Default:</b> `true`

<b>Example:</b> ``` "cluster" : {"selection": {"add_subclusters": true}} ```

<a id="tag_max_clusters_per_order"></a>
## cluster.selection.max_per_order
<b>Description:</b><br>
Number of clusters for each CCE order. <br>
This is useful for CCE simulations with order ≥ 3, leading tons of clusters <br>
The list has exactly `order + 1` entries and the index is the cluster order: entry n limits the clusters of order n
(entry 0 is the order-0 slot). A 0 means no limit for that order, so all zeros include every cluster. A list of
another length is refused. Only entries 3 and up are used (with `pcce`, orders of partitions), with `cluster.mode`
standard and pcce: entries 0, 1 and 2 are never applied, and `dsj`, `itb`, `dsjitb` apply none. A nonzero entry that
is not applied is ignored with one warning line (`Warning: cluster.selection.max_per_order entry 1 is ignored ...`).
The clusters of order n are grown from the kept clusters of order n − 1, ranked by their strength, the sum of
\|J<sub>zz</sub>\| of the point-dipole coupling over every spin pair of the cluster (with `pcce`, of the partitions:
see [connections](#tag_pcce_connections)), and the strongest `max_per_order[n]` are kept. `cluster.mode` bundle
requires all 0. <br>
<b>Data Type:</b> Int List (length `order + 1`)

<b>Default:</b> all entries 0 (no limit)

<b>Example:</b> ``` "cluster" : {"selection": {"max_per_order": [0, 0, 0, 30000]}} ``` (order 3: at most 30000 clusters of three spins)

<a id="tag_pcce"></a>
## cluster.pcce
<b>Description:</b><br>
The settings of `cluster.mode` = `pcce`. Refused with any
other `cluster.mode`. `example/tutorial/07_pCCE_and_bundle/pCCE2_Hahn_states/` runs `partition_size` 1 to 3 (4 on request).
How the partitions are found and connected: [pCCE](../methods/pcce.md). <br>

<b>Data Type:</b> Object

<b>Default:</b> - (required with `cluster.mode` pcce, for its `partition_size`; every other key has its default)

<b>Example:</b> ``` "cluster" : {"mode": "pcce", "pcce": {"partition_size": 2, "trials": 1000, "max_iterations": 1000, "connections": "spins"}} ```

<a id="tag_pcce_partition_size"></a>
## cluster.pcce.partition_size
<b>Description:</b><br>
The number of bath spins in one partition. The bath is divided into partitions of
`partition_size` spins by a k-means search, and a pCCE cluster of order n is n partitions, n × `partition_size`
spins. Required with `cluster.mode` = `pcce`. <br>
**The rest spins are dropped:** with N bath spins (within `rbath`) there are N / `partition_size` partitions (rounded
down), and the N mod `partition_size` spins that are left over, the **last ones in the order of the bath file**, take no
part in the calculation: they are in no cluster. Which spins those are depends on how the file is written, not on their
distance from the qubit, and the run log does not list them; choose `rbath` or the bath file so that N is a multiple of
`partition_size` to keep every spin ([pCCE](../methods/pcce.md#pcce_rest)). <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 3}} ```

<a id="tag_pcce_trials"></a>
## cluster.pcce.trials
<b>Description:</b><br>
The number of partition trials, each a k-means run from its own starting centres;
the assignment with the best score is kept. With [pcce.partition](#tag_pcce_partition) =
`deterministic` the trials are split over the ranks in a fixed way, so the partition follows
`seed`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "trials": 1000}} ```

<a id="tag_pcce_max_iterations"></a>
## cluster.pcce.max_iterations
<b>Description:</b><br>
The largest number of k-means iterations in one partition trial. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "max_iterations": 1000}} ```

<a id="tag_pcce_kmeans_pp"></a>
## cluster.pcce.kmeans_pp
<b>Description:</b><br>
Start every k-means trial from k-means++ centres instead of random ones. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "kmeans_pp": true}} ```

<a id="tag_pcce_print_trials"></a>
## cluster.pcce.print_trials
<b>Description:</b><br>
The number of partition trials rank 0 prints the details of; `0` prints none. Before 2.0.0 it was a switch beside
`print_limit`, the most printed (two keys for one setting): `true` with `print_limit` L is L now, `false` is 0, and the
switch and `print_limit` are refused with this form. <br>

<b>Data Type:</b> Integer (nonnegative)

<b>Default:</b> `0`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "print_trials": 8}} ```

<a id="tag_pcce_partition"></a>
## cluster.pcce.partition
<b>Description:</b><br>
How the partition of the bath is searched. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `deterministic` | the partition follows `seed`; trials are distributed deterministically over the ranks |
| `legacy` | the original random-number and trial distribution |

<b>Default:</b> `deterministic`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "partition": "deterministic"}} ```

<a id="tag_pcce_connections"></a>
## cluster.pcce.connections
<b>Description:</b><br>
Which distance decides that two partitions are connected, so that they can be
in one cluster of order 2 or more. The inner and the outer `rdip` bound that
distance under either rule. (`cluster.mode` bundle connects its bundles by [bundle.connections](#tag_bundle_connections).) <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `centers` | the distance between the two k-means centres, the centroids of the partitions. The strength that orders the clusters is the point-dipole T<sub>zz</sub> between the centroids. |
| `spins` | any spin of the one partition between the inner and the outer `rdip` of any spin of the other, the rule of the bundles. The strength is the sum of \|T<sub>zz</sub>\| over the spin pairs of the two partitions. |

With `partition_size` = 1 the two rules are the same. Under `centers` two partitions whose spins
are close can stay unconnected when their centroids are far apart, and the more so
the larger `partition_size`; `spins` connects them. The rank-0 log says which rule ran and how
many partition pairs it connected. On the P1 bath of tutorial `03_P1` the two rules were compared: at that tutorial's `rdip` (620 Å) the two rules
differ less than pCCE2 differs from CCE4, and `spins` costs 27-65 % more; at an `rdip` comparable
to the partition size, `spins` is the closer of the two to CCE4.

<b>Default:</b> `centers` (the rule of every earlier version)

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "connections": "spins"}} ```

<a id="tag_pcce_cache"></a>
## cluster.pcce.cache
<b>Description:</b><br>
Reuse a stored partition instead of searching again. Needs
[pcce.partition](#tag_pcce_partition) = `deterministic` and method `cCCE` (refused with gCCE, meCCE and megCCE). `mode`: <br>

| mode | Description |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition; stop if there is none |

`file`: the partition cache file, relative to the working directory or absolute; its directory must exist. Required
with `mode` other than `off`. With [configuration.k_range](InputFile_General.md#tag_configuration_k_range) `{k}`
in it is the configuration's k, and a `file` without `{k}` is refused (one cache file per configuration, e.g.
`"partition_{k}.cache"`). <br>

<b>Data Type:</b> Object `{"mode": ..., "file": ...}`

<b>Default:</b> `{"mode": "off"}`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "cache": {"mode": "readwrite", "file": "./results/partition.cache"}}} ```

<a id="tag_bundle"></a>
## cluster.bundle
<b>Description:</b><br>
The settings of `cluster.mode` = `bundle`. Refused
with any other `cluster.mode`. `max_dim` is required. How the bundles are made and connected:
[bundle clustering](../methods/bundle-adaptive.md). <br>

<b>Data Type:</b> Object

<b>Default:</b> - (required with `cluster.mode` bundle, for its `max_dim`; every other key has its default)

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8, "merge_cutoff": 1.0}} ```

<a id="tag_bundle_max_dim"></a>
## cluster.bundle.max_dim
<b>Description:</b><br>
The largest Hilbert dimension one bundle may have. Required with `cluster.mode` =
`bundle`. The dimension of a cluster follows as the product of its bundle dimensions (gCCE: times the
qubit dimension); a cluster that cannot be handled is refused during the run. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8}} ```

<a id="tag_bundle_merge_cutoff"></a>
## cluster.bundle.merge_cutoff
<b>Description:</b><br>
The smallest coupling strength for which two sites are put into one bundle. The sites are taken in pairs, strongest
first; a pair whose strength is at least `merge_cutoff` (and above 0) joins the two bundles holding its sites when the
product of their dimensions is at most [max_dim](#tag_bundle_max_dim). The strength of a site pair is
Σ 3 I<sub>i</sub> I<sub>j</sub> ‖J<sub>ij</sub>‖<sub>F</sub> over the spins of the two sites (a defect with its on-site
spins is one site), ‖J‖<sub>F</sub> the Frobenius norm of their point-dipole tensor: a bound on the size of the
coupling, independent of the field direction ([the coupling strength](../methods/bundle-adaptive.md#bundle_coupling)).
With `0` every candidate pair may merge and `max_dim` alone limits the bundles; a large value keeps every site alone. <br>

<b>Data Type:</b> Double (nonnegative), or `{"value": ..., "unit": ...}` with the unit `rad/s`, `rad/ms` or `rad/us`
([Units](InputFile.md#Units))

<b>Unit:</b> rad/ms

<b>Default:</b> `0`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "merge_cutoff": 1.0}} ```

<a id="tag_bundle_connections"></a>
## cluster.bundle.connections
<b>Description:</b><br>
How bundles are connected into clusters. `checked`, the coherence-checked connections, is what an input without `rdip`
and without this key runs, and can be written (it was a mode only the absence of the two keys chose): every link is kept at first; then up to 8 of the weakest bundle pairs (while the
coupling score of `connection_budget` below stays under about 3e-3) are proposed for omission, the domains around them
and the clusters of their local expansion are computed directly for every bath state and time point, and the omissions
are made only if the expansion with them stays within 1e-7 + 1e-5 \|L\| of the direct curves (and the expansion
without them does too); otherwise every link is kept. The cost is one more pass over the bath states with at most 64
domain curves per state. The outcome is the log line `[adaptive-check] status=...` and the report
`<output file>.adaptive.json`; it is a local check, not a convergence test of the whole bath
([checked connections](../methods/bundle-adaptive.md#bundle_connections)). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `distance` | pairs within `rdip`; `rdip` is required |
| `auto` | the legacy budget mode: weak pairs are left out up to [connection_budget](#tag_bundle_connection_budget), which is required |
| `checked` | the coherence-checked connections above; `rdip` is refused with it |

The mode used, with the keys written (the run log names it, `cluster.bundle.connections : ...`):

| `rdip` | `connections` | mode |
|---|---|---|
| given | absent | `distance` |
| absent | absent | `checked` |
| given or absent | `auto` | `auto` (without `rdip`: no distance limit) |
| given | `distance` | `distance` |
| absent | `checked` | `checked` |

<b>Default:</b> `distance` with `rdip`, `checked` without

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}} ```

<a id="tag_bundle_connection_budget"></a>
## cluster.bundle.connection_budget
<b>Description:</b><br>
With `connections` = `auto`, the allowed score of the left-out bundle pairs:
t<sub>max</sub> × Σ w<sub>AB</sub> over the pairs left out must stay at most this value, with t<sub>max</sub> =
`deltat` × (`nstep` − 1) in ms and w<sub>AB</sub> the [coupling strength](#tag_bundle_merge_cutoff) of the two bundles
in rad/ms, so the score is in rad. Pairs outside `rdip` (when given) are left out first; the weakest of the rest are left
out while the score allows. A score that the required omissions already exceed is refused. It is a coupling score, not
an error bound. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> -

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}} ```

<a id="tag_bundle_max_clusters"></a>
## cluster.bundle.max_clusters
<b>Description:</b><br>
The largest number of bundle sets the expansion of a bundle partition may list: the connected sets of up to `order`
bundles and their subsets, single bundles included. A partition whose list would be longer stops the run before any
dynamics ("expansion_max_clusters exceeded"); the list is never cut short. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `200000`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "max_clusters": 200000}} ```

<a id="tag_bundle_max_work"></a>
## cluster.bundle.max_work
<b>Description:</b><br>
The largest number of elementary steps the expansion of a bundle partition may take: one per link read, per neighbour
tried while a connected set is grown, per subset visited and per power updated. Past it the run stops before any
dynamics ("expansion_max_work exceeded"). <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `20000000`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "max_work": 20000000}} ```

<a id="tag_adaptive"></a>
## cluster.adaptive
<b>Description:</b><br>
Turns on the adaptive refinement: instead of one calculation the run repeats calculating,
measuring the error and correcting the bundles until `atol` and `rtol` are met or `max_rounds`
is reached. Belongs to `cluster.mode` `pcce` or `bundle` (refused with another mode); the refinement driver itself runs
with `cluster.mode` = `bundle`, `bundle.max_dim` and an explicit [seed](InputFile_General.md#tag_seed). The keys are the
driver's options with `_` for `-` (and `workdir` for the driver's `output`); `workdir`, `atol` and `rtol` are required,
and every step goes into `workdir` (one folder per step and `summary.json`). All keys with their defaults are in [the keys of cluster.adaptive](#tag_adaptive_keys),
the results in [the results](#tag_adaptive_results), and the method (patches, probes, the local error, the merges,
regroup) in [adaptive refinement](../methods/bundle-adaptive.md#adaptive_round).
A top-level `adaptive` block, or the block under an older name, is refused (with the place it belongs). Its
[probe_memory_gb](#tag_adaptive_probe_memory_gb) is the memory of the probes, read by every step. <br>

<b>Data Type:</b> Object

<b>Options:</b><br>
| Option | Description |
|---|---|
| `workdir` | folder for the steps and the summary (required; `output` before 2.0.0) |
| `atol`, `rtol` | absolute and relative tolerance of the coherence (required) |
| `max_rounds` | the largest number of rounds, default 8 |
| `probe_memory_gb` | the memory of the probes running at once on a node ([its entry](#tag_adaptive_probe_memory_gb)) |

<b>Default:</b> - (absent: one calculation)

<b>Example:</b> ``` "cluster" : {"adaptive": { "workdir" : "./results/adaptive", "atol" : 1e-3, "rtol" : 0, "max_rounds" : 4 }} ```

<a id="tag_adaptive_keys"></a>
## The keys of cluster.adaptive
<b>Description:</b><br>
Every key of the block. The names are the refinement driver's options with `_` for `-` (the stand-alone
`ccex-adaptive` takes the same options as `--max-rounds` etc.), read without regard to case; in the input a key is
written with `_` only (`"max-rounds"` is refused with `"max_rounds"`), and the driver's `--output` is `workdir`
(`output` read like the top-level output block; it is refused with the new name). A `true/false` key is on
only when `true`. An unknown key, a value of the wrong type or a missing required key stops the run before anything is
computed (`Error: cluster.adaptive.<key> ...`). `input`, `binary`, `cwd`, `mpiexec` and `ranks` cannot be given: the
input is this file, the binary is the running process and the ranks are those of `mpirun`. The steps go into `workdir`:
one folder per step (`0001-run`, `0002-probe`, `0003-plan`, ...) and `summary.json`; the bath is read once, since the steps run in the
same process. The `input.json` a step writes has no `adaptive` block. `configuration.seed` is required (every round
recomputes the same bath states). <br>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `workdir` | path | required | folder of the steps and of `summary.json` (`output` before 2.0.0) |
| `atol` | real | required | absolute tolerance of the coherence: a patch passes when \|L<sub>expansion</sub> − L<sub>direct</sub>\| ≤ `atol` + `rtol` \|L<sub>direct</sub>\| at every bath state and every time point |
| `rtol` | real | required | relative tolerance of the coherence (`atol` + `rtol` > 0) |
| `policy` | string | `merge` | not needed: merge (merge the bundles whose merge reduces the local error most) is the only policy; a written one is reported as not used, the removed policies `error-gated`, `budgeted`, `auto-partition`, `cost-preserving` and `auto-order` are refused |
| `max_rounds` | integer | `8` | the largest number of rounds |
| `max_patches` | integer | `64` | the largest number of patches checked in one round (a patch: a few neighbouring bundles, at order 2 a bundle and two of its most strongly linked neighbours, whose expansion is compared with their direct calculation) |
| `max_probes` | integer | `4096` | the largest number of probe curves held (a probe: one bundle set of a patch computed directly as one cluster, all states and time points) |
| `max_candidates` | integer | `16` | the largest number of merge candidates per round |
| `neighbors` | integer | `8` | the number of neighbouring bundles a patch takes |
| `cache_mb` | integer | `256` | memory (MB) for the cached probe curves |
| `timeout` | real | `600` | time limit (s) of one probe run |
| `raw_cache` | path | - | file in which probe results are kept and reused across runs |
| `probe_batch_clusters` | integer | `0` | clusters per probe batch; `0` derives it from `cache_mb` |
| `merge_priority` | string | `gain-per-cost` | order of the merges: `gain-per-cost` or `worst-error` |
| `max_total_work` | real | `8000000000000.0` | work limit of the whole run (D³ × states × time points) |
| `prefer_physical` | true/false | off | prefer candidates whose curve stays physical (magnitude of L at most 1) |
| `hub_candidates` | true/false | off | add candidates built around a hub bundle (needs `prefer_physical`) |
| `regroup` | true/false | off | instead of merge's selection, regroup the bundles where the magnitude of L exceeds 1 (cCCE order 2 only; not with `prefer_physical` or `hub_candidates`) |
| `regroup_max_regions` | integer | `64` | with `regroup`: regions checked per round |
| `regroup_max_candidates` | integer | `256` | with `regroup`: candidates per region |
| `regroup_factors` | string | `probes` | with `regroup`: `probes` or `factors` (one `factors` step that returns the values at the worst state and time only, instead of probing the culprit pairs over all states) |
| `regroup_whole` | string | `probes` | with `regroup`: `probes` or `except` (one `except` step that returns the product without the terms a candidate changes, instead of probing every term; needs an x87 80-bit long double build) |
| `probe_memory_gb` | real | `0` (automatic) | the memory of the probes running at once on a node ([its entry](InputFile_Backend.md#tag_probe_memory_gb)) |

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8}, "adaptive": {"workdir": "./adaptive_out", "atol": 1e-8, "rtol": 2e-6}} ```

<a id="tag_adaptive_results"></a>
## The results of cluster.adaptive
<b>Description:</b><br>
What a refinement leaves in `workdir`, and how it ends ([details](../methods/bundle-adaptive.md#adaptive_stop)). <br>
- **The result is the last run**: `summary.json` names its folder in `final_directory`, and that folder's
  `coherence_wiDiv` (`coherence_state<s>_wiDiv` with sampled states) is the curve. The first run, without refinement,
  is `0001-run`.
- `summary.json`: `status`, the settings, one record per step (`runs`) and per round (`rounds`: the patches with their
  scaled errors and worst state and time step, the largest local error, the change since the previous round, the
  accepted merge with the local error before and after), `final_directory`, `final_bundles`, `limitations`.
  `final-input.json` (with `final-request.json`) reruns the chosen partition alone.
- **Status and exit status of `main.out`:** `validated_full_bath_model` and `passed_convergence_checks` (merge) and
  `physical` (regroup) exit 0; `round_limit`, `no_feasible_improving_merge`, `probe_budget_limited`, and for regroup
  `no_positive_culprit`, `no_candidate`, `budget_limited`, `change_limit` exit 2, the last run being the result;
  `failed` exits 1. With `configuration.k_range` the next configuration is not started after a non-zero status.
- **Cost:** every round computes the whole bath again, plus the probes of the patches checked; use it on baths whose
  single run is short. The local errors are checks of the patches, not error bounds of the whole bath, and in the cases
  studied so far the refinement changed the curves only a little; compare with a run without it.

`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_adaptive/` (merge; it ends with `no_feasible_improving_merge`, exit 2)
and `cCCE2_Hahn_regroup/` (regroup) show both. <br>

<a id="tag_adaptive_probe_memory_gb"></a>
## cluster.adaptive.probe_memory_gb
<b>Description:</b><br>
The memory, in GB per node, that the probes of the adaptive refinement running at once may hold (the driver copies it
into every step's input). The probes of a bath state go to the ranks in turn, and a probe waits while the probes
already running on its node would exceed this budget; a probe larger than the budget runs alone. `0` takes 80 % of the
memory the job may use on the node and, when clusters may go to the GPU, of the smallest GPU. The run prints the budget
(`[probe-budget]`). How a probe's memory is estimated: [Execution tags](InputFile_Backend.md#tag_probe_memory_gb). <br>

<b>Data Type:</b> Double (nonnegative)

<b>Unit:</b> GB per node

<b>Default:</b> `0` (automatic)

<b>Example:</b> ``` "cluster" : {"adaptive": {..., "probe_memory_gb": 60}} ```

<a id="tag_adaptive_work_file"></a>
## cluster.adaptive.work_file
<b>Description:</b><br>
The plan/probe protocol file the adaptive refinement driver writes for each step it runs. The driver sets it itself, in a step's input, in the `cluster.adaptive` block (beside the user's
`probe_memory_gb` only); a hand-written input does not need it, and it is refused beside the refinement options. <br>

<b>Data Type:</b> String (path)

<b>Default:</b> -

<b>Example:</b> ``` "cluster" : {"adaptive": {"work_file": "./work.json"}} ```
