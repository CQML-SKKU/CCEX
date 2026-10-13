# Cluster-related tags
How the bath is cut into clusters: the [cluster_mode](#tag_cluster_mode), the keys every mode reads, and the blocks of
the `pcce` and `bundle` modes and of the adaptive refinement. The keys before 2026-10-13 are in the table of
[InputFile.md](InputFile.md#OldNames).

- [order](InputFile_General.md#tag_order)        
- [method](InputFile_General.md#tag_method)       
- [cluster_mode](#tag_cluster_mode)
- [add_subclusters](#tag_add_subclusters)   
- [max_clusters_per_order](#tag_max_clusters_per_order)           
- [pcce](#tag_pcce): [partition_size](#tag_pcce_partition_size), [trials](#tag_pcce_trials),
  [max_iterations](#tag_pcce_max_iterations), [kmeans_pp](#tag_pcce_kmeans_pp), [print_trials](#tag_pcce_print_trials),
  [print_limit](#tag_pcce_print_limit), [partition](#tag_pcce_partition), [connections](#tag_pcce_connections),
  [cache](#tag_pcce_cache) (`mode`, `file`)
- [bundle](#tag_bundle): [max_dim](#tag_bundle_max_dim), [merge_cutoff](#tag_bundle_merge_cutoff),
  [connections](#tag_bundle_connections), [connection_budget](#tag_bundle_connection_budget),
  [max_clusters](#tag_bundle_max_clusters), [max_work](#tag_bundle_max_work)
- [adaptive](#tag_adaptive) (and its internal [work_file](#tag_adaptive_work_file))

<a id="tag_cluster_mode"></a>
## cluster_mode
<b>Description:</b><br>
How the bath is divided into clusters (new values on 2026-10-13: the cluster shapes that were values of `method`, and
`bundle` for what was `adaptive`). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description | method |
|---|---|---|
| `standard` | clusters of up to `order` spins within `rdip` | cCCE, gCCE, meCCE, megCCE |
| `pcce` | partition CCE: the bath is cut into partitions of [pcce.partition_size](#tag_pcce_partition_size) spins, a cluster is up to `order` partitions (`method: "pCCE"` before) | cCCE |
| `bundle` | bundle clustering: strongly coupled spins are tied into bundles first, a cluster is up to `order` bundles (`cluster_mode: "adaptive"` before) | cCCE, gCCE |
| `dsj` | disjoint clusters (`method: "dsj"` before) | cCCE |
| `itb` | inter-bath clusters, the pair interaction between disjoint clusters | cCCE |
| `dsjitb` | disjoint + inter-bath (order 2 gives the same result as standard CCE2) | cCCE |

gCCE with `pcce`, `dsj`, `itb` or `dsjitb` is refused (not implemented yet), also when `-m gCCE` makes the run gCCE;
meCCE and megCCE with any mode but `standard` likewise.
The [pcce](#tag_pcce) block belongs to `pcce`, the [bundle](#tag_bundle) block to `bundle`, the [adaptive](#tag_adaptive)
block to `pcce` and `bundle`; each is refused with another mode. `"adaptive"` as a value is refused with its new name.

With `bundle`, [order](InputFile_General.md#tag_order) counts bundles per cluster and
[bundle.max_dim](#tag_bundle_max_dim) (required) caps the Hilbert dimension of one bundle. With
[rdip](InputFile_General.md#tag_rdip) given, bundles are connected by distance; with `rdip` omitted,
connections are chosen by a coherence check and the weak ones are dropped. Adding the [adaptive](#tag_adaptive) block turns
the adaptive refinement on top of it. The words are defined in
[Clustering terminology](../methods/clustering-terminology.md), the method in [Adaptive cluster mode](../methods/adaptive-cluster-mode.md). <br>

<b>Default:</b> `standard`

<b>Example:</b> ``` "cluster_mode" : "bundle", "bundle" : {"max_dim": 8} ```

<b>Related variable:</b> [Config](#Config)::cluster_mode (inside: `adaptive` for `bundle`), [Config](#Config)::method

<a id="tag_add_subclusters"></a>
## add_subclusters
<b>Description:</b><br>
For CCE simulations with order ≥ 3, enabling this option allows all sub-clusters of the highest order to be included
(`addsubclus` until 2026-10-13, which also took 0 or 1). <br>
This leads to more stable and robust results, and may also reduce computation time. <br>
<b>Data Type:</b> Boolean (`true` or `false`)

<b>Default:</b> `true`

<b>Example:</b> ``` "add_subclusters" : true ```

<b>Related variable:</b> [Cluster](#Cluster)::addsubclus

<a id="tag_max_clusters_per_order"></a>
## max_clusters_per_order
<b>Description:</b><br>
Number of clusters for each CCE order (`nk` until 2026-10-13). <br>
This is useful for CCE simulations with order ≥ 3, leading tons of clusters <br>
The list has exactly `order + 1` entries and the index is the cluster order: entry n limits the clusters of order n
(entry 0 is the order-0 slot and is kept, usually 0). A 0 means no limit for that order, so all zeros include every
cluster. A list of another length is refused. <br>
<b>Data Type:</b> Int List (length `order + 1`)

<b>Default:</b> all entries 0 (no limit)

<b>Example:</b> ``` "max_clusters_per_order" : [0, 0, 0, 30000] ``` (order 3: at most 30000 clusters of three spins)

<b>Related variable:</b> [Cluster](#Cluster)::nk

<a id="tag_pcce"></a>
## pcce
<b>Description:</b><br>
The settings of `cluster_mode` = `pcce` (the top-level keys `sK`, `max_trial`, ... until 2026-10-13). Refused with any
other `cluster_mode`. `example/tutorial/04_Diamond_NV_P1Bath/SpinEcho_pCCE2_1ppm/` runs `partition_size` 1 to 4. <br>

<b>Data Type:</b> Object

<b>Example:</b> ``` "cluster_mode" : "pcce", "pcce" : {"partition_size": 2, "trials": 1000, "max_iterations": 1000, "connections": "spins"} ```

<a id="tag_pcce_partition_size"></a>
## pcce.partition_size
<b>Description:</b><br>
The number of bath spins in one partition (`sK` until 2026-10-13). The bath is divided into partitions of
`partition_size` spins by a k-means search, and a pCCE cluster of order n is n partitions, n × `partition_size`
spins. Required with `cluster_mode` = `pcce`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "pcce" : {"partition_size": 3} ```

<b>Related variable:</b> [Cluster](#Cluster)::sK

<a id="tag_pcce_trials"></a>
## pcce.trials
<b>Description:</b><br>
The number of partition trials (`max_trial` until 2026-10-13), each a k-means run from its own starting centres;
the assignment with the best score is kept. With [pcce.partition](#tag_pcce_partition) =
`deterministic` the trials are split over the ranks in a fixed way, so the partition follows
`seed`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "trials": 1000} ```

<b>Related variable:</b> [Cluster](#Cluster)::max_trial

<a id="tag_pcce_max_iterations"></a>
## pcce.max_iterations
<b>Description:</b><br>
The largest number of k-means iterations in one partition trial (`max_iter` until 2026-10-13). <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "max_iterations": 1000} ```

<b>Related variable:</b> [Cluster](#Cluster)::max_iter

<a id="tag_pcce_kmeans_pp"></a>
## pcce.kmeans_pp
<b>Description:</b><br>
Start every k-means trial from k-means++ centres instead of random ones. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "kmeans_pp": true} ```

<b>Related variable:</b> [Cluster](#Cluster)::kmeans_pp

<a id="tag_pcce_print_trials"></a>
## pcce.print_trials
<b>Description:</b><br>
Rank 0 prints the details of the partition trials, at most [print_limit](#tag_pcce_print_limit) of them
(`iter_detail` until 2026-10-13). <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "print_trials": true} ```

<b>Related variable:</b> [Cluster](#Cluster)::iter_detail

<a id="tag_pcce_print_limit"></a>
## pcce.print_limit
<b>Description:</b><br>
With [print_trials](#tag_pcce_print_trials): rank 0 prints the details of at most this many
trials; `0` prints none (`pcce_detail_limit` until 2026-10-13). <br>

<b>Data Type:</b> Integer (nonnegative)

<b>Default:</b> `32`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "print_trials": true, "print_limit": 8} ```

<b>Related variable:</b> [Config](#Config)::pcce_detail_limit

<a id="tag_pcce_partition"></a>
## pcce.partition
<b>Description:</b><br>
How the partition of the bath is searched (`pcce_partition` until 2026-10-13). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `deterministic` | the partition follows `seed`; trials are distributed deterministically over the ranks |
| `legacy` | the original random-number and trial distribution |

<b>Default:</b> `deterministic`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "partition": "deterministic"} ```

<b>Related variable:</b> [Config](#Config)::pcce_partition

<a id="tag_pcce_connections"></a>
## pcce.connections
<b>Description:</b><br>
Which distance decides that two partitions are connected, so that they can be
in one cluster of order 2 or more (`pcce_connections` until 2026-10-13). The inner and the outer `rdip` bound that
distance under either rule. (`cluster_mode` bundle connects its bundles by [bundle.connections](#tag_bundle_connections).) <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `centers` | the distance between the two k-means centres, the centroids of the partitions. The strength that orders the clusters is the point-dipole T<sub>zz</sub> between the centroids. |
| `spins` | any spin of the one partition between the inner and the outer `rdip` of any spin of the other, the rule of the bundles. The strength is the sum of \|T<sub>zz</sub>\| over the spin pairs of the two partitions. |

With `partition_size` = 1 the two rules are the same. Under `centers` two partitions whose spins
are close can stay unconnected when their centroids are far apart, and the more so
the larger `partition_size`; `spins` connects them. The rank-0 log says which rule ran and how
many partition pairs it connected. The comparison of the two rules on the P1 bath of tutorial 04
is (internal notes): at that tutorial's `rdip` (620 Å) the two rules
differ less than pCCE2 differs from CCE4, and `spins` costs 27-65 % more; at an `rdip` comparable
to the partition size, `spins` is the closer of the two to CCE4.

<b>Default:</b> `centers` (the rule of every earlier version)

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "connections": "spins"} ```

<b>Related variable:</b> [Config](#Config)::pcce_connections

<a id="tag_pcce_cache"></a>
## pcce.cache
<b>Description:</b><br>
Reuse a stored partition instead of searching again (`pcce_cache_mode` and `pcce_cache_file` until 2026-10-13). Needs
[pcce.partition](#tag_pcce_partition) = `deterministic`. `mode`: <br>

| mode | Description |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition; stop if there is none |

`file`: the partition cache file, relative to the working directory or absolute; its directory must exist. Required
with `mode` other than `off`. <br>

<b>Data Type:</b> Object `{"mode": ..., "file": ...}`

<b>Default:</b> `{"mode": "off"}`

<b>Example:</b> ``` "pcce" : {"partition_size": 2, "cache": {"mode": "readwrite", "file": "./results/partition.cache"}} ```

<b>Related variable:</b> [Config](#Config)::pcce_cache_mode, [Config](#Config)::pcce_cache_file

<a id="tag_bundle"></a>
## bundle
<b>Description:</b><br>
The settings of `cluster_mode` = `bundle` (the top-level keys `bundle_*` and `expansion_*` until 2026-10-13). Refused
with any other `cluster_mode`. `max_dim` is required. <br>

<b>Data Type:</b> Object

<b>Example:</b> ``` "cluster_mode" : "bundle", "bundle" : {"max_dim": 8, "merge_cutoff": 1.0} ```

<a id="tag_bundle_max_dim"></a>
## bundle.max_dim
<b>Description:</b><br>
The largest Hilbert dimension one bundle may have (`bundle_max_dim` until 2026-10-13). Required with `cluster_mode` =
`bundle`. The dimension of a cluster follows as the product of its bundle dimensions (gCCE: times the
qubit dimension); a cluster that cannot be handled is refused during the run. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "bundle" : {"max_dim": 8} ```

<b>Related variable:</b> [Config](#Config)::bundle_max_dim

<a id="tag_bundle_merge_cutoff"></a>
## bundle.merge_cutoff
<b>Description:</b><br>
The smallest coupling strength for which two groups are merged into one bundle (`bundle_merge_cutoff` until 2026-10-13). <br>

<b>Data Type:</b> Double (nonnegative)

<b>Unit:</b> rad/ms

<b>Default:</b> `0`

<b>Example:</b> ``` "bundle" : {"max_dim": 8, "merge_cutoff": 1.0} ```

<b>Related variable:</b> [Config](#Config)::bundle_merge_cutoff

<a id="tag_bundle_connections"></a>
## bundle.connections
<b>Description:</b><br>
How bundles are connected into clusters (`bundle_connections` until 2026-10-13). Omit it together
with `rdip` to get the coherence-checked connections (the run log then says `checked`; that
word is not a value to write). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `distance` | pairs within `rdip`; `rdip` is required |
| `auto` | the legacy budget mode: weak pairs are left out up to [connection_budget](#tag_bundle_connection_budget), which is required |

<b>Default:</b> `distance`

<b>Example:</b> ``` "bundle" : {"max_dim": 8, "connections": "auto", "connection_budget": 0.01} ```

<b>Related variable:</b> [Config](#Config)::bundle_connections

<a id="tag_bundle_connection_budget"></a>
## bundle.connection_budget
<b>Description:</b><br>
With `connections` = `auto`, the allowed score of the left-out pairs: the largest
time times the sum of their coupling weights must stay below it (`bundle_connection_budget` until 2026-10-13). <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> -

<b>Example:</b> ``` "bundle" : {"max_dim": 8, "connections": "auto", "connection_budget": 0.01} ```

<b>Related variable:</b> [Config](#Config)::bundle_connection_budget

<a id="tag_bundle_max_clusters"></a>
## bundle.max_clusters
<b>Description:</b><br>
The largest number of clusters the expansion of a bundle partition may list (`expansion_max_clusters` until 2026-10-13).
`order` may not exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `200000`

<b>Example:</b> ``` "bundle" : {"max_dim": 8, "max_clusters": 200000} ```

<b>Related variable:</b> [Config](#Config)::expansion_max_clusters

<a id="tag_bundle_max_work"></a>
## bundle.max_work
<b>Description:</b><br>
The work limit of the expansion and of the dimension search of a bundle partition (`expansion_max_work` until 2026-10-13). <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `20000000`

<b>Example:</b> ``` "bundle" : {"max_dim": 8, "max_work": 20000000} ```

<b>Related variable:</b> [Config](#Config)::expansion_max_work

<a id="tag_adaptive"></a>
## adaptive
<b>Description:</b><br>
Turns on the adaptive refinement: instead of one calculation the run repeats calculating,
measuring the error and correcting the bundles until `atol` and `rtol` are met or `max_rounds`
is reached. Belongs to `cluster_mode` `pcce` or `bundle` (refused with another mode); the refinement driver itself runs
with `cluster_mode` = `bundle`, `bundle.max_dim` and an explicit [seed](InputFile_General.md#tag_seed). The keys are the
driver's options with `_` for `-`; `output`, `atol` and `rtol` are required, and every step goes into `output` (one folder
per step and `summary.json`). The full key table with defaults is in [the `adaptive` block](../methods/adaptive-input-block.md).
The block was written `Adaptive` until 2026-10-13 (keys are read without regard to case, so that spelling still reads)
and `Refine` until 2026-09-23 (that name is refused). <br>

<b>Data Type:</b> Object

<b>Options:</b><br>
| Option | Description |
|---|---|
| `policy` | `merge`, the only policy and the default. `error-gated`, `budgeted`, `auto-partition` and `cost-preserving` were removed on 2026-10-10 and `auto-order` on 2026-09-25; naming one is refused |
| `output` | folder for the steps and the summary (required) |
| `atol`, `rtol` | absolute and relative tolerance of the coherence (required) |
| `max_rounds` | the largest number of rounds, default 8 |

<b>Default:</b> - (absent: one calculation)

<b>Example:</b> ``` "adaptive" : { "policy" : "merge", "output" : "./results/adaptive", "atol" : 1e-3, "rtol" : 0, "max_rounds" : 4 } ```

<a id="tag_adaptive_work_file"></a>
## adaptive.work_file
<b>Description:</b><br>
The plan/probe protocol file the adaptive refinement driver writes for each step it runs (`adaptive_work_file` until
2026-10-13). The driver sets it itself, in a step's input, as the only key of the `adaptive` block; a hand-written input
does not need it, and it is refused beside the refinement options. <br>

<b>Data Type:</b> String (path)

<b>Default:</b> -

<b>Example:</b> ``` "adaptive" : {"work_file": "./work.json"} ```

<b>Related variable:</b> [Config](#Config)::adaptive_work_file
