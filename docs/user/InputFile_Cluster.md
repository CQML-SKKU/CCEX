# Cluster tags

The `cluster` block says how the bath is cut into clusters. `order` and `method` stay at the top level
([General tags](InputFile_General.md#tag_order)). Old key names are in [Upgrading old inputs](InputFile_OldNames.md). How the
modes work is in Methods: [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundles and the adaptive refinement](../methods/bundle-adaptive.md).

1. [cluster.mode](#clu_mode): the six ways of cutting the bath.
2. [cluster.selection](#clu_selection): `add_subclusters` and `max_per_order`.
3. [cluster.pcce](#clu_pcce): the partitions of pCCE.
4. [cluster.bundle](#clu_bundle): the bundles of bundle clustering.
5. [cluster.adaptive](#clu_adaptive): the adaptive refinement, its keys and its results.

<a id="clu_mode"></a>
## 1. cluster.mode

<a id="tag_cluster_mode"></a>
| Type | Default |
|---|---|
| string | `"standard"` |

How the bath is divided into clusters. The `pcce` and `bundle` blocks belong to their modes, and `adaptive` to
`bundle`.

| Value | Clusters | Methods |
|---|---|---|
| `standard` | up to `order` spins connected within `rdip` | cCCE, gCCE, meCCE, megCCE |
| `pcce` | up to `order` partitions of [partition_size](#tag_pcce_partition_size) spins | cCCE, gCCE, meCCE, megCCE |
| `bundle` | up to `order` bundles of strongly coupled spins (it needs [bundle.max_dim](#tag_bundle_max_dim)) | cCCE, gCCE |
| `dsj` | disjoint groups of at most `order` spins, each computed exactly as one cluster, and no other cluster | cCCE |
| `itb` | only the factor that the pairs between different groups add | cCCE |
| `dsjitb` | both, the groups and the pairs between them ([dsj, itb and dsjitb](../methods/cce.md#cce_dsj)) | cCCE |

With `bundle`, a bundle is a set of sites (a bath spin, or a defect with its on-site spins) that always enters a
cluster whole. Without `rdip`, a coherence check chooses which bundles are connected
([bundle.connections](#tag_bundle_connections)).

```json
"cluster": {"mode": "bundle", "bundle": {"max_dim": 8}}
```

<a id="clu_selection"></a>
## 2. cluster.selection

<a id="tag_add_subclusters"></a>
<a id="tag_max_clusters_per_order"></a>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `add_subclusters` | boolean | `true` | From order 3, every sub-cluster is listed first. `false` adds them when needed, with the same result and sometimes faster. Not read with bundle. |
| `max_per_order` | list of `order` + 1 integers | all 0 (no limit) | The most clusters kept at each order (entries 3 and up, standard and pcce). The strongest are kept. |

```json
"cluster": {"selection": {"add_subclusters": false, "max_per_order": [0, 0, 0, 30000]}}
```

<a id="clu_pcce"></a>
## 3. cluster.pcce

<a id="tag_pcce"></a>
| Type | Default |
|---|---|
| object | required with `cluster.mode` pcce (for `partition_size`) |

The settings of `cluster.mode` pcce. A k-means search divides the bath into partitions ([pCCE](../methods/pcce.md)).
[Advanced 10](../../example/tutorial/advanced/10_cluster_modes/README.md) runs `partition_size` 1 to 3.

```json
"cluster": {"mode": "pcce", "pcce": {"partition_size": 2, "trials": 1000, "connections": "spins"}}
```

<a id="tag_pcce_partition_size"></a>
### cluster.pcce.partition_size

| Type | Default |
|---|---|
| positive integer | required |

The number of bath spins in one partition. A cluster of order n is n partitions. The N mod `partition_size` spins
left over (the last in the bath file) take no part in the calculation, so choose `rbath` or the bath so that N is a
multiple of it ([pCCE](../methods/pcce.md#pcce_rest)).

```json
"cluster": {"pcce": {"partition_size": 3}}
```

<a id="tag_pcce_trials"></a>
<a id="tag_pcce_max_iterations"></a>
<a id="tag_pcce_kmeans_pp"></a>
<a id="tag_pcce_print_trials"></a>
<a id="tag_pcce_partition"></a>
<a id="tag_pcce_connections"></a>
### The search keys of cluster.pcce

| Key | Type | Default | Meaning |
|---|---|---|---|
| `trials` | positive integer | `30000` | The k-means trials, each from its own starting centres. The best partition is kept. |
| `max_iterations` | positive integer | `30000` | The most k-means iterations in one trial. |
| `kmeans_pp` | boolean | `true` | Start every trial from k-means++ centres instead of random ones. |
| `print_trials` | integer, 0 or more | `0` | The number of trials whose details are printed. |
| `partition` | string | `"deterministic"` | With `deterministic` the partition follows `seed`. `legacy` uses the random numbers of earlier versions. |
| `connections` | string | `"centers"` | Two partitions are connected when their centroids (`centers`) or any of their spins (`spins`, more costly) are within `rdip`. |

```json
"cluster": {"pcce": {"partition_size": 2, "trials": 1000, "connections": "spins"}}
```

<a id="tag_pcce_cache"></a>
### cluster.pcce.cache

| Type | Default |
|---|---|
| object `{"mode", "file"}` | `{"mode": "off"}` |

Reuse a stored partition instead of searching again. It needs [partition](#tag_pcce_partition) `deterministic` and
cCCE. `file` is required with a mode other than `off`, and needs `{k}` with
[configuration.k_range](InputFile_General.md#tag_configuration_k_range).

| mode | Does |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition |

```json
"cluster": {"pcce": {"partition_size": 2, "cache": {"mode": "readwrite", "file": "results/partition_{k}.cache"}}}
```

<a id="clu_bundle"></a>
## 4. cluster.bundle

<a id="tag_bundle"></a>
| Type | Default |
|---|---|
| object | required with `cluster.mode` bundle (for `max_dim`) |

The settings of `cluster.mode` bundle ([bundle clustering](../methods/bundle-adaptive.md)).

```json
"cluster": {"mode": "bundle", "bundle": {"max_dim": 8, "merge_cutoff": 1.0}}
```

<a id="tag_bundle_max_dim"></a>
### cluster.bundle.max_dim

| Type | Default |
|---|---|
| positive integer | required |

The largest Hilbert dimension of one bundle. A cluster's dimension is the product of its bundles' dimensions.

```json
"cluster": {"bundle": {"max_dim": 8}}
```

<a id="tag_bundle_connections"></a>
### cluster.bundle.connections

| Type | Default |
|---|---|
| string | `"distance"` with `rdip`, `"checked"` without |

How bundles are connected into clusters ([checked connections](../methods/bundle-adaptive.md#bundle_connections)).

| Value | Connections |
|---|---|
| `distance` | the pairs within `rdip` (`rdip` required) |
| `checked` | every link, except the weakest that a local coherence check allows to leave out (no `rdip`) |
| `auto` | the weakest pairs are left out up to [connection_budget](#tag_bundle_connection_budget) |

```json
"cluster": {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}}
```

<a id="tag_bundle_merge_cutoff"></a>
<a id="tag_bundle_connection_budget"></a>
<a id="tag_bundle_max_clusters"></a>
<a id="tag_bundle_max_work"></a>
### The other keys of cluster.bundle

| Key | Type | Default | Meaning |
|---|---|---|---|
| `merge_cutoff` | number, 0 or more (rad/ms, or `{"value", "unit"}`) | `0` | The smallest [coupling strength](../methods/bundle-adaptive.md#bundle_coupling) that puts two sites into one bundle within `max_dim`. With `0` only `max_dim` limits. |
| `connection_budget` | number, 0 or more (rad) | required with `connections` `auto` | The coupling the left-out pairs may carry, the last time times the sum of their coupling strengths. It is not an error bound. |
| `max_clusters` | positive integer | `200000` | The most bundle sets the expansion may list. More stops the run. |
| `max_work` | positive integer | `20000000` | The most steps the expansion may take to list them. More stops the run. |

```json
"cluster": {"bundle": {"max_dim": 8, "merge_cutoff": 1.0, "connections": "auto", "connection_budget": 0.01}}
```

<a id="clu_adaptive"></a>
## 5. cluster.adaptive

<a id="tag_adaptive"></a>

| Type | Default |
|---|---|
| object | none (one calculation) |

Turns on the adaptive refinement. The run calculates, checks a few patches of the bath and merges or regroups the
bundles there, until `atol` and `rtol` are met or `max_rounds` is reached
([adaptive refinement](../methods/bundle-adaptive.md#adaptive_round)). It needs `cluster.mode` bundle with
`bundle.max_dim`, `rdip` (or an explicit `bundle.connections`), method cCCE or gCCE with `quantity` coherence,
`nstep` 2 or more, a given [configuration.seed](InputFile_General.md#tag_seed), and one configuration (no `k_range`).
Every step goes into `workdir`, not into `-o` or `output.path`.

```json
"cluster": {"mode": "bundle", "bundle": {"max_dim": 8},
            "adaptive": {"workdir": "results/adaptive", "atol": 1e-3, "rtol": 0, "max_rounds": 4}}
```

<a id="tag_adaptive_keys"></a>
### The main keys

| Key | Type | Default | Meaning |
|---|---|---|---|
| `workdir` | path | required | The folder of the steps and of `summary.json`. |
| `atol`, `rtol` | numbers | required | The tolerance of a patch. A patch passes when \|L<sub>expansion</sub> − L<sub>direct</sub>\| ≤ `atol` + `rtol` \|L<sub>direct</sub>\| at every state and time. `atol` + `rtol` must be above 0. |
| `max_rounds` | integer | `8` | The largest number of rounds. |

### The limits of a round

| Key | Type | Default | Meaning |
|---|---|---|---|
| `max_patches` | integer | `64` | The most patches checked in a round, the most at risk first. |
| `neighbors` | integer | `8` | How many of a bundle's most strongly linked bundles are tried when its patch is formed. |
| `max_candidates` | integer | `16` | The most merge candidates tried in a round. |
| `max_probes` | integer | `4096` | The most probe curves (direct calculations of a patch's subsets) kept. |
| `cache_mb` | integer | `256` | The memory in MB for the cached probe curves. |
| `probe_batch_clusters` | integer | `0` | The probes computed in one batch. With `0` all the missing probes are computed in one batch. |
| `max_total_work` | number | `8e12` | A patch whose missing probes would need more work than this is not checked. The work is the sum of D³ × states × time points over those probes, with D the dimension including the qubit. |
| `timeout` | number | `600` | The time limit in seconds of one probe run. |

<a id="tag_adaptive_probe_memory_gb"></a>
<a id="tag_adaptive_work_file"></a>
### The other keys

| Key | Type | Default | Meaning |
|---|---|---|---|
| `merge_priority` | `"gain-per-cost"` or `"worst-error"` | `"gain-per-cost"` | The order of the merge candidates. `gain-per-cost` divides the gain by the cost D³, and `worst-error` puts the largest local error first. |
| `prefer_physical` | boolean | `false` | Prefer candidates whose curve stays finite and bounded. |
| `hub_candidates` | boolean | `false` | Also try growing a bundle around a strongly linked hub first. It needs `prefer_physical`. |
| `raw_cache` | path | none | A file that keeps the probe curves across runs. |
| `probe_memory_gb` | number (GB per node) | `0` (80 % of what the job may use) | The memory of the probes running at once on a node ([Execution tags](InputFile_Backend.md#tag_probe_memory_gb)). |
| `work_file` | path | none | Internal. The refinement writes it into the input of each step. |

### Regroup

| Key | Type | Default | Meaning |
|---|---|---|---|
| `regroup` | boolean | `false` | Regroup the bundles where \|L\| exceeds 1 instead of merging. It needs cCCE, order 2 and `connections` distance with `rdip`, and does not go with `prefer_physical` or `hub_candidates`. |
| `regroup_max_regions` | integer | `64` | The most regions regrouped in a round. |
| `regroup_max_candidates` | integer | `256` | The most candidates tried in a region. |
| `regroup_factors`, `regroup_whole` | `"probes"`, or `"factors"` and `"except"` | `"probes"` | How the culprit values and the whole curves are obtained. `"except"` needs a build whose long double is the x87 80-bit type. |

```json
"adaptive": {"workdir": "results/adaptive", "atol": 1e-8, "rtol": 2e-6, "regroup": true}
```

<a id="tag_adaptive_results"></a>
### The results

The final curve is `coherence_wiDiv` in the folder that `summary.json` names as `final_directory`
([Output files](OutputFiles.md#out_adaptive), [exit codes](BuildingAndRunning.md#ExitCodes)). The patch checks are not
an error bound of the whole bath, so compare with a run without the refinement.
[Advanced 10](../../example/tutorial/advanced/10_cluster_modes/README.md) runs the merge and the regroup.
