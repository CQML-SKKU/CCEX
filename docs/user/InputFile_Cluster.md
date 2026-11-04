# Cluster tags

The `cluster` block says how the bath is cut into clusters: the mode, the `selection` every mode reads, the blocks of
the `pcce` and `bundle` modes, and the adaptive refinement. `order` and `method` stay at the top level
([General tags](InputFile_General.md#tag_order)). Each key has a small table (its type, default and unit), what it
does, and an example. Old key names are in [Old names](InputFile_OldNames.md). How the modes work is in Methods:
[cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundles and the adaptive refinement](../methods/bundle-adaptive.md).

1. [cluster.mode](#clu_mode): the six ways of cutting the bath.
2. [cluster.selection](#clu_selection): `add_subclusters` and `max_per_order`.
3. [cluster.pcce](#clu_pcce): the partitions of pCCE.
4. [cluster.bundle](#clu_bundle): the bundles of bundle clustering.
5. [cluster.adaptive](#clu_adaptive): the adaptive refinement, its keys and its results.

```json
"cluster": {"mode": "pcce", "selection": {"add_subclusters": true}, "pcce": {"partition_size": 2}}
```

<a id="clu_mode"></a>
## 1. cluster.mode

<a id="tag_cluster_mode"></a>
### cluster.mode

| Type | Default |
|---|---|
| string | `"standard"` |

How the bath is divided into clusters.

| Value | Clusters | Methods |
|---|---|---|
| `standard` | up to `order` spins connected within `rdip` | cCCE, gCCE, meCCE, megCCE |
| `pcce` | the bath is cut into partitions of [partition_size](#tag_pcce_partition_size) spins; a cluster is up to `order` partitions | cCCE, gCCE, meCCE, megCCE |
| `bundle` | strongly coupled spins are tied into bundles first; a cluster is up to `order` bundles | cCCE, gCCE |
| `dsj` | disjoint groups of at most `order` spins (strongest pairs within `rdip` first), each computed whole, no other cluster | cCCE |
| `itb` | only the pair corrections between connected spins of different groups (what `dsjitb` adds to `dsj`) | cCCE |
| `dsjitb` | disjoint groups plus the pairs between them; order 2 gives the result of standard CCE2 | cCCE |

Other combinations of method and mode are refused. gCCE and the master-equation methods run `pcce` on the CPU only.
The [pcce](#tag_pcce) block belongs to `pcce`, the [bundle](#tag_bundle) block to `bundle`, and the
[adaptive](#tag_adaptive) block to `pcce` and `bundle`; each is refused with another mode.

With `bundle`:

- [order](InputFile_General.md#tag_order) counts bundles, and [bundle.max_dim](#tag_bundle_max_dim) is required.
- With [rdip](InputFile_General.md#tag_rdip), bundles are connected by distance. Without `rdip`, a coherence check
  chooses the connections ([bundle.connections](#tag_bundle_connections)).
- With every bundle a single site, `bundle` is standard CCE.

A *site* is one bath spin, or a defect with its on-site spins; a *bundle* is a set of sites that always enters a
cluster whole; a *cluster* is up to `order` connected bundles.

```json
"cluster": {"mode": "bundle", "bundle": {"max_dim": 8}}
```

<a id="clu_selection"></a>
## 2. cluster.selection

<a id="tag_add_subclusters"></a>
### cluster.selection.add_subclusters

| Type | Default |
|---|---|
| boolean | `true` |

For order 3 and above: whether every sub-cluster of the kept clusters (each subset of two or more of its spins, or
partitions with `pcce`) is listed before the expansion is computed. That includes sub-clusters whose spins are not
all within `rdip` of each other, and those that [max_per_order](#tag_max_clusters_per_order) cut at a lower order. With
`false`, a sub-cluster the expansion needs is added when it is found missing. Both give the same result; `false` can
be a little faster. Not read with `cluster.mode` bundle, whose expansion always holds every subset.

```json
"cluster": {"selection": {"add_subclusters": false}}
```

<a id="tag_max_clusters_per_order"></a>
### cluster.selection.max_per_order

| Type | Default |
|---|---|
| list of `order` + 1 integers | all 0 (no limit) |

The largest number of clusters kept at each order, for order 3 and above, where the number of clusters grows fast.
Entry n limits the clusters of order n, and 0 means no limit. The list must have exactly `order` + 1 entries.

- Only entries 3 and up are applied, and only with `cluster.mode` standard and pcce (for pcce, orders of partitions).
  A nonzero entry that is not applied is ignored with a warning.
- `cluster.mode` bundle needs all entries 0.

The clusters of order n are grown from the kept clusters of order n − 1, ranked by the sum of \|J<sub>zz</sub>\| of
the point-dipole coupling over their spin pairs, and the strongest are kept.

```json
"cluster": {"selection": {"max_per_order": [0, 0, 0, 30000]}}
```

<a id="clu_pcce"></a>
## 3. cluster.pcce

<a id="tag_pcce"></a>
### cluster.pcce

| Type | Default |
|---|---|
| object | required with `cluster.mode` pcce (for `partition_size`) |

The settings of `cluster.mode` pcce. A k-means search divides the bath into partitions
([pCCE](../methods/pcce.md)). [Advanced 07](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) runs
`partition_size` 1 to 3.

```json
"cluster": {"mode": "pcce", "pcce": {"partition_size": 2, "trials": 1000, "connections": "spins"}}
```

<a id="tag_pcce_partition_size"></a>
### cluster.pcce.partition_size

| Type | Default |
|---|---|
| positive integer | required |

The number of bath spins in one partition. A pCCE cluster of order n is n partitions.

The spins left over are dropped. With N bath spins within `rbath` there are N / `partition_size` partitions, rounded
down, and the N mod `partition_size` spins left over are the last ones in the order of the bath file. They take no part
in the calculation, whatever their distance from the qubit, and the log does not list them. To keep every spin, choose
`rbath` or the bath file so that N is a multiple of `partition_size` ([pCCE](../methods/pcce.md#pcce_rest)).

```json
"cluster": {"pcce": {"partition_size": 3}}
```

<a id="tag_pcce_trials"></a>
### cluster.pcce.trials

| Type | Default |
|---|---|
| positive integer | `30000` |

The number of partition trials. Each trial is a k-means run from its own starting centres, and the best assignment is
kept.

```json
"cluster": {"pcce": {"partition_size": 2, "trials": 1000}}
```

<a id="tag_pcce_max_iterations"></a>
### cluster.pcce.max_iterations

| Type | Default |
|---|---|
| positive integer | `30000` |

The largest number of k-means iterations in one trial.

```json
"cluster": {"pcce": {"partition_size": 2, "max_iterations": 1000}}
```

<a id="tag_pcce_kmeans_pp"></a>
### cluster.pcce.kmeans_pp

| Type | Default |
|---|---|
| boolean | `true` |

Start every trial from k-means++ centres instead of random ones.

```json
"cluster": {"pcce": {"partition_size": 2, "kmeans_pp": false}}
```

<a id="tag_pcce_print_trials"></a>
### cluster.pcce.print_trials

| Type | Default |
|---|---|
| integer, 0 or more | `0` |

The number of trials whose details are printed; `0` prints none.

```json
"cluster": {"pcce": {"partition_size": 2, "print_trials": 8}}
```

<a id="tag_pcce_partition"></a>
### cluster.pcce.partition

| Type | Default |
|---|---|
| string | `"deterministic"` |

How the partition is searched.

| Value | Search |
|---|---|
| `deterministic` | the trials are spread over the ranks in a fixed way, so the partition follows `seed` |
| `legacy` | the random numbers and trial distribution of earlier versions |

```json
"cluster": {"pcce": {"partition_size": 2, "partition": "legacy"}}
```

<a id="tag_pcce_connections"></a>
### cluster.pcce.connections

| Type | Default |
|---|---|
| string | `"centers"` |

Which distance decides that two partitions are connected, and so can be in one cluster. Both rules use the inner and
outer `rdip`.

| Value | Connected when |
|---|---|
| `centers` | the centroids of the two partitions are within `rdip`; the clusters are ranked by the point-dipole T<sub>zz</sub> between the centroids |
| `spins` | any spin of one partition is within `rdip` of any spin of the other; the clusters are ranked by the sum of \|T<sub>zz</sub>\| over the spin pairs |

With `partition_size` 1 the two are the same. With `centers`, two partitions with close spins can stay unconnected
when their centroids are far apart; `spins` connects them, at a higher cost. The log says which rule ran and how many
pairs it connected.

```json
"cluster": {"pcce": {"partition_size": 2, "connections": "spins"}}
```

<a id="tag_pcce_cache"></a>
### cluster.pcce.cache

| Type | Default |
|---|---|
| object `{"mode", "file"}` | `{"mode": "off"}` |

Reuse a stored partition instead of searching again. It needs [partition](#tag_pcce_partition) `deterministic` and
method cCCE.

| mode | Does |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition, and stop if there is none |

`file` is the cache file, relative to the working directory or absolute, in a folder that exists; it is required with a
mode other than `off`. With [configuration.k_range](InputFile_General.md#tag_configuration_k_range) the name needs
`{k}`, one cache file per configuration.

```json
"cluster": {"pcce": {"partition_size": 2, "cache": {"mode": "readwrite", "file": "results/partition_{k}.cache"}}}
```

<a id="clu_bundle"></a>
## 4. cluster.bundle

<a id="tag_bundle"></a>
### cluster.bundle

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

The largest Hilbert dimension of one bundle. The dimension of a cluster is the product of its bundles' dimensions
(times the qubit dimension for gCCE).

```json
"cluster": {"bundle": {"max_dim": 8}}
```

<a id="tag_bundle_merge_cutoff"></a>
### cluster.bundle.merge_cutoff

| Type | Default | Unit |
|---|---|---|
| number, 0 or more | `0` | rad/ms; or `{"value", "unit"}` with `rad/s`, `rad/ms`, `rad/us` |

The smallest coupling strength at which two sites go into one bundle. The site pairs are taken strongest first; a pair
of strength at least `merge_cutoff` (and above 0) joins the bundles of its two sites when the product of their
dimensions is at most [max_dim](#tag_bundle_max_dim). With `0` every pair may merge and `max_dim` alone limits the
bundles; a large value keeps every site alone.

The strength of a site pair is Σ 3 I<sub>i</sub> I<sub>j</sub> ‖J<sub>ij</sub>‖<sub>F</sub> over the spins of the two
sites, with ‖J‖<sub>F</sub> the Frobenius norm of their point-dipole tensor
([the coupling strength](../methods/bundle-adaptive.md#bundle_coupling)).

```json
"cluster": {"bundle": {"max_dim": 8, "merge_cutoff": 1.0}}
```

<a id="tag_bundle_connections"></a>
### cluster.bundle.connections

| Type | Default |
|---|---|
| string | `"distance"` with `rdip`, `"checked"` without |

How bundles are connected into clusters. The log names the rule that ran.

| Value | Connections |
|---|---|
| `distance` | the pairs within `rdip` (`rdip` required) |
| `checked` | every link at first; the weakest are left out only if a coherence check allows it (`rdip` refused) |
| `auto` | the weakest pairs are left out up to [connection_budget](#tag_bundle_connection_budget), which is required |

`checked` proposes up to 8 of the weakest bundle pairs for omission and computes the domains around them directly, for
every bath state and time point. The omissions are made only if the expansion with them stays within
1e-7 + 1e-5 \|L\| of the direct curves; otherwise every link is kept. It costs one more pass over the bath states and
writes `[adaptive-check]` to the log and `<output.file>.adaptive.json`. It is a local check, not a convergence test of
the whole bath ([checked connections](../methods/bundle-adaptive.md#bundle_connections)).

```json
"cluster": {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}}
```

<a id="tag_bundle_connection_budget"></a>
### cluster.bundle.connection_budget

| Type | Default | Unit |
|---|---|---|
| number, 0 or more | required with `connections` `auto` | rad |

The allowed score of the bundle pairs left out with `connections` `auto`: t<sub>max</sub> × Σ w<sub>AB</sub> over the
pairs left out, with t<sub>max</sub> = `deltat` × (`nstep` − 1) and w<sub>AB</sub> the
[coupling strength](#tag_bundle_merge_cutoff) of the two bundles. Pairs outside `rdip` (when given) are left out
first, then the weakest of the rest while the score allows. The score measures coupling; it does not bound the error.

```json
"cluster": {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}}
```

<a id="tag_bundle_max_clusters"></a>
### cluster.bundle.max_clusters

| Type | Default |
|---|---|
| positive integer | `200000` |

The largest number of bundle sets the expansion may list (the connected sets of up to `order` bundles and their
subsets). A longer list stops the run before the calculation; it is never cut short.

```json
"cluster": {"bundle": {"max_dim": 8, "max_clusters": 500000}}
```

<a id="tag_bundle_max_work"></a>
### cluster.bundle.max_work

| Type | Default |
|---|---|
| positive integer | `20000000` |

The largest number of elementary steps the expansion may take while it lists the bundle sets. Past the limit the run
stops before the calculation.

```json
"cluster": {"bundle": {"max_dim": 8, "max_work": 50000000}}
```

<a id="clu_adaptive"></a>
## 5. cluster.adaptive

<a id="tag_adaptive"></a>
### cluster.adaptive

| Type | Default |
|---|---|
| object | no block: one calculation |

Turns on the adaptive refinement. The run repeats three steps until `atol` and `rtol` are met or `max_rounds` is
reached: calculate, measure the error of a few patches of the bath, and merge or regroup the bundles there. The method
is described in [adaptive refinement](../methods/bundle-adaptive.md#adaptive_round).

- The block belongs to `cluster.mode` pcce or bundle. The refinement runs with `cluster.mode` bundle,
  `bundle.max_dim` and a given [configuration.seed](InputFile_General.md#tag_seed), since every round recomputes the
  same bath states.
- `workdir`, `atol` and `rtol` are required. `-o` and `output.path` are refused with the block; every step goes into
  `workdir`.
- It does not run with [configuration.k_range](InputFile_General.md#tag_configuration_k_range); refine one
  configuration at a time.

```json
"cluster": {"mode": "bundle", "bundle": {"max_dim": 8},
            "adaptive": {"workdir": "results/adaptive", "atol": 1e-3, "rtol": 0, "max_rounds": 4}}
```

<a id="tag_adaptive_keys"></a>
### The keys of cluster.adaptive

Keys are read without regard to case and written with `_` (`"max-rounds"` is refused). An unknown key, a wrong type or
a missing required key stops the run before anything is computed.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `workdir` | path | required | the folder of the steps and of `summary.json` |
| `atol` | number | required | absolute tolerance: a patch passes when \|L<sub>expansion</sub> − L<sub>direct</sub>\| ≤ `atol` + `rtol` \|L<sub>direct</sub>\| at every bath state and time point |
| `rtol` | number | required | relative tolerance (`atol` + `rtol` > 0) |
| `max_rounds` | integer | `8` | the largest number of rounds |
| `max_patches` | integer | `64` | the largest number of patches checked in one round (a patch is a bundle and a few of its most strongly linked neighbours) |
| `max_probes` | integer | `4096` | the largest number of probe curves held (a probe is one bundle set of a patch computed directly) |
| `max_candidates` | integer | `16` | the largest number of merge candidates per round |
| `neighbors` | integer | `8` | the number of neighbouring bundles a patch takes |
| `cache_mb` | integer | `256` | the memory (MB) for cached probe curves |
| `timeout` | number | `600` | the time limit (s) of one probe run |
| `raw_cache` | path | none | a file in which probe results are kept and reused across runs |
| `probe_batch_clusters` | integer | `0` | clusters per probe batch; `0` derives it from `cache_mb` |
| `merge_priority` | string | `gain-per-cost` | the order of the merges: `gain-per-cost` or `worst-error` |
| `max_total_work` | number | `8e12` | the work limit of the whole run (D³ × states × time points) |
| `prefer_physical` | boolean | `false` | prefer candidates whose curve stays physical (\|L\| ≤ 1) |
| `hub_candidates` | boolean | `false` | also try candidates built around a strongly linked bundle (needs `prefer_physical`) |
| `regroup` | boolean | `false` | instead of merging, regroup the bundles where \|L\| exceeds 1 (cCCE order 2 only; not with `prefer_physical` or `hub_candidates`) |
| `regroup_max_regions` | integer | `64` | with `regroup`: the regions checked per round |
| `regroup_max_candidates` | integer | `256` | with `regroup`: the candidates per region |
| `regroup_factors` | string | `probes` | with `regroup`: `probes`, or `factors` to take the values at the worst state and time only |
| `regroup_whole` | string | `probes` | with `regroup`: `probes`, or `except` to compute the product without the terms a candidate changes (needs a build with an 80-bit long double) |
| `probe_memory_gb` | number | `0` (automatic) | the memory of the probes running at once ([probe_memory_gb](#tag_adaptive_probe_memory_gb)) |

`policy` is not needed: `merge` is the only policy, and a written one is reported as not used.

```json
"adaptive": {"workdir": "results/adaptive", "atol": 1e-8, "rtol": 2e-6, "regroup": true}
```

<a id="tag_adaptive_results"></a>
### The results of cluster.adaptive

The refinement writes one folder per step into `workdir` (`0001-run`, `0002-probe`, ...) and `summary.json`; the
first step, `0001-run`, is the calculation without refinement. The result is the last run: `summary.json` names its
folder in `final_directory`, whose `coherence_wiDiv` (`coherence_state<s>_wiDiv` with sampled states) is the curve.
`final-input.json` reruns the chosen bundles as one calculation. The files are listed in
[Output files](OutputFiles.md#out_adaptive).

| Exit code of `main.out` | Status |
|---|---|
| 0 | `validated_full_bath_model`, `passed_convergence_checks` (merge); `physical` (regroup) |
| 2 | stopped, the last run being the result: `round_limit`, `no_feasible_improving_merge`, `probe_budget_limited` (merge); `no_positive_culprit`, `no_candidate`, `budget_limited`, `change_limit` (regroup) |
| 1 | `failed` |

Every round computes the whole bath again, plus the probes, so use it on baths whose single run is short. The local
errors are checks of the patches, not error bounds of the whole bath; compare with a run without the refinement.
[Advanced 07](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) runs both the merge (it ends with
`no_feasible_improving_merge`, exit 2) and the regroup.

<a id="tag_adaptive_probe_memory_gb"></a>
### cluster.adaptive.probe_memory_gb

| Type | Default | Unit |
|---|---|---|
| number, 0 or more | `0` (automatic) | GB per node |

The memory the probes running at once on a node may hold. A probe waits while the probes already running on its node
would exceed it; a probe larger than the budget runs alone. `0` takes 80 % of the memory the job may use on the node
(and, when clusters may go to the GPU, of the smallest GPU). The run prints the budget (`[probe-budget]`). How a
probe's memory is estimated is under [Execution tags](InputFile_Backend.md#tag_probe_memory_gb).

```json
"adaptive": {"workdir": "results/adaptive", "atol": 1e-3, "rtol": 0, "probe_memory_gb": 60}
```

<a id="tag_adaptive_work_file"></a>
### cluster.adaptive.work_file (internal)

A file the refinement writes into the input of each of its steps. A hand-written input does not need it, and it is
refused beside the refinement's keys.
