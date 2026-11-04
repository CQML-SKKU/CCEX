# Partition CCE (pCCE)

pCCE is `cluster.mode` pcce, and its keys are in [cluster.pcce](../user/InputFile_Cluster.md#tag_pcce). The bath is
cut into *partitions* of exactly [partition_size](../user/InputFile_Cluster.md#tag_pcce_partition_size) spins that are close
in space. The [cluster expansion](cce.md#cce_expansion) then runs over partitions instead of single spins:

- A cluster of order n is n connected partitions, that is n × `partition_size` spins.
- A partition always enters a cluster whole, so the couplings inside it are treated exactly.
- With `partition_size` 1, pCCE is standard CCE.

pCCE runs with cCCE, gCCE, meCCE and megCCE. The gCCE and master-equation forms run on the CPU only, and cCCE also runs
on the GPU ([GPU execution](gpu.md)). The partition search and the cluster list are always made on the CPU.

1. [The spins that do not fill a partition](#pcce_rest): which spins are dropped.
2. [The partition search](#pcce_search): the k-means search and its random numbers.
3. [Which partitions are connected](#pcce_connections): the two rules, centers and spins.
4. [The partition cache](#pcce_cache): when a stored partition is reused.

<a id="pcce_rest"></a>
## 1. The spins that do not fill a partition

With N bath spins (the spins within `rbath`) there are $\lfloor N / K \rfloor$ partitions of K = `partition_size`
spins. The remaining N mod K spins are dropped before the search and take no part in the calculation. They are the
last spins of the bath file among those within `rbath`, whatever their position. A bath file sorted by distance from
the qubit therefore loses its outermost spins, and an unsorted one may lose a near one.

The run log names them in one line, `pCCE spins left out : ...`, with the bath index, the species and the distance to
the nearest qubit of the first eight and the count of the rest. To keep every spin, choose `rbath` (or the bath file)
so that N is a multiple of `partition_size`.

<a id="pcce_search"></a>
## 2. The partition search

The partitions are found by a k-means search on the spin positions with equal-size clusters:

1. A trial shuffles the spins and picks $\lfloor N/K \rfloor$ starting centres, by k-means++
   ([kmeans_pp](../user/InputFile_Cluster.md#tag_pcce_kmeans_pp), the default) or at random.
2. Assignment: each spin in turn goes to the nearest centre that still has room (every centre takes exactly K spins).
   Update: each centre moves to the mean position of its spins. The two steps repeat until no centre moves by more than
   1e-4 Å in any coordinate, until a cycle of assignments repeats, or for
   [max_iterations](../user/InputFile_Cluster.md#tag_pcce_max_iterations).
3. Of all [trials](../user/InputFile_Cluster.md#tag_pcce_trials), the partition with the smallest inertia (the sum of
   squared distances of the spins from their centre) is kept. A tie goes to the larger silhouette coefficient.

[partition](../user/InputFile_Cluster.md#tag_pcce_partition) chooses how the trials draw their random numbers:

| `partition` | Random numbers | Trials per rank | Result |
|---|---|---|---|
| `deterministic` (default) | Trial t has its own stream from `seed` and t. | The trials are spread over the ranks in turn. | It follows `seed`. |
| `legacy` | Each rank's C library random stream, as earlier versions did. | trials / ranks each. | It can change with the number of ranks. |

[print_trials](../user/InputFile_Cluster.md#tag_pcce_print_trials) prints the course of the trials.

<a id="pcce_connections"></a>
## 3. Which partitions are connected

Two partitions can be in one cluster of order 2 or more when they are connected. There are two rules
([connections](../user/InputFile_Cluster.md#tag_pcce_connections)), and in both the inner and the outer `rdip` bound
the distance:

| Rule | Connected when | Strength (ranks the clusters for `max_per_order`) |
|---|---|---|
| `centers` (default) | The distance between the two centroids is within `rdip`. | The $\lvert J_{zz}\rvert$ of a point dipole between the centroids, with a fixed, equal γ. It is a relative measure. |
| `spins` | Any spin of the one partition is within `rdip` of any spin of the other. | The sum of $\lvert J_{zz}\rvert$ over the spin pairs of the two partitions. |

With `centers`, two partitions whose spins are close can stay unconnected when their centroids are far apart. This
happens more often for a larger `partition_size`, and `spins` connects them. `spins` looks at every spin pair, a cost
that grows as N², and usually connects more partition pairs, so it gives more clusters.

The rank-0 log prints the rule and how many partition pairs it connected (`pCCE connections : ...`). The clusters are
then made from this graph as in standard CCE, with
[max_per_order](../user/InputFile_Cluster.md#tag_max_clusters_per_order) and
[add_subclusters](../user/InputFile_Cluster.md#tag_add_subclusters) counting partitions.

<a id="pcce_cache"></a>
## 4. The partition cache

The search with many trials can take longer than the dynamics. [cache](../user/InputFile_Cluster.md#tag_pcce_cache)
stores the partition in a file and reuses it (method cCCE, partition `deterministic`). A stored partition is reused only
when everything that decides it is the same:

- The search algorithm's version.
- The compiler and its floating-point settings.
- `seed`, `partition_size`, `trials`, `max_iterations` and `kmeans_pp`.
- Every bath spin (name, spin and position), including the dropped ones.

With `configuration.k_range` the file name holds `{k}`, so there is one cache per configuration.

`example/tutorial/advanced/07_pCCE_and_bundle/pCCE2_Hahn_states/` runs `partition_size` 1 to 3 on a P1 bath, and
`example/tutorial/advanced/11_GPU/pCCE2_Hahn_backends/` runs pCCE on the GPU.

References: to be added by the CCEX authors.
