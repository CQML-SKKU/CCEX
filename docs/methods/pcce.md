# Partition CCE (pCCE)

`cluster.mode` pcce (keys: [cluster.pcce](../user/InputFile_Cluster.md#tag_pcce)). The bath is cut into
*partitions* of exactly [partition_size](../user/InputFile_Cluster.md#tag_pcce_partition_size) spins that are close
in space, and the [cluster expansion](cce.md#cce_expansion) runs over partitions instead of single spins: a cluster of
order n is n connected partitions, n × `partition_size` spins, and a partition always enters a cluster whole, so the
couplings inside it are treated exactly. With `partition_size` 1 pCCE is standard CCE.

pCCE runs with cCCE, gCCE, meCCE and megCCE; the gCCE and master-equation forms on the CPU only, cCCE also on the GPU
([GPU execution](gpu.md)). The partition search and the cluster list are made on the CPU in every case.

<a id="pcce_rest"></a>
## The spins that do not fill a partition

With N bath spins (the spins within `rbath`) there are $\lfloor N / K \rfloor$ partitions of K = `partition_size`
spins. The remaining **N mod K spins are dropped** before the search: the last ones of the bath in the order of the
bath file (among the spins within `rbath`). They are in no cluster and take no part in the calculation. For a large
bath and a small `partition_size` this is a few spins at the edge of the file; but which spins are last depends only
on how the file is written, not on where they are, so a bath file sorted by distance from the qubit loses its outermost
spins, while an unsorted one may lose a near one. The run log does not name the dropped spins; to keep every spin,
choose `rbath` (or the bath file) so that N is a multiple of `partition_size`.

<a id="pcce_search"></a>
## The partition search

The partitions are found by a k-means search on the spin positions with equal-size clusters:

1. A trial shuffles the spins and picks $\lfloor N/K \rfloor$ starting centres, by k-means++
   ([kmeans_pp](../user/InputFile_Cluster.md#tag_pcce_kmeans_pp), the default) or at random.
2. Assignment: each spin in turn goes to the nearest centre that still has room (every centre takes exactly K spins).
   Update: each centre moves to the mean position of its spins. The two steps repeat until no centre moves by more than
   1e-4 Å in any coordinate, until a cycle of assignments repeats, or for
   [max_iterations](../user/InputFile_Cluster.md#tag_pcce_max_iterations).
3. Of all [trials](../user/InputFile_Cluster.md#tag_pcce_trials), the partition with the smallest inertia (the sum of
   squared distances of the spins from their centre) is kept; a tie goes to the larger silhouette coefficient.

[partition](../user/InputFile_Cluster.md#tag_pcce_partition) `deterministic` (the default) gives trial t its own
random stream from `seed` and t, and spreads the trials over the ranks in turn, so the partition follows `seed`. `legacy` splits the trials evenly over the ranks (trials / ranks each)
and draws from each rank's C library random stream, as earlier versions did, so the partition can change with the
number of ranks. [print_trials](../user/InputFile_Cluster.md#tag_pcce_print_trials) prints the course of the trials.

<a id="pcce_connections"></a>
## Which partitions are connected

Two partitions can be in one cluster of order 2 or more when they are connected, by one of two rules
([connections](../user/InputFile_Cluster.md#tag_pcce_connections)); the inner and the outer `rdip` bound the distance
in both:

| Rule | Connected when | Strength (ranks the clusters for `max_per_order`) |
|---|---|---|
| `centers` (default) | the distance between the two centroids is within `rdip` | $\lvert J_{zz}\rvert$ of a point dipole between the centroids (with a fixed, equal γ: a relative measure) |
| `spins` | any spin of the one partition is within `rdip` of any spin of the other | the sum of $\lvert J_{zz}\rvert$ over the spin pairs of the two partitions |

With `centers` two partitions whose spins are close can stay unconnected when their centroids are far apart, more so
for a larger `partition_size`; `spins` connects them. `spins` looks at every spin pair (a cost that grows as N²) and
usually connects more partition pairs, so it gives more clusters. The rank-0 log prints the rule and how many partition pairs it
connected (`pCCE connections : ...`). The clusters themselves are then made from this graph as
in standard CCE, with [max_per_order](../user/InputFile_Cluster.md#tag_max_clusters_per_order) and
[add_subclusters](../user/InputFile_Cluster.md#tag_add_subclusters) counting partitions.

<a id="pcce_cache"></a>
## The partition cache

A search with many trials can take longer than the dynamics. [cache](../user/InputFile_Cluster.md#tag_pcce_cache)
stores the partition in a file and reuses it (method cCCE, partition `deterministic`). A stored partition is reused only
when everything that decides it is the same: the search algorithm's version, the compiler and its floating-point
settings, `seed`, `partition_size`, `trials`, `max_iterations`, `kmeans_pp`, and every bath spin (name, spin and
position) including the dropped ones. With `configuration.k_range` the file name holds `{k}`: one cache per
configuration.

`example/tutorial/07_pCCE_and_bundle/pCCE2_Hahn_states/` runs `partition_size` 1 to 3 on a P1 bath;
`example/tutorial/08_GPU/pCCE2_Hahn_backends/` runs pCCE on the GPU.

References: to be added by the CCEX authors.
