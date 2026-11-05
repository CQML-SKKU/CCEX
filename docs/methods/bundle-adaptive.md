# Bundle clustering and adaptive refinement

Bundle clustering is `cluster.mode` bundle, with its keys in [cluster.bundle](../user/InputFile_Cluster.md#tag_bundle).
It puts strongly coupled bath spins into *bundles* that always enter a cluster whole, and makes the bundles once, from
the couplings. The adaptive refinement, the [cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive) block, runs
on top of it. It repeats the calculation, measures the error of the expansion locally and changes the bundles.

1. [Bundle clustering](#bundle_part): the coupling strength, the bundles, their links and the expansion over them.
2. [Adaptive refinement](#adaptive_part): a round, when it stops, the cost and limits, and regroup.

The words used on this page:

| Word | Meaning |
|---|---|
| site | One bath spin, or one defect with its on-site spins. A site cannot be split. |
| bundle | A set of sites that is treated exactly inside every cluster it is in. |
| cluster | The union of up to `order` connected bundles, the unit whose Hamiltonian is built and propagated. |
| partition | The assignment of every site to exactly one bundle. |

With every bundle one site, bundle clustering is standard CCE ([CCE](cce.md)). [pCCE](pcce.md) has the same layer with
fixed-size groups made from the positions, while bundles have variable size and are made from the couplings.

Bundle clustering runs with cCCE and gCCE, `quantity` coherence, ensemble and sampled states, on the CPU and the GPU
(`execution.backend`). `cluster.selection.max_per_order` must stay all 0 and `cluster.pcce.cache` off.

<a id="bundle_part"></a>
## 1. Bundle clustering

<a id="bundle_coupling"></a>
### 1.1 The coupling strength

For two sites $a$ and $b$ the coupling strength is a bound on the size of their point-dipole interaction,

$$w_{ab} = \sum_{i \in a,\; j \in b} 3\, I_i I_j\, \lVert \mathsf{J}_{ij} \rVert_F, \qquad \lVert \mathsf{J}_{ij}\rVert_F = \sqrt{6}\,\frac{\mu_0}{4\pi}\frac{\hbar\,\lvert\gamma_i\gamma_j\rvert}{r_{ij}^3},$$

summed over every spin of the two sites (a defect's on-site spins included), in rad/ms. It does not depend on the field
direction. It is not an estimate of the coherence error, and it only orders the merges and the links. Between two
bundles the strength is the sum over their site pairs.

<a id="bundle_partition"></a>
### 1.2 Making the bundles

1. Every site starts as its own bundle, whose dimension is the product of its spins' $2I+1$. A site larger than
   [max_dim](../user/InputFile_Cluster.md#tag_bundle_max_dim) is refused.
2. The candidate site pairs are taken in order of decreasing $w_{ab}$ (ties by site index). With
   `connections` distance, a pair is a candidate only when some spin of the one site is between the inner and the outer
   `rdip` of some spin of the other.
3. A pair with $w_{ab}$ > 0 and at least [merge_cutoff](../user/InputFile_Cluster.md#tag_bundle_merge_cutoff) joins the
   two bundles that hold its sites, if the product of their dimensions is at most `max_dim`. Otherwise the pair is
   counted as blocked. `merge_cutoff` 0 (the default) lets every candidate pair merge as long as `max_dim` allows, so
   `max_dim` alone sets how large the bundles grow. A large `merge_cutoff` keeps the sites apart.
4. Two bundles are *linked*, and can be in one cluster, when a candidate pair joins them. The cluster dimension is the
   product of its bundles' dimensions, times the qubit dimension with gCCE. It is reported, not capped.

The log line `[adaptive] ... sites= bundles= pairs= merges= blocked_bundle=` says what was made, and one
`[adaptive-bundle]` line per bundle gives its dimension and the 0-based bath indices of its sites. The partition does
not use the random numbers of the bath states and does not depend on the number of ranks or the backend.

<a id="bundle_connections"></a>
### 1.3 Which bundles are linked

[connections](../user/InputFile_Cluster.md#tag_bundle_connections) chooses the links:

| `connections` | Needs | The links |
|---|---|---|
| `distance` | `rdip` | The links of step 4, by distance. |
| `auto` | [connection_budget](../user/InputFile_Cluster.md#tag_bundle_connection_budget) | The weakest links are left out while their coupling score stays within the budget. |
| `checked` | `rdip` and `connections` both left out | All links are kept at first, and a coherence check on a few domains decides which weak ones may go. |

<a id="bundle_auto"></a>
**auto.** Every site pair is evaluated. With `rdip` given, pairs outside it can only be left out. After the merges,
links are left out weakest first while the score of what is left out, $t_{\max} \sum_{\text{left out}} w_{AB}$ with
$t_{\max}$ = `deltat` × (`nstep` − 1), stays at most the budget. The score is in rad (rad/ms × ms). Links that must be
left out, those outside `rdip` and every link at `order` 1, count first, and a run whose required omissions already
exceed the budget is refused. The budget is a coupling score, not an error bound. The log line is
`[adaptive-connections]`.

<a id="bundle_checked"></a>
**checked.**

1. Candidates. The bundle pairs are taken in order of increasing strength while $t_{\max}\sum w_{AB}$ stays below
   $\sqrt{10^{-7} + 10^{-5}}$, at most 8 links.
2. Domains. Each candidate gets a domain, its two bundles and the most strongly coupled third bundle that fits. When all
   the bundles together fit (the dimension stays within the integer range), the whole bath is the one domain.
3. Comparison. Each domain and every cluster of its local expansion, with and without the candidate links, is computed
   directly for every bath state and every time point, within 64 such curves per bath state and 64 MiB for the stored
   states and curves. The errors are scaled as $|L_a - L_b| / (10^{-7} + 10^{-5}\,|L_b|)$ for three comparisons. These
   are the expansion before the omission against the direct domain, after against the direct domain, and after against
   before.
4. Decision. All candidates are left out only if every scaled error is at most 1. Otherwise all are kept, and a budget
   that is too small also keeps every link.

The result and every domain's errors are written to `<output file>.adaptive.json` (with `"accuracy_certified":
false`). The line `[adaptive-check] status=...` gives the outcome, one of `local_checks_passed`,
`pruning_failed_links_restored`, `base_unconverged_links_restored`, `probe_budget_unverified_links_kept`,
`memory_budget_unverified_links_kept`, `single_bundle_no_links` and `order_one_no_link_expansion_unchecked`.

The check covers the selected domains only, not the convergence of the whole bath. It costs one extra pass over the
bath states with the domain curves.

<a id="bundle_expansion"></a>
### 1.4 The expansion over bundles

The clusters are the connected sets of up to `order` linked bundles (paths and stars as well as cliques), closed under
subsets. Each gets the integer power of the [inclusion-exclusion](cce.md#cce_expansion) product, and the sets with a
non-zero power are computed. `order` counts bundles.

Two limits stop a run before any dynamics when the expansion is too large. They never truncate it.

- [max_clusters](../user/InputFile_Cluster.md#tag_bundle_max_clusters) limits the number of bundle sets listed, single
  bundles included.
- [max_work](../user/InputFile_Cluster.md#tag_bundle_max_work) limits the number of elementary steps of the
  enumeration, one per link read, per neighbour tried while growing a set, per subset visited and per power update.

The line `[adaptive-expansion] ... terms=` gives how many clusters are propagated.

<a id="adaptive_part"></a>
## 2. Adaptive refinement

<a id="adaptive_round"></a>
### 2.1 A round of the merge policy

With the `cluster.adaptive` block the run is a sequence of *steps*, each a full CCEX calculation in the same process,
in the folder [workdir](../user/InputFile_Cluster.md#tag_adaptive_keys). A round goes through:

1. Run. The current partition is run (step folder `NNNN-run`), giving the coherence of every bath state.
2. Patches. A patch is a small set of neighbouring bundles where the error of the expansion is measured.
   - At `order` 2 it is a bundle with two of its [neighbors](../user/InputFile_Cluster.md#tag_adaptive_keys) most
     strongly linked bundles (a triple). At `order` 1 it is a pair, and above `order` 2 a bounded beam of larger sets.
   - With at most 8 bundles in all and the work affordable, the whole bath is also a patch, an exact reference for the
     whole model.
   - At most `max_patches` are checked, the most at risk first (bundles whose curve is small and that have many links),
     then the most strongly coupled.
3. Probe. Every non-empty subset of the bundles of a patch is computed directly, as one cluster, for every bath state
   and time point (step folders `NNNN-probe`). Probe curves are cached (`cache_mb`, `raw_cache`) and limited in number
   (`max_probes`) and in memory per node ([probe_memory_gb](../user/InputFile_Cluster.md#tag_adaptive_probe_memory_gb)).
4. Local error. For each patch, the expansion restricted to it (its bundles and links, to `order`) is compared with
   the direct curve of the whole patch, at every bath state and every time point, scaled as
   $|L_{\rm exp} - L_{\rm direct}| / (\text{atol} + \text{rtol}\,|L_{\rm direct}|)$. A patch passes when the largest
   value is at most 1. This is the error at the bundle boundaries, what the links between the bundles of a patch miss
   against treating the patch exactly.
5. Merge. Every pair of bundles in a patch whose merged dimension is at most `max_dim` is a candidate. Its gain is
   how much the merge lowers the patch's local error, and its score is the gain divided by $D^3$, with $D$ the merged
   dimension. With `merge_priority` worst-error the local error comes first. Up to `max_candidates` candidates are
   tried in that order, and a `plan` step for each checks that CCEX can make the partition. The first that works is
   adopted, and the next round runs it.

`prefer_physical` prefers candidates whose curve stays finite and bounded, and `hub_candidates` (with it) also tries
growing a bundle around a strongly linked hub first.

<a id="adaptive_stop"></a>
### 2.2 When it stops

| Status (`summary.json`) | When | Exit status |
|---|---|---|
| `validated_full_bath_model` | The whole bath is a patch, the expansion matches its direct curve (scaled error ≤ 1), and every curve is finite (cCCE: $\lvert L\rvert \le$ 1 + atol). | 0 |
| `passed_convergence_checks` | Every checked patch passes, the curves changed by at most the tolerance since the previous round, every curve is finite and bounded as above, and no patch was left out by the limits. | 0 |
| `round_limit` | `max_rounds` rounds passed without passing the checks. | 2 |
| `no_feasible_improving_merge` | No candidate improves a patch, or none can be made. | 2 |
| `probe_budget_limited` | The probe curves needed would exceed `max_probes` or the memory of `cache_mb`. | 2 |
| `failed` | An error during a step. The summary keeps what was done and the error. | 1 |

Whatever the status, the result is the last run. `summary.json` names its folder (`final_directory`), whose
`coherence_wiDiv` is the curve, one file per state with sampled states (`_state<s>`). The first run is in `0001-run`.
[Output files](../user/OutputFiles.md#out_adaptive) describes the folder. Besides the status, `summary.json` keeps one
record per round, with the bundles, the patches and their errors, the change since the previous round and the accepted
merge.

`main.out` exits with the status above, so a script should not treat 2 as a crash. With `configuration.k_range` the
next configuration is not started after a non-zero status.
`example/tutorial/advanced/11_cluster_shapes/cCCE2_Hahn_adaptive/run.sh` shows how to read it.

<a id="adaptive_outputs"></a>
<a id="adaptive_cost"></a>
### 2.3 Cost and limits

Every round computes the whole bath again, plus the probes of its patches. A patch whose missing probes would need
more than `max_total_work` (in units of $D^3$ × states × time points) is not checked and is counted in the round's
coverage. The refinement needs an explicit `configuration.seed`, since every round must see the same bath states.

Limits of the method:

- The local errors and the change between rounds are checks, not error bounds of the whole bath.
- A greedy merge cannot undo a poor first partition.
- One partition is used for all time points and states.

Check its effect on your system against a calculation without it.

<a id="adaptive_regroup"></a>
### 2.4 regroup

`"regroup": true` replaces the merge search by a search aimed at curves with $|L| > 1$. It needs cCCE, `order` 2 and
`connections` distance with `rdip`. It is refused otherwise, and with `prefer_physical` or `hub_candidates`.

- The run stops as soon as every state and time point is finite with $|L| \le 1.01$ (status `physical`, exit 0).
- Otherwise, at the worst point it measures $g = \log\lvert L_{ij}/(L_i L_j)\rvert$ for every linked bundle pair.
- Around the pairs with $g > 0$ it looks for one change of the bundles (a merge, a split, moving or exchanging a site)
  that lowers the local error without raising that of the neighbouring bundles or the number of unphysical states.
- An accepted change is checked by running it, and undone if the check fails.

The other statuses are `no_positive_culprit`, `no_candidate`, `budget_limited`, `change_limit` (exit 2) and `failed`
(exit 1). `physical` means bounded, not accurate.

`example/tutorial/advanced/11_cluster_shapes/cCCE2_Hahn_regroup/` runs it. `regroup_factors` and `regroup_whole`
choose how the culprit values and the whole curves are obtained. The second needs a build whose long double is the x87
80-bit type.

Bundle clustering and the adaptive refinement were developed for CCEX and have not been published yet.
