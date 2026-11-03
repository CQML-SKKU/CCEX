# Bundle clustering and adaptive refinement

`cluster.mode` bundle (keys: [cluster.bundle](../user/InputFile_Cluster.md#tag_bundle)) and the
[cluster.adaptive](../user/InputFile_Cluster.md#tag_adaptive) block on top of it. Bundle clustering puts strongly
coupled bath spins into *bundles* that always enter a cluster whole, once, from the couplings; the adaptive refinement
repeats the calculation, measures the error of the expansion locally and changes the bundles.

The words: a **site** is one bath spin, or one defect with its on-site spins (it cannot be split); a **bundle** is a
set of sites that is treated exactly inside every cluster it is in; a **cluster** is the union of up to `order`
connected bundles, the unit whose Hamiltonian is built and propagated; a **partition** is the assignment of every site
to exactly one bundle. With every bundle one site, bundle clustering is standard CCE ([CCE](cce.md)). pCCE
([pCCE](pcce.md)) has the same layer with fixed-size groups made from the positions; bundles have variable size and are
made from the couplings.

Bundle clustering runs with cCCE and gCCE, `quantity` coherence, ensemble and sampled states, on the CPU and the GPU
(`execution.backend`). `cluster.selection.max_per_order` must stay all 0 and `cluster.pcce.cache` off.

<a id="bundle_coupling"></a>
## The coupling strength

For two sites $a$ and $b$ the coupling strength is a bound on the size of their point-dipole interaction,

$$w_{ab} = \sum_{i \in a,\; j \in b} 3\, I_i I_j\, \lVert \mathsf{J}_{ij} \rVert_F, \qquad \lVert \mathsf{J}_{ij}\rVert_F = \sqrt{6}\,\frac{\mu_0}{4\pi}\frac{\hbar\,\lvert\gamma_i\gamma_j\rvert}{r_{ij}^3},$$

summed over every spin of the two sites (a defect's on-site spins included), in rad/ms. It does not depend on the field
direction and is not an estimate of the coherence error; it only orders the merges and the connections. Between two
bundles the strength is the sum over their site pairs.

<a id="bundle_partition"></a>
## Making the bundles

1. Every site starts as its own bundle; its dimension is the product of its spins' $2I+1$. A site larger than
   [max_dim](../user/InputFile_Cluster.md#tag_bundle_max_dim) is refused.
2. The candidate site pairs are taken in order of decreasing $w_{ab}$ (ties by site index). With
   `connections` distance, a pair is a candidate only when some spin of the one site is between the inner and the outer
   `rdip` of some spin of the other.
3. A pair with $w_{ab}$ > 0 and at least [merge_cutoff](../user/InputFile_Cluster.md#tag_bundle_merge_cutoff) joins the
   two bundles that hold its sites, if the product of their dimensions is at most `max_dim`; otherwise the pair is
   counted as blocked. `merge_cutoff` 0 (the default) lets every candidate pair merge as long as `max_dim` allows, so
   `max_dim` alone sets how large the bundles grow; a large `merge_cutoff` keeps the sites apart.
4. Two bundles are *linked* (can be in one cluster) when a candidate pair joins them. The cluster dimension is the
   product of its bundles' dimensions (gCCE: times the qubit dimension); it is reported, not capped.

The log line `[adaptive] ... sites= bundles= pairs= merges= blocked_bundle=` and one `[adaptive-bundle]` line per bundle
(its dimension and the 0-based bath indices of its sites) say what was made. The partition does not use the random
numbers of the bath states and does not depend on the number of ranks or the backend.

<a id="bundle_connections"></a>
## Which bundles are linked

[connections](../user/InputFile_Cluster.md#tag_bundle_connections):

- `distance` (with `rdip`): the links of step 4, by distance (the log names it `coupling-bound-v1`).
- `auto` with [connection_budget](../user/InputFile_Cluster.md#tag_bundle_connection_budget): every site pair is
  evaluated (with `rdip` given, pairs outside it can only be left out). After the merges, links are left out weakest
  first while the score of what is left out,
  $t_{\max} \sum_{\text{left out}} w_{AB}$ with $t_{\max}$ = `deltat` × (`nstep` − 1), stays at most the budget; the
  score is in rad (rad/ms × ms). Links that must be left out (outside `rdip`; every link at `order` 1) count first, and
  a run whose required omissions already exceed the budget is refused. The budget is a coupling score, not an error
  bound (log name `coupling-budget-v1`, line `[adaptive-connections]`).
- `checked`: with `rdip` and `connections` both left out, the links are first all kept, and a bounded coherence check
  decides which weak ones may go (log name `coherence-checked-links-v1`):
  1. Candidates: the bundle pairs in order of increasing strength, while $t_{\max}\sum w_{AB}$ stays below
     $\sqrt{10^{-7} + 10^{-5}}$; at most 8 links.
  2. For each candidate a *domain*: its two bundles and the most strongly coupled third bundle that fits; when all the
     bundles together fit (the dimension stays within the integer range), the whole bath is the one domain.
  3. Each domain and every cluster of its local expansion, with and without the candidate links, is computed directly
     (at most 64 such curves per bath state, and 64 MiB for the stored states and curves), for every bath state and
     every time point. The errors are scaled as $|L_a - L_b| / (10^{-7} + 10^{-5}\,|L_b|)$: the expansion before the
     omission against the direct domain, after against the direct domain, and after against before.
  4. All candidates are left out only if every scaled error is at most 1; otherwise all are kept. A budget that is too
     small keeps every link.

  The result and every domain's errors are written to `<output file>.adaptive.json` (with `"accuracy_certified":
  false`), and the line `[adaptive-check] status=...` gives the outcome: `local_checks_passed`,
  `pruning_failed_links_restored`, `base_unconverged_links_restored`, `probe_budget_unverified_links_kept`,
  `memory_budget_unverified_links_kept`, `single_bundle_no_links`, `order_one_no_link_expansion_unchecked`. The check
  covers the selected domains only, not the convergence of the whole bath; it costs one extra pass over the bath states
  with the domain curves.

<a id="bundle_expansion"></a>
## The expansion over bundles

The clusters are the connected sets of up to `order` linked bundles (paths and stars as well as cliques), closed under
subsets; each gets the integer power of the [inclusion-exclusion](cce.md#cce_expansion) product, and the sets with a
non-zero power are computed. `order` counts bundles. Two limits stop a run before any dynamics when the expansion is
too large, never by truncating it: [max_clusters](../user/InputFile_Cluster.md#tag_bundle_max_clusters), the number
of bundle sets listed (single bundles included), and [max_work](../user/InputFile_Cluster.md#tag_bundle_max_work), the
number of elementary steps of the enumeration (one per link read, per neighbour tried while growing a set, per subset
visited and per power update). The line `[adaptive-expansion] ... terms=` gives how many clusters are propagated.

<a id="adaptive_round"></a>
## Adaptive refinement (the merge policy)

With the `cluster.adaptive` block the run is not one calculation but a sequence of *steps*, each a full CCEX
calculation in the same process, in the folder [workdir](../user/InputFile_Cluster.md#tag_adaptive_keys):

1. **Run** the current partition (step folder `NNNN-run`): the coherence of every bath state.
2. **Patches.** A patch is a small set of neighbouring bundles where the error of the expansion is measured. At
   `order` 2 it is a bundle with two of its [neighbors](../user/InputFile_Cluster.md#tag_adaptive_keys) most strongly
   linked bundles (a triple); at `order` 1 a pair; above `order` 2 a bounded beam of larger sets. With at most 8 bundles
   in all and the work affordable, the whole bath is also a patch (an exact reference for the whole model). At most
   `max_patches` are checked, the most at risk first (bundles whose curve is small and that have many links), then the
   most strongly coupled.
3. **Probe** (step folders `NNNN-probe`): every non-empty subset of the bundles of a patch is computed directly, as one
   cluster, for every bath state and time point. Probe curves are cached (`cache_mb`, `raw_cache`) and limited in number
   (`max_probes`) and in memory per node ([probe_memory_gb](../user/InputFile_Cluster.md#tag_adaptive_probe_memory_gb)).
4. **Local error.** For each patch, the expansion restricted to it (its bundles and links, to `order`) is compared with
   the direct curve of the whole patch, at every bath state and every time point, scaled as
   $|L_{\rm exp} - L_{\rm direct}| / (\text{atol} + \text{rtol}\,|L_{\rm direct}|)$; a patch passes when the largest
   value is at most 1. This is the "error at the bundle boundaries": what the links between the bundles of a patch miss
   against treating the patch exactly.
5. **Merge.** Every pair of bundles in a patch whose merged dimension is at most `max_dim` is a candidate; its gain is
   how much the merge lowers the patch's local error, and its score the gain divided by $D^3$ ($D$ the merged
   dimension), or the local error first with `merge_priority` worst-error. Up to `max_candidates` candidates are tried
   in that order (a `plan` step each checks that CCEX can make the partition); the first that works is adopted, and
   the next round runs it.

<a id="adaptive_stop"></a>
### When it stops

| Status (`summary.json`) | When | Exit status |
|---|---|---|
| `validated_full_bath_model` | the whole bath is a patch, the expansion matches its direct curve (scaled error ≤ 1) and every curve is finite (cCCE: $\lvert L\rvert \le$ 1 + atol) | 0 |
| `passed_convergence_checks` | every checked patch passes, the curves changed by at most the tolerance since the previous round, every curve is finite and bounded as above, and no patch was left out by the limits | 0 |
| `round_limit` | `max_rounds` rounds without passing | 2 |
| `no_feasible_improving_merge` | no candidate improves a patch, or none can be made | 2 |
| `probe_budget_limited` | the probe curves needed would exceed `max_probes` or the memory of `cache_mb` | 2 |
| `failed` | an error during a step (the summary keeps what was done and the error) | 1 |

Whatever the status, the result is the **last run**: `summary.json` names its folder (`final_directory`), whose
`coherence_wiDiv` (with sampled states one file per state, `_state<s>`) is the curve; the first run is in `0001-run`.
`main.out` exits with the status above, so a script should not treat 2 as a crash; with `configuration.k_range` the
next configuration is not started after a non-zero status
(`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_adaptive/run.sh` shows how to read it).

<a id="adaptive_outputs"></a>
### What the output folder holds

- one folder per step, `0001-run`, `0002-probe`, `0003-plan`, ..., each with the step's `input.json` (the user's input,
  its refinement options replaced by the step's request), `request.json`, `run.log` and the step's results;
- `summary.json`: `status`, the input and the settings, one record per step (`runs`: folder, mode, wall time, return
  code), one record per round (`rounds`: the number of bundles, the patches with their scaled and absolute errors and
  the worst state and time step, the largest local error, whether the curves are finite and bounded, the change since
  the previous round, the coverage of the patches, the accepted merge with its local error before and after),
  `final_directory`, `final_bundles` (the sites of each bundle), the number of cached probe curves and of non-finite or
  zero probe values, the checksums of the input files, and `limitations`;
- `final-bundles.json`, and `final-input.json` with `final-request.json` to run the chosen partition again on its own.

<a id="adaptive_cost"></a>
### Cost and limits

Every round computes the whole bath again, plus the probes of its patches; a patch whose missing probes would need
more than `max_total_work` (in units of $D^3$ × states × time points) is not checked and is counted in the round's
coverage. The refinement needs an explicit `configuration.seed`, since
every round must see the same bath states. The local errors and the change between rounds are checks, not error bounds
of the whole bath; a greedy merge cannot undo a poor first partition; and one partition is used for all time points and
states. In the studies made with it so far the refinement changed the curves only a little on the baths tried
(`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_adaptive/` is one such case), so check its effect on your system
against a calculation without it.

Options: `prefer_physical` prefers candidates whose curve stays finite and bounded; `hub_candidates` (with it) also
tries growing a bundle around a strongly linked hub first.

<a id="adaptive_regroup"></a>
### regroup

`"regroup": true` replaces the merge search by a search aimed at curves with $|L| > 1$ (cCCE, `order` 2,
`connections` distance with `rdip`; refused otherwise and with `prefer_physical` or `hub_candidates`). The run stops as
soon as every state and time point is finite with $|L| \le 1.01$ (status `physical`, exit 0). Otherwise, at the worst
point it measures, for every linked bundle pair, $g = \log\lvert L_{ij}/(L_i L_j)\rvert$ and, around the pairs with
$g > 0$, looks for one change of the bundles (a merge, a split, moving or exchanging a site) that lowers the local error
without raising that of the neighbouring bundles or the number of unphysical states; an accepted change is checked by
running it, and undone if the check fails. Other statuses: `no_positive_culprit`, `no_candidate`, `budget_limited`,
`change_limit` (exit 2), `failed` (exit 1). `physical` means bounded, not accurate.
`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_regroup/` runs it. `regroup_factors` and `regroup_whole` choose how
the culprit values and the whole curves are obtained (the second needs a build whose long double is the x87 80-bit
type).

Bundle clustering and the adaptive refinement were developed for CCEX and have not been published yet.
