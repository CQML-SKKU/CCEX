# Output files
<a id="OutputFiles"></a>

Every file a run writes, its name and its columns. The keys that choose them are in [Output tags](InputFile_Output.md);
the folder is `-o DIR`, else [output.path](InputFile_Output.md#tag_output_path), else the working directory, and the
name `<file>` below is [output.file](InputFile_Output.md#tag_outfile) (by default the input file's name without `.json`
or `.jsonc`).

| File | Written when | Section |
|---|---|---|
| `<file>_wiDiv`, `<file>_noDiv` | every calculation | [The coherence files](#out_coherence) |
| `<file>_state<s>_wiDiv`, `_noDiv` | `configuration.s_range`, one pair per sampled state | [Names](#out_names) |
| `<file>_k<k>_...` | `configuration.k_range`, one set per configuration | [Names](#out_names) |
| the same files with blocks of rows | `quantity` `dm` | [The density matrix](#out_dm) |
| three more columns per row | `output.savemode` `coherence_tagged` | [Tagged rows](#out_tagged) |
| `<file>_<cluster>_wiDiv`, `<file>_<cluster>_full_...` | `output.savemode` `cluster_coherence`, `cluster_density_matrix` | [One file per cluster](#out_clusters) |
| `<file>.clusters.json` | `output.cluster_only` | [The cluster plan](#out_plan) |
| `<file>.contributions.state-<s>.tsv`, `<file>.contributions.sites.tsv` | `output.cluster_contributions` | [Cluster contributions](#out_contributions) |
| `<file>_rdip<r>_...` | `rdip_ladder` (written by the automatic convergence into its runs) | [Names](#out_names) |
| the `cluster.adaptive.workdir` folder | a `cluster.adaptive` block | [Adaptive refinement](#out_adaptive) |
| the `auto_convergence.workdir` folder | an `auto_convergence` block | [Automatic convergence](#out_autoconv) |
| `<-fo>` files of the analyzer | `example/analysis/Analyzer` | [Analyzer outputs](#out_analyzer) |

The run log goes to standard output (and the refusals to standard error); [The run log](GettingStarted.md#RunLog)
describes it.

<a id="out_coherence"></a>
## The coherence files: `_wiDiv` and `_noDiv`

One line per time step, `nstep` lines: the time in ms (`%15g`), a tab, the complex value as `+re+imj` with ten decimals
(real part `%+15.10f` followed by the imaginary part and `j`), and a tab. The first line is t = 0; the step is `deltat`.

```
              0	  +1.0000000000+0.0000000000j	
          0.001	  +0.9243627460+0.0000003519j	
          0.002	  +0.8232539603+0.0000033616j	
```

(the first lines of `example/tutorial/00_quickstart`'s `ccein_wiDiv`). A value that is not a number is written `nan`.

- **`_wiDiv`** is the CCE result: the product of the irreducible cluster contributions, each cluster's coherence divided
  by those of its sub-clusters,

  $$L(t) = \prod_{C} \tilde L_C(t), \qquad \tilde L_C = \frac{L_C}{\prod_{C' \subset C} \tilde L_{C'}} .$$

  This is the coherence to read and plot. Where a sub-cluster's coherence passes near zero the division can make it
  large; |L| > 1 is a sign of that (a divergence), not a physical value.
- **`_noDiv`** is the product of the cluster coherences with positive powers only, without the division. It does not
  diverge but it is **not** a coherence: the analyzer's error correction compares `_wiDiv` with it to find the divergent
  points ([Analysis](Analysis.md)). `output.savemode` `{"nodiv": false}` leaves it out.

Read them with `read_coherence` of `example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py` (time and a complex
array), or in any language: split each line at the tab and parse the second field as a complex number (Python's
`complex()` reads `+re+imj` after the spaces are removed).

<a id="out_names"></a>
## Names

The name is built in this order, each part only when it applies:

```
<file>[_k<k>][_rdip<r>][_<cluster>][_full][_state<s>]_wiDiv     (and _noDiv)
```

| Part | When |
|---|---|
| `_k<k>` | with [configuration.k_range](InputFile_General.md#tag_configuration_k_range): one set of outputs per configuration k |
| `_rdip<r>` | the extra rdip values of `rdip_ladder` (an internal key of the [automatic convergence](InputFile_AutoConvergence.md#tag_ac_rdip_ladder)) |
| `_<cluster>` | `output.savemode` `cluster_coherence` or `cluster_density_matrix`: the bath indices of the cluster joined by `_` |
| `_full` | `cluster_density_matrix`: the density matrix of the whole qubit + cluster space |
| `_state<s>` | [configuration.s_range](InputFile_General.md#tag_s_range): the sampled state s = first ... last; absent for the ensemble (no `s_range`) |

So the ensemble writes `ccein_wiDiv`, `s_range` [1, 20] writes `ccein_state1_wiDiv` ... `ccein_state20_wiDiv` (and [11, 20]
the last ten of them), and `k_range`
[1, 4] in the ensemble writes `ccein_k1_wiDiv` ... `ccein_k4_wiDiv`. The states are not averaged by CCEX: average them
(and the configurations) as complex numbers afterwards ([Analysis](Analysis.md)).

<a id="out_dm"></a>
## The density matrix (`quantity` dm)

With [quantity](InputFile_General.md#tag_quantity) `dm` (gCCE) the same files hold, per time step, the reduced density
matrix of the qubits, D × D with D the qubit dimension (3 for one NV, 3ⁿ for n spin-1 qubits). One block per time step:
its first line starts with the time, the elements follow row by row (ρ₀₀, ρ₀₁, ..., ρ₀,D−1, ρ₁₀, ...), three to a line,
and each continuation line starts with an empty time field (15 spaces) and a tab. The basis is the $S_z$ basis of each
qubit, m = +S first down to −S (for an NV: ms +1, 0, −1), and for several qubits their product with the first qubit the
slowest index (kron(q₀, q₁, ...), the last qubit varies fastest).

The cluster expansion is applied element by element. An element the free evolution keeps at zero is left at zero; when
clusters would put a value into it the run prints a `[dm]` line with the largest such value.

`example/tutorial/plot_coherence.py` reads these blocks (its `read` function); the analyzer's `read_coherence` reads only
the one-value files.

<a id="out_tagged"></a>
## Tagged rows (`coherence_tagged`)

The names are those of `coherence`. Every time row ends, after the coherence, with three more tab-separated fields: the
zx and the zz component of the first bath spin's hyperfine tensor (in rad/ms, the unit CCEX holds the tensors in, ten
decimals) and the name of the first bath file. They tell the results of many bath files apart after they are collected.

<a id="out_clusters"></a>
## One file per cluster (`cluster_coherence`, `cluster_density_matrix`)

`cluster_coherence` writes the result as `coherence` does and, for every cluster and state, that cluster's own coherence
$L_C$ (not the divided $\tilde L_C$): `<file>_<cluster>[_state<s>]_wiDiv` and `_noDiv`, where `<cluster>` is the cluster's bath
indices joined by `_` (`ccein_curves_519_557_wiDiv`); the two files of one cluster hold the same curve.
`cluster_density_matrix` writes all of that and, with gCCE, for every cluster the density matrix of the whole qubit +
cluster space at every time step, as blocks of rows like [the density matrix](#out_dm) with the full dimension, to
`<file>_<cluster>_full[_state<s>]_wiDiv` (`<file>_full_...` for the qubits alone). cCCE has no qubit + bath density
matrix: there it writes what `cluster_coherence` writes. Both are for small systems: one file per cluster.

<a id="out_plan"></a>
## The cluster plan (`<file>.clusters.json`)

[output.cluster_only](InputFile_Output.md#tag_cluster_only) writes the clusters and stops before the dynamics; the log
prints one `[cluster-plan]` line with the number of sites and clusters and the path. The file is one JSON object:

| Key | Value |
|---|---|
| `schema` | 1 |
| `dynamics_performed` | `false` |
| `method`, `cluster_mode` | as the run read them |
| `input_order` | `order` of the input |
| `max_site_order` | the largest number of sites in a cluster |
| `loaded_sites`, `retained_sites` | the bath spins read, and those kept after the bath was cut (pCCE drops the rest of a partition) |
| `qubit_dimension` | the qubit dimension in a cluster (1 for cCCE) |
| `zero_power` | the power of the zeroth cluster (the qubit alone) |
| `sites` | one object per site: `id` (its index), `name`, `xyz` (Å), `spin`, `component_dimensions` (its own dimension and those of its on-site spins), `defect_axis` when the bath has Defect species |
| `bundles` | the groups of sites that enter a cluster whole (one site each in `standard` mode); `null` for pCCE |
| `clusters` | one object per cluster: `power` (its power in the product) and `sites` (site indices) |

<a id="out_contributions"></a>
## Cluster contributions (`.contributions.state-<s>.tsv`)

[output.cluster_contributions](InputFile_Output.md#tag_cluster_contributions) (cCCE and gCCE, `quantity` coherence)
writes, for every state, `<file>.contributions.state-<s>.tsv` (s = 0 for the ensemble, first ... last of `s_range` for sampled states)
and once `<file>.contributions.sites.tsv`. The table's first line is a header comment beginning
`# schema=ccex-contributions-v1`, with the number of time points (`nstep=`), the time step in ms (`deltat_ms=`), the
method, and `normalization=` (`factor_t0` for gCCE, `none` otherwise); then one line of column names and one row per
cluster, tab-separated:

| Column | Meaning |
|---|---|
| `sites` | the bath indices of the cluster, comma-separated |
| `bundles` | the bundles it is made of (with `cluster.mode` standard, each site is its own bundle) |
| `production_power` | its power in the product of the result; 0 for a sub-cluster added only for this table |
| `dimension` | the Hilbert dimension of the cluster (with the qubit for gCCE) |
| `max_abs_log_tilde` | the largest \|log L̃_C(t)\| over the time points: the logarithm of the cluster's coherence relative to the qubit alone, minus those of all its sub-clusters (the irreducible contribution), with the imaginary part folded into (−π, π]; for gCCE each log is first referred to its value at t = 0. `nan` when any point is undefined |
| `undefined_points` | the time points where it is not finite (a zero or an overflow) |
| `worst_step` | the time index (0 ... nstep − 1) of the largest value |

The largest values are the clusters that dominate the decay (or diverge): `sort -t "$(printf '\t')" -k5 -g -r` on the table lists
them first. `.contributions.sites.tsv` has one row per site: `site`, `name`, `spin`, `gyro` (rad/ms/G), `x_A`, `y_A`,
`z_A` (Å). The log prints `[cluster-contributions]` lines with the rows added and the time it took.

<a id="out_adaptive"></a>
## Adaptive refinement (`cluster.adaptive.workdir`)

With a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block every step of the refinement is an ordinary run in
its own folder under `workdir`, numbered in order with its kind: `0001-run`, `0002-probe`, ... Each holds `request.json`
(what the step was asked), `input.json` (the run's input), `run.log`, the step's coherence `coherence_wiDiv` (and
`_noDiv`, `_state<s>_` with sampled states) and its response files. In `output` itself:

| File | Content |
|---|---|
| `summary.json` | the report, rewritten after every step |
| `final-bundles.json` | the bundles of the final step |
| `final-input.json`, `final-request.json` | an input that reruns the final partition as one calculation |

`summary.json` holds `status`, `input` (the input as read), `settings`, `binary_sha256`, `input_file_sha256`, `runs`
(one record per step: `directory`, `mode`, `wall_s`, `returncode`, `timeout`), `rounds` (one record per round of
the policy), `limitations`, and at the end `final_directory`, `final_bundles`, `cached_probe_clusters`,
`raw_nonfinite_values`, `raw_zero_values` (and `raw_cache` when a probe cache file is used). **The final curve is
`<final_directory>/coherence_wiDiv`.**

The status and the exit code of `main.out`:

| Status | Exit code |
|---|---|
| `validated_full_bath_model`, `passed_convergence_checks` (merge), `physical` (regroup) | 0 |
| `round_limit`, `no_feasible_improving_merge`, `probe_budget_limited` (merge); `no_candidate`, `no_positive_culprit`, `budget_limited`, `change_limit` (regroup) | 2 |
| an error during the refinement | 1 |

A status with exit code 2 is not a failure of the run: the refinement stopped there and its last step is the result
(`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_adaptive` ends with `no_feasible_improving_merge`). The run prints the
status and the path of `summary.json` as its last line before the wall time.

<a id="out_autoconv"></a>
## Automatic convergence (`auto_convergence.workdir`)

With an [auto_convergence](InputFile_AutoConvergence.md) block everything goes under `workdir`:

| Path | Content |
|---|---|
| `runs/<step>/<k>-<key>/` | one folder per CCEX run (step, configuration k, the first characters of a hash of what the run depends on): its `input.json`, the outputs `out_wiDiv` / `out_state<s>_wiDiv` (and `_noDiv`), `done.json` once it finished |
| `decisions.jsonl` | one line per decision of the current plan (earlier plans: `decisions.<n>.jsonl`) |
| `split.json` | the test and production configurations, when they were found from the files |
| `pcce_cache/` | the partition cache files of pCCE runs |
| `result/summary.txt`, `result/json/summary.json` | the status, the values chosen, the production (A with its SE, T2, n), the core-hours and every warning |
| `result/progress.txt`, `result/json/progress.json` | where the run is, rewritten at every step |
| `result/coherence_mean.txt` | the production mean: columns t (ms), Re, Im, \|L\| |
| `result/1_window.png` ... `result/6_converged_coherence.png` | the figures of the steps and of the converged coherence |
| `lock` | held while a runner works in the folder |

Its exit codes (0 converged, 3 unresolved, 2 refused, 4 stopped, 1 error) and how to continue a cut job are in
[Results and exit codes](InputFile_AutoConvergence.md#tag_ac_results).

<a id="out_analyzer"></a>
## Analyzer outputs

The analyzer (`example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py`, [Analysis](Analysis.md)) writes files named by
its `-fo` option, in the same two-column format (time, then `+re+imj`, ten decimals):

| Option | File |
|---|---|
| `-err` | the error-corrected curve per input file: the `-fo` name with the variable names replaced by their values (`-fo results/ccein_kCONF_ErrCorr` gives `results/ccein_k1_ErrCorr`, ...) |
| `-ea` | the average over all files read: `<-fo>_conv<N>`, N the number of curves averaged (`results/ccein_EnsAvg_conv4`) |
| `-fit_w` | `<-fo>_params`: one row per variable value, columns value, T2 (in the `-uc` time unit), stretching power n (with `-ea`, one row for the average, value N) |
| `-pl`, `-pl_m`, `-pl_i`, `-pl_T`, `-pl_p` | figures, `<-fo>.png` (and `_T2.png`, `_Power.png`) |
