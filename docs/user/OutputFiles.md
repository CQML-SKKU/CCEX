# Output files
<a id="OutputFiles"></a>

This page lists the files a run writes and what is in them. The keys that choose the files are in
[Output tags](InputFile_Output.md).

Below, `<file>` is the start of every name: [output.file](InputFile_Output.md#tag_outfile), or by default the input
file's name without `.json` or `.jsonc`. The folder is `-o DIR`, else [output.path](InputFile_Output.md#tag_output_path),
else the working directory. The run log is not a file; it goes to standard output
([The run log](BuildingAndRunning.md#RunLog)).

Sections 1 and 2 apply to every run. The other sections cover files that only some keys write.

1. [The coherence files](#out_coherence): `<file>_wiDiv` and `<file>_noDiv`, the result of every run.
2. [Names](#out_names): how the names change with several configurations (`k_range`) or sampled states (`s_range`).
3. [The density matrix](#out_dm): the same files with `quantity` dm, one block of rows per time.
4. [Tagged rows](#out_tagged): three more columns per line, with `output.savemode` `coherence_tagged`.
5. [One file per cluster](#out_clusters): the coherence of every cluster, with `output.savemode` `cluster_coherence`
   or `cluster_density_matrix`.
6. [The cluster plan](#out_plan): the list of clusters, `<file>.clusters.json`, with `output.cluster_only`.
7. [Cluster contributions](#out_contributions): a table of what each cluster contributes, with
   `output.cluster_contributions`.
8. [Adaptive refinement](#out_adaptive): the folder of a `cluster.adaptive` run.
9. [Automatic convergence](#out_autoconv): the folder of an `auto_convergence` run.
10. [Analyzer outputs](#out_analyzer): the files that `example/analysis/Analyzer` writes.

<a id="out_coherence"></a>
## 1. The coherence files: _wiDiv and _noDiv

Each file has one line per time step, `nstep` lines in all. A line holds the time in ms, a tab and the complex value
as `+re+imj` with ten decimals. The first line is t = 0, and the step is `deltat`. The first lines of
`results/ccein_wiDiv` in `example/tutorial/basic/01_first_run` are:

```
              0	  +1.0000000000+0.0000000000j	
          0.001	  +0.9482416411-0.0000013139j	
          0.002	  +0.8296882434-0.0000000881j	
```

A value that is not a number is written `+nan` or `-nan`.

`_wiDiv` is the CCE result, the coherence to read and plot. It is the product of the irreducible contributions of the
clusters. Each contribution is the coherence of a cluster divided by those of its sub-clusters:

$$L(t) = \prod_{C} \tilde L_C(t), \qquad \tilde L_C = \frac{L_C}{\prod_{C' \subset C} \tilde L_{C'}} .$$

Where a sub-cluster's coherence passes near zero, the division can blow up. A value with |L| > 1 is the sign of
this; it is not physical.

`_noDiv` is the same product without the division. It cannot blow up, but it is not a coherence. The analyzer's error
correction compares `_wiDiv` with it to find the points that blew up ([Analysis](Analysis.md)). `output.savemode`
`{"nodiv": false}` leaves this file out.

In Python, `read_coherence` of `example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py` returns the time and a
complex array. Elsewhere, split each line at the tab and read the second field as a complex number (Python's
`complex()` reads `+re+imj` once the spaces are removed).

<a id="out_names"></a>
## 2. Names

A name is built from these parts, in this order; each part appears only when it applies:

```
<file>[_k<k>][_<cluster>][_full][_state<s>]_wiDiv     (and _noDiv)
```

| Part | When |
|---|---|
| `_k<k>` | with [configuration.k_range](InputFile_General.md#tag_configuration_k_range): one set of files per configuration k |
| `_<cluster>` | with `output.savemode` `cluster_coherence` or `cluster_density_matrix`: the bath indices of the cluster, joined by `_` |
| `_full` | with `cluster_density_matrix`: the density matrix of the whole qubit + cluster space |
| `_state<s>` | with [configuration.s_range](InputFile_General.md#tag_s_range): the sampled state s; absent for the ensemble |

For an input `ccein.jsonc`:

- the ensemble writes `ccein_wiDiv`;
- `s_range` [1, 20] writes `ccein_state1_wiDiv` ... `ccein_state20_wiDiv`;
- `k_range` [1, 10] writes `ccein_k1_wiDiv` ... `ccein_k10_wiDiv`;
- both write `ccein_k1_state1_wiDiv` and so on.

CCEX does not average over states or configurations. Average the files afterwards, as complex numbers
([Analysis](Analysis.md)).

<a id="out_dm"></a>
## 3. The density matrix (quantity dm)

With [quantity](InputFile_General.md#tag_quantity) `dm` (gCCE), the same files hold the reduced density matrix of the
qubits at each time step. It is D × D, where D is the qubit dimension: 3 for one NV, $3^n$ for n spin-1 qubits.

Each time step is one block of lines. The first line of a block starts with the time, and the following lines with an
empty time field. The elements follow row by row ($\rho_{00}, \rho_{01}, \ldots$), three to a line. The basis is the
$S_z$ basis of each qubit from m = +S down to −S; for an NV that is $m_s$ = +1, 0, −1. For several qubits the basis is
their product, with the last qubit changing fastest.

The cluster expansion is applied element by element. An element that the free evolution keeps at zero stays zero; if
clusters would put a value into it, the run prints a `[dm]` line with the largest such value.

`example/tutorial/plot_coherence.py` reads these blocks with its `read` function. The analyzer's `read_coherence`
reads only the one-value files.

<a id="out_tagged"></a>
## 4. Tagged rows (coherence_tagged)

The files are those of `coherence`, with three more tab-separated fields at the end of each line: the zx and zz
components of the first bath spin's hyperfine tensor (rad/ms) and the name of the first bath file. They identify the
bath a file came from when the results of many bath files are collected in one place.

<a id="out_clusters"></a>
## 5. One file per cluster (cluster_coherence, cluster_density_matrix)

Both modes write one file per cluster, so they are for small systems.

`cluster_coherence` writes the usual result and, for every cluster and state, the factor that cluster contributes to
the product: $(L_C/L_0)^p$, with $p$ its power, so a sub-cluster with a negative power appears inverted. The files are `<file>_<cluster>[_state<s>]_wiDiv` and `_noDiv`, with the cluster's bath
indices joined by `_`, for example `ccein_curves_0_10_wiDiv`. The two files of one cluster hold the same curve.

`cluster_density_matrix` writes the same files. With gCCE it also writes, for every cluster, the density matrix of the
whole qubit + cluster space at every time step, in blocks as in [the density matrix](#out_dm). These files are
`<file>_<cluster>_full[_state<s>]_wiDiv`, and `<file>_full[_state<s>]_wiDiv` for the qubits alone. cCCE has no qubit + bath density
matrix, so there this mode writes what `cluster_coherence` writes.

<a id="out_plan"></a>
## 6. The cluster plan (<file>.clusters.json)

[output.cluster_only](InputFile_Output.md#tag_cluster_only) writes the list of clusters and stops before the dynamics.
The log prints a `[cluster-plan]` line with the number of sites and clusters and the path of the file. The file is one
JSON object. Its main keys are:

| Key | Value |
|---|---|
| `sites` | one object per bath spin: `id`, `name`, `xyz` (Å), `spin` |
| `bundles` | the groups of sites that enter a cluster together (one site each in `standard` mode; `null` for pCCE) |
| `clusters` | one object per cluster: `sites` (site indices) and `power` (its power in the product) |

The other keys record how the run was set up: the method and cluster mode, the order, the number of sites read and
kept, and the qubit dimension.

<a id="out_contributions"></a>
## 7. Cluster contributions (.contributions.state-<s>.tsv)

[output.cluster_contributions](InputFile_Output.md#tag_cluster_contributions) works with cCCE and gCCE and `quantity`
coherence. It writes a table per state, `<file>.contributions.state-<s>.tsv`, with s = 0 for the ensemble and the
states of `s_range` otherwise. It also writes `<file>.contributions.sites.tsv` once, with one row per site: `site`,
`name`, `spin`, `gyro` (rad/ms/G) and the position `x_A`, `y_A`, `z_A` (Å).

The table starts with a comment line `# schema=ccex-contributions-v1 ...` that gives `nstep`, `deltat_ms` and the
method. Then come the column names and one tab-separated row per cluster:

| Column | Meaning |
|---|---|
| `sites` | the bath indices of the cluster, comma-separated |
| `bundles` | the bundles it is made of (in `standard` mode, each site is its own bundle) |
| `production_power` | its power in the result; 0 for a sub-cluster added only for this table |
| `dimension` | the Hilbert dimension of the cluster (with the qubit for gCCE) |
| `max_abs_log_tilde` | the largest $\lvert\log \tilde L_C(t)\rvert$ over time: how much this cluster's own contribution changes the coherence; `nan` if any point is undefined |
| `undefined_points` | the number of time points where it is not finite |
| `worst_step` | the time index of the largest value |

The clusters with the largest `max_abs_log_tilde` dominate the decay, or blow up. `sort -k5 -g -r` on the table lists
them first (the columns are tab-separated and hold no spaces).

<a id="out_adaptive"></a>
## 8. Adaptive refinement (cluster.adaptive.workdir)

With a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block, every step of the refinement is an ordinary run in
its own folder under `workdir`. The folders are numbered in order and named by their kind: `0001-run`, `0002-probe`,
and so on. Each holds the step's input (`input.json`, the user's input with the step's request in place of the refinement
options), the request (`request.json`), its log (`run.log`) and its coherence (`coherence_wiDiv`).

`workdir` itself holds:

| File | Content |
|---|---|
| `summary.json` | the report, rewritten after every step: the status, the input and the settings, one record per step and per round, `final_directory`, `final_bundles`, the probe cache counts, the checksums of the input files and `limitations` |
| `final-bundles.json` | the bundles of the final step |
| `final-input.json`, `final-request.json` | an input and its request that rerun the final bundles as one calculation |

The final curve is `<final_directory>/coherence_wiDiv`. The last line of the run log, before the wall time, prints the
status and the path of `summary.json`. The status sets the exit code of `main.out`:

| Status | Exit code |
|---|---|
| `validated_full_bath_model`, `passed_convergence_checks` (merge), `physical` (regroup) | 0 |
| `round_limit`, `no_feasible_improving_merge`, `probe_budget_limited` (merge); `no_candidate`, `no_positive_culprit`, `budget_limited`, `change_limit` (regroup) | 2 |
| an error during the refinement | 1 |

With exit code 2 the refinement stopped at that status, and its last step is still the result.
`example/tutorial/advanced/07_pCCE_and_bundle/cCCE2_Hahn_adaptive` ends this way, with `no_feasible_improving_merge`.

<a id="out_autoconv"></a>
## 9. Automatic convergence (auto_convergence.workdir)

With an [auto_convergence](InputFile_AutoConvergence.md) block, everything goes under `workdir`:

| Path | Content |
|---|---|
| `result/summary.txt`, `result/json/summary.json` | the status, the values chosen, the production result (A with its standard error, T2, n), the core-hours and every warning |
| `result/coherence_mean.txt` | the production mean: columns t (ms), Re, Im, \|L\| |
| `result/1_window.png` ... `result/6_converged_coherence.png` | the figures of the steps and of the converged coherence |
| `result/progress.txt`, `result/json/progress.json` | where the run is, rewritten at every step |
| `runs/<step>/<k>-<key>/` | one folder per CCEX run, with its `input.json`, its outputs and `done.json` once finished |
| `decisions.jsonl` | one line per decision |
| `split.json` | the test and production configurations, when they were found from the files |
| `pcce_cache/` | the partition cache of pCCE runs |
| `lock` | held while a runner works in the folder |

The exit code says how the run ended ([Exit codes](BuildingAndRunning.md#ExitCodes)).
[Results and exit codes](../methods/auto-convergence.md#tag_ac_results) explains each status and how to continue a
stopped job.

<a id="out_analyzer"></a>
## 10. Analyzer outputs

The analyzer (`example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py`, see [Analysis](Analysis.md)) names its files
after its `-fo` option. They have the two-column format of the coherence files, with one difference: the time is in
the unit of `-uc`, not always in ms.

| Option | File |
|---|---|
| `-err` | one error-corrected curve per input file: the `-fo` name with the variable names replaced by their values (`-fo results/ccein_kCONF_ErrCorr` gives `results/ccein_k1_ErrCorr`, ...) |
| `-ea` | the average of all files read: `<-fo>_conv<N>`, where N is the number of curves averaged |
| `-fit_w` | `<-fo>_params`: one row per variable value with the value, T2 (in the `-uc` unit) and the power n; with `-ea`, one row for the average |
| `-pl` | one figure per value, the `-fo` name with the variable replaced by the value, `.png`; with `-ea`, `<-fo>_conv<N>.png` |
| `-pl_m`, `-pl_i` | `<-fo>.png` |
| `-pl_T`, `-pl_p` | with `-fit_w`: `<-fo>_T2.png` and `<-fo>_Power.png` |
