# Advanced 02: Analysis of the output

Folder: `example/tutorial/advanced/02_analysis/`

A run writes one coherence curve for each bath configuration and sampled state, so a run of many states ends with many
files rather than one T2. This page shows what to do with these files.

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

This page has two cases:

1. [Turning many curves into one T2](#a02_analyzer) (`cCCE2_Hahn_analyzer/`). A P1 bath is computed for 20 sampled
   states, and the analyzer of `example/analysis` ([Analysis](../../../../docs/user/Analysis.md)) repairs, averages
   and fits the 20 curves.
2. [Finding which clusters matter most](#a02_report) (`cCCE2_Hahn_clusterReport/`). CCEX reports how much each
   cluster of bath spins contributes, which shows the clusters that change the coherence most.

<a id="a02_analyzer"></a>
## 1. Turning many curves into one T2

Folder: `cCCE2_Hahn_analyzer/`

<a id="a02_analyzer_input"></a>
### 1.1 The input

In a P1 bath every P1 points up or down. CCEX computes the bath in a few such states drawn at random (sampled states,
`s_range`) instead of the ensemble, and each state gives its own curve. The coherence of the bath is their average.

The input is `cCCE2_Hahn_states` of
[Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md) cut to two bath configurations of ten
sampled states each, so the run writes 20 curves. It also uses `deltat` 0.01 and `nstep` 300, the seed 1 for both
configurations and `alpha_state` ms 1, so its T2 differs from that of Basic 06:

```jsonc
"configuration" : {"gyrofile": "../../../basic/06_defect_bath/bath/DiaP1_gyro",
                   "bathfile": ["../../../basic/06_defect_bath/bath/14N_1ppm/bath_DiaP1_1ppm_{k}"],
                   "k_range": [1, 2], "s_range": [1, 10], "seed": [1, 1]}
```

<a id="a02_analyzer_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

`run.sh` runs CCEX, then the analyzer `example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py`, and prints the number
of curves, T2 and n. It takes 3 s and 0.05 GB per rank, and needs numpy, scipy and matplotlib.

The analyzer works in three steps:

1. **Repair** (`-err`). CCE divides the coherence of each cluster by those of its smaller clusters, and the division
   blows up where one of them is near zero. The analyzer finds the points where the divided curve (`_wiDiv`) and the
   undivided one (`_noDiv`) differ by more than a tolerance (`-err_tol`, default 2) and repairs them.
2. **Average** (`-ea`). The 20 repaired curves are averaged as complex numbers, into `results/ccein_Avg_conv20`.
3. **Fit** (`-fit_w`). $\exp(-(t/T_2)^n)$ is fitted to the average. `results/ccein_Avg_params` holds the number of
   curves, T2 in ms and n.

`run.sh` makes two calls, the repair and then the average with the fit:

```sh
python3 $A -d results/ -fi ccein_kCONF_stateST_wiDiv -err_d results/ -err_fi ccein_kCONF_stateST_noDiv \
    -vn CONF ST -v $K -v $S -uc ms -err -ntn -fo results/ccein_kCONF_stateST_ErrCorr
python3 $A -d results/ -fi ccein_kCONF_stateST_ErrCorr -vn CONF ST -v $K -v $S -uc ms -ea -ntn \
    -fo results/ccein_Avg -fit_w -fit_space 5
```

`-fi` is the file name with two variables, `CONF` and `ST`. `-vn` names them, and each `-v` gives the values of one
(`K` is 1 2 and `S` is 1 to 10). Keep `-uc ms` in both calls, because the files of the first call have their times in
the `-uc` unit and the second call reads them as ms. [Analysis](../../../../docs/user/Analysis.md) has every option.

<a id="a02_analyzer_result"></a>
### 1.3 The result

![cCCE2_Hahn_analyzer: expected result](cCCE2_Hahn_analyzer/expected.png)

The fit gives T2 348 µs and n 0.82. `plot_coherence.py` prints 307 µs for the same average, because the two define T2
differently:

- The analyzer fits the whole curve and takes T2 from the fit.
- `plot_coherence.py` takes the last time where |L| ≥ 1/e.

No point was repaired here, since no division blew up. `run.sh` prints the number of points the repair changed, 0.

Next, `-pl -fit` draws the analyzer's own figure, and a smaller `-err_tol` shows the correction at work.

Input: [`ccein.jsonc`](cCCE2_Hahn_analyzer/ccein.jsonc)

<a id="a02_report"></a>
## 2. Finding which clusters matter most

Folder: `cCCE2_Hahn_clusterReport/`

<a id="a02_report_input"></a>
### 2.1 The input

CCE splits the bath into small clusters of spins, computes each, and multiplies their contributions into the curve. A run can therefore report what each cluster contributes. The measure is
$\max_t \lvert\log \tilde L_C(t)\rvert$, how much the cluster's own contribution changes the coherence.

Three inputs on the system of Basic 01 ask for three reports. They differ only in the `output` block, and the last one
also in `rbath`:

```jsonc
// ccein_contributions.jsonc
"output" : {"cluster_contributions": {"enabled": true}, "savemode": {"mode": "coherence", "nodiv": false}}
// ccein_plan.jsonc
"output" : {"cluster_only": true}
// ccein_curves.jsonc
"rbath"  : 15,
"output" : {"savemode": "cluster_coherence"}
```

- `ccein_contributions.jsonc` writes a table with one row per cluster and its contribution. Use it to find the
  clusters that change the coherence most.
- `ccein_plan.jsonc` writes the list of clusters and computes no dynamics. Use it to see the size of a run before
  running it.
- `ccein_curves.jsonc` writes the contribution of every cluster as a curve, one file per cluster. It writes many
  files, so it is for small baths, here `rbath` 15 Å.

The files and their columns are in [Output files](../../../../docs/user/OutputFiles.md#out_contributions).

<a id="a02_report_run"></a>
### 2.2 Running it

```sh
NP=4 bash run.sh
```

It runs the three inputs, in 1.3 s and 0.05 GB per rank. The contribution table is
`results/ccein_contributions.contributions.state-0.tsv`. Its first two lines are a comment and the column names, and
the fifth column is the contribution, so the largest come first with:

```sh
tail -n +3 results/ccein_contributions.contributions.state-0.tsv | sort -k5 -g -r | head -3
```

```
621,655    621,655    1    4    0.57992775548967024    0    1950
341,342    341,342    1    4    0.46308521343074627    0    1999
615,646    615,646    1    4    0.44953573569736255    0    1999
```

The first column is the cluster's sites (the three rows above are pairs of ¹³C), and `results/ccein_contributions.contributions.sites.tsv`
gives the position of each site. The plan is `results/ccein_plan.clusters.json`, and the curves are
`results/ccein_curves_<sites>_wiDiv`, 97 clusters in all.

<a id="a02_report_result"></a>
### 2.3 The result

![cCCE2_Hahn_clusterReport: expected result](cCCE2_Hahn_clusterReport/expected.png)

The figure shows the coherence of the contributions run (`rbath` 50 Å) and of the curves run (`rbath` 15 Å). T2 is
1.084 ms with rbath 50 Å, the curve of Basic 01, and |L| is 0.73 at 2 ms with rbath 15 Å. In the table the pair of
sites 621 and 655 contributes most, with $\lvert\log \tilde L\rvert$ 0.58.

Next, `savemode` `cluster_density_matrix` with gCCE writes the density matrix of each cluster, for a few spins.

Inputs: [`ccein_contributions.jsonc`](cCCE2_Hahn_clusterReport/ccein_contributions.jsonc), [`ccein_curves.jsonc`](cCCE2_Hahn_clusterReport/ccein_curves.jsonc), [`ccein_plan.jsonc`](cCCE2_Hahn_clusterReport/ccein_plan.jsonc)
