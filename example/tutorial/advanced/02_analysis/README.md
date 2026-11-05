# Advanced 02: Analysis of the output

Folder: `example/tutorial/advanced/02_analysis/`

This topic is about what to do with the output of a run. The analyzer of `example/analysis`
([Analysis](../../../../docs/user/Analysis.md)) corrects, averages and fits the curves, and the run itself can report
which clusters matter.

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

This page has two cases:

1. [The error correction, the average and the T2 fit](#a02_analyzer) (`cCCE2_Hahn_analyzer/`).
2. [What a run reports about its clusters](#a02_report) (`cCCE2_Hahn_clusterReport/`).

<a id="a02_analyzer"></a>
## 1. The error correction, the average and the T2 fit

Folder: `cCCE2_Hahn_analyzer/`

<a id="a02_analyzer_input"></a>
### 1.1 The input

The CCEX input is `cCCE2_Hahn_states` of
[Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md) cut to two configurations of ten sampled
states each, 20 curves:

```jsonc
"configuration" : {"gyrofile": "../../../basic/06_defect_bath/bath/DiaP1_gyro",
                   "bathfile": ["../../../basic/06_defect_bath/bath/14N_1ppm/bath_DiaP1_1ppm_{k}"],
                   "k_range": [1, 2], "s_range": [1, 10], "seed": [1, 1]}
```

The analyzer, `example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py`, then works on the 20 output files in three
steps:

- `-err` repairs each `_wiDiv` against its `_noDiv` where the two differ by more than the tolerance.
- `-ea` averages the 20 corrected curves as complex numbers (`results/ccein_Avg_conv20`).
- `-fit_w` fits $\exp(-(t/T_2)^n)$ to the average and writes `results/ccein_Avg_params` (the number of curves, T2 in
  ms and n).

The options that matter are `-fi` with the variables in the file name (`ccein_kCONF_stateST_wiDiv`), `-vn CONF ST` and
one `-v` per variable. `-uc ms` is given in both calls, because the files the analyzer writes keep their time in the
`-uc` unit and the next call reads them as ms.

<a id="a02_analyzer_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

`run.sh` runs CCEX, then the two analyzer calls, and prints the number of curves, T2 and n. It takes 3 s and 0.05 GB
per rank, and needs numpy, scipy and matplotlib.

<a id="a02_analyzer_result"></a>
### 1.3 The result

![cCCE2_Hahn_analyzer: expected result](cCCE2_Hahn_analyzer/expected.png)

The fit gives T2 348 µs and n 0.82. `plot_coherence.py` prints 307 µs for the same average, because its T2 is the last
time |L| ≥ 1/e, not a fit. The error correction repairs no point here, since no division blew up.

Next, `-pl -fit` draws the analyzer's own figure, and a smaller `-err_tol` shows the correction at work.

Input: [`ccein.jsonc`](cCCE2_Hahn_analyzer/ccein.jsonc)

<a id="a02_report"></a>
## 2. What a run reports about its clusters

Folder: `cCCE2_Hahn_clusterReport/`

<a id="a02_report_input"></a>
### 2.1 The input

Three inputs on the system of Basic 01 ask for three reports. They differ only in the `output` block, and the last
one also in `rbath`:

```jsonc
// ccein_contributions.jsonc
"output" : {"cluster_contributions": {"enabled": true}, "savemode": {"mode": "coherence", "nodiv": false}}
// ccein_plan.jsonc
"output" : {"cluster_only": true}
// ccein_curves.jsonc
"rbath"  : 15,
"output" : {"savemode": "cluster_coherence"}
```

- `ccein_contributions.jsonc` writes `results/ccein_contributions.contributions.state-0.tsv`, one row per cluster with
  its sites, its dimension, the largest |log L̃| over the time points and the step where it is reached.
  `.contributions.sites.tsv` gives each site's name, spin, γ and position. `savemode` with `nodiv` false writes
  `_wiDiv` alone.
- `ccein_plan.jsonc` writes the cluster list, `results/ccein_plan.clusters.json`, and computes no dynamics.
- `ccein_curves.jsonc` writes, beside the result, the coherence of each of its 97 clusters,
  `results/ccein_curves_<bath indices>_wiDiv` and `_noDiv`.

<a id="a02_report_run"></a>
### 2.2 Running it

```sh
NP=4 bash run.sh
```

It runs the three inputs, in 1.3 s and 0.05 GB per rank. Sorting the table by its fifth column (`sort -k5 -g -r`, since
the columns are tab-separated and hold no spaces) lists the clusters that matter most, here the pair 621,655
(|log L̃| 0.58).

<a id="a02_report_result"></a>
### 2.3 The result

![cCCE2_Hahn_clusterReport: expected result](cCCE2_Hahn_clusterReport/expected.png)

T2 is 1.084 ms with rbath 50 Å, the curve of Basic 01, and |L| is 0.73 at 2 ms with rbath 15 Å.

Next, `savemode` `cluster_density_matrix` with gCCE writes the density matrix of each cluster, for a few spins.

Inputs: [`ccein_contributions.jsonc`](cCCE2_Hahn_clusterReport/ccein_contributions.jsonc), [`ccein_curves.jsonc`](cCCE2_Hahn_clusterReport/ccein_curves.jsonc), [`ccein_plan.jsonc`](cCCE2_Hahn_clusterReport/ccein_plan.jsonc)
