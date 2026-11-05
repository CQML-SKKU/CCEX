# Analysis

CCEX writes the coherence as text files ([Output files](OutputFiles.md)). This page covers the tools that turn them
into numbers and figures. They read the files and do not run CCEX. Without a display, set `MPLBACKEND=Agg`.

1. [The plot scripts](#an_plot): `example/tutorial/plot_coherence.py` and `plot_convergence.py`, which the tutorials
   use to average curves, print and fit T2 and draw figures. They need numpy and matplotlib.
2. [The analyzer](#an_analyzer): `example/analysis/Analyzer/CoherenceAnalyzer.v1/`, which corrects divisions that blew
   up, averages and fits T2. It needs numpy, scipy and matplotlib.
3. [Density matrices](#an_dm): how to read the files of `quantity` dm.
4. [Convergence](#an_convergence): where the convergence checks are.

<a id="an_plot"></a>
## 1. The plot scripts

`plot_coherence.py OUT.png LABEL=GLOB ...` draws one curve per `LABEL=GLOB`. The files that GLOB matches, for example
all configurations or all states, are averaged as complex numbers, and |L| is plotted. For each curve it prints the
number of files and T2, taken as the last time |L| is above 1/e. From the second curve on it also prints the largest
difference from the first curve.

| Option | Does |
|---|---|
| `--unit ms`, `us`, `ns` | the time unit of the figure and of the numbers printed |
| `--fit` | fits $\exp(-(t/T_2)^n)$ to the tops of the curve and draws the fit ([Basic 01](../../example/tutorial/basic/01_first_run/README.md)) |
| `--members` | draws every averaged file as a thin grey line under the mean ([Basic 02](../../example/tutorial/basic/02_bath_and_states/README.md)) |
| `--spread` | prints the T2 of every averaged file, their deviation and the standard error of their mean ([Basic 03](../../example/tutorial/basic/03_convergence/README.md)) |
| `--dashed` | draws the first curve thick and the others dashed on top, so that curves on top of each other stay visible |
| `--bold-first` | draws the first curve thick and the others faint |
| `--log-time` | a logarithmic time axis, for decays of very different lengths |

`plot_convergence.py` draws a convergence scan: the curves of one input key at several values, and how far each is
from the curve of the last value. [Basic 03: Checking convergence](../../example/tutorial/basic/03_convergence/README.md)
uses it for `rbath`, `rdip`, `order` and the number of configurations.

A density-matrix file (`quantity` dm) is read too. `plot_coherence.py` then takes the element between the two levels
that are coherent at t = 0, divided by its value there. When the qubit starts in one level, it takes that level's
population instead.

<a id="an_analyzer"></a>
## 2. The analyzer

The analyzer is the tool the group has used for its published results. It works in three steps.

1. **Error correction** (`-err`). CCE divides each cluster's coherence by those of its sub-clusters, and the division
   can blow up where a sub-cluster's coherence is near zero. The analyzer repairs the divided curve (`_wiDiv`) against
   the undivided one (`_noDiv`). Only points where the two differ by at least a tolerance (`-err_tol`, default 2.0)
   change. A bad first point becomes 0. A bad last point is moved to just inside the tolerance. Any other bad point
   becomes the mean of the corrected point before it and the next good point. The results in this repository were made
   with this rule, so call it rather than rewrite it.
2. **Average** (`-ea`). The mean of all files named by the pattern and the variable values, as complex numbers. It
   writes `<-fo>_conv<N>`, where N is the number of curves.
3. **Plots and fits.** `-pl` draws one figure per file, or with `-ea` the average alone. `-pl_m` draws all curves in
   one figure, and `-pl_i` the coherence against time and a variable as an intensity plot. These three take one
   variable (`-vn`) unless `-ea` is given. `-fit` adds a fit of $\exp(-(t/T_2)^n)$ to the `-pl` figure, and `-fit_w`
   writes T2 and n to `<-fo>_params`.

The file names carry variables. `-fi` gives the name with variable names in it, `-vn` the variable names and `-v` the
values. For example, `-fi ccein_kCONF_wiDiv -vn CONF -v 1 2 3` reads `ccein_k1_wiDiv`, `ccein_k2_wiDiv` and
`ccein_k3_wiDiv`. `-uc` sets the time unit, `-ntn` turns `nan` into 0, and `-fo` gives the output name.
`python3 run.py -h` lists every option.

Keep `-uc ms` when one call reads the files that another call wrote. The analyzer writes the time of its files in the
`-uc` unit, but reads every file as ms. With `-uc us` in both calls the times are converted twice, and T2 comes out
1000 times too large.

For the ten bath configurations of
[Basic 02](../../example/tutorial/basic/02_bath_and_states/README.md) (`cCCE2_Hahn_ensemble`, files
`results/ccein_k<k>_wiDiv` and `_noDiv`):

```sh
A=../../../../analysis/Analyzer/CoherenceAnalyzer.v1/run.py
python3 $A -d results/ -fi ccein_kCONF_wiDiv -err_d results/ -err_fi ccein_kCONF_noDiv -vn CONF -v $(seq 1 10) \
  -uc ms -err -ntn -fo results/ccein_kCONF_ErrCorr
python3 $A -d results/ -fi ccein_kCONF_ErrCorr -vn CONF -v $(seq 1 10) -uc ms -ea -ntn -fo results/ccein_EnsAvg -fit_w -fit_space 5
cat results/ccein_EnsAvg_params
```

The first call writes `results/ccein_k1_ErrCorr` ... `ccein_k10_ErrCorr`. The second writes the average,
`results/ccein_EnsAvg_conv10`, and its fit, `results/ccein_EnsAvg_params`: one row with N, T2 in ms and n.

With sampled states, the files are `ccein_state<s>_wiDiv`, or `ccein_k<k>_state<s>_wiDiv` together with `k_range`.
Configurations and states can be averaged in one call with two variables, every curve with the same weight:
`-fi ccein_kCONF_stateST_ErrCorr -vn CONF ST -v $(seq 1 10) -v $(seq 1 10)`.
[Advanced 02: Analysis of the output](../../example/tutorial/advanced/02_analysis/README.md) runs the three steps on a
P1 bath.

The fit is made by least squares on the real part of the curve, at the points of its upper envelope. `-fit_space` sets
how many points ahead the envelope looks (default 1). The fit starts where the curve is closest to 1/e, with n =
`-fit_initp` (default 2). A fit with n outside (0.1, 4], or with T2 beyond the last time, is rejected, and the next
start is tried.

To make a figure of your own, read the curves with `coherence.py`:

```python
import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
from coherence import read_coherence
t, L = read_coherence('results/ccein_EnsAvg_conv10', time_scale=1e3)  # t in us, L complex
```

<a id="an_dm"></a>
## 3. Density matrices

With `quantity` dm the files hold blocks of rows ([The density matrix](OutputFiles.md#out_dm)). The analyzer does not
read them. Read them with the `read` function of `example/tutorial/plot_coherence.py`, as described above, or parse
them yourself. The first field of a block's first line is the time, and the D² complex values that follow, over the
lines of the block, are ρ row by row. A population is the real part of a diagonal element. Average states and
configurations element by element.

<a id="an_convergence"></a>
## 4. Convergence

[Basic 03: Checking convergence](../../example/tutorial/basic/03_convergence/README.md) checks `rbath`, `rdip`, `order`,
the number of configurations and the time axis by hand, with `plot_convergence.py`.
[Advanced 01: Automatic convergence](../../example/tutorial/advanced/01_auto_convergence/README.md) does the search by
itself with the `auto_convergence` block.
