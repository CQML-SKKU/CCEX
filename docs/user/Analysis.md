# Analysis

The analyzer in `example/analysis/` works on the output of CCEX; it does not run CCEX. It needs numpy. The plots and
fits also need scipy and matplotlib. On WSL use `PATH=~/anaconda3/bin:$PATH`, and `MPLBACKEND=Agg` without a display.

| Tool | Answers | Entry point |
|---|---|---|
| `Analyzer/` | make the raw `_noDiv` / `_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (command line) |

## The files

Each file has one line per time step: the time in ms, a tab, and the complex coherence as `+re+imj`. Two real columns
are read too.

- CCEX writes `_noDiv`, the product of the cluster coherences before the division by the sub-clusters, and `_wiDiv`,
  after it.
- The analyzer writes `_ErrCorr` and `<name>_EnsAvg_conv<N>`.

Every file and its columns are in [Output files](OutputFiles.md).

## The analyzer

The analyzer works in three steps. The tutorials' `run.sh` do not call it: they plot the raw `_wiDiv` with
`example/tutorial/plot_coherence.py`.

1. **Error correction** (`-err`). CCE divides a cluster's coherence by those of its sub-clusters, and the division
   blows up near zeros. The divided curve (`_wiDiv`) is therefore repaired against the undivided one (`_noDiv`) with a
   tolerance (`-err_tol`, default 2.0).
   - Only a point where |wiDiv − noDiv| reaches the tolerance is changed.
   - A bad last point is first moved to noDiv + tolerance − 1e-6, so it is no longer bad and stays there.
   - A bad first point is set to 0.
   - Any other bad point becomes the mean of the previous corrected value and the next good one.

   This rule is the group's convention, and every recorded result was made with it. Call it; do not rewrite it.
2. **Ensemble average** (`-ea`) over the values of the variables (`-vn`, `-v`). It writes one `<-fo>_conv<N>` file,
   the mean of the N curves read, as complex numbers.
3. **Plots and fits.**
   - `-pl`: one figure per file.
   - `-pl_m`: all curves in one figure.
   - `-fit`: a stretched-exponential fit with T2 and the power (`-fit_initp`, `-fit_space`).
   - `-pl_i`: an intensity plot of coherence against time and a variable.

The file names carry the variables. `-fi` gives the pattern with the variable names in it, and `-v` the values, so
`-fi ccein_kCONF_wiDiv -vn CONF -v 1 2 3` is three files. Other common options:

| Option | Meaning |
|---|---|
| `-uc` | the time unit of the plots (`ms`, `us`, `ns`, `s`) |
| `-ntn` | turns a `nan` into 0 |
| `-fo` | the output name |

A complete call, for the four bath configurations of `example/tutorial/basic/02_bath_and_states/cCCE2_Hahn_ensemble`
(`configuration.k_range` [1, 4], written as `results/ccein_k<k>_wiDiv` and `_noDiv`):

```
ANALYZER=example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi ccein_kCONF_wiDiv -err_d results/ -err_fi ccein_kCONF_noDiv -vn "CONF" -v 1 2 3 4 -uc "ms" -err -ntn -fo results/ccein_kCONF_ErrCorr
# ensemble average, then a plot with a fit
python3 $ANALYZER -d results/ -fi ccein_kCONF_ErrCorr -vn "CONF" -v 1 2 3 4 -uc "ms" -ea -ntn -fo results/ccein_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
```

**Sampled bath states.** With `configuration.s_range`, CCEX writes one file per state: `ccein_state<s>_wiDiv` and
`_noDiv`, for s = first ... last of `s_range`. With `k_range` the names are `ccein_k<k>_state<s>_...`. Average the
states like the configurations, as complex numbers. For the 20 states of
`example/tutorial/basic/02_bath_and_states/cCCE2_Hahn_states`:

```
python3 $ANALYZER -d results/ -fi ccein_stateST_wiDiv -err_d results/ -err_fi ccein_stateST_noDiv -vn "ST" -v $(seq 1 20) -uc "ms" -err -ntn -fo results/ccein_stateST_ErrCorr
python3 $ANALYZER -d results/ -fi ccein_stateST_ErrCorr -vn "ST" -v $(seq 1 20) -uc "ms" -ea -ntn -fo results/ccein_StateAvg
```

The second call writes `results/ccein_StateAvg_conv20`. `-ea` averages every file the pattern and the values name.
Configurations and states together are therefore one call with two variables, each curve with the same weight:
`-fi ccein_kCONF_stateST_ErrCorr -vn CONF ST -v 1 2 3 4 -v $(seq 1 20)`.

**Density matrices** (`quantity` dm) are written as blocks of rows ([The density matrix](OutputFiles.md#out_dm)).
`read_coherence` and `run.py` do not read them. There are two ways to read them:

- The `read` function of `example/tutorial/plot_coherence.py`. It returns the element between the two levels that are
  coherent at t = 0, divided by its value there. When the qubit starts in one level, it returns that level's
  population instead (as in `example/tutorial/advanced/01_P1/gCCE2_Ramsey_population`).
- Parse the blocks yourself. The first field of a block's first line is the time. The $D^2$ complex values that follow,
  over the block's lines, are ρ row by row.

A population is the real part of a diagonal element. Average states and configurations element by element.

**The fit.** The fitted function is $\exp(-(t/T_2)^n)$.

- It is fitted by least squares to the real part of the curve, at the points of its upper envelope. `-fit_space` sets
  how many points ahead the envelope looks (default 1).
- The fit starts from the time where the curve is closest to 1/e, and from n = `-fit_initp` (default 2).
- A fit with n outside (0.1, 4] or T2 beyond the last time is rejected, and the next start is tried. A value it cannot
  fit is reported and left out of `_params`.

The options that show or write the fit:

| Option | Does |
|---|---|
| `-fit` | draws the fit on the `-pl` figure |
| `-fit_w` | writes it, for one variable (`-vn` with one name), to `<-fo>_params`: one row per value of the variable, sorted by it, with the columns value, T2 and the power n (ten decimals, tab-separated; T2 in the `-uc` unit) |
| `-pl_T`, `-pl_p` | plot T2 and n against the variable (`<-fo>_T2.png`, `<-fo>_Power.png`) |

With `-ea`, `-fit_w` fits the average alone, with any number of variables. It writes one row whose value is N, the
number of curves averaged. The same row comes from a call of its own on the average file, with one value and a
variable name that does not occur in the file name:

```
python3 $ANALYZER -d results/ -fi ccein_StateAvg_conv20 -vn "XX" -v 20 -uc "ms" -fo results/ccein_StateAvg_fit -fit_w -fit_space 5
```

writes `results/ccein_StateAvg_fit_params` (the value 20, T2 in ms, n).

`python3 run.py -h` lists every option. For a new figure, read the curves with `coherence.py` and draw them yourself.
That is simpler than the plot classes:

```python
import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
from coherence import read_coherence, error_correct, ensemble_average, write_coherence
t, L = read_coherence('results/ccein_EnsAvg_conv4', time_scale=1e3)  # t in us
```

`example/tutorial/plot_coherence.py` (the tutorials' figures) is an example of that.

## The convergence test

There is no separate convergence-test script. The automatic convergence (`auto_convergence` block,
`example/tutorial/advanced/12_auto_convergence/`) chooses the time window, rdip, rbath and the number of configurations. To test
one variable by hand, run its values and compare the curves with `read_coherence`.
