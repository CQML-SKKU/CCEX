# example/analysis/

One Python tool that post-processes the coherence files CCEX writes (it does not
run CCEX). This page is written so that a reader can decide quickly what to
call and what to write themselves.

| Tool | Answers | Entry point | Needs |
|---|---|---|---|
| `Analyzer/` | make the raw `_noDiv`/`_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (CLI) | numpy; the CLI and plots also scipy and matplotlib |

On WSL it needs `PATH=~/anaconda3/bin:$PATH` (the system python3 has no numpy)
and, without a display, `MPLBACKEND=Agg`.

Removed before 2.0.0 (git history has them): `ConvTest/` (the convergence over
one variable, `do_convTest.py`; the automatic convergence of `tools/AutoConvergence/`
and `configuration.k_range` took its place) and `AdaptiveStudy/` (the scripts of the
adaptive-refinement study).

## Use as-is, do not rewrite

**The error-correction rule** (`coherence.error_correct`). CCE divides a cluster's
coherence by its sub-clusters' and the division blows up near zeros; the
analyzer repairs the divided curve against the undivided one with a tolerance,
by a rule that is the group's convention (`Analyzer/README.md` states it). Every
recorded result and every published curve was made with it, so a
re-implementation that differs even slightly changes what "the same result"
means. Call the function, or the CLI with `-err`.

**Reading and writing the files** (`coherence.read_coherence`,
`coherence.write_coherence`). The format is simple (time in ms, tab, `+re+imj`),
but the reader also handles the two-column real/imag form and a `nan` in the
coherence column, and the writer's 10-decimal format is what every recorded
result was written with.

**The ensemble average** (`coherence.ensemble_average`), so that the `EnsAvg`
files it writes stay comparable with the recorded ones.

## Fine to do yourself

**Plots and fits.** `DataSpace` and `CoherencePlot_v2` draw and fit (stretched
exponential, T2 and power) the way the 2023 figures were made. For a new figure,
reading the curves with `coherence.py` and drawing them directly is simpler than
learning the class; `example/tutorial/plot_coherence.py` is such a script (|L| of the
raw `_wiDiv` files, averaged over files, and a T2).

## How the tutorials use it

A tutorial's `run.sh` runs CCEX and then
`example/tutorial/plot_coherence.py` on the raw `_wiDiv` files; it does not call the
analyzer. `test/code_verification/check_input_forms` (`errcorr.sh`) calls `run.py -err`;
`test/code_verification/check_numerics/` compares raw outputs and does not call it.
