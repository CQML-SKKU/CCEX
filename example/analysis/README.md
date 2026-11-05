# example/analysis/

One Python tool that post-processes the coherence files CCEX writes; it does not run CCEX. This page says which parts
to call as they are and which are fine to write yourself.

| Tool | Answers | Entry point | Needs |
|---|---|---|---|
| `Analyzer/` | make the raw `_noDiv`/`_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (CLI) | numpy; the CLI and plots also scipy and matplotlib |

Without a display, set `MPLBACKEND=Agg` for the plots.

## Use as-is, do not rewrite

**The error-correction rule** (`coherence.error_correct`). CCE divides a cluster's coherence by its sub-clusters',
and the division blows up near zeros. The analyzer repairs the divided curve against the undivided one with a
tolerance, by a rule that is the group's convention ([Analyzer/README.md](Analyzer/README.md) states it). Recorded
and published curves were made with it, so a re-implementation that differs even slightly gives different numbers.
Call the function, or the CLI with `-err`.

**Reading and writing the files** (`coherence.read_coherence`,
`coherence.write_coherence`). The format is simple (time in ms, tab, `+re+imj`),
but the reader also handles the two-column real/imag form and a `nan` in the
coherence column, and the writer's 10-decimal format is what every recorded
result was written with.

**The ensemble average** (`coherence.ensemble_average`), so that its `EnsAvg` files stay comparable with
recorded ones.

## Fine to do yourself

**Plots and fits.** `DataSpace` and `CoherencePlot_v2` draw and fit (a stretched exponential: T2 and the power) the
way the group's 2023 figures were made. For a new figure it is simpler to read the curves with `coherence.py` and draw
them directly; `example/tutorial/plot_coherence.py` is such a script.

## Who calls it

[Advanced 02: Analysis of the output](../tutorial/advanced/02_analysis/README.md) and
[Advanced 03: Runs on a cluster](../tutorial/advanced/03_server_runs/README.md) run it on their results. The other
tutorials plot the raw `_wiDiv` files with `example/tutorial/plot_coherence.py`. The check
`test/code_verification/check_input_forms/errcorr.sh` calls `run.py -err`.
