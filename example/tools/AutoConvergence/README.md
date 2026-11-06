# AutoConvergence: automatic convergence of a CCEX calculation

An ordinary CCEX input with an `auto_convergence` block makes `main.out -f input.json` choose the time window, rdip,
rbath and the number of configurations (and states) of a coherence calculation by itself, then run the production.
Supported: cCCE with `cluster.mode` standard or pcce (order ≥ 2), gCCE (order ≥ 1), fixed bundles (`cluster.mode`
bundle without a `cluster.adaptive` block), CPU or GPU backends. You prepare the bath files (and, if you like, state
files with `configuration.s_range`). Almost everything else is found. Every calculation is an ordinary CCEX run in its
own folder under `workdir`.

The full reference (the steps and their rules, every key and default, what is refused, the analysis commands) is
[Automatic convergence tags](../../../docs/user/InputFile_AutoConvergence.md). This page is the short version.

## The smallest input

```json
{
  "method": "cCCE", "quantity": "coherence", "order": 2,
  "Qubit": {"qubit": [{"file": "data/qubit"}]},
  "configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"]},
  "bfield": 500, "pulse": {"count": 1},
  "auto_convergence": {"workdir": "auto"}
}
```

`bathfile` is one name with `{k}`. Every k whose bath file exists is a configuration (or those of
`configuration.k_range`). `rbath`, `rdip`, `nstep` and `deltat` are found, and `nstep` of the input is the number of time
points of every run (150 when absent).

## Running it

```sh
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8 (rank 0 alone runs the runner)
main.out -f input.json               # every CCEX run on one rank
```

- `main.out` starts the Python runner, which needs Python 3.10 or newer with numpy, scipy and matplotlib.
  `CCEX_PYTHON` picks the interpreter.
- The automatic run takes only `-f FILE`. Other options are refused (set them in the input).
- `"enabled": false` in the block runs the ordinary calculation instead.
- To resume after a cut job, call the same command again. Finished runs whose inputs, files and binary are unchanged
  are reused. `"resume": "fresh"` starts over (the earlier workdir is kept as `<workdir>.old<n>`).
- Where a run is (nothing computed, no mpirun): `python -m autoconv status input.json`, from this folder or with it on
  `PYTHONPATH`.

## Results

Everything is under `workdir`: `runs/` (one folder per CCEX run), `decisions.jsonl` and `result/`.
`result/summary.txt` (and `result/json/summary.json`) gives the status, the chosen window, rbath, rdip and
configurations, the production's area A with its SE, T2 and n, and every warning. `result/progress.txt` shows where a
run is, and `result/` holds the figures of every step.

| Exit code | Meaning |
|---|---|
| 0 | converged |
| 3 | unresolved: `window-unresolved`, `truncation-unresolved` (make larger baths), `sampling-unresolved`, `numerical-failure` |
| 2 | refused |
| 4 | stopped (walltime or a signal), so call again to continue |
| 1 | error |

The runner does not make baths. When they are too small or too few, the message says what to add.

## Tests

In this folder:

```sh
python -m pytest -q -n 8 tests            # in parallel (pytest-xdist), about 2.5 min
python -m unittest discover -s tests -v   # one at a time, about 13 min
```
