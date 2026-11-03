# Automatic convergence tags

An `auto_convergence` block in an ordinary input makes `main.out -f input.json` choose the time
window, rdip, rbath and the number of configurations (and states) of a coherence calculation by
itself, then run the production. Every calculation it makes is an ordinary CCEX run in its own
folder under `workdir`. Supported: cCCE with `cluster.mode` standard or pcce (order ≥ 2), gCCE (order ≥ 1), fixed bundles
(`cluster.mode` bundle without a `cluster.adaptive` block), CPU or GPU backends. The outputs and exit codes are [below](#tag_ac_results), what to do when a job is cut or the baths do not
suffice in [When a job is cut](#tag_ac_resume) and [When the baths are too small or too few](#tag_ac_baths);
`example/tutorial/09_auto_convergence/` is a worked input. The runner reads the keys (of the block and of the rest of the
input) without regard to case, as CCEX does, and refuses two keys of one block that differ only in case.

```
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8 (rank 0 alone runs the runner)
main.out -f input.json               # every CCEX run on one rank
```

(The block has no `launcher` key, and one is refused; the runs use the `mpirun -n N` main.out was started
with.)

Where a run is (nothing computed) is the tool's own command ([Commands](#tag_ac_commands)), run from `example/tools/AutoConvergence` or with that
folder on `PYTHONPATH` (`main.out` itself has no status option):

```
PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status input.json
```

The automatic run takes only `-f FILE`; other command-line options are refused (set them in the
input). `configuration.s_range` of the input chooses the approach: without it the ensemble (only the number of
configurations is chosen), with it sampled states. Its numbers are not used here (the start says so): `sampling.test_states`,
`states` and `max_states` set how many states. `configuration.seed` is one integer, 1 when absent (the start says so); it
fixes the split of the configurations and the state and axis files the runner makes (a list, one per configuration of
`k_range`, is CCEX's form for one run of several configurations and is refused here). `rbath`, `rdip`, `nstep` and `deltat`
of the input are not needed; `nstep` is the number of time points of every run (150 when absent).

The block is grouped by step. A section left out uses its defaults; a key given wins over what
would be found; `null` or an empty list is refused (leave the key out to have it found).
Unknown keys and the runner's earlier flat names are refused with a message that says where each
went.

<a id="tag_ac_configurations"></a>
## The bath configurations

The configurations are the CCEX [configuration](InputFile_General.md#tag_configuration) block's:
`bathfile` (one name with `{k}`), and when given `bath_state_file`, `defect_axis_file` and `defect_state_file` (one name
each, `{k}` too), all in `configuration.path` when it is given. With `configuration.k_range` [first, last] the
configurations are those k (every one must have its files); without it, every k whose bath file exists (`config_7`, not
`config_07`). A state file name holds `{s}` for the state number, or the number is appended (as CCEX reads them).
`bath_state_file` left out with `s_range`: the states are made per run and removed (with the on-site states of a Defect
species with `n_onsite_spins`); `defect_axis_file` left out with a Defect species that has axis states: the axes are made
per run. The rule of the made files is `example/tools/StateFiles`' (`make_state_files.py`, [Tools](Tools.md#tag_tool_statefiles)): m uniform over
the 2S+1 levels from numpy `default_rng([seed, k, s])`, axes from `default_rng([seed, k, 0])`, on-site states from
`default_rng([seed, k, s, 1])`; that tool writes the same files for an input beforehand. Each run's input holds one
configuration: the names with k in place of `{k}`, absolute, without `k_range` and `path`, its states as
`s_range` [1, S].

The block has no section of its own for these file names: a former `auto_convergence.configurations` section is
refused ([the old names](InputFile.md#OldNames) say where each key went), and
`example/tools/InputKeys/convert_input_keys.py` rewrites such an input.

<a id="tag_ac_fixed"></a>
## Fixed values instead of convergence

A section is found unless its value is given (before 2.0.0 `"converge": false` said so beside the value, and `converge`
is now refused with this rule):

| Section with | used | effect |
|---|---|---|
| `time_window.value` | a number of ms, or `{"value": ..., "unit": ...}` ([Units](InputFile.md#Units)) | the window is used as it is; $T_d$ is only measured |
| `rdip.value` | Å | the rdip step and its check at the final rbath are skipped (`candidates`, `num_candidates` are then refused) |
| `rbath.value` | Å | the rbath step is skipped; a plateau stops instead of raising rbath |
| `sampling.configurations` | that many configurations (and `sampling.states` states each with `s_range`) | the sampling step is skipped; the production uses exactly that many |

All four given is refused (run CCEX itself).

<a id="tag_ac_keys"></a>
## Keys

| Section | Key | Meaning | Default |
|---|---|---|---|
| (top) | `enabled` | false: `main.out -f` runs the ordinary calculation and the block is not read | true |
| | `resume` | `continue`: reuse what the workdir holds; `fresh`: move it to `<workdir>.old<n>` and start over | `continue` |
| | `workdir` | folder for the runs (`runs/`), the decisions (`decisions.jsonl`) and the results (`result/`) | required |
| `time_window` | `value` | see above: a number of ms or `{"value": ..., "unit": ...}` (`s`, `ms`, `us`, `ns`; [Units](InputFile.md#Units)) | found |
| | `decay_level` | $T_d$ is where the T2 fit reaches this \|L\| level | 0.01 |
| | `revival_check_factor` | the test runs go to this × $T_d$ to see a revival | 2.5 |
| | `production_factor` | production window = this × $T_d$ | 1.5 |
| `rdip` | `value` | see above | found |
| | `candidates` | rdip candidates (Å); the largest is the reference | found |
| | `num_candidates` | how many to start with when found | 7 |
| | `threshold` | relative difference of the fitted area A = T2 Γ(1+1/n) at which two rdip agree | 0.02 |
| `rbath` | `value` | see above | found |
| | `candidates` | rbath candidates (Å) | found |
| | `num_candidates` | how many to start with when found | 4 |
| | `threshold` | relative difference of A at which two rbath agree | 0.02 |
| `t2_fit` | `points` | points of every T2 fit: `peaks` (t = 0 and the local maxima; all points when fewer than 3) or `all` | `peaks` |
| | `floor` | only the points before \|mean\| falls below this are fitted | 0.05 |
| | `n_range` | allowed stretching exponent n of exp(−(t/T2)^n); outside it the fit fails | [0.3, 4.0] |
| `sampling` | `test`, `production` | configuration IDs for the convergence tests and for the production (never the same; within `configuration.k_range` when it is given) | found |
| | `configurations`, `states` | see above | found |
| | `threshold` | SE(A)/A target of the sampling | 0.05 |
| | `test_configurations`, `test_states` | configurations and states of the test runs (states with `s_range`) | 10, 5 |
| | `states_threshold` | the states' share of the sampling error at which the state test stops | 0.05 |
| | `max_states` | the most states the state test computes (with state files, at most the files there) | 30 |
| | `min_configurations` | the fewest production configurations | 10 |
| | `max_configurations` | the most production configurations used | every production file |

The earlier names are refused with the new one: `converge` (every section; a value given fixes it),
`time_window.unit` (write the value as `{"value", "unit"}`), `time_window.threshold` (now `decay_level`: an \|L\| level,
not the relative thresholds of the other sections), `sampling.number` (now `configurations`) and `sampling.test_number`
(now `test_configurations`); `example/tools/InputKeys/convert_input_keys.py` rewrites an old input.

<a id="tag_ac_refused"></a>
## Refused together with the block

Other cluster modes (`dsj`, `itb`, `dsjitb`, the adaptive refinement) or quantities than coherence, `rdip_ladder` in the
input, a `cluster.adaptive` block, `cluster.pcce.cache` (the runner gives its own), `cluster.pcce.partition` legacy,
`output.cluster_contributions`, `output.cluster_only`, a nonzero `cluster.selection.max_per_order`, a nonzero
`configuration.bath_shift`, more than one name in `configuration.bathfile` (one bath file per configuration) or a name
without `{k}`, keys in their former places (CCEX refuses them too, see [the old names](InputFile.md#OldNames)), the block's `launcher` and
`configurations`, more than one qubit, a list together with its count, `rdip.value` with candidates, `sampling.states`
without `sampling.configurations`, `configuration.nstate` (now `s_range`), a `configuration.seed` list, and options out of
range. `output.path` and `output.file` of the input are refused (they were dropped without a word): every run writes
into its own folder under `workdir`, as `-o` is refused.

<a id="tag_ac_rdip_ladder"></a>
## rdip_ladder (an internal key)

`rdip_ladder` is a CCEX input key the runner writes into its runs' inputs: a list of 1 to 64 increasing rdip values
computed from one order-2 run (cCCE, gCCE, cluster.mode pcce), each written as `<output.file>_rdip<r>_...` beside the
run's own rdip. It serves the rdip step, which gets every smaller candidate from the run at the largest one; it is not
part of the documented input (left out of the key lists and of `main.out -e`), and an input given to the automatic
convergence may not hold it.

<a id="tag_ac_results"></a>
## Results and exit codes

Everything is under `workdir`: `runs/` (one folder per CCEX run, each an ordinary run with its own `input.json`),
`decisions.jsonl` (one line per decision) and `result/`. `result/summary.txt` (and `result/json/summary.json`) gives the
status, which keys were found and which sections converged, the window, the rbath and rdip chosen, the production
(configurations, states, the area A with its standard error SE, T2, n), the core-hours (reported only: seconds × launcher
ranks × `OMP_NUM_THREADS`; "unknown" when either cannot be read) and every warning: the start, the combination check,
|L| > 1.01 points (a CCE division by a near-zero cluster term; the data are left as they are), a given window the curve
does not decay in. While a run is going on it says `running`; a refusal, a stop or an error is written there too.
`result/` also holds the figures of the steps and the converged coherence (`example/tutorial/09_auto_convergence/`
lists them).

Decisions belong to one plan (the input, what was found, the source files, the binary and the decision code): when any of
them changes, the earlier `decisions.jsonl` is moved to `decisions.<n>.jsonl` and the choices are made again (the CCEX
runs are reused).

| Exit code | Meaning |
|---|---|
| 0 | converged |
| 3 | unresolved: `window-unresolved`, `truncation-unresolved` (with "make larger baths" when rbath reaches the bath files), `sampling-unresolved`, `numerical-failure` |
| 2 | refused |
| 4 | stopped (walltime or a signal); call again to continue |
| 1 | error |

<a id="tag_ac_resume"></a>
## When a job is cut, and where a run is

Calling the same command again continues: finished CCEX runs (whose inputs, files and binary are unchanged) and the
decisions are reused, only the rest is computed. Under Slurm a run is not started with less walltime left than 1.5 × the
longest run of its step + 300 s (exit 4). A job cut from outside (walltime, a node, a kill) leaves
`result/json/summary.json` saying `running` with nobody holding the workdir's lock: the next call says the run before was
interrupted and continues it. `"resume": "fresh"` starts over instead (the earlier workdir is kept as
`<workdir>.old<n>`).

`result/progress.txt` (and `result/json/progress.json`) is rewritten at every step and after every finished CCEX run: the
status, the step now, the steps done and left, finished CCEX runs per step, the values chosen so far, an earlier
interruption and the command to continue. The `status` command ([Commands](#tag_ac_commands)) prints `running` (a process
holds the lock), `interrupted`, `stopped` or how it ended, then `progress.txt` (and `summary.txt` once it ended).

<a id="tag_ac_baths"></a>
## When the baths are too small or too few

The runner does not make baths; it says what to make:

- **too small** (a plateau, or rbath not converging at the radius the files cover): `truncation-unresolved` (or
  `window-unresolved` for a plateau) with "make larger baths, covering at least R Å around the qubit" (R = the largest
  rbath × 1.25). Larger files under the same names are new data: the next call starts a new plan (finished runs of
  unchanged files are still reused).
- **too few** (the production used every configuration and the SE is still above the target, or fewer than
  max(20, `min_configurations`) exist): `sampling-unresolved` with the measured SE, about how many are needed and how many
  bath (and state) files to add. When the test and production configurations were found from the files, the split is
  kept in `<workdir>/split.json`: the files added are put after the production, and the next call goes on from the last
  batch (the configurations computed are not computed again). A production list given by hand is not extended: add the
  new IDs to `sampling.production` (a new plan).

<a id="tag_ac_limits"></a>
## What a result does not say

The numbers hold for the chosen finite setting, sample, time grid and fit rule. The chosen rdip and rbath are each within
their threshold (default 2 %) of the largest computed, not the combination (the combination check reports it) and not an
infinite bath; the SE is an estimate from the configurations computed (with 10 configurations the true value may be
0.7–1.8 × the estimate); the CCE order and whether the time step is fine enough are not checked; a curve that stays at a
level and oscillates (a physical plateau) is not handled. The start estimates are starts only: the steps converge from
them. Curves that are not one stretched exponential (a collapse and revival, spiky or noisy curves) may stop at the time
window or fail to converge in rbath.

<a id="tag_ac_commands"></a>
## Commands

The runner's other commands are run from `example/tools/AutoConvergence` (or with that folder on `PYTHONPATH`); none of
them starts CCEX:

```
python3 -m autoconv status INPUT.json      # where an automatic convergence is (the workdir is read from the input)
python3 -m autoconv analyze RUN_DIR [RUN_DIR ...] [--config INPUT] [--report-dir DIR]
python3 -m autoconv analyze --manifest manifest.json [--config INPUT]
python3 -m autoconv judge --setting RUN_DIR ... --reference RUN_DIR ... [--report-dir DIR]
```

`analyze` judges existing runs: one `RUN_DIR` per bath configuration, each with the CCEX outputs
`out_state<s>_wiDiv`. It reports, as JSON on standard output, the decay time $T_d$ and whether the window is long enough
(`resolved`, `needs-extension`, `window-unresolved`, or `numerical-failure` for a non-finite state or an overflow), and the
sampling verdict (`t2-met` or `sampling-unresolved`: the standard error of the fitted area against its target) with the
cheapest allocation of configurations and states. `--report-dir DIR` also writes figures and text summaries. The same
directory twice, or a run with fewer states or time points than expected, is refused; a manifest
(`{"configurations": [{"id", "run", "states"}], "available_configurations": N, "grid": {"deltat", "nstep"}}`) gives each
run its ID and its own state count. `judge` compares runs of one setting (a smaller rdip or rbath) with reference runs of
the same configurations, as the rdip and rbath steps do; the paired folders' `input.json` must agree on every key except
`rdip`, `rbath`, `rdip_ladder` and `output.file`. `--config INPUT` (both commands) is the input with the
`auto_convergence` block, read as the runner reads it (the driver file with keys of its own was removed): the time
window is judged unless `time_window.value` is given, the sampling unless `sampling.configurations` is, with
`time_window.decay_level`, `production_factor`, `t2_fit.points` and `sampling.threshold`; the expected states from
`configuration.s_range`, the time points from `nstep`, the bath files available from `configuration.bathfile`. Without it
their defaults are used.

<a id="tag_ac_example"></a>
## Example

```
"configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"], "k_range": [1, 60]},   // the ensemble
"auto_convergence": {
  "workdir": "auto",
  "time_window": {"value": {"value": 15, "unit": "us"}},
  "rbath": {"value": 288},
  "sampling": {"threshold": 0.03}
}
```
