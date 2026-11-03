# AutoConvergence — automatic convergence of a CCEX calculation

An ordinary CCEX input with an `auto_convergence` block makes `main.out -f input.json` choose the time window, rdip, rbath
and the number of configurations (and states) of a coherence calculation by itself - cCCE with `cluster.mode` standard or pcce (order ≥ 2) or gCCE (order ≥ 1), fixed
bundles (`cluster.mode` bundle without a `cluster.adaptive` block), CPU or GPU backends, DFT tensor files - then run the
production (the order is the input's; the adaptive refinement is excluded). The rdip ladder (one run gives every smaller rdip)
exists for order 2 of cCCE, gCCE and pCCE with standard clusters; otherwise (order 3 and above, fixed bundles) each rdip
candidate is its own run, which costs more but is judged alike. You
prepare the bath files (and, if you like, state files with `configuration.s_range`); almost everything else is found. Every calculation is an
ordinary CCEX run in its own folder.

```sh
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8 (rank 0 alone runs the runner)
main.out -f input.json               # every CCEX run on one rank
```

`"enabled": false` in the block runs the ordinary calculation instead (the block is not read). The automatic run takes
only `-f FILE`: another option (`-o`, `-v`, ...) would not reach the runner, so it is refused - set it in the input. (The
command-line options `--auto_conv` and `--auto_status` are refused: the block starts the run, and the
status is this tool's `status` command below.) Under `mpirun -n N`
every rank but 0 leaves at once; rank 0 drops the outer launch's MPI variables (Open MPI refuses an mpirun started from inside
one otherwise) and the runner starts every CCEX run with `mpirun -n N` (the `mpirun` on PATH); without mpirun every run is
on one rank. (A `launcher` key in the block is refused.) `summary.txt` section 0 says which
launcher was used. `main.out` starts the Python runner
(`python -m autoconv run`; Python 3.10 or newer with numpy, scipy and matplotlib, `CCEX_PYTHON` picks the interpreter). Calling it again
continues: finished runs are reused when their inputs, files and binary are unchanged.

## When a job is cut, and where a run is

Calling the same command again continues: finished CCEX runs and the decisions are reused, only the rest is computed. Under
Slurm a run is not started with too little walltime left (exit 4). A job cut from outside (walltime, a node, a kill) leaves
`result/json/summary.json` saying `running` with nobody holding the workdir's lock: the next call says the run before was interrupted
and continues it. `"resume": "fresh"` in the block starts over instead (the earlier workdir is kept as `<workdir>.old<n>`).

- `result/progress.txt` (and `result/json/progress.json`), rewritten at every step and after every finished CCEX run: the status, the step now,
  the steps done and left, finished CCEX runs per step, the values chosen so far, an earlier interruption, and the command to
  continue.
- `python -m autoconv status input.json` (run in this folder, or with it on `PYTHONPATH`; no mpirun, nothing computed; the
  `workdir` is read from the input, relative to the input's folder): `running` (a process holds the lock), `interrupted`,
  `stopped` or how it ended, then `progress.txt` (and `summary.txt` once it ended).

## When the baths are too small or too few

The runner does not make baths (a `bath_generation` section is refused); it says what to make:

- **too small** — a plateau or rbath not converging at the radius the files cover — `truncation-unresolved` (or
  `window-unresolved` for a plateau) with "make larger baths, covering at least R Å around the qubit" (R = the largest rbath
  × 1.25). Larger files under the same names are new data: the next call starts a new plan (finished runs of unchanged files
  are still reused).
- **too few** — the production used every configuration and the SE is still above the target, or fewer than
  max(20, `min_configurations`) exist — `sampling-unresolved` with the measured SE, about how many are needed and how many
  bath (and state) files to add. When the test and production configurations were found from the files, the split is kept in
  `<workdir>/split.json`: the files added are put after the production and the next call goes on from the last batch (the
  configurations computed are not computed again). A production list given by hand is not extended: add the new IDs to
  `sampling.production` (a new plan).

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

The bath configurations are the `configuration` block's: `bathfile` (one name with `{k}`), and when
given `bath_state_file`, `defect_axis_file` and `defect_state_file` (one name each, `{k}` too), all in `configuration.path`
when it is given. With `configuration.k_range` [first, last] the configurations are those k (each must have its files);
without it, every k whose bath file exists. Each run's input holds one configuration: the names with k in place of `{k}`,
absolute, without `k_range` and `path`. (Before 2.0.0 the patterns were `auto_convergence.configurations.bathfile`,
`statefile`, `avaaxfile` and `exstatefile`, and the rest of that section is now `auto_convergence.sampling`; both old forms
are refused, and `example/tools/InputKeys/convert_input_keys.py` rewrites an old input.)

Divergent points: wherever a curve is read, a point with |wiDiv − noDiv| ≥ 2 or a nan/inf wiDiv takes the noDiv
value (not the rule of the group's analyzer: its `error_correct` replaces such a point by the mean of the previous corrected
point and the next good wiDiv, which gives other numbers); a configuration that still has |L| > 10 or nan/inf (noDiv bad too,
or no noDiv with `"output": {"savemode": {"mode": "coherence", "nodiv": false}}`) is left out of every average and judgement, and more than 5 % of the production left out
is a numerical failure. `summary.txt` lists the corrected points and the configurations left out, with the nearest bath spin.

`configuration.s_range` chooses the approach (`configuration.nstate` before 2.0.0, now refused): **without
it** the ensemble, one curve per configuration, only the number of configurations is chosen; **with it** states sampled (its
numbers are not used here - the start says so -: `sampling.test_configurations`, default 10 configurations, and
`sampling.test_states`, default 5 states, set the test), and `configuration.bath_state_file` (`data/cfg{k}_state{s}` reads
`data/cfg{k}_state1`, `..._state2`, ...; without `{s}` the number is appended, `data/cfg{k}_state` reads the same files) is
optional: when given, the state files of each
configuration are counted (every configuration needs at least `test_states` of them; the state test below uses more
when the test configurations hold them); when left out, the runner makes them for each run in the run folder and removes them after (no stored state
files; m_s drawn evenly from each species' sublevels with `default_rng([seed, k, s])`, one draw per bath-file line, so a spin
keeps its state at every rbath and rdip and the paired comparisons hold; up to `sampling.max_states` states, default 30;
with them the on-site states of a Defect species with `n_onsite_spins`; `configuration.seed` of the input - one
integer -, else 1, which the start says). The rule is `example/tools/StateFiles`' (`make_state_files.py`), which writes the
same files for an input beforehand (one rule; `examples/nv13c/make_states.py` uses it too). `output.path` and
`output.file` of the input are refused: the runs go to `workdir`.
CCEX's own random states are not used: they change with rbath. `rbath`, `rdip`, `nstep` and `deltat` of the
input are not needed: they are found. `nstep` of the input is the number of time points of every run (150 when absent); a
`deltat` given only sets the first window ((nstep − 1) × deltat): nstep stays fixed and the time step follows the window.
`examples/nv13c/input_auto.jsonc` is this form; `input.jsonc` and `input_ensemble.jsonc` give the candidates.

## Steps

| Step | What it does | Evidence |
|---|---|---|
| 0. start | Finds the configurations: every existing file of the `bathfile` pattern (`{k}` a non-negative integer); a configuration with a missing state file is refused, not dropped. Splits them: `test_configurations` (default 10) for the convergence tests (`test`; fewer when `production` would hold less than needed), chosen evenly over the distance from the qubit to each configuration's nearest bath spin (the closer that spin, the faster the decay), and the rest for the production (`production`) in a seeded order (input `seed`, default 1). Counts the spins of every species around the qubit in the test configurations and estimates rough starts (the steps converge from them): the first rbath candidates (`rbath.num_candidates` values 1.25 apart, the largest where 300 spins are counted - bulk, a layer, a qubit above a layer, clusters alike), a decay time from the couplings T_seed = 51.4 / Σ_s ρ_s (D_qs² D_ss)^(1/3) (ρ_s the species' count in that sphere over its volume, D_ab = 1.0546 \|γ_a γ_b\| rad·Å³/ms; uncalibrated for layers and mixtures), and the first rdip candidates (`rdip.num_candidates` values 0.8 apart below 1.25² × the largest species' (D_ss T_seed / 0.25)^(1/3)). A species absent from the start sphere is named. The radius the bath files cover is the distance from the qubit to the nearest face of the bath's box; a thin bath (a layer: `opt2D`, or a shortest side below 0.2 of the longest) is no limit in its thin direction - √(a² + h²) from the in-plane edge distance a and the qubit height h - and when the in-plane size is too small the run asks for larger baths | `summary.txt` section 0 |
| 1. time window | The test configurations run to `revival_check_factor` × T_seed on `nstep` points of the input (150 when it has none) (deltat = window / (points − 1)). T2 fit exp(−(t/T2)^n) on the upper envelope of their mean; T_d = T2 (ln 100)^(1/n) (where the fit reaches `time_window.decay_level`). The window grows (× 1.5, at most 8 times) until the test runs reach `revival_check_factor` × T_d, or is cut once when far too long; when the mean is already below 0.05 within the first 10 points the grid is too coarse and the window is cut to the decay instead (at most 3 times). A revival (the mean comes back above 0.01) extends it. A plateau (decayed, then stays above 0.05, also in a 1.5 × longer window) raises rbath and the rdip top by 1.25 (full steps within the bath files) **Noise floor:** the fit floor and the threshold are raised to 2 × the median SE of the test mean over the last third of the window when that is higher (a mean of single-state curves stays near it), and the later fits of the run stop there too; a plateau must be above twice it. **Fit points:** every T2 fit of the run (window, rdip, rbath, sampling, production) uses the peaks of its curve: t = 0 and the points higher than both neighbours (y[i] > y[i−1], y[i] ≥ y[i+1]), up to the first one below the fit floor; a curve with fewer than 3 peaks (a smooth decay) is fitted on all its points. The T2 is a convergence measure, compared between settings, not the exact T2 of the system (summary.txt says so). `1_window_fit.png` shows the test mean, the fit points and their fit, with the fit on all points for reference. nstep is never changed (150 when the input has none) and whether it suffices is not checked. | `1_window.png`, `1_window_fit.png`, `summary.txt` section 1 |
| 2. rdip | One run at the largest rdip gives every smaller one (`rdip_ladder`; order 2 of CCE, gCCE, pCCE, otherwise one run per value). Two values agree when the area under their T2 fits, A = T2 Γ(1+1/n), differs by at most `rdip.threshold`. **When either T2 fit fails**, that comparison (here, in rbath and in the combination check) judges both curves on the raw area ∫\|L̄\|dt over [0, T_d] (trapezoid on the configuration mean) with the same threshold - a fail there is an actual fail; `summary.txt`, the figure (open circles) and the `truncation` decision in `decisions.jsonl` (`raw_area`) say which comparisons were; such a comparison gives no trend estimate. Fitted comparisons are unchanged. The largest is raised by 1.25 (up to 2 × rbath) until the last comparisons agree — two in a row when the candidates were found — and the smallest rdip from which every larger one agrees with the largest is chosen | `2_rdip.png` |
| 3. rbath | The same rule with `rbath.threshold`, one run per value, up to the radius the bath files cover. Then rdip is judged again at the final rbath (the same consecutive passes at its top); when it does not hold there, rdip is converged again at that rbath (raised if needed) and rbath once more | `3_rbath.png`, `summary.txt` sections 2 and 3 |
| 3b. window again | When the final reference (largest rbath, rdip) is 1.25 × or more above the one the window was found at, T_d is measured there; outside 0.8–1.25 × the earlier T_d, steps 1–3 run again from it (at most 3 rounds). The fitted area cannot see a shape change: a smaller rbath may still be chosen whose production curve then does not decay in the new window, and the run stops there as `window-unresolved` | `summary.txt` section 1 ("found again") |
| 4. sampling | From the test configurations at the chosen rbath and rdip. **States (a check):** the states' share of the sampling error of A (the test configurations kept, only their states drawn again) at 2, 3, ..., `test_states`; above `sampling.states_threshold` (default 0.05) at `test_states`, 5 more states (each new block computed alone, at most `max_states` (default 30) or the state files the test configurations hold) and measured again; an S above `test_states` or none within the threshold is a warning only (`4_states.png`: the mean curves at 1, 2, ... states, the errors and the threshold). **Configurations:** the states stay `test_states`; the variance of A between configurations and between states gives the configurations needed for SE(A)/A ≤ `sampling.threshold`, never fewer than `min_configurations` | `4_states.png`, `summary.txt` sections 4 and 5 |
| 5. production | The chosen setting on `production` configurations, on the same number of points over `production_factor` × T_d; generously first, then more when short: 1.5 × the pilot's estimate of the configurations needed, at least max(20, `min_configurations`), rounded up to a multiple of 5, at most the production holds; judged on everything computed; when short, 1.5 × the shortage the measured SE implies (at least 20) more, until the SE meets the target or the files are used up (then reported with how many to add); `5_configurations.png`: the production's configurations growing (2, 5, 10, ..., all; 2, 3, ... up to 10), left their mean curves, right the SE(A)/A measured at each count and the threshold; the production curve must decay inside the window; the chosen setting is compared once with the largest rbath and rdip on the first 3 production configurations (a warning when it differs by more than the larger threshold); **Production window:** when the production curve has not decayed inside its window, the production is computed again on a longer window (production_factor × its T_d, at least 1.5 ×), the same configurations first, at most 3 times; `summary.txt` lists the windows. | `5_configurations.png`, `6_converged_coherence.png`, `summary.txt`, `coherence_mean.txt` |

## The `auto_convergence` block

Grouped by step. A section left out uses its defaults; a key given wins over what would be found; `null` (also `deltat: null`) or an
empty list is refused - leave the key out to have it found. A section is found unless its value is given (the earlier
`converge` switch said the same twice and is refused now):

| Section with | used | effect |
|---|---|---|
| `time_window.value` | ms, or `{"value": ..., "unit": ...}` | the window is used as it is; T_d is only measured (a warning when the curve does not decay in it) |
| `rdip.value` | Å | step 2 and the rdip check at the final rbath are skipped |
| `rbath.value` | Å | step 3 is skipped; a plateau stops instead of raising rbath |
| `sampling.configurations` | that many (and `sampling.states` states with `s_range`) | step 4 is skipped; the production uses exactly that many configurations, its SE is reported |

Step 3b runs only when the time window, rdip and rbath are all found. All four given is refused (run CCEX itself). Example,
rdip and the configurations only:

```json
"auto_convergence": {
  "workdir": "auto",
  "time_window": {"value": {"value": 15, "unit": "us"}},
  "rbath": {"value": 288}
}
```

| Section | Key | Meaning | Default |
|---|---|---|---|
| (top) | `enabled` | false: `main.out -f` runs the ordinary calculation | true |
| | `resume` | `continue`: reuse what the workdir holds; `fresh`: move it to `<workdir>.old<n>` and start over | continue |
| | `workdir` | folder for the runs (`runs/`), decisions (`decisions.jsonl`) and results (`result/`) | required |
| `t2_fit` | `points` | points of every T2 fit (window, rdip, rbath, sampling, production): `peaks` (t = 0 and the points higher than both neighbours; all points when fewer than 3) or `all` | peaks |
| | `floor` | only the points before \|mean\| falls below this are fitted | 0.05 |
| | `n_range` | allowed stretching exponent n of exp(−(t/T2)^n); outside it the fit fails | [0.3, 4.0] |
| `sampling` | `test`, `production` | configuration IDs for the convergence tests and for the production (never the same) | found |
| | `configurations`, `states` | see above | found |
| | `threshold` | SE(A)/A target of the sampling | 0.05 |
| | `min_configurations` | the fewest production configurations | 10 |
| | `max_configurations` | the most production configurations used (raising it later continues the same run) | every production file |
| | `test_configurations`, `test_states` | configurations and states of the test runs (states with `s_range`) | 10, 5 |
| | `states_threshold` | the states' share of the sampling error at which the state test stops | 0.05 |
| | `max_states` | the most states the state test computes (with `s_range`; with state files, at most the files there) | 30 |
| `time_window` | `value` | see above (units `s`, `ms`, `us`, `ns` in `{"value", "unit"}`) | found |
| | `decay_level` | T_d is where the fit reaches this \|L\| level | 0.01 |
| | `revival_check_factor` | the test runs go to this × T_d to see a revival | 2.5 |
| | `production_factor` | production window = this × T_d | 1.5 |
| `rdip` | `value` | see above | found |
| | `candidates` | rdip candidates (Å); the largest is the reference | found |
| | `num_candidates` | how many to start with when found | 7 |
| | `threshold` | relative difference of A at which two rdip agree | 0.02 |
| `rbath` | `value` | see above | found |
| | `candidates` | rbath candidates (Å) | found |
| | `num_candidates` | how many to start with when found | 4 |
| | `threshold` | relative difference of A at which two rbath agree | 0.02 |

Renamed or removed before 2.0.0, refused with the new name: `converge` (every section), `time_window.unit`,
`time_window.threshold` (now `decay_level`), `sampling.number` (now `configurations`), `sampling.test_number` (now
`test_configurations`).

pCCE partition cache: the runner gives every pCCE run `pcce.cache.mode` `readwrite` and a file
`<workdir>/pcce_cache/<key>.cache`, the key a hash of the child input without `rdip`, `rdip_ladder`, `configuration.s_range`, `configuration.bath_state_file`, `configuration.defect_state_file`,
`nstep` and `deltat`, of the bath (and other) file contents and of the binary. Runs that differ only in rdip (or states and
time grid) therefore search the k-means partition once and read it afterwards; CCEX itself compares the stored partition's
inputs byte for byte and searches again on any difference. The key is part of the run key; `done.json` records whether the
run read (`hit`) or stored the file and its hash, and a cache file is read only when its hash is the one recorded when it was
stored. Outputs are the same bytes as without the cache. Changing to this version recomputes the pCCE runs of an earlier
workdir once (their keys changed); other methods' runs are reused.

Refused: other methods, cluster modes (`dsj`, `itb`, `dsjitb`) or quantities, `rdip_ladder` in the input, a `cluster.adaptive`
block, `cluster.pcce.cache` (the runner gives its own), `cluster.pcce.partition` legacy (its partitions depend on the MPI size),
`output.cluster_contributions`, `output.cluster_only`, a nonzero `cluster.selection.max_per_order` (the files `Qubit.qubit[0].file`,
`configuration.gyrofile`, `hamiltonian.common.hyperfine.file`, `quadrupole.file` and `quadrupole.file_without_defect` are
hashed and allowed), `pilot_states` (use `sampling.test_states`), `configuration.nstate` (now `s_range`), a `configuration.seed`
list, `output.path` and `output.file`, `configuration.bath_shift`, more than one name in
`configuration.bathfile` (one bath file per configuration), a `bathfile` without `{k}`, the keys in their places before
2.0.0 (the configuration keys at the top; `cluster_mode`, `outfile`, `savemode`, `backend`, ... at the top and
`n_nuclear_spins`/`n_configurations` - CCEX refuses them too), the block's `launcher` and `configurations`,
more than one qubit, a list together with its count, `rdip.value` with candidates, `sampling.states` without
`sampling.configurations`,
options out of range, unknown keys, and the earlier flat names (`threshold_rdip`, `candi_rdip`, `rdip_rungs`, `se_target`,
`time_window` as a number, ...) or removed keys (`max_spins`, `max_pairs`, `max_core_hours`, `cost_config`, `cost_state`,
`margin_s`, `run_timeout_s`, `truncation_max_tol`, `truncation_area_tol`, `version`, `source_radius`, `max_nstep`,
`states_available`, `recheck`, `default_nstep`, `time_window.method` (always the T2 fit), `time_window.smoothing_points`
(always 9: a running median that ignores spikes of up to 4 points), `sampling.pilot`/`pool` - now `test`/`production`, `fit_points` - now `t2_fit.points`, `t2_fit.step_threshold`, `t2_fit.quality`, `t2_fit.step_limit` and `time_window.max_points`) — the message says
where each went. rbath never goes above the radius the bath files cover: the distance from the qubit to the nearest face of
the box their spins span, the smallest over the configurations.

## Results and exit codes

`result/summary.txt` (and `result/json/summary.json`) gives the status, which keys were found and which sections converged, the window, the
rbath and rdip chosen, the production (configurations, states, A with its SE, T2, n), the core-hours (reported only:
seconds × launcher ranks × `OMP_NUM_THREADS`; "unknown" when either cannot be read) and every warning: the start, the
combination check, |L| > 1.01 points (a CCE division by a near-zero cluster term; the data are left as they are), a given
window the curve does not decay in. While a run is going on it says `running`; a refusal, a stop or an error is written there
too. Decisions belong to one plan (the input, what was found, the source files, the binary and the decision code): when any
of them changes, the earlier `decisions.jsonl` is moved to `decisions.<n>.jsonl` and the choices are made again (the CCEX
runs are reused). Under Slurm a run is not started with less walltime left than 1.5 × the longest run of its step + 300 s
(exit 4; call again to continue).

Exit codes: 0 converged; 3 unresolved (`window-unresolved`, `truncation-unresolved` with "make larger baths" when rbath
reaches the bath files, `sampling-unresolved`, `numerical-failure`); 2 refused; 4 stopped (walltime or a signal; call again
to continue); 1 error.

## What a result does not say

The numbers hold for the chosen finite setting, sample, time grid and fit rule. The chosen rdip and rbath are each within
their threshold (default 2 %) of the largest computed, not the combination (the check reports it) and not an infinite bath;
the SE is an estimate from the configurations computed (with 10 configurations the true value may be 0.7–1.8 × the
estimate); the CCE order and whether the time step is fine enough are not checked; a curve that stays at a level and
oscillates (a physical plateau) is not handled. The start estimates (300 spins, the clock constant 51.4 from the NV-¹³C
example and P1 100 ppm 500 G) are starts only: the steps converge from them. Tested end to end on the NV-¹³C example with
the candidates given, and from the automatic starts on dense baths: P1 cCCE2 and gCCE2 (ensemble and
single) and an h-BN gCCE2 ensemble at 3 T converged; P1 pCCE2 does not converge in rbath (spiky curves: the fitted areas
disagree, and where the fit fails the raw area disagrees too), and h-BN at 500 G (collapse and revival, not one stretched exponential) and noisy curves stop at the time window.

## Other commands (analysis of existing runs)


```sh
cd example/tools/AutoConvergence
python -m autoconv status INPUT.json     # where an automatic convergence is, nothing computed
python -m autoconv analyze RUN_DIR [RUN_DIR ...] [--config INPUT] [--cost-config A --cost-state B]
                            [--expect-states N] [--expect-points N]
python -m autoconv analyze --manifest manifest.json [--config INPUT] [--cost-config A --cost-state B]
python -m autoconv judge --setting RUN_DIR ... --reference RUN_DIR ... [--setting-prefix out_rdip251.359] [--reference-prefix out]
                          [--decay-time T] [--name NAME] [--report-dir DIR]
python -m pytest -q -n 8 tests            # the checks in parallel (pytest-xdist), about 2.5 min
python -m unittest discover -s tests -v   # the same one at a time, about 13 min
```

`--report-dir DIR` writes figures and text summaries (`autoconv/report.py`): `sampling.png`/`.txt` (\|mean\| with its fit and
T_d; SE(A)/A, the relative standard error of the fitted area A = T2 Γ(1+1/n), against C with the target, the bootstrap point and
the allocation) and `truncation_<name>.png`/`.txt` (the two
\|mean\| with their fits; fit-line and smoothed differences with the tolerance and the jackknife band; the residuals). With an
`rdip_ladder` run (CCEX writes `out_rdip<r>_state<s>_wiDiv` beside `out_state<s>_wiDiv`) each rung is judged against the
main output of the same runs with `--setting-prefix out_rdip<r>`. The setting and reference runs are paired by order but must be the same configuration: each run folder's
`input.json` must agree on every key except `rdip`, `rbath`, `rdip_ladder` and `output.file` (`outfile` in an input before
2.0.0) (bath, state and axis files, seed, physics; a key of the configuration or output block is named `configuration.<key>`,
`output.<key>`), otherwise `judge` refuses; `--unverified-pairs` accepts folders without `input.json` with a warning. A T2 fit with
T2 <= 0, n outside [0.3, 4] or non-finite values is a failure (`t2-unresolved`), not a number.

One `RUN_DIR` per bath configuration, each holding the CCEX outputs `out_state<s>_wiDiv` (the coherence; `_noDiv` is
not read). `--config INPUT` is a CCEX input with an `auto_convergence` block, read as the runner reads it:
the time window is judged unless `time_window.value` is given, the sampling unless `sampling.configurations` is;
`time_window.decay_level` (0.01), `time_window.production_factor` (1.5), `t2_fit.points` (peaks) and `sampling.threshold`
(the SE(A)/A target, 0.05) are taken from it, the expected states from `configuration.s_range`, the points from `nstep`, the
bath files available from `configuration.bathfile`. Without it those defaults are used. The error budget (tolerance 0.01
split into rdip/rbath/sampling), the area norm and judge's tolerances (0.01, 0.005) are the analysis's own, not input keys. The sampling error (proxy, bootstrap,
allocation) is averaged over [0, T_d] when a decay time is found, so a longer window with a near-zero tail does not make it
look smaller; the report says which points were evaluated (`evaluated_on`). The report is JSON on standard output.

Every `RUN_DIR` is one independent configuration: the same directory is refused twice, and so is a run with fewer states or
time points than expected or another time grid (`configuration.s_range`/`nstep`/`deltat` of the input with `--config`, or the options) —
a partial or stale run would otherwise look like a smaller complete sample. Without a manifest, identical outputs under two
names are refused and the report says `provenance: unverified`. A manifest (`{"configurations": [{"id", "run", "states"}],
"available_configurations": N, "grid": {"deltat", "nstep"}}`, run paths relative to it) gives each run its ID and its own
state count (a mixed 10/30-state pilot), turns identical bytes under different IDs into a warning, and bounds the allocation
by the unused configurations it lists. Every decision (window, sampling) uses the same balanced sample (the first
S_min states of each configuration); extra states only feed `within_sd_mean_all_states`. The report lists per configuration the non-finite
states, the unphysical states (non-finite or |L| > 1.01, U) and the largest |L|.

Status words: `numerical-failure` (a non-finite state, or a pooled mean, variance or bootstrap that overflows; the report
names the failed stage; nothing is estimated and no configuration is dropped), `resolved`, `needs-extension` (decay found but the window is shorter than
`window_margin × T_d`), `window-unresolved` (no decay inside the window). The sampling verdict `sampling.status` is
`t2-met` or `sampling-unresolved`: the larger of the linearised and bootstrap SE(A)/A against `sampling.threshold` (default
0.05; the report's names say t2, the quantity is the area A); `sampling-unresolved` also when the T2 fit fails. Beside it,
`sampling.proxy_status` (`proxy-met` / `sampling-unresolved`, `sampling_rms_proxy` of the present sample against
`sampling_tolerance`) is a diagnostic only. `total` is always `unresolved` in this version
and lists what is not measured (rdip, rbath, and sampling when it is off). None of them certifies accuracy: the envelope only sees the
computed window, the proxy is an expected-error scale and not a confidence limit (the bootstrap quantile is reported
next to it), and a few configurations give rough variance estimates. The allocation starts at S = 2 (the analysis needs two
states per configuration) and, for a fixed list of bath files, does not exceed their number: an unreachable target is
reported as infeasible with the best feasible allocation and its shortfall; fewer than 2 available configurations give
`source-limited` and no allocation; costs must be finite.

### `analyze --config` before 2.0.0

Before 2.0.0 `--config` was a driver file of its own (`{"base_input", "auto_convergence": {"tolerance", "steps", ...},
"bath_source", "state_source", "server"}`) with keys the runner refuses (`decay_threshold`, `window_margin`,
`t2_fit_points`, ...). It was removed: `--config` is the CCEX input itself.

### Example of `analyze` (existing CR0 runs of the study, 10 configurations × 30 states, 500 G, all-pairs rdip)

`needs-extension` (T_d 11.3 µs, recommended 17 µs, window 15 µs); `sampling-unresolved` (proxy 0.0153 for target 0.0067); with equal configuration and
state costs the cheapest allocation is S = 3, C = 192.
