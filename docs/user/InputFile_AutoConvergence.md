# Automatic convergence tags

An `auto_convergence` block makes `main.out -f input.json` choose the time window, rdip, rbath and the number of
configurations (and states) by itself, and then run the production. Every calculation it makes is an ordinary CCEX run
in its own folder under `workdir`. You prepare the bath files; almost everything else is found.
[Advanced 12: Automatic convergence](../../example/tutorial/advanced/12_auto_convergence/README.md) is a worked input.

It supports cCCE with `cluster.mode` standard or pcce (order 2 or more), gCCE (order 1 or more) and fixed bundles
(`cluster.mode` bundle without `cluster.adaptive`), on the CPU or the GPU, with DFT tensor files. The order is the
input's.

1. [Running it](#tag_ac_running): the command, and the smallest input.
2. [The bath configurations](#tag_ac_configurations): where the configurations and states come from.
3. [The steps](#tag_ac_steps): what the runner does, how it fits T2, and fixed values.
4. [Keys](#tag_ac_keys): every key of the block, by section.
5. [Results](#tag_ac_results): the outputs, the exit codes, continuing a job, and what a result does not say.
6. [Refused](#tag_ac_refused): what cannot go with the block.
7. [Commands](#tag_ac_commands): the runner's own commands.

<a id="tag_ac_running"></a>
## 1. Running it

```sh
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8
main.out -f input.json               # every CCEX run on one rank
```

Under `mpirun -n N` only rank 0 stays; it runs the Python runner, which starts every CCEX run with the `mpirun -n N`
on `PATH`. The runner needs Python 3.10 or newer with numpy, scipy and matplotlib; `CCEX_PYTHON` picks the
interpreter. The automatic run takes no option but `-f FILE`. `"enabled": false` in the block runs the ordinary
calculation instead.

Calling the same command again continues: finished runs are reused when their inputs, files and binary are unchanged
([When a job is cut](#tag_ac_resume)).

<a id="tag_ac_smallest"></a>
### The smallest input

```json
{
  "method": "cCCE", "quantity": "coherence", "order": 2,
  "Qubit": {"qubit": [{"file": "data/qubit"}]},
  "configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"]},
  "bfield": 500, "pulse": {"count": 1},
  "auto_convergence": {"workdir": "auto"}
}
```

`rbath`, `rdip` and `deltat` are found. `nstep` is the number of time points of every run (150 when absent); a
`deltat` given only sets the first window, (nstep − 1) × deltat.

<a id="tag_ac_configurations"></a>
## 2. The bath configurations

The configurations come from the [configuration](InputFile_General.md#tag_configuration) block:

- `bathfile`: one name with `{k}`; `bath_state_file`, `defect_axis_file` and `defect_state_file`, when given, one name
  each with `{k}` too; all in `configuration.path` when it is given.
- With `configuration.k_range` [first, last], those k are used, and each must have its files. Without it, every k whose
  bath file exists (`config_7`, not `config_07`).
- A state file name holds `{s}` for the state, or the number is appended: `data/cfg{k}_state{s}` and `data/cfg{k}_state`
  both read `data/cfg{k}_state1`, `..._state2`, .... A configuration with a missing state file is refused, not dropped.

`configuration.s_range` chooses the approach:

- **Without it:** the ensemble, one curve per configuration. Only the number of configurations is chosen.
- **With it:** sampled states. Its numbers are not used; `sampling.test_states` and `sampling.max_states` set the
  states.

Where the states come from:

- With `bath_state_file`, the state files of each configuration are counted; every configuration needs at least
  `test_states` of them.
- Without it, the runner makes the states of each run in the run folder, by the rule of
  [StateFiles](Tools.md#tag_tool_statefiles): each spin's m drawn uniformly from numpy `default_rng([seed, k, s])`, so a
  spin keeps its state at every rbath and rdip. The on-site states of a Defect species and the axes are made the same
  way. CCEX's own random states are not used, since they change with rbath.

`configuration.seed` is one integer, 1 when absent. It fixes the split of the configurations and the states the runner
makes. A seed list is refused.

<a id="tag_ac_steps"></a>
## 3. The steps

| Step | What it does | Figures and report |
|---|---|---|
| 0. start | finds and splits the configurations, estimates rough starts for rbath, the decay time and rdip | `summary.txt` section 0 |
| 1. time window | finds the decay time $T_d$ and the window from a T2 fit of the test runs | `1_window.png`, `1_window_fit.png` |
| 2. rdip | chooses rdip by comparing the fitted areas with the largest rdip | `2_rdip.png` |
| 3. rbath | chooses rbath the same way, then checks rdip again at the final rbath | `3_rbath.png` |
| 3b. window again | measures $T_d$ again at the final setting and repeats steps 1–3 when it moved | `summary.txt` section 1 |
| 4. sampling | checks the number of states and finds the number of configurations | `4_states.png` |
| 5. production | runs the chosen setting on the production configurations | `5_configurations.png`, `6_converged_coherence.png`, `coherence_mean.txt` |

**0. Start.** `sampling.test_configurations` (default 10) configurations go to the convergence tests, chosen evenly
over the distance from the qubit to each configuration's nearest bath spin; the rest are the production, in a seeded
order. From the spins counted around the qubit the runner estimates first rbath and rdip candidates and a decay time.
rbath never goes above the radius the bath files cover (the distance from the qubit to the nearest face of their box;
for a layer, the in-plane distance).

**1. Time window.** The test runs go to `revival_check_factor` × the decay time. A T2 fit $\exp[-(t/T_2)^n]$ of their
mean gives $T_d$, where the fit reaches `decay_level`. The window grows until it covers the decay, or is cut when far
too long. A revival extends the window; a plateau (the curve decays, then stays above 0.05) raises rbath and the
largest rdip. nstep is never changed, and whether it suffices is not checked.

**2. rdip.** One run at the largest rdip gives every smaller one ([rdip_ladder](#tag_ac_rdip_ladder), order 2 of
cCCE, gCCE and pCCE; otherwise one run per value). Two values agree when the area under their T2 fits,
$A = T_2 \Gamma(1+1/n)$, differs by at most `rdip.threshold`. The largest is raised by × 1.25 (up to 2 × rbath) until
the last comparisons agree; the chosen rdip is the smallest from which every larger one agrees with the largest.

**3. rbath.** The same rule with `rbath.threshold`, one run per value. Then rdip is checked again at the final rbath,
and converged again there if needed.

**3b. Window again.** Only when the window, rdip and rbath are all found: when the final setting is well above the one
the window was found at, $T_d$ is measured there, and steps 1–3 run again when it moved by more than × 1.25 (at most 3
rounds).

**4. Sampling.** With `s_range`, the states' share of the sampling error of A is measured at 2 ... `test_states`
states; above `sampling.states_threshold`, more states are computed, up to `max_states`. Then the variance of A
between configurations and between states gives the number of configurations for SE(A)/A ≤ `sampling.threshold`, at
least `min_configurations`.

**5. Production.** The chosen setting runs on the production configurations over `production_factor` × $T_d$: first
1.5 × the estimate of the configurations needed (at least max(20, `min_configurations`)), then more while the SE is short of the target or until
the files run out. The production curve must decay inside its window; otherwise it is computed again on a longer
window (at most 3 times). The chosen setting is compared once with the largest rbath and rdip on 3 production
configurations, with a warning when they differ by more than the threshold.

<a id="tag_ac_fits"></a>
### T2 fits and divergent points

With `t2_fit.points` `peaks` (the default) every T2 fit uses t = 0 and the local maxima of its curve, up to the first
point below `t2_fit.floor`; a curve with fewer than 3 peaks is fitted on all its points. The T2 is a measure to compare
settings, not the exact T2 of the system. When either fit of a comparison fails, the comparison uses the raw area
$\int |\bar L|\,dt$ over [0, $T_d$] with the same threshold; `summary.txt` marks those comparisons.

Wherever a curve is read, a point with |wiDiv − noDiv| ≥ 2 or a nan/inf wiDiv takes the noDiv value. (The analyzer of
[Analysis](Analysis.md) uses another rule, so its numbers differ.) A configuration that still has |L| > 10 or nan/inf
is left out of every average; more than 5 % of the production left out is a numerical failure.

<a id="tag_ac_fixed"></a>
### Fixed values

A section is found unless its value is given:

| Given | Effect |
|---|---|
| `time_window.value` | the window is used as it is; $T_d$ is only measured |
| `rdip.value` | the rdip step and its check at the final rbath are skipped |
| `rbath.value` | the rbath step is skipped; a plateau stops instead of raising rbath |
| `sampling.configurations` | the sampling step is skipped; the production uses exactly that many configurations (and `sampling.states` states each) |

Step 3b runs only when the window, rdip and rbath are all found. All four given is refused: run CCEX itself.

<a id="tag_ac_keys"></a>
## 4. Keys

A section left out uses its defaults, and a key given wins over what would be found. `null` or an empty list is
refused; leave the key out to have it found. Unknown keys and old names are refused with a message that names the
current key ([Old names](InputFile_OldNames.md)). Keys are read without regard to case.

### The block itself

| Key | Type | Default | Meaning |
|---|---|---|---|
| `workdir` | string | required | the folder of the runs (`runs/`), the decisions and the results (`result/`) |
| `enabled` | boolean | `true` | `false`: `main.out -f` runs the ordinary calculation and the block is not read |
| `resume` | string | `"continue"` | `"continue"` reuses what the workdir holds; `"fresh"` moves it to `<workdir>.old<n>` and starts over |

### time_window

| Key | Type | Default | Meaning |
|---|---|---|---|
| `value` | number (ms), or `{"value", "unit"}` with `s`, `ms`, `us`, `ns` | found | a fixed window |
| `decay_level` | number | 0.01 | $T_d$ is where the T2 fit reaches this \|L\| |
| `revival_check_factor` | number | 2.5 | the test runs go to this × $T_d$, to see a revival |
| `production_factor` | number (≥ 1) | 1.5 | the production window is this × $T_d$ |

### rdip and rbath

| Key | Type | Default (rdip / rbath) | Meaning |
|---|---|---|---|
| `value` | number (Å) | found | a fixed value |
| `candidates` | list of numbers (Å) | found | the values to try; the largest is the reference |
| `num_candidates` | integer | 7 / 4 | how many to start with when found |
| `threshold` | number | 0.02 / 0.02 | the relative difference of the fitted area A at which two values agree |

### t2_fit

| Key | Type | Default | Meaning |
|---|---|---|---|
| `points` | string | `"peaks"` | the points of every T2 fit: `"peaks"` or `"all"` |
| `floor` | number | 0.05 | only the points before \|mean\| falls below this are fitted |
| `n_range` | `[min, max]` | [0.3, 4.0] | the allowed exponent n; outside it the fit fails |

### sampling

| Key | Type | Default | Meaning |
|---|---|---|---|
| `test`, `production` | lists of k | found | the configurations of the tests and of the production (never the same) |
| `configurations`, `states` | integers | found | a fixed production size ([Fixed values](#tag_ac_fixed)) |
| `threshold` | number | 0.05 | the SE(A)/A target |
| `test_configurations`, `test_states` | integers | 10, 5 | the configurations and states of the test runs |
| `states_threshold` | number | 0.05 | the states' share of the sampling error at which the state test stops |
| `max_states` | integer | 30 | the most states the state test computes (with state files, at most the files there) |
| `min_configurations` | integer | 10 | the fewest production configurations |
| `max_configurations` | integer | every production file | the most production configurations; raising it later continues the same run |

<a id="tag_ac_example"></a>
### Example

The window and rbath fixed, rdip and the configurations found, a tighter target:

```json
"configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"], "k_range": [1, 60]},
"auto_convergence": {
  "workdir": "auto",
  "time_window": {"value": {"value": 15, "unit": "us"}},
  "rbath": {"value": 288},
  "sampling": {"threshold": 0.03}
}
```

<a id="tag_ac_results"></a>
## 5. Results

Everything is under `workdir`: `runs/` (one folder per CCEX run), `decisions.jsonl` (one line per decision) and
`result/`. The files are listed in [Output files](OutputFiles.md#out_autoconv).

`result/summary.txt` (and `result/json/summary.json`) gives the status, which keys were found, the window, rbath and
rdip, the production (configurations, states, the area A with its standard error, T2, n), the core-hours and every
warning. It says `running` while the run goes on.

| Exit code | Meaning |
|---|---|
| 0 | converged |
| 1 | error |
| 2 | refused |
| 3 | unresolved: `window-unresolved`, `truncation-unresolved`, `sampling-unresolved` or `numerical-failure` |
| 4 | stopped (walltime or a signal); call again to continue |

When the input, the files or the binary change, the earlier decisions are set aside (`decisions.<n>.jsonl`) and made
again; the CCEX runs are reused.

<a id="tag_ac_resume"></a>
### When a job is cut

Calling the same command again continues; finished runs and decisions are reused. Under Slurm a run is not started
with less walltime left than 1.5 × the longest run of its step + 300 s (exit 4). A job cut from outside is noticed on
the next call and continued. `"resume": "fresh"` starts over and keeps the earlier workdir as `<workdir>.old<n>`.

`result/progress.txt` is rewritten after every run: the step now, the steps done and left, and the values chosen so
far. The [status command](#tag_ac_commands) prints it without computing anything.

<a id="tag_ac_baths"></a>
### When the baths are too small or too few

The runner does not make baths; it says what to make.

- **Too small** (a plateau, or rbath not converged at the radius the files cover): `truncation-unresolved` or
  `window-unresolved`, with "make larger baths, covering at least R Å around the qubit". Larger files under the same
  names start a new plan.
- **Too few** (the SE is above the target with every configuration used, or fewer than max(20, `min_configurations`)
  exist): `sampling-unresolved` with how many configurations to add. Files added under the same pattern are put after
  the production, and the next call goes on from there. A `sampling.production` list given by hand is not extended.

<a id="tag_ac_limits"></a>
### What a result does not say

- The chosen rdip and rbath are each within their threshold of the largest computed; the combination is checked once,
  and an infinite bath not at all.
- The SE is an estimate from the configurations computed; with 10 configurations the true value may be 0.7 to 1.8 ×
  the estimate.
- The CCE order and whether the time step is fine enough are not checked.
- Curves that are not one stretched exponential (a collapse and revival, spiky or noisy curves) may stop at the time
  window or not converge in rbath. P1 with pCCE2 and h-BN at 500 G are such cases.

<a id="tag_ac_refused"></a>
## 6. Refused

Together with the block:

- other methods, `cluster.mode` dsj, itb or dsjitb, a `cluster.adaptive` block, a quantity other than coherence, more
  than one qubit;
- `rdip_ladder`, `cluster.pcce.cache` (the runner sets its own), `cluster.pcce.partition` legacy,
  `output.cluster_contributions`, `output.cluster_only`, a nonzero `cluster.selection.max_per_order`;
- `output.path` and `output.file` (every run writes into its own folder);
- a nonzero `configuration.bath_shift`, more than one bath file per configuration, a `bathfile` without `{k}`, a seed
  list;
- `rdip.value` with candidates, `sampling.states` without `sampling.configurations`, values out of range, unknown keys
  and old names.

<a id="tag_ac_pcce_cache"></a>
### The pCCE partition cache

Every pCCE run gets `cluster.pcce.cache` `readwrite` with a file in `<workdir>/pcce_cache/`, keyed on the input without
rdip, the states and the time grid, on the files and on the binary. Runs that differ only in those search the k-means
partition once. The outputs are the same as without the cache.

<a id="tag_ac_rdip_ladder"></a>
### rdip_ladder

An internal CCEX key the runner writes into its runs: a list of up to 64 smaller rdip values computed from one order-2
run, written as `<output.file>_rdip<r>_...` ([rdip_ladder](InputFile_General.md#tag_rdip_ladder)). An input given to
the automatic convergence may not hold it.

<a id="tag_ac_commands"></a>
## 7. Commands

The runner's own commands compute nothing. Run them from `example/tools/AutoConvergence`, or with that folder on
`PYTHONPATH`:

```sh
python3 -m autoconv status INPUT.json                       # where a run is (the workdir is read from the input)
python3 -m autoconv analyze RUN_DIR [RUN_DIR ...] [--config INPUT] [--report-dir DIR]
python3 -m autoconv judge --setting RUN_DIR ... --reference RUN_DIR ... [--report-dir DIR]
```

- `status` prints `running`, `interrupted`, `stopped` or how the run ended, then `progress.txt` (and `summary.txt`
  once it ended). It needs no mpirun.
- `analyze` judges existing runs, one `RUN_DIR` per configuration with the outputs `out_state<s>_wiDiv`. It reports, as
  JSON, the decay time and whether the window is long enough, and the sampling verdict with the cheapest allocation of
  configurations and states. A directory given twice, or a run with fewer states or time points than expected, is
  refused.
- `judge` compares runs of one setting (a smaller rdip or rbath) with reference runs of the same configurations, as
  the rdip and rbath steps do. The pairs must agree in every input key but `rdip`, `rbath`, `rdip_ladder` and
  `output.file`.
- `--config INPUT` is a CCEX input with an `auto_convergence` block, read as the runner reads it; its `decay_level`,
  `production_factor`, `t2_fit.points` and `sampling.threshold` are used. `--report-dir DIR` writes figures and text
  summaries. `-h` lists every option.
