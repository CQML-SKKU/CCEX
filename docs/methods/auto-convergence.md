# Automatic convergence
<a id="AutoConvergenceMethod"></a>

With an `auto_convergence` block, `main.out -f input.json` does not run one calculation. Rank 0 starts the runner,
`example/tools/AutoConvergence`, which needs Python 3.10 or newer with numpy, scipy and matplotlib (`CCEX_PYTHON`
picks the interpreter). The runner chooses the time window, rdip, rbath and the number of configurations (and states), and then
runs the production. Every calculation it makes is an ordinary CCEX run in its own folder under `workdir`, with the
`mpirun -n N` that `main.out` was started with. The keys are in
[Automatic convergence tags](../user/InputFile_AutoConvergence.md), and
[Advanced 12](../../example/tutorial/advanced/12_auto_convergence/README.md) runs it.

1. [The steps](#tag_ac_steps): what the runner does, step by step.
2. [When two settings agree](#tag_ac_fits): the comparison every step uses.
3. [Results](#tag_ac_results): the summary, the exit codes, continuing a job, and the limits.

<a id="tag_ac_steps"></a>
## 1. The steps

| Step | What it does | Figure |
|---|---|---|
| 0. start | Splits the configurations into test and production, and estimates first values. | |
| 1. time window | Finds the decay time from a T2 fit of the test runs. | `1_window.png` |
| 2. rdip | Takes the smallest rdip from which every larger one agrees with the largest. When the two largest do not agree, the largest is raised first. | `2_rdip.png` |
| 3. rbath | Does the same for rbath, then checks rdip again. | `3_rbath.png` |
| 4. sampling | Chooses the number of configurations for the target error. With sampled states it also checks whether `test_states` states are enough, and warns when they are not. The production keeps `test_states` states. | `4_states.png` |
| 5. production | Runs the chosen setting on the production configurations, and compares the chosen rbath and rdip with the largest on 3 of them. | `5_configurations.png`, `6_converged_coherence.png` |

<a id="tag_ac_fixed"></a>
A section whose value is given is not searched. These values are `time_window.value`, `rdip.value`, `rbath.value` and
`sampling.configurations` (with `sampling.states`). Giving all four is refused, since then CCEX itself does the run.

<a id="tag_ac_fits"></a>
## 2. When two settings agree

Two settings agree when the areas under their T2 fits, $A = T_2 \Gamma(1+1/n)$, differ by at most the section's
`threshold`. When a fit fails, the raw areas under the curves are compared instead. The T2 is a measure to compare
settings, not the exact T2 of the system.

<a id="tag_ac_results"></a>
## 3. Results

Everything is under `workdir` ([Output files](../user/OutputFiles.md#out_autoconv)). `result/summary.txt` gives the
status, the values found, the production result (the area A with its standard error, T2, n), the core-hours and every
warning.

| Exit code | Meaning |
|---|---|
| 0 | Converged. |
| 1 | Error. |
| 2 | Refused. |
| 3 | Unresolved (`window-unresolved`, `truncation-unresolved`, `sampling-unresolved` or `numerical-failure`). |
| 4 | Stopped (walltime or a signal). Call again to continue. |

<a id="tag_ac_resume"></a>
**Continuing.** Calling the same command again reuses the finished runs and decisions. A job cut by walltime or from
outside goes on from there. `result/progress.txt` says where it is.

<a id="tag_ac_baths"></a>
**Too small or too few baths.** The runner does not make baths. It says what to make, either larger baths ("covering at
least R Å around the qubit") or more configurations. Files added under the same names are used on the next call.

<a id="tag_ac_limits"></a>
**What it does not check.** The runner does not check the CCE order, the time step or the limit of an infinite bath.
Curves that are not one stretched exponential (a collapse and revival, spiky curves) may not converge.
