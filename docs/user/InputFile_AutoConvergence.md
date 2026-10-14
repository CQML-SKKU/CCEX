# Automatic convergence tags

An `auto_convergence` block in an ordinary input makes `main.out -f input.json` choose the time
window, rdip, rbath and the number of configurations (and states) of a coherence calculation by
itself, then run the production. Every calculation it makes is an ordinary CCEX run in its own
folder under `workdir`. Supported: cCCE with `cluster_mode` standard or pcce (order ≥ 2), gCCE (order ≥ 1), fixed bundles
(`cluster_mode` bundle without an `adaptive` block), CPU or GPU backends. How each step decides,
the outputs and the exit codes are in `example/tools/AutoConvergence/README.md`; tutorial 11
([Examples_11.md](Examples_11.md)) is a worked input. The runner reads the keys (of the block and of the rest of the
input) without regard to case, as CCEX does, and refuses two keys of one block that differ only in case (since review
20261013-180049 R4; `"Auto_Convergence"` was refused as missing and `"Method": "gCCE"` ran as cCCE).

```
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8 (rank 0 alone runs the runner)
main.out -f input.json               # with the block's launcher, else one rank
main.out --auto_status -f input.json # where a run is (nothing computed)
```

The automatic run takes only `-f FILE`; other command-line options are refused (set them in the
input). `nstate` of the input chooses the approach: `0` is the ensemble (only the number of
configurations is chosen), `> 0` samples states (`configurations.test_states` and `max_states`
set how many). `rbath`, `rdip`, `nstep` and `deltat` of the input are not needed; `nstep` is the
number of time points of every run (150 when absent).

The block is grouped by step. A section left out uses its defaults; a key given wins over what
would be found; `null` or an empty list is refused (leave the key out to have it found).
Unknown keys and the earlier flat names (`threshold_rdip`, `candi_rdip`, `se_target`, ...) are
refused with a message that says where each went.

<a id="tag_ac_fixed"></a>
## Fixed values instead of convergence

Every section converges by default. `"converge": false` uses the value given instead, and that
value is then required:

| Section with `"converge": false` | value used | effect |
|---|---|---|
| `time_window` | `value` (+ `unit`) | the window is used as it is; T_d is only measured |
| `rdip` | `value` (Å) | the rdip step and its check at the final rbath are skipped |
| `rbath` | `value` (Å) | the rbath step is skipped; a plateau stops instead of raising rbath |
| `configurations` | `number` (and `states` when nstate > 0) | the sampling step is skipped; the production uses exactly that many |

All four `false` is refused (run CCEX itself).

<a id="tag_ac_keys"></a>
## Keys

| Section | Key | Meaning | Default |
|---|---|---|---|
| (top) | `enabled` | false: `main.out -f` runs the ordinary calculation and the block is not read | true |
| | `resume` | `continue`: reuse what the workdir holds; `fresh`: move it to `<workdir>.old<n>` and start over | `continue` |
| | `workdir` | folder for the runs (`runs/`), the decisions (`decisions.jsonl`) and the results (`result/`) | required |
| | `launcher` | how each CCEX run is started, e.g. `["mpirun", "-n", "4"]` or `["srun"]` | `mpirun -n N` when main.out was started with `mpirun -n N`, else one rank |
| `time_window` | `converge`, `value`, `unit` | see above; units `s`, `ms`, `us`, `ns` | true, —, `ms` |
| | `threshold` | T_d is where the T2 fit reaches this level | 0.01 |
| | `revival_check_factor` | the test runs go to this × T_d to see a revival | 2.5 |
| | `production_factor` | production window = this × T_d | 1.5 |
| `rdip` | `converge`, `value` | see above | true |
| | `candidates` | rdip candidates (Å); the largest is the reference | found |
| | `num_candidates` | how many to start with when found | 7 |
| | `threshold` | relative difference of the fitted area A = T2 Γ(1+1/n) at which two rdip agree | 0.02 |
| `rbath` | `converge`, `value` | see above | true |
| | `candidates` | rbath candidates (Å) | found |
| | `num_candidates` | how many to start with when found | 4 |
| | `threshold` | relative difference of A at which two rbath agree | 0.02 |
| `t2_fit` | `points` | points of every T2 fit: `peaks` (t = 0 and the local maxima; all points when fewer than 3) or `all` | `peaks` |
| | `floor` | only the points before \|mean\| falls below this are fitted | 0.05 |
| | `n_range` | allowed stretching exponent n of exp(−(t/T2)^n); outside it the fit fails | [0.3, 4.0] |
| `configurations` | `bathfile` | bath file pattern with `{k}` (a non-negative integer; `config_7`, not `config_07`) | required |
| | `statefile` | state file prefix pattern (nstate > 0), written into each run's input as `bath_state_file`; left out, the states are made per run and removed | — |
| | `avaaxfile`, `exstatefile` | optional per-configuration file patterns, written into each run's input as `Defect.axis_file` and `Defect.state_file` | — |

The `configurations` keys keep these names (the tool's own pattern names); they were not renamed with the input keys on
2026-10-13.
| | `test`, `production` | configuration IDs for the convergence tests and for the production (never the same) | found |
| | `converge`, `number`, `states` | see above | true, —, — |
| | `threshold` | SE(A)/A target of the sampling | 0.05 |
| | `test_number`, `test_states` | configurations and states of the test runs (states for nstate > 0) | 10, 5 |
| | `states_threshold` | the states' share of the sampling error at which the state test stops | 0.05 |
| | `max_states` | the most states the state test computes (with state files, at most the files there) | 30 |
| | `min_configurations` | the fewest production configurations | 10 |
| | `max_configurations` | the most production configurations used | every production file |

<a id="tag_ac_refused"></a>
## Refused together with the block

Other cluster modes (`dsj`, `itb`, `dsjitb`, the adaptive refinement) or quantities than coherence, `rdip_ladder` in the
input, an `adaptive` block, `pcce.cache` (the runner gives its own), `pcce.partition` legacy, `cluster_contributions`,
a nonzero `max_clusters_per_order`, `Defect.axis_file`/`Defect.state_file` outside `configurations`, a nonzero `bath_shift`
and more than one qubit, a list together with its count, a fixed value with `converge` true, and
options out of range.

<a id="tag_ac_example"></a>
## Example

```
"auto_convergence": {
  "workdir": "auto",
  "configurations": {"bathfile": "data/config_{k}"},
  "time_window": {"converge": false, "value": 15, "unit": "us"},
  "rbath": {"converge": false, "value": 288}
}
```
