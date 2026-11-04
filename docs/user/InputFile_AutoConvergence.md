# Automatic convergence tags

An `auto_convergence` block makes `main.out -f input.json` choose the time window, rdip, rbath and the number of
configurations (and states) by itself, and then run the production. You prepare the bath files, and the rest is found.
How the runner works is in [Automatic convergence](../methods/auto-convergence.md), and
[Advanced 12: Automatic convergence](../../example/tutorial/advanced/12_auto_convergence/README.md) runs it.

It supports cCCE with `cluster.mode` standard or pcce (order 2 or more), gCCE (order 1 or more) and fixed bundles, on
the CPU or the GPU. The order is the input's. The run needs one qubit, `quantity` coherence, and one bath file per
configuration, named with `{k}` and without a `configuration.bath_shift`. Other combinations are refused with a
message.

1. [The smallest input](#tag_ac_running): what the block needs.
2. [The bath configurations](#tag_ac_configurations): where the configurations and states come from.
3. [Keys](#tag_ac_keys): every key of the block, by section.

<a id="tag_ac_running"></a>
<a id="tag_ac_smallest"></a>
## 1. The smallest input

```json
{
  "method": "cCCE", "quantity": "coherence", "order": 2,
  "Qubit": {"qubit": [{"file": "data/qubit"}]},
  "configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"]},
  "bfield": 500, "pulse": {"count": 1},
  "auto_convergence": {"workdir": "auto"}
}
```

`rbath`, `rdip` and `deltat` are found. `nstep` is the number of time points of every run (150 when absent).

<a id="tag_ac_refused"></a>
<a id="tag_ac_pcce_cache"></a>
<a id="tag_ac_rdip_ladder"></a>
The runner sets `cluster.pcce.cache` and [rdip_ladder](InputFile_General.md#tag_rdip_ladder) in its runs itself, so an
input given to it may not hold them, nor `output.path` or `output.file`.

<a id="tag_ac_configurations"></a>
## 2. The bath configurations

The configurations come from the [configuration](InputFile_General.md#tag_configuration) block. The name in
`bathfile`, and any state or axis file name, holds `{k}`. With `k_range` those k are used, and without it every k whose
bath file exists. `configuration.seed` is one integer, 1 when it is left out.

- **Without `s_range`:** each configuration gives one curve of the ensemble, and only the number of configurations is
  chosen.
- **With `s_range`:** each configuration gives sampled states, and `sampling.test_states` and `sampling.max_states`
  set how many. The states are read from `bath_state_file` when it is given. Otherwise they are made by the rule of
  [StateFiles](Tools.md#tag_tool_statefiles), so a spin keeps its state at every rbath and rdip.

<a id="tag_ac_keys"></a>
## 3. Keys

The block has five sections. A section that is left out uses its defaults, and a value that is left out is found by
the runner. `null` is refused.

A value that is given is used as it is, and its convergence is not tested. These values are `time_window.value`,
`rdip.value`, `rbath.value`, and `sampling.configurations` with `sampling.states`. For example, with `rbath.value`
given, the rbath step is skipped and every run uses that rbath. Giving all four is refused, since nothing is then left
to converge. Run CCEX itself instead.

<a id="ac_keys_block"></a>
### The block itself

| Key | Type | Default | Meaning |
|---|---|---|---|
| `workdir` | string | required | The folder of the runs and the results. |
| `enabled` | boolean | `true` | With `false` the input runs as an ordinary calculation. |
| `resume` | `"continue"` or `"fresh"` | `"continue"` | `"continue"` goes on from the runs already in `workdir`. `"fresh"` keeps the old folder as `<workdir>.old<n>` and starts over. |

<a id="ac_keys_window"></a>
### time_window

| Key | Type | Default | Meaning |
|---|---|---|---|
| `value` | number (ms), or `{"value", "unit"}` | found | Fixes the window instead of finding it. |
| `decay_level` | number | 0.01 | The decay time is where the T2 fit reaches this \|L\|. |
| `revival_check_factor` | number, 0 or more | 2.5 | The test runs reach this many times the decay time (or `production_factor` times, if larger), so that a revival can be seen. `0` skips the revival check. |
| `production_factor` | number, 1 or more | 1.5 | The production window is this many times the decay time. |

<a id="ac_keys_radii"></a>
### rdip and rbath

The two sections take the same keys.

| Key | Type | Default (rdip / rbath) | Meaning |
|---|---|---|---|
| `value` | number (Å) | found | Fixes the value instead of finding it. |
| `candidates` | list of numbers (Å) | found | The values to try. The largest is the reference. |
| `num_candidates` | integer | 7 / 4 | How many candidates are tried at first, when `candidates` is not given. |
| `threshold` | number | 0.02 / 0.02 | Two values agree when their areas A differ by at most this fraction. |

<a id="ac_keys_fit"></a>
### t2_fit

How the T2 fit that every step compares is made.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `points` | `"peaks"` or `"all"` | `"peaks"` | Fit the local maxima, or every point. |
| `floor` | number | 0.05 | Only the points before \|L\| falls below this are fitted. |
| `n_range` | `[min, max]` | `[0.3, 4.0]` | The range allowed for the exponent n. |

<a id="ac_keys_sampling"></a>
### sampling

How many configurations (and states) the production uses.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `test`, `production` | lists of k | found | The configurations of the test runs and of the production. |
| `configurations`, `states` | integers | found | Fix the size of the production instead of finding it. |
| `threshold` | number | 0.05 | The target relative error SE(A)/A of the production. |
| `test_configurations`, `test_states` | integers | 10, 5 | The number of configurations and states of the test runs. |
| `states_threshold` | number | 0.05 | The state test adds states until their share of the error falls below this. The test only warns, and the production keeps `test_states` states. |
| `max_states` | integer | 30 | The most states computed per configuration. |
| `min_configurations` | integer | 10 | The fewest production configurations. |
| `max_configurations` | integer | all production configurations | The most production configurations. |

<a id="tag_ac_example"></a>
An example that fixes the window and rbath, lets the runner find rdip and the configurations, and asks for a tighter
target:

```json
"configuration": {"gyrofile": "data/gyro_13C", "bathfile": ["data/config_{k}"], "k_range": [1, 60]},
"auto_convergence": {
  "workdir": "auto",
  "time_window": {"value": {"value": 15, "unit": "us"}},
  "rbath": {"value": 288},
  "sampling": {"threshold": 0.03}
}
```
