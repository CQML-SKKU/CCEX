# Pulse tags

The `pulse` block sets the pulse sequence of the qubit and, if needed, DEER pulses on bath spins. The block and its
`count` are required. Old key names are in [Upgrading old inputs](InputFile_OldNames.md).

1. [The sequence](#pul_qubit): `count`, `name` and `sequence`.
2. [Finite and per-qubit pulses](#pul_finite): `duration`, `detuning_factor` and `per_qubit`.
3. [DEER bath pulses](#tag_pulse_bath): the `pulse.bath` block.

<a id="pul_qubit"></a>
## 1. The sequence

<a id="tag_pulse"></a>
<a id="tag_pulse_count"></a>
<a id="tag_pulse_name"></a>
Most runs need only `count`. A sequence is chosen like this:

| Sequence | Written as | Pulses |
|---|---|---|
| Ramsey (free induction) | `"count": 0` | none |
| Hahn echo | `"count": 1` | one π pulse at the middle |
| CPMG | `"count": n`, with n > 1 | n equally spaced π pulses |
| WAHUHA | `"count": 4, "name": "WAHUHA"` | X 90°, Y 270°, Y 90°, X 270° (gCCE or megCCE) |
| any other | `"name": "Manual"` and [sequence](#tag_sequence) | as written, `count` of them |

| Key | Type | Default | Meaning |
|---|---|---|---|
| `count` | integer | required | The number of pulses. |
| `name` | `"Ramsey"`, `"HahnEcho"`, `"CPMG"`, `"WAHUHA"` or `"Manual"` | chosen from `count` | The name of the sequence. A written name must fit `count`. |

cCCE takes ideal π flips about X or Y only, and other angles need gCCE. meCCE takes its sequence from `count` only.

```json
"pulse": {"count": 2}
"pulse": {"count": 4, "name": "WAHUHA"}
```

<a id="tag_sequence"></a>
### pulse.sequence

| Type | Default |
|---|---|
| list of `[fraction, axis, angle]` | the sequence built from `count` |

The sequence written out, with `name` `Manual`. It has `count` + 1 entries. Each entry is the free evolution before a
pulse, as a fraction of the total time (`"1/4"` or a number), then the axis of the pulse (`X`, `Y`, `Z`, or `I` for
none) and its angle in degrees. The fractions add up to 1, and the last entry is the evolution after the last pulse.
[Basic 04](../../example/tutorial/basic/04_pulses/README.md) writes a CPMG-2 this way.

```json
"pulse": {"count": 2, "name": "Manual", "sequence": [["1/4", "X", 180], ["1/2", "X", 180], ["1/4", "I", 0]]}
```

<a id="pul_finite"></a>
## 2. Finite and per-qubit pulses

<a id="tag_pulse_duration"></a>
<a id="tag_detuning_factor"></a>
<a id="tag_pulse_per_qubit"></a>
By default every pulse is ideal and instantaneous. gCCE can also take finite pulses, and pulses that address each
qubit separately. cCCE, meCCE and megCCE take ideal pulses only.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `detuning_factor` | number | `1.0` | `1.0` gives ideal pulses. Another value drives each pulse for `duration` at this factor times the qubit frequency. |
| `duration` | number in ns, or `{"value", "unit"}` | `0` | The length of a finite pulse, read with `detuning_factor` other than 1. |
| `per_qubit` | boolean | `false` | With `true` each qubit gets its own pulse between its two levels of [alpha_state and beta_state](InputFile_Qubit.md#tag_alpha_state), which must be given as `{"ms": [...]}`. Only gCCE reads it, and meCCE, megCCE, bath pulses, noise and pseudospin refuse it. |

```json
"pulse": {"count": 1, "duration": 10.0, "detuning_factor": 0.98}
```

<a id="tag_pulse_bath"></a>
## 3. DEER bath pulses

| Type | Default |
|---|---|
| object `pulse.bath` | none (no bath pulse) |

Pump pulses on the bath spins of one species during the qubit's sequence, for gCCE or megCCE with `cluster.mode`
standard or pcce, on the CPU. How they work is in [DEER bath pulses](../methods/deer.md).
[Advanced 05](../../example/tutorial/advanced/05_P1/README.md) pumps the central P1 line.

```json
"pulse": {"count": 1, "bath": {"count": 1, "defect": "P1", "energy_shift": [0.0],
                               "sequence": [["0.50", "X", 180.0], ["0.50", "I", 0.0]]}}
```

<a id="tag_pulse_bath_count"></a>
<a id="tag_pulse_bath_sequence"></a>
<a id="tag_pulse_bath_defect"></a>
<a id="tag_pulse_bath_energy_shift"></a>
<a id="tag_pulse_bath_bandwidth"></a>
<a id="tag_pulse_bath_levels"></a>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `count` | integer, 0 or more | `0` | The number of bath pulses. `0` turns the block off. |
| `sequence` | list of `[fraction, axis, angle]` | required when `count` > 0 | The bath pulses, written like [pulse.sequence](#tag_sequence) on the same time axis. It has `count` + 1 entries, the last with no pulse (`I`). |
| `defect` | string (case is not significant) | required when `count` > 0 | The name of the bath spins the pump can drive, as in the bath file, or the name of a Defect sub-spin. |
| `energy_shift` | number or list in MHz, or `{"value", "unit"}` | required when `count` > 0 | The pump tones. A spin is driven when a tone is within `bandwidth` of its [Defect species detuning](InputFile_Defect.md#tag_df_detuning). |
| `bandwidth` | positive number in MHz, or `{"value", "unit"}` | 10 rad/ms (about 1.59 kHz) | The half width of the window of a tone. |
| `levels` | list of two numbers | `[0.5, -0.5]` | The two levels `[m_a, m_b]` the pulse rotates between. All the driven spins have the same S. |

```json
"bath": {"count": 1, "defect": "P1", "energy_shift": [0.0], "bandwidth": {"value": 5, "unit": "kHz"}}
```
