# Pulse tags

The `pulse` block sets the qubit's pulse sequence and, optionally, DEER pulses on the bath. The block and its `count`
are required. Each key has a small table (its type, default and unit), what it does, and an example. Old key names
are in [Old names](InputFile_OldNames.md).

1. [The qubit sequence](#pul_qubit): `count`, `name`, `sequence`, `per_qubit`, `duration`, `detuning_factor`.
2. [DEER bath pulses](#tag_pulse_bath): the `pulse.bath` block.

<a id="pul_qubit"></a>
## 1. The qubit sequence

<a id="tag_pulse"></a>
### pulse

| Type | Default |
|---|---|
| object | required |

The pulse sequence of the run. Without the block the run stops (`pulse is not found in the input file`). Limits by
method:

- cCCE: `detuning_factor` 1 and `duration` 0 (ideal instantaneous pulses), and every pulse is an ideal π flip about X
  or Y; other angles need gCCE.
- meCCE and megCCE: `per_qubit` false, `detuning_factor` 1 and `duration` 0. meCCE takes its pulses from `count` only;
  a `Manual` sequence is refused.

```json
"pulse": {"count": 2, "name": "Manual", "sequence": [["1/4", "X", 180], ["1/2", "X", 180], ["1/4", "I", 0]]}
```

This is a CPMG-2: wait t/4, π about X, wait t/2, π about X, wait t/4.

<a id="tag_pulse_count"></a>
### pulse.count

| Type | Default |
|---|---|
| integer | required |

The number of pulses. Without [sequence](#tag_sequence) the run builds the sequence from it:

| Value | Sequence |
|---|---|
| `0` | Ramsey (free induction) |
| `1` | Hahn echo |
| n > 1 | CPMG with n equally spaced π pulses |

With [name](#tag_pulse_name) `Manual` it is the number of entries of `sequence` minus one.

```json
"pulse": {"count": 1}
```

<a id="tag_pulse_name"></a>
### pulse.name

| Type | Default |
|---|---|
| string | chosen from `count` |

The name of the sequence. A written name must fit `count`, otherwise the run stops with
`pulse.name X needs pulse.count Y`. Without it the run takes `Ramsey`, `HahnEcho` or `CPMG` from `count` and says so
in the log, which also prints the sequence as a table of intervals, axes and angles.

| Value | Sequence | count |
|---|---|---|
| `Ramsey` | free induction | 0 |
| `HahnEcho` | one π pulse | 1 |
| `CPMG` | n equally spaced π pulses | 2 or more |
| `WAHUHA` | four 90° pulses (gCCE) | 4 |
| `Manual` | the sequence written out in [sequence](#tag_sequence) | any |

`WAHUHA` is, as a `Manual` sequence, `[["1/6", "X", 90], ["1/6", "Y", 270], ["2/6", "Y", 90], ["1/6", "X", 270],
["1/6", "I", 0]]`.

```json
"pulse": {"count": 4, "name": "WAHUHA"}
```

<a id="tag_sequence"></a>
### pulse.sequence

| Type | Default |
|---|---|
| list of `[fraction, axis, angle]` | the sequence built from `count` |

The sequence written out, one entry per pulse:

- `fraction`: the free evolution before the pulse, as a fraction of the total time, a string such as `"1/6"` or a
  number. The fractions must add up to 1.
- `axis`: `X`, `Y`, `Z` or `I` (no pulse).
- `angle`: the rotation in degrees.

It is required, with `count` + 1 entries, when [name](#tag_pulse_name) is `Manual`, and refused otherwise. The last
entry is the free evolution after the last pulse.
[Basic 04](../../example/tutorial/basic/04_pulses/README.md) writes a CPMG-2 this way and compares it with the built-in
one.

```json
"pulse": {"count": 2, "name": "Manual", "sequence": [["1/4", "X", 180], ["1/2", "X", 180], ["1/4", "I", 0]]}
```

<a id="tag_pulse_per_qubit"></a>
### pulse.per_qubit

| Type | Default |
|---|---|
| boolean | `false` |

With `false` a pulse rotates the qubits as a whole. With `true` it is applied to each qubit separately, between that
qubit's own two levels, the m<sub>s</sub> of [Qubit.alpha_state and beta_state](InputFile_Qubit.md#tag_alpha_state)
written as `{"ms": [...]}`. `true` is refused with `{"index"}` or `{"vector"}` states, with a finite pulse
([detuning_factor](#tag_detuning_factor) other than 1), with bath pulses, and with meCCE, megCCE, noise and pseudospin.

Known problem: with gCCE, `true` stops with an internal matrix error for Hahn-echo, CPMG and `Manual` sequences
(Ramsey is not affected). Keep the default `false` with gCCE.

```json
"pulse": {"count": 1, "per_qubit": true}
```

<a id="tag_pulse_duration"></a>
### pulse.duration

| Type | Default | Unit |
|---|---|---|
| number | `0` | ns; or `{"value", "unit"}` with `s`, `ms`, `us`, `ns` ([Units](InputFile.md#Units)) |

The length of a pulse, used only by the finite pulses of gCCE ([detuning_factor](#tag_detuning_factor) other than 1).
With gCCE and `detuning_factor` 1 it is reported as not used. cCCE, meCCE and megCCE refuse a value other than 0.

```json
"pulse": {"count": 1, "duration": 10.0, "detuning_factor": 0.98}
```

<a id="tag_detuning_factor"></a>
### pulse.detuning_factor

| Type | Default |
|---|---|
| number | `1.0` |

`1.0` applies every pulse as an ideal instantaneous rotation. Another value, with gCCE only, propagates each pulse as a
finite drive of [duration](#tag_pulse_duration): the drive frequency is `detuning_factor` times the qubit transition
frequency, the detuning is their difference, and the Rabi frequency follows from the angle and the duration. cCCE,
meCCE and megCCE refuse a value other than 1.

```json
"pulse": {"count": 1, "duration": 10.0, "detuning_factor": 0.98}
```

<a id="tag_pulse_bath"></a>
## 2. DEER bath pulses

### pulse.bath

| Type | Default |
|---|---|
| object | no block: no bath pulse |

Pump pulses on the bath spins of one species during the qubit's sequence, to measure the qubit's coupling to one group
of spins (a pumped line). They run with gCCE or megCCE, `cluster.mode` standard or pcce, on the CPU: `execution.backend`
`auto` is set to `cpu`, and `cuda` is refused. [DEER bath pulses](../methods/deer.md) explains what is supported, and
[Advanced 01](../../example/tutorial/advanced/01_P1/README.md) pumps the central P1 line against the plain Hahn echo.

Without the block, or with `count` 0, no bath pulse is applied. The block takes the six keys below; any other key is
refused.

```json
"pulse": {"count": 1, "bath": {"count": 1, "defect": "P1", "energy_shift": [79.975539],
                               "sequence": [["0.50", "X", 180.0], ["0.50", "I", 0.0]]}}
```

<a id="tag_pulse_bath_count"></a>
### pulse.bath.count

| Type | Default |
|---|---|
| integer, 0 or more | `0` |

The number of bath pulses. `0` turns the block off.

```json
"bath": {"count": 1}
```

<a id="tag_pulse_bath_sequence"></a>
### pulse.bath.sequence

| Type | Default |
|---|---|
| list of `[fraction, axis, angle]` | required when `count` > 0 |

The bath pulses written out like [pulse.sequence](#tag_sequence), on the same time axis: `count` + 1 rows, the
fractions adding up to 1 (within 1e-9). The last row is the interval after the last pulse, and its pulse must be the
identity (`I` or angle 0). The run merges it with the qubit sequence and prints the merged sequence; pulses within
1e-9 of the sequence act at the same time.

```json
"sequence": [["0.505", "X", 180.0], ["0.495", "I", 0.0]]
```

This is a π pump 0.005 of the sequence after the Hahn echo's qubit π.

<a id="tag_pulse_bath_defect"></a>
### pulse.bath.defect

| Type | Default |
|---|---|
| string (case is not significant) | required when `count` > 0 |

The name of the bath spins the pump can drive, as in the bath file, or a Defect sub-spin's name. A name that no spin
of the run has is refused.

```json
"defect": "P1"
```

<a id="tag_pulse_bath_energy_shift"></a>
### pulse.bath.energy_shift

| Type | Default | Unit |
|---|---|---|
| number or list of numbers | required when `count` > 0 | MHz; or `{"value", "unit"}` with `Hz`, `kHz`, `MHz`, `GHz` |

The pump tones. A spin of [defect](#tag_pulse_bath_defect) is driven when a tone is within
[bandwidth](#tag_pulse_bath_bandwidth) of its detuning, the `detuning` of its
[Defect species](InputFile_Defect.md#tag_df_detuning). When no spin is within reach of a tone, the run warns once and
the pump does nothing.

```json
"energy_shift": [79.975539, 200.0]
```

<a id="tag_pulse_bath_bandwidth"></a>
### pulse.bath.bandwidth

| Type | Default | Unit |
|---|---|---|
| positive number | 10 rad·kHz (0.0015915494309189536 MHz) | MHz; or `{"value", "unit"}` with `Hz`, `kHz`, `MHz`, `GHz` |

The half width of the frequency window of a tone: a spin is driven when its detuning is closer than this to a tone.
The default is a fixed provisional value, not derived from the pulse length.

```json
"bandwidth": {"value": 5, "unit": "kHz"}
```

<a id="tag_pulse_bath_levels"></a>
### pulse.bath.levels

| Type | Default |
|---|---|
| list of two numbers | `[0.5, -0.5]` |

The two levels `[m_a, m_b]` the bath pulse rotates between (the \|a⟩ and \|b⟩ of [DEER bath pulses](../methods/deer.md)).
They are two different levels of the driven spins' S, and all spins named by [defect](#tag_pulse_bath_defect) must
have the same S.

```json
"levels": [0.5, -0.5]
```
