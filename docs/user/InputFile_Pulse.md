# Pulse tags
The `pulse` block. The block and its `count` are required. Older inputs had its keys at the top level under other
names; see the table of [InputFile.md](InputFile.md#OldNames).

``` "pulse" : {"count": 1} ```

- [pulse](#tag_pulse)
    - [count](#tag_pulse_count)
    - [name](#tag_pulse_name)
    - [sequence](#tag_sequence)
    - [per_qubit](#tag_pulse_per_qubit)
    - [duration](#tag_pulse_duration)
    - [detuning_factor](#tag_detuning_factor)
    - [bath](#tag_pulse_bath) (DEER bath pulses)
        - [count](#tag_pulse_bath_count)
        - [sequence](#tag_pulse_bath_sequence)
        - [defect](#tag_pulse_bath_defect)
        - [energy_shift](#tag_pulse_bath_energy_shift)
        - [bandwidth](#tag_pulse_bath_bandwidth)
        - [levels](#tag_pulse_bath_levels)

<a id="tag_pulse"></a>
## pulse
<b>Description:</b><br>
The pulse sequence of the run. Without it the run stops (`pulse is not found in the input file`). <br>
Limits by method:

- cCCE: `detuning_factor` 1 and `duration` 0 (ideal instantaneous pulses).
- meCCE and megCCE: `per_qubit` false, `detuning_factor` 1 and `duration` 0.
- meCCE takes its pulses from `count` only; a `Manual` sequence is refused.

<b>Data Type:</b> Object

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [["1/4", "X", 180], ["1/2", "X", 180], ["1/4", "I", 0]]} ``` (CPMG-2: wait t/4, π about X, wait t/2, π about X, wait t/4)

<a id="tag_pulse_count"></a>
## pulse.count
<b>Description:</b><br>
The number of pulses. Required. <br>
Without [sequence](#tag_sequence) the run builds the default sequence from it:

- `0`: a Ramsey experiment (free induction);
- `1`: a Hahn echo;
- `n` > 1: CPMG with n equally spaced π pulses.

With [name](#tag_pulse_name) = `Manual` it is the number of entries of `sequence` minus one. <br>
`method` = `cCCE` applies every pulse as an ideal π flip about X or Y. Other angles need gCCE. <br>

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "pulse" : {"count": 1} ```

<a id="tag_pulse_name"></a>
## pulse.name
<b>Description:</b><br>
The name of the pulse sequence. A written name is one of the five below and has to fit `count`; otherwise the run
stops with `pulse.name X needs pulse.count Y`. <br>
Without it the run derives the sequence from [count](#tag_pulse_count) and reports the name it chose (`Ramsey`,
`HahnEcho` or `CPMG`). The run log prints the sequence as a table of intervals, axes and angles. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `Ramsey` | free induction: `count` 0 |
| `HahnEcho` | one π pulse: `count` 1 |
| `CPMG` | n equally spaced π pulses: `count` 2 or more |
| `WAHUHA` | the four-pulse WAHUHA sequence (90° pulses; gCCE): `count` must be 4 |
| `Manual` | the sequence written out in [sequence](#tag_sequence), any `count` |

`WAHUHA` is a fixed sequence of four pulses and five free intervals, in fractions of the total time: 1/6, X 90°, 1/6,
Y 270°, 2/6, Y 90°, 1/6, X 270°, 1/6. As a `Manual` sequence it is [["1/6", "X", 90], ["1/6", "Y", 270], ["2/6", "Y", 90],
["1/6", "X", 270], ["1/6", "I", 0]]. The run refuses `WAHUHA` with a `count` other than 4.

<b>Default:</b> - (derived from `count`)

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [...]} ```

<a id="tag_sequence"></a>
## pulse.sequence
<b>Description:</b><br>
The pulse sequence written out, one entry `[fraction, axis, angle]` per pulse:

- `fraction`: the free evolution before the pulse, as a fraction of the total time (a string such as `"1/6"` or a
  number). The fractions must sum to 1.
- `axis`: `X`, `Y`, `Z` or `I` (no pulse).
- `angle`: the rotation in degrees.

Required, with `count` + 1 entries, when [name](#tag_pulse_name) = `Manual`; refused with any other name or none. The
last entry is the free evolution after the last pulse. <br>
Without it the run builds the default sequence of `count` equally spaced π pulses (Hahn echo, CPMG).
`example/tutorial/basic/04_pulses/cCCE2_Manual_ensemble/` writes a CPMG-2 this way and compares it with the built-in one. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (the default sequence)

<b>Example:</b> ``` "pulse" : {"count": 4, "name": "Manual", "sequence" : [["1/6", "X", 90.0], ["1/6", "Y", 270.0], ["2/6", "Y", 90.0], ["1/6", "X", 270.0], ["1/6", "I", 0.0]]} ```

<a id="tag_pulse_per_qubit"></a>
## pulse.per_qubit
<b>Description:</b><br>
`false`: a pulse rotates the qubit array as a whole. `true`: the pulse is applied to each qubit separately and the
rotations are combined. <br>
With `true` each qubit is turned between its own two levels, the m<sub>s</sub> of
[Qubit.alpha_state and beta_state](InputFile_Qubit.md#tag_alpha_state) `{"ms": [...]}`. Not available with `true`:

- a finite pulse ([detuning_factor](#tag_detuning_factor) other than 1): not implemented;
- `{"index"}` or `{"vector"}` states: refused.

**Known defect (not fixed):** with gCCE, `true` stops with an internal matrix error for ideal Hahn-echo and CPMG pulses
(the per-qubit X operator is left empty) and for a `Manual` sequence (a matrix of the wrong size). Ramsey (no pulse) is
not affected. Use the default `false` until it is fixed. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "pulse" : {"count": 1, "per_qubit": true} ```

<a id="tag_pulse_duration"></a>
## pulse.duration
<b>Description:</b><br>
The duration of a pulse. Used only by the finite-pulse propagation of gCCE, that is with
[detuning_factor](#tag_detuning_factor) other than 1. With gCCE and `detuning_factor` 1 it is reported as not used.
cCCE, meCCE and megCCE refuse a value other than 0. <br>

<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the unit `s`, `ms`, `us` or `ns` ([Units](InputFile.md#Units))

<b>Unit:</b> ns

<b>Default:</b> `0.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<a id="tag_detuning_factor"></a>
## pulse.detuning_factor
<b>Description:</b><br>
gCCE only; cCCE, meCCE and megCCE refuse a value other than 1. <br>
`1.0` applies every pulse as an ideal instantaneous rotation. Any other value propagates each pulse as a finite drive of
[duration](#tag_pulse_duration):

- the drive frequency is `detuning_factor` times the qubit transition frequency;
- the detuning is their difference;
- the Rabi frequency follows from the angle and the duration.

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<a id="tag_pulse_bath"></a>
## pulse.bath
<b>Description:</b><br>
DEER bath pulses: pump pulses on the bath spins of one species during the qubit's sequence. Use it to measure the
coupling of the qubit to one group of bath spins (a pumped line). `example/tutorial/advanced/01_P1/gCCE2_Hahn_DEER/` runs a pump
on the central P1 line against the plain Hahn echo. <br>
For `method` = `gCCE` or `megCCE` (`cluster.mode` standard or pcce) on the CPU: `execution.backend` `auto` is set to
`cpu`, and `cuda` is refused. What it does, supports and refuses is in [deer.md](../methods/deer.md). <br>
Without the block, or with `count` 0, no bath pulse is applied. Its keys are the six below; any other key is refused
(`pulse.bath has no key "..."`). Old key names are refused with their new names ([the old names](InputFile.md#OldNames)). <br>

<b>Data Type:</b> Object

<b>Default:</b> no block: no bath pulse.

<b>Example:</b> ``` "pulse" : {"count": 1, "bath": {"count": 1, "defect": "P1", "energy_shift": [79.975539], "sequence": [["0.50", "X", 180.0], ["0.50", "I", 0.0]]}} ```

<a id="tag_pulse_bath_count"></a>
## pulse.bath.count
<b>Description:</b><br>
The number of bath pulses. `0` turns the block off; a negative number is refused. <br>

<b>Data Type:</b> Integer

<b>Default:</b> `0`

<b>Example:</b> ``` "bath" : {"count": 1, ...} ```

<a id="tag_pulse_bath_sequence"></a>
## pulse.bath.sequence
<b>Description:</b><br>
The bath pulses written out like [pulse.sequence](#tag_sequence), on the same time axis. Required when `count` > 0. <br>
`count` + 1 rows `[fraction, axis, angle]`:

- `fraction`: the interval before a pulse, as a fraction of the total time (a number or a string such as `"1/4"`). The
  fractions sum to 1 within 1e-9.
- `axis`: `X`, `Y`, `Z` or `I`.
- `angle`: in degrees.

The last row is the interval after the last pulse; its pulse must be identity (`I` or angle 0). <br>
The run merges it with the qubit sequence. Pulses closer than 1e-9 of the sequence act at the same time. Two pulses of
one sequence that would land on one edge are refused, so each pulse is applied exactly once. The run log prints the
merged sequence. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "sequence" : [["0.505", "X", 180.0], ["0.495", "I", 0.0]] ``` (a π pump 0.005 of the sequence after the Hahn echo's qubit π)

<a id="tag_pulse_bath_defect"></a>
## pulse.bath.defect
<b>Description:</b><br>
The name of the bath spins the pump can drive, as in the bath file, or a Defect sub-spin's name. The comparison
ignores case. A name no spin of the run has is refused. Required when `count` > 0. <br>

<b>Data Type:</b> String

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "defect" : "P1" ```

<a id="tag_pulse_bath_energy_shift"></a>
## pulse.bath.energy_shift
<b>Description:</b><br>
The pump tone(s). Required when `count` > 0, at least one number. <br>
A spin of [defect](#tag_pulse_bath_defect) is driven when one tone is within [bandwidth](#tag_pulse_bath_bandwidth)
(default 10 rad·kHz, 1.59 kHz) of its detuning, the `detuning` of its Defect configuration
([Defect](InputFile_Defect.md#tag_df_detuning)). When no spin of that name is within reach of a tone, the run warns once
and goes on; the pump is then the identity. <br>

<b>Data Type:</b> Double or list of Double, or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> MHz

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "energy_shift" : [79.975539, 200.0] ```

<a id="tag_pulse_bath_bandwidth"></a>
## pulse.bath.bandwidth
<b>Description:</b><br>
The half width of the frequency window of the bath pulse. A spin of [defect](#tag_pulse_bath_defect) is driven when
its detuning is closer than this to one of the [tones](#tag_pulse_bath_energy_shift). A number that is not positive and
finite is refused. <br>
The default is a fixed window of 10 rad·kHz, a provisional value not derived from the pulse length. Written out it is
`0.0015915494309189536` MHz, $10/(2\pi \cdot 1000)$, which gives back 10 rad·kHz exactly. <br>

<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> MHz

<b>Default:</b> 10 rad·kHz = 0.0015915494309189536 MHz

<b>Example:</b> ``` "bandwidth" : {"value": 5, "unit": "kHz"} ```

<a id="tag_pulse_bath_levels"></a>
## pulse.bath.levels
<b>Description:</b><br>
The two levels m the bath pulse rotates between, `[m_a, m_b]`: the |a> and |b> of the generalized Pauli operators in
[deer.md](../methods/deer.md). They are two different levels of the driven spins' own S (m = -S, -S+1, ..., S). The
spins named by [defect](#tag_pulse_bath_defect) must all have the same S. <br>
The old keys `alphams`, `betams` and `spin` of this block are refused; `example/tools/InputKeys/convert_input_keys.py`
rewrites them. <br>

<b>Data Type:</b> List of two Double

<b>Default:</b> `[0.5, -0.5]`

<b>Example:</b> ``` "levels" : [0.5, -0.5] ```
