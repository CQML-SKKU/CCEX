# Pulse tags
The `pulse` block (older inputs had its keys at the top level under other names, see the table of
[InputFile.md](InputFile.md#OldNames)). The block and its `count` are required.

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
    - [bath_drive](#tag_pulse_bath_drive) (continuous drive of bath spins)
        - [tones](#tag_pulse_bath_drive_tones)
        - [shifts](#tag_pulse_bath_drive_shifts)
        - [qubit_transverse](#tag_pulse_bath_drive_qubit_transverse)
        - [line_tolerance](#tag_pulse_bath_drive_line_tolerance)

<a id="tag_pulse"></a>
## pulse
<b>Description:</b><br>
The pulse sequence of the run. Without it the run stops (`pulse is not found in the input file`). `method` meCCE and
megCCE need `per_qubit` false, `detuning_factor` 1 and `duration` 0, and meCCE takes its pulses from `count` only (a
`Manual` sequence is refused). cCCE needs `detuning_factor` 1 and `duration` 0 (ideal instantaneous pulses; the two
keys were ignored there without a word). <br>
<b>Data Type:</b> Object

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [["1/4", "I", 0], ["1/2", "X", 180], ["1/4", "X", 180]]} ```

<a id="tag_pulse_count"></a>
## pulse.count
<b>Description:</b><br>
The number of pulses. Without [sequence](#tag_sequence) the run builds the default sequence
from it: `0` is a Ramsey experiment (free induction), `1` a Hahn echo, `n` > 1 CPMG with n
equally spaced π pulses. With [name](#tag_pulse_name) = `Manual` it is the number of
entries of `sequence` minus one. `method` = `cCCE` applies every pulse as an ideal π flip about
X or Y; other angles need gCCE. Required. <br>

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "pulse" : {"count": 1} ```

<a id="tag_pulse_name"></a>
## pulse.name
<b>Description:</b><br>
The name of the pulse sequence. Without it the run derives the sequence from
[count](#tag_pulse_count) and reports the name it chose (`Ramsey`, `HahnEcho` or `CPMG`); the
run log prints the sequence as a table of intervals, axes and angles. A written name is one of the five below and has to
fit `count` (an unknown name was kept and `count` decided without a word, so `"Hahn"` ran a Hahn echo and
`{"count": 0, "name": "CPMG"}` a Ramsey): otherwise the run stops with `pulse.name X needs pulse.count Y`. <br>

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
Y 270°, 2/6, Y 90°, 1/6, X 270°, 1/6 (as a `Manual` sequence: [["1/6", "X", 90], ["1/6", "Y", 270], ["2/6", "Y", 90],
["1/6", "X", 270], ["1/6", "I", 0]]). `count` must be 4: the run refuses `WAHUHA` with any other count.

<b>Default:</b> - (derived from `count`)

<b>Example:</b> ``` "pulse" : {"count": 2, "name": "Manual", "sequence": [...]} ```

<a id="tag_sequence"></a>
## pulse.sequence
<b>Description:</b><br>
The pulse sequence written out, one entry per pulse: `[fraction, axis, angle]`, where
`fraction` is the free evolution before the pulse as a fraction of the total time (a string
such as `"1/6"` or a number; the fractions must sum to 1), `axis` is `X`, `Y`, `Z` or `I` (no
pulse) and `angle` is the rotation in degrees. Required, with `count` + 1 entries, when
[name](#tag_pulse_name) = `Manual`, and refused with any other name (or none); the last entry is the free evolution after the last pulse.
Without it the run builds the default sequence of `count` equally spaced π pulses (Hahn echo,
CPMG). `example/tutorial/02_pulses/cCCE2_Manual_ensemble/` writes a CPMG-2 this way and compares it with the built-in one. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (the default sequence)

<b>Example:</b> ``` "pulse" : {"count": 4, "name": "Manual", "sequence" : [["1/6", "X", 90.0], ["1/6", "Y", 270.0], ["2/6", "Y", 90.0], ["1/6", "X", 270.0], ["1/6", "I", 0.0]]} ```

<a id="tag_pulse_per_qubit"></a>
## pulse.per_qubit
<b>Description:</b><br>
`false`: a pulse rotates the qubit array as a whole. `true`: the pulse is applied to each qubit
separately and the rotations are combined. With a finite pulse ([detuning_factor](#tag_detuning_factor)
other than 1) `true` is not implemented. Each qubit is turned between its own two levels, the m<sub>s</sub> of
[Qubit.alpha_state and beta_state](InputFile_Qubit.md#tag_alpha_state) `{"ms": [...]}`: `true` with `{"index"}` or
`{"vector"}` is refused. <br>
**Known defect (not fixed):** with gCCE, `true` also stops with an internal
matrix error for ideal Hahn-echo and CPMG pulses (the per-qubit X operator is left empty) and for a `Manual` sequence (a
matrix of the wrong size); Ramsey (no pulse) is not affected. Use the default `false` until it is fixed. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "pulse" : {"count": 1, "per_qubit": true} ```

<a id="tag_pulse_duration"></a>
## pulse.duration
<b>Description:</b><br>
The duration of a pulse. Used only by the finite-pulse propagation of gCCE, i.e. with
[detuning_factor](#tag_detuning_factor) other than 1: with gCCE and `detuning_factor` 1 it is reported as not used,
and cCCE, meCCE and megCCE refuse a value other than 0. <br>

<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the unit `s`, `ms`, `us` or `ns` ([Units](InputFile.md#Units))

<b>Unit:</b> ns

<b>Default:</b> `0.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<a id="tag_detuning_factor"></a>
## pulse.detuning_factor
<b>Description:</b><br>
gCCE (cCCE, meCCE and megCCE refuse a value other than 1). `1.0` applies every pulse as an ideal instantaneous
rotation. Any other value propagates
each pulse as a finite drive of [duration](#tag_pulse_duration): the drive frequency is
`detuning_factor` times the qubit transition frequency, the detuning is their difference, and the
Rabi frequency follows from the angle and the duration. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "pulse" : {"count": 1, "duration": 10.0, "detuning_factor": 0.98} ```

<a id="tag_pulse_bath"></a>
## pulse.bath
<b>Description:</b><br>
DEER bath pulses: pump pulses on the bath spins of one species during the qubit's sequence, for `method` = `gCCE` or
`megCCE` (`cluster.mode` standard or pcce) on the CPU (`execution.backend` `auto` is set to `cpu`, `cuda` is refused). What it does, what it supports
and what it refuses is in [deer.md](../methods/deer.md). Without the block, or with `count` 0, no bath pulse is applied. Its
keys are the six below; any other key is refused (`pulse.bath has no key "..."`). The earlier key
names are refused with their new names ([the old names](InputFile.md#OldNames)). Use it to measure the coupling of the qubit
to one group of bath spins (a pumped line); `example/tutorial/03_P1/gCCE2_Hahn_DEER/` runs a pump on the central P1 line
against the plain Hahn echo. <br>

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
The bath pulses written out like [pulse.sequence](#tag_sequence), on the same time axis: `count` + 1 rows
`[fraction, axis, angle]`, each the interval before a pulse as a fraction of the total time (a number or a string such
as `"1/4"`; the fractions sum to 1 within 1e-9), the axis `X`, `Y`, `Z` or `I`, the angle in degrees. The last row is
the interval after the last pulse; its pulse must be identity (`I` or angle 0). The run merges it with the qubit
sequence: pulses closer than 1e-9 of the sequence act at the same time, and two pulses of one sequence that would land
on one edge are refused (each pulse is applied exactly once). The run log prints the merged sequence. Required when
`count` > 0. <br>

<b>Data Type:</b> List of [String or Double, String, Double]

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "sequence" : [["0.505", "X", 180.0], ["0.495", "I", 0.0]] ``` (a π pump 0.005 of the sequence after the Hahn echo's qubit π)

<a id="tag_pulse_bath_defect"></a>
## pulse.bath.defect
<b>Description:</b><br>
The name of the bath spins the pump can drive, as in the bath file (or a Defect sub-spin's name); the comparison
ignores case. A name no spin of the run has is refused. Required when `count` > 0. <br>

<b>Data Type:</b> String

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "defect" : "P1" ```

<a id="tag_pulse_bath_energy_shift"></a>
## pulse.bath.energy_shift
<b>Description:</b><br>
The pump tone(s). A spin of [defect](#tag_pulse_bath_defect) is driven when one tone is within
[bandwidth](#tag_pulse_bath_bandwidth) (default 10 rad·kHz, 1.59 kHz) of its detuning, the `detuning` of its Defect
configuration ([Defect](InputFile_Defect.md#tag_df_detuning)). When no spin of that name is within reach of a
tone, the run warns once and goes on (the pump is then the identity). Required when `count` > 0, at least one number. <br>

<b>Data Type:</b> Double or list of Double, or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> MHz

<b>Default:</b> - (required when `count` > 0)

<b>Example:</b> ``` "energy_shift" : [79.975539, 200.0] ```

<a id="tag_pulse_bath_bandwidth"></a>
## pulse.bath.bandwidth
<b>Description:</b><br>
The half width of the frequency window of the bath pulse: a spin of [defect](#tag_pulse_bath_defect)
is driven when its detuning is closer than this to one of the [tones](#tag_pulse_bath_energy_shift). The default is a
fixed window of 10 rad·kHz (a provisional value, not derived from the pulse length); written out it is
`0.0015915494309189536` MHz (10 / (2π·1000), which gives back 10 rad·kHz exactly). A
number that is not positive and finite is refused. <br>

<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> MHz

<b>Default:</b> 10 rad·kHz = 0.0015915494309189536 MHz

<b>Example:</b> ``` "bandwidth" : {"value": 5, "unit": "kHz"} ```

<a id="tag_pulse_bath_levels"></a>
## pulse.bath.levels
<b>Description:</b><br>
The two levels m the bath pulse rotates between, `[m_a, m_b]` (the |a> and |b> of the generalized Pauli operators in
[deer.md](../methods/deer.md)): two different levels of the driven spins' S (m = -S, -S+1, ..., S). S is the driven
spins' own; the spins named by [defect](#tag_pulse_bath_defect) all have to have one S. Before 2.0.0 `alphams` and `betams`
were the two levels and `spin` repeated S (refused when it differed); the three are refused with this key. <br>

<b>Data Type:</b> List of two Double

<b>Default:</b> `[0.5, -0.5]`

<b>Example:</b> ``` "levels" : [0.5, -0.5] ```

<a id="tag_pulse_bath_drive"></a>
## pulse.bath_drive
<b>Description:</b><br>
Continuous driving of bath spins during the qubit's sequence, in the rotating frame of each drive tone with the
rotating-wave approximation ([bath-drive.md](../methods/bath-drive.md)): the setting of García Oliván et al.,
arXiv 2512.06948. The qubit gets the instantaneous π pulses of [count](#tag_pulse_count) (Ramsey, Hahn echo, CPMG);
electron-like spin-1/2 bath spins (γ < 0) on the lines of a [Defect](InputFile_Defect.md) table are driven by
[tones](#tag_pulse_bath_drive_tones). Every bath spin is put in a frame: the tone's when it drives the spin's line,
the spin's own line otherwise; only what is static there is kept (pairs in different frames keep T<sub>zz</sub> J<sub>z</sub>J<sub>z</sub>
only). Runs with `method` cCCE, `cluster.mode` standard or pcce, one qubit, `quantity` coherence, `bfield` along +z and
the ensemble (no `configuration.s_range`; sampled states with a mean field are not implemented yet). Not with
[pulse.bath](#tag_pulse_bath), a `Manual` sequence, `duration` or `detuning_factor` other than 1. Every cluster goes
through the drive kernel (no CPU batch, no GPU path). A carrier within 3 W of a line it does not target, or two
frames within 3 W of each other, is reported once (the per-line rotating-wave approximation drops a near-resonant
process there). `tones` `[]` puts every spin in the frame of its own line without a field: the undriven reference of
the same model. <br>
<b>Data Type:</b> Object

<b>Default:</b> no block: no drive.

<b>Example:</b> ``` "pulse" : {"count": 1, "bath_drive": {"tones": [{"lines": [59.91], "rabi": 7.0}, {"lines": [-59.91], "rabi": 7.0}], "qubit_transverse": true}} ```

<a id="tag_pulse_bath_drive_tones"></a>
## pulse.bath_drive.tones
<b>Description:</b><br>
The drive tones, each `{"lines": [...], "rabi": W, "offset": o, "phase": p}`: a field of Rabi frequency W (H = W J<sub>x</sub>
gives Rabi oscillation at W) on the bath spins whose detuning is within [line_tolerance](#tag_pulse_bath_drive_line_tolerance)
of one of `lines`, at the frequency line + `offset`, with the phase `p` (degrees) in the frame. `lines`, `rabi` and
`offset` are MHz, as [Defect.species[].detuning](InputFile_Defect.md#tag_df_detuning); every key has the default 0
(a tone with `rabi` 0 drives nothing). `rabi` < 0 is refused (a sign goes into `phase`), a spin on the lines of two tones
is refused, and two tones may not share a carrier (combine them into one tone). <br>
<b>Data Type:</b> List of objects

<b>Default:</b> - (required in the block; `[]` allowed)

<b>Example:</b> ``` "tones" : [{"lines": [-59.91], "rabi": 7.0, "offset": 4.95}, {"lines": [59.91], "rabi": 7.0}] ```

<a id="tag_pulse_bath_drive_shifts"></a>
## pulse.bath_drive.shifts
<b>Description:</b><br>
Add the leading lab-frame corrections: the Bloch-Siegert and AC Stark shifts of the driven spins from every tone
(an AC Stark term only outside 3 W of the tone), the Van Vleck shift m<sub>s</sub><sup>2</sup>|A<sub>⊥</sub>|<sup>2</sup>/2ω of the
qubit's transverse coupling (with `qubit_transverse` false; with true it is exact already) and the drive-assisted
cross-line flip-flop −(T<sub>zz</sub>/2)(W/D) J<sub>⊥</sub>J<sub>z</sub>. <br>
<b>Data Type:</b> Boolean

<b>Default:</b> `false` (the paper's model)

<a id="tag_pulse_bath_drive_qubit_transverse"></a>
## pulse.bath_drive.qubit_transverse
<b>Description:</b><br>
Keep the qubit coupling A<sub>zx</sub>, A<sub>zy</sub> of every bath spin exactly: in branch m<sub>s</sub> spin k precesses about
ω<sub>k</sub>ẑ + m<sub>s</sub>A<sub>k</sub>, so each branch has its own tilted axes, the rotating-wave approximation is taken about them and the
frame changes between the branches are applied at every π pulse and at the end. `false` keeps A<sub>zz</sub> only (García Oliván
Eq. S.17). <br>
<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<a id="tag_pulse_bath_drive_line_tolerance"></a>
## pulse.bath_drive.line_tolerance
<b>Description:</b><br>
A bath spin is on a line of a tone when its detuning is closer than this. Two frames are the same only when they come
from the same line and offset. <br>
<b>Data Type:</b> Double

<b>Unit:</b> MHz

<b>Default:</b> 0.001

