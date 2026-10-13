# DEER bath pulses

What the `pulse.bath` block of the input file does: pump pulses on bath spins (double electron-electron resonance,
DEER) during the qubit's own sequence, for gCCE and megCCE (ME-gCCE) on the CPU, with `cluster_mode` standard or pcce.
Ported from the adding_DEER / adding_MECCE branches (`origin/adding_MECCE` 2c63eac52), whose flat keys are refused with
the new names. The keys are listed in [InputFile_Pulse.md](../user/InputFile_Pulse.md#tag_pulse_bath); the test is
`test/code_verification/check_deer/`.

## The physics

**Which spins a pump drives.** Bath spin $j$ is driven when its name is `pulse.bath.defect` (case ignored) and one of the
pump tones $\omega_k$ (`pulse.bath.energy_shift`) lies within a fixed window of its detuning $\Delta_j$:

$$\text{name}_j = \texttt{defect} \quad\text{and}\quad \exists k:\ \vert\omega_k-\Delta_j\vert < 10\ \text{rad·kHz}.$$

$\Delta_j$ is the `detuning` of the spin's Defect configuration (`Defect.species[].detuning`, MHz, the line of its
principal axis), stored in rad·kHz like the tones ($2\pi\cdot10^3$ per MHz). The window, 10 rad·kHz (1.59 kHz), is an
**arbitrary constant** taken over from the branch (`adding_MECCE:src/simulator_gcce.cpp:546-553`), not derived from a
pulse bandwidth; it is not an input key. A tone either hits a line or misses it: there is no partial excitation and no
pulse shape.

**The pulse on a driven spin.** In the two levels $\vert a\rangle,\vert b\rangle$ given by `pulse.bath.alphams` and
`pulse.bath.betams` of a spin of size `pulse.bath.spin`, with the generalized Pauli operators

$$\sigma_x = \vert a\rangle\langle b\vert + \vert b\rangle\langle a\vert,\quad
\sigma_y = -i\vert a\rangle\langle b\vert + i\vert b\rangle\langle a\vert,\quad
\sigma_z = \vert a\rangle\langle a\vert - \vert b\rangle\langle b\vert,$$

a row `[fraction, axis, angle]` of `pulse.bath.sequence` applies $U_j = \exp(-i\theta\,\sigma_{axis}/2)$ ($\theta$ the angle
in radians), the same form as an ideal qubit pulse. A spin that is not driven gets the identity, and the pulse on the
cluster is $U_b = \bigotimes_j U_j$ over the cluster's spins.

**When.** `pulse.bath.sequence` is written like the qubit's `pulse.sequence`: each row is an interval (a fraction of the
total time) with its pulse at the end, and the last row's pulse is never applied (it must be identity). The two
sequences are merged into one: every pulse edge of either becomes an edge, and edges closer than $10^{-9}$ of the
sequence are one edge, where both pulses act at once. With the merged free evolutions $U_{\text{free},i}$,

$$U = U_{\text{free},N}\prod_{i=N-1}^{0}\big(U_{q,i}\otimes U_{b,i}\big)\,U_{\text{free},i},$$

with $U_{q,i}$ or $U_{b,i}$ the identity at an edge that has only the other pulse. Every pulse is applied exactly once:
two pulses of one channel that the merge would put on one edge are refused.

**megCCE.** The density matrix of the cluster (qubit and bath, Liouville space, row-major) is propagated by
$e^{\mathcal L\tau_i}$ over the merged intervals, and at an edge by the pulse superoperator
$\mathcal U_i=U_i\otimes U_i^{*}$ with the same $U_i=U_{q,i}\otimes U_{b,i}$ as above ($\mathcal U_i[\rho]=U_i\rho U_i^\dagger$).
The 0th cluster has no bath spin: $U_b=1$ there. At rate 0 (no jump operators) this is the gCCE result (`check_deer`).

**What a pump does at order 1.** The bath spins do not interact within a cluster of one spin, and a bath $\pi$ at the
fraction $f$ of a Hahn echo leaves the qubit with the coupling $a_j$ of the driven spins for a net time
$2\min(f,1-f)\,t$ instead of refocusing it. Averaged over the bath state,

$$V(t) = \prod_{j\ \text{driven}} \cos\!\big(a_j\, t\,\min(f,1-f)\big),\qquad
a_j = 2\pi\cdot 52.04\ \text{MHz nm}^3\,\frac{1-3\cos^2\theta_j}{r_j^3},$$

which `check_deer` compares with the run (within $2\times10^{-4}$; at $f=0.5$ the argument is $a_j t/2$).

## Supported and refused

| | supported | refused (the message names the reason) |
|---|---|---|
| method | `gCCE`, `megCCE` | cCCE and meCCE (no bath spin operator in their matrices: the bath is propagated under the qubit-projected Hamiltonians) |
| `cluster_mode` | `standard`, `pcce` (inside pgCCE, pmegCCE) | `bundle` (gCCE); `dsj`, `itb`, `dsjitb` run with cCCE only |
| backend | `cpu`; `auto` is set to `cpu` (logged) | `cuda` |
| `propagator` | gCCE: `eigen`, `expm`; megCCE: `eigen` (the name; its step is a matrix exponential) | `sparse` (it applies pulses as $U_q\otimes\mathbb 1$); megCCE refuses `expm` as without bath pulses |
| `evolution` × `evolution_update` | gCCE: `matrix`, `vector` × `exact`, `recursive`; megCCE: its own Liouville-space step (`exact`) | megCCE: `vector`, `recursive` (as without bath pulses) |
| GPU warp kernel | — | `gpu_options.warp_max_dim` > 0 |
| pulses | ideal (`pulse.detuning_factor` 1), `pulse.per_qubit` false, any qubit sequence | `detuning_factor` ≠ 1, `per_qubit` true |
| `pulse.bath` | its seven keys | an unknown key; `count` < 0; a sequence with other than `count`+1 rows, fractions not summing to 1, a last row that is not identity, two pulses of one channel on one edge; `count` > 0 without `defect` or `energy_shift` |
| driven spins | | a `defect` no spin has (bath file and Defect sub-spins); a `spin` other than that spin's S; `alphams`/`betams` equal or not levels of S |
| old keys | | the branch's `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift`, `bspin`, `balphams`, `bbetams` (each named with its new key) |

When spins of the `defect` exist but no tone reaches any of them, the run goes on with one warning
(`Warning: pulse.bath drives no spin ...`): the pump is then the identity, which is a legitimate control run
(`offres` in the test). `pulse.bath` absent or `count` 0 gives the output bytes of a run without the block.

## Example

The test's `onres` case: one NV (S = 1, levels −1 and 0), twelve P1 at 6.5–12 nm, a Hahn echo with a P1 $\pi$ pump at
the echo centre on the line of configuration 1 (79.975539 MHz):

```
"pulse": {
    "count": 1,
    "bath": {
        "count": 1,
        "defect": "P1",
        "energy_shift": [79.975539],
        "spin": 0.5, "alphams": 0.5, "betams": -0.5,
        "sequence": [["0.50", "X", 180.0], ["0.50", "I", 0.0]]
    }
},
"Defect": {"axis_file": "avaax_DEERtest",
           "species": [{"name": "P1", "n_nuclear_spins": 0, "n_configurations": 8,
                        "detuning": [[1, "e", 79.975539], [2, "e", 59.579552], ...]}]}
```

The run log prints the merged sequence (qubit and bath axis and angle per interval) under the pulse report.

## Where it is in the code

| piece | file |
|---|---|
| reading `pulse.bath`, refusing the branch keys | `src/input/option.cpp` (`cJSON_readOptionPulse`, `refuseOldInputKeys`) |
| merging the sequences, installing the merged one | `src/physical/pulse_bath.cpp` (`mergePulseSequences`, `Pulse_resolveBathSequence`) |
| method/backend refusals | `src/main.cpp` (`validateBathPulses`) |
| driven-spin checks, $U_b$, $U_q\otimes U_b$ in the expanded pulses | `src/calculation/simulator_gcce_bath.cpp` (`validateBathPulseTargets`, `bathPulseOperator`, `applyBathPulsesGcce`) |
| $U_q\otimes U_b$ in the megCCE pulse superoperator | `src/calculation/simulator_megcce.cpp` (`calSuperPropagatorMegcce`) |

`applyBathPulsesGcce` replaces the cached $U_q\otimes\mathbb 1$ of a row by $U_q\otimes U_b$ where a spin of the cluster is
driven, and then makes the vector path take the dense product (its block product assumes $U_q\otimes\mathbb 1$), so
matrix and vector, exact and recursive read the same operators.
