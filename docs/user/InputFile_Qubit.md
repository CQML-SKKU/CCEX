# Qubit tags
<a id="QubitTags"></a>

The `Qubit` block describes the qubits: a list with one entry per qubit, the interactions between them, and the two
states whose coherence is computed. Each key has a small table (its type, default and unit), what it does, and an
example. Old key names are in [Old names](InputFile_OldNames.md).

1. [The block](#qb_block): `Qubit`, and a whole example.
2. [The qubits](#qb_qubits): `qubit` and its keys `name`, `spin`, `gyro`, `xyz` or `file`, `detuning`.
3. [The interactions](#qb_intmap): `intmap` and its keys `between`, `tensor`, `unit`, `frame`.
4. [The states](#qb_states): `alpha_state`, `beta_state` and `psi0`.

The Overhauser field of the bath on the qubits is a key of the `hamiltonian` block,
[qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field).

<a id="qb_block"></a>
## 1. The block

<a id="tag_Qubit"></a>
### Qubit

| Type | Default |
|---|---|
| object | required |

The qubits. A run without a `Qubit` block, or with an empty `qubit` list, stops. One NV read from a defect file:

```json
"Qubit": {"qubit": [{"file": "bath/bath_Diamond_NV_natab_defect"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}}
```

Two qubits with their interactions written out:

```json
"Qubit": {
    "qubit": [
        {"name": "Q1", "spin": 0.5, "gyro": -17608.597050, "xyz": [1.0, 2.0, 3.0]},
        {"name": "Q2", "spin": 1.0, "gyro": -17608.597050, "xyz": [4.0, 5.0, 6.0]}
    ],
    "intmap": [
        {"between": ["Q2", "Q2"], "tensor": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]},
        {"between": ["Q1", "Q2"], "tensor": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]}
    ],
    "alpha_state": {"ms": [0.5, 0]},
    "beta_state":  {"vector": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]},
    "psi0": [0.0, 0.70710678, 0.0, 0.70710678, 0.0, 0.0]
}
```

- The first `intmap` entry is the zero-field splitting of Q2 (`between` names the same qubit twice); the second is the
  Q1–Q2 coupling.
- `alpha_state` gives one m<sub>s</sub> per qubit, in the order of `qubit`: \|+1/2, 0⟩.
- The `vector` of `beta_state` and `psi0` span the whole qubit space, here 2 × 3 = 6, ordered as kron(Q1, Q2) with the
  last qubit varying fastest. This `beta_state` is \|−1/2, +1⟩, the same as `{"ms": [-0.5, 1]}`.

<a id="qb_qubits"></a>
## 2. The qubits

<a id="tag_qubit"></a>
### qubit

| Type | Default |
|---|---|
| list of objects | required (at least one entry) |

One entry per qubit: [name](#tag_name), [spin](#tag_spin), [gyro](#tag_gyro), its position ([xyz](#tag_xyz) or
[file](#tag_Qubit_file)) and [detuning](#tag_detuning). The number of qubits is the length of the list, and
[alpha_state](#tag_alpha_state) and [beta_state](#tag_beta_state) give their levels in this order.

```json
"qubit": [{"name": "Q1", "spin": 0.5, "xyz": [1.0, 2.0, 3.0]}, {"name": "Q2", "spin": 1.0, "xyz": [4.0, 5.0, 6.0]}]
```

<a id="tag_name"></a>
### name

| Type | Default |
|---|---|
| string | `"q0"`, `"q1"`, ... (the position in the list, from 0) |

The qubit's name, used by [between](#tag_between) and by a qubit jump operator.

```json
"name": "NV"
```

<a id="tag_spin"></a>
### spin

| Type | Default |
|---|---|
| number | `1.0` |

The spin of the qubit.

```json
"spin": 0.5
```

<a id="tag_gyro"></a>
### gyro

| Type | Default | Unit |
|---|---|---|
| number | `-17608.597050` (the electron) | rad kHz/G (= rad/ms/G); or `{"value", "unit"}` with the units of [Defect gyros](InputFile_Defect.md#tag_df_gyros) |

The gyromagnetic ratio of the qubit.

```json
"gyro": -17608.597050
```

<a id="tag_xyz"></a>
### xyz

| Type | Default | Unit |
|---|---|---|
| `[x, y, z]` | required, unless `file` is given | Å; or `{"value", "unit"}` with `A`, `nm` |

The position of the qubit, in the frame of the bath file.

```json
"xyz": [1.0, 2.0, 3.0]
```

<a id="tag_Qubit_file"></a>
### file

| Type | Default |
|---|---|
| string (file path) | none: the entry gives `xyz` or `file` |

In place of [xyz](#tag_xyz), the file the position is read from: its first line holds x y z in Å, in the frame of the
bath file (the defect file that `example/tools/BathGenerator` writes, for example).

- The path is relative to the working directory; [configuration.path](InputFile_General.md#tag_configuration_path)
  does not apply to it.
- A file that does not open stops the run, and an entry with both `file` and `xyz` is refused.

```
10.0 12.0 13.0
```

```json
"qubit": [{"file": "defect", "name": "NV", "spin": 1.0}]
```

<a id="tag_detuning"></a>
### detuning

| Type | Default | Unit |
|---|---|---|
| number | `0` | **kHz**; or `{"value", "unit"}` with `Hz`, `kHz`, `MHz`, `GHz` |

Adds detuning × S<sub>z</sub> to the qubit Hamiltonian, a static shift of the qubit levels. The phase of the coherence
turns faster while \|L\| stays the same; in gCCE the resonances of the qubit with the bath spins move.

A plain number is kHz here but MHz in [Defect.species[].detuning](InputFile_Defect.md#tag_df_detuning), so the log
notes a plain number other than 0. A list `[x]` is refused: write `x`.

```json
"detuning": {"value": 2.14, "unit": "MHz"}
```

<a id="qb_intmap"></a>
## 3. The interactions

<a id="tag_intmap"></a>
### intmap

| Type | Default |
|---|---|
| list of objects | zero self terms; point-dipole tensors between different qubits |

The interactions between the qubits, and of a qubit with itself (its zero-field splitting). Each entry names a pair
([between](#tag_between)) and gives its [tensor](#tag_tensor), with an optional [unit](#tag_intmap_unit) and, with
the rotation on, its [frame](#tag_intmap_frame).

- A self entry names the qubit twice (`["q0", "q0"]` for an unnamed first qubit). Not given, it is zero.
- A pair of different qubits that is not given gets the point-dipole tensor from their `xyz` and `gyro`.

```json
"intmap": [{"between": ["q0", "q0"], "unit": "MHz", "tensor": [[-960, 0, 0], [0, -960, 0], [0, 0, 1920]]}]
```

<a id="tag_between"></a>
### between

| Type | Default |
|---|---|
| two names | required in every entry |

The two qubits of the entry; the same name twice for a zero-field splitting.

```json
"between": ["Q1", "Q2"]
```

<a id="tag_tensor"></a>
### tensor

| Type | Default | Unit |
|---|---|---|
| 3 × 3 list | required in every entry | the entry's [unit](#tag_intmap_unit), kHz without it |

The interaction tensor of the pair. CCEX supplies the factor 2π when it converts it to angular frequency.

```json
"tensor": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
```

<a id="tag_intmap_unit"></a>
### unit (intmap)

| Type | Default |
|---|---|
| `"Hz"`, `"kHz"`, `"MHz"` or `"GHz"` | `"kHz"` |

The unit of this entry's tensor; each entry may have its own. Units are case-sensitive (`"mhz"` is refused). A Defect
`zfs` is MHz, so the log notes a nonzero tensor written without `unit`.

```json
"unit": "MHz"
```

<a id="tag_intmap_frame"></a>
### frame (intmap)

| Type | Default |
|---|---|
| `"file"` or `"qubit"` | none: required in every written entry with [frame.rotation](InputFile_General.md#tag_frame_rotation) on |

The frame the tensor components are written in. With the rotation off it is not used (reported).

| Value | Meaning |
|---|---|
| `"file"` | the bath file's frame: the tensor is turned with the positions, T → R T R<sup>T</sup> |
| `"qubit"` | the computational frame, whose z is `qubit_axis`: used as written |

`"qubit"` is the one computational frame of the run, not a qubit's own local frame. A second qubit with another
orientation is described by its own ZFS tensor in that common frame, off-diagonal elements included. An entry that is
not written (zero, or the point-dipole tensor of a pair) is always turned.

```json
"frame": "file"
```

<a id="qb_states"></a>
## 4. The states

<a id="tag_alpha_state"></a>
### alpha_state

| Type | Default |
|---|---|
| `{"ms": [...]}`, `{"index": n}` or `{"vector": [...]}` | `{"ms": [1, 1, ...]}` (m<sub>s</sub> = +1 on every qubit) |

The state \|α⟩ of the qubits, written in exactly one of three ways:

| Form | Meaning |
|---|---|
| `"ms": [m1, m2, ...]` | one m<sub>s</sub> per qubit, in the order of [qubit](#tag_qubit) (a list also for one qubit); \|α⟩ is their product |
| `"index": n` | the n-th eigenstate of the qubit Hamiltonian (zero-field splitting, Zeeman, intmap) by **descending** energy, from 0 |
| `"vector": [...]` | the state over the whole qubit space (kron of the qubits, the last varying fastest); each component a number or `[re, im]`; normalized when read |

- The order of `"index"` changes with the field: an NV at 500 G along its axis orders m<sub>s</sub> = +1, −1, 0, so
  `{"index": 2}` and `{"index": 1}` is the 0 ↔ −1 transition. To name a sublevel, use `"ms"`. The log names the
  m<sub>s</sub> each state mostly is (`psia is mostly`, `psib is mostly`).
- [pulse.per_qubit](InputFile_Pulse.md#tag_pulse_per_qubit) needs both states as `"ms"`.
- Each state is decided on its own: `beta_state` may be written another way.

```json
"alpha_state": {"ms": [1]}
"alpha_state": {"index": 0}
"alpha_state": {"vector": [[0.7071, 0], [0, 0.7071], 0]}
```

<a id="tag_beta_state"></a>
### beta_state

| Type | Default |
|---|---|
| as [alpha_state](#tag_alpha_state) | `{"ms": [0, 0, ...]}` (m<sub>s</sub> = 0 on every qubit) |

The state \|β⟩, written as [alpha_state](#tag_alpha_state).

```json
"beta_state": {"ms": [0]}
```

<a id="tag_psi0"></a>
### psi0

| Type | Default |
|---|---|
| list of numbers or `[re, im]` | (\|α⟩ + \|β⟩), normalized |

The initial state of the qubits, over the whole qubit space, for another superposition than the default. gCCE and
megCCE read it; cCCE and meCCE work from \|α⟩ and \|β⟩ alone and report it as not used. It is used as written, not
normalized.

```json
"psi0": [0.70710678, 0.70710678, 0.0]
```
