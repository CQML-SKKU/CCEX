# Qubit tags
<a id="QubitTags"></a>

The `Qubit` block describes the qubits, the interactions between them, and the two states whose coherence is
computed. The Overhauser field of the bath on the qubit is a key of the `hamiltonian` block,
[qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field). Old key names are in
[Upgrading old inputs](InputFile_OldNames.md).

1. [The qubits](#qb_qubits): `qubit` and its keys `name`, `spin`, `gyro`, `xyz` or `file`, `detuning`.
2. [The interactions](#qb_intmap): `intmap` and its keys `between`, `tensor`, `unit`, `frame`.
3. [The states](#qb_states): `alpha_state`, `beta_state` and `psi0`.

<a id="qb_block"></a>
<a id="tag_Qubit"></a>
The block is required. The first example is one NV read from a defect file.

```json
"Qubit": {"qubit": [{"file": "bath/bath_Diamond_NV_natab_defect"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}}
```

The second has two qubits with their interactions written out. The first `intmap` entry is the zero-field splitting of Q2, the second
the Q1–Q2 coupling. The `vector` and `psi0` span the 2 × 3 = 6 states of kron(Q1, Q2), the last qubit varying fastest.

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

<a id="qb_qubits"></a>
## 1. The qubits

<a id="tag_qubit"></a>
### qubit

| Type | Default |
|---|---|
| list of objects | required (at least one entry) |

One entry per qubit, with its name, spin, gyromagnetic ratio, position (`xyz` or `file`) and detuning (the table
below). The states give the qubits' levels in this order.

```json
"qubit": [{"name": "Q1", "spin": 0.5, "xyz": [1.0, 2.0, 3.0]}, {"name": "Q2", "spin": 1.0, "xyz": [4.0, 5.0, 6.0]}]
```

<a id="tag_name"></a>
<a id="tag_spin"></a>
<a id="tag_gyro"></a>
<a id="tag_xyz"></a>
<a id="tag_Qubit_file"></a>
<a id="tag_detuning"></a>
The keys of an entry:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | string | `"q0"`, `"q1"`, ... | The qubit's name, used by `between` and by a qubit jump operator. |
| `spin` | number | `1.0` | The spin of the qubit. |
| `gyro` | number in rad kHz/G, or `{"value", "unit"}` with the units of [Defect gyros](InputFile_Defect.md#tag_df_gyros) | `-17608.597050` (the electron) | The gyromagnetic ratio. |
| `xyz` | `[x, y, z]` in Å, or `{"value", "unit"}` with `A` or `nm` | required unless `file` is given | The position, in the frame of the bath file. |
| `file` | string (a file path) | none (give `xyz` or `file`, not both) | A file whose first line of three numbers is x y z in Å, such as the defect file of BathGenerator. The path is relative to the working directory. |
| `detuning` | number in **kHz** (MHz in the Defect block), or `{"value", "unit"}` with `Hz`, `kHz`, `MHz` or `GHz` | `0` | Shifts the qubit levels by adding detuning × S<sub>z</sub> to the qubit Hamiltonian. The phase of L turns faster and \|L\| stays the same. With gCCE the resonances of the qubit with the bath spins move. |

<a id="qb_intmap"></a>
## 2. The interactions

<a id="tag_intmap"></a>
### intmap

| Type | Default |
|---|---|
| list of objects | no self terms, and point-dipole tensors between different qubits |

The interactions between the qubits, and of a qubit with itself (its zero-field splitting, an entry that names the
qubit twice). A pair of different qubits that is not given gets the point-dipole tensor from their positions.

```json
"intmap": [{"between": ["q0", "q0"], "unit": "MHz", "tensor": [[-960, 0, 0], [0, -960, 0], [0, 0, 1920]]}]
```

<a id="tag_tensor"></a>
<a id="tag_between"></a>
<a id="tag_intmap_unit"></a>
<a id="tag_intmap_frame"></a>
The keys of an entry:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `between` | two names | required | The two qubits, the first not after the second in the order of `qubit`. The same name twice gives a zero-field splitting. |
| `tensor` | 3 × 3 list | zero, which also removes the point-dipole default of that pair | The interaction tensor of the pair in frequency units. CCEX supplies the 2π. |
| `unit` | `"Hz"`, `"kHz"`, `"MHz"` or `"GHz"` (case-sensitive) | `"kHz"` | The unit of this entry's tensor. |
| `frame` | `"file"` or `"qubit"` | required with [frame.rotation](InputFile_General.md#tag_frame_rotation) on | The frame of the tensor. `"file"` is the bath file's frame, turned with the positions, and `"qubit"` is the computational frame, used as written. |

<a id="qb_states"></a>
## 3. The states

<a id="tag_alpha_state"></a>
<a id="tag_beta_state"></a>
### alpha_state and beta_state

| Key | Type | Default |
|---|---|---|
| `alpha_state` | `{"ms": [...]}`, `{"index": n}` or `{"vector": [...]}` | `{"ms": [1, 1, ...]}` |
| `beta_state` | the same forms | `{"ms": [0, 0, ...]}` |

The two qubit states \|α⟩ and \|β⟩ whose coherence the run computes. Both keys take the same three forms.

| Form | Meaning |
|---|---|
| `"ms": [m1, m2, ...]` | One m<sub>s</sub> per qubit, in the order of [qubit](#tag_qubit). |
| `"index": n` | The n-th eigenstate of the qubit Hamiltonian, counted from 0 by descending energy. |
| `"vector": [...]` | The state over the whole qubit space. Each component is a number or `[re, im]`. |

The defaults hold for integer spins only, so give both states for a spin-1/2 qubit. The order of `"index"` changes with
the field, so use `"ms"` to name a sublevel. [pulse.per_qubit](InputFile_Pulse.md#tag_pulse_per_qubit) needs both
states as `"ms"`.

```json
"alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}
"alpha_state": {"vector": [[0.7071, 0], [0, 0.7071], 0]}
```

<a id="tag_psi0"></a>
For gCCE and megCCE the initial state can also be given.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `psi0` | list of numbers or `[re, im]` | (\|α⟩ + \|β⟩), normalised | The initial state over the whole qubit space, for another superposition. |
