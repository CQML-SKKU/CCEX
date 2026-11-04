# Defect tags
<a id="DefectTags"></a>

The `Defect` block describes the paramagnetic defects of the bath, such as P1 centres: a list of defect species, the
model. The files of the defects' axes and states are keys of the configuration,
[defect_axis_file](InputFile_General.md#tag_defect_axis_file) and
[defect_state_file](InputFile_General.md#tag_defect_state_file). Each key has a small table (its type, default and
unit), what it does, and an example. Old key names are in [Old names](InputFile_OldNames.md).

1. [The bath model](#DefectModel): what a defect species is, and the two forms of its values.
2. [The block](#df_block): `Defect` and `species`.
3. [A species](#df_species_keys): `name`, `n_onsite_spins`, `n_axis_states`, `frame`.
4. [The on-site nuclear spins](#df_onsite): `types`, `spins`, `gyros`, `quadrupole_moments`, `nuclear_positions`,
   `hyperfine`, `efg`.
5. [The electron spin](#df_electron): `zfs` (and its `axis`) and `detuning`.

<a id="DefectModel"></a>
## 1. The bath model

A bath spin whose name in the bath file is a species' `name` is a paramagnetic defect: an electron spin (its spin and
γ from the gyrofile) with a set of indexed *axis states*, `1 .. n_axis_states`. An axis state stands for one state of
the electron's surroundings, for example a Jahn-Teller axis of a P1 centre, or a frozen state of its nucleus. Each
defect spin is given one axis state, from [defect_axis_file](InputFile_General.md#tag_defect_axis_file) or at random
from `configuration.seed`, and its tensors are those of that state.

A nucleus on the defect site enters in one of two ways:

- through the electron's [detuning](#tag_df_detuning) only, with `n_onsite_spins` 0: the hyperfine coupling is a
  static shift A·m<sub>I</sub> of the electron, with one axis state per axis and m<sub>I</sub> (the twelve detunings of
  a P1 with <sup>14</sup>N);
- as an on-site spin of its own (`n_onsite_spins` > 0): it joins its electron in every cluster the defect is in.

With [configuration.s_range](InputFile_General.md#tag_s_range) the bath is sampled in product states. The electrons'
states come from `bath_state_file` or the seed, the on-site spins' from `defect_state_file` or the seed, and the spins
outside a cluster act as a mean field ([bath_bath_mean_field](InputFile_General.md#tag_bath_bath_mean_field)).
[Basic 06](../../example/tutorial/basic/06_defect_bath/README.md) runs a P1 bath.

<a id="tag_df_forms"></a>
**Two forms of the values.**

- The list form `[[axis state, name, value], ...]` gives a value per axis state.
- The object form gives one value for every axis state (`zfs`), or wraps the list with a unit:
  `{"values": [...], "unit": ...}` (`"value"` is taken too, but not both).

Units are case-sensitive (`"mhz"` is refused). A plain frequency here is MHz, while a qubit's
([Qubit.qubit[].detuning](InputFile_Qubit.md#tag_detuning), [Qubit.intmap](InputFile_Qubit.md#tag_intmap)) is kHz.

<a id="df_block"></a>
## 2. The block

<a id="tag_Defect"></a>
### Defect

| Type | Default |
|---|---|
| `{"species": [...]}` | no block: no bath spin is a defect |

The paramagnetic defects of the bath. A whole P1-like species with one on-site <sup>15</sup>N in four axis states:

```json
"configuration": {
    "bathfile":          ["./bath_DiaP1_1ppm_1"],
    "defect_axis_file":  ["./bathJT_DiaP1_1ppm_1"],
    "defect_state_file": ["./stateEx_DiaP1_1ppm_"]
},
"Defect": {"species": [{
    "name": "NV0", "n_onsite_spins": 1, "n_axis_states": 4,
    "types": ["15N"], "spins": [0.5], "gyros": [-2.71261804],
    "nuclear_positions": [[1, "15N", [ 1.36404, -0.78843, -0.56598]],
                          [2, "15N", [-1.36801, -0.79025, -0.56701]],
                          [3, "15N", [ 0.00000,  1.57603, -0.56511]],
                          [4, "15N", [ 0.00000, -0.00733,  1.67433]]],
    "hyperfine": [[1, "15N", [-10.33219,  1.01563,  1.50162,  1.01563, -5.18317, -0.16410,  1.50162, -0.16410,  -6.50270]],
                  [2, "15N", [ -7.34999, -2.73741, -0.89293, -2.73741, -8.16537, -1.21839, -0.89293, -1.21839,  -6.50270]],
                  [3, "15N", [ -5.59086,  1.72177, -0.60870,  1.72177, -9.92450,  1.38249, -0.60870,  1.38249,  -6.50270]],
                  [4, "15N", [ -6.01395,  0.00000,  0.00000,  0.00000, -5.59086, -1.82620,  0.00000, -1.82620, -10.41325]]]
}]}
```

<a id="tag_df_species"></a>
### species

| Type | Default |
|---|---|
| list of objects | required in the block |

The defect species, one entry each. A bath spin whose name is a species' `name` is that defect.

```json
"species": [{"name": "NV_0", "n_axis_states": 3, "frame": "file",
             "zfs": {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"},
             "detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}}]
```

<a id="df_species_keys"></a>
## 3. A species

<a id="tag_df_name"></a>
### name

| Type | Default |
|---|---|
| string | required |

The name of the defect, as in the bath file and the gyrofile.

```json
"name": "P1"
```

<a id="tag_df_n_onsite_spins"></a>
### n_onsite_spins

| Type | Default |
|---|---|
| integer | `0` |

The number of on-site nuclear spins. Each joins its electron in every cluster the defect is in, so a cluster of one
P1 is its electron and its <sup>14</sup>N. With 0 the nucleus enters only through the electron's
[detuning](#tag_df_detuning). [Advanced 01](../../example/tutorial/advanced/01_P1/README.md) compares the two.

```json
"n_onsite_spins": 1
```

<a id="tag_df_n_axis_states"></a>
### n_axis_states

| Type | Default |
|---|---|
| integer | `0` |

The number of axis states. The indexed entries of `nuclear_positions`, `hyperfine`, `efg`, `zfs` and `detuning` must
name an axis state in 1 .. `n_axis_states`; a larger index is refused. An object-form `zfs` holds for every axis
state.

```json
"n_axis_states": 4
```

<a id="tag_df_frame"></a>
### frame (species)

| Type | Default |
|---|---|
| `"file"` or `"qubit"` | none: required with [frame.rotation](InputFile_General.md#tag_frame_rotation) on, for a species that writes `nuclear_positions`, `hyperfine`, `efg` or `zfs` |

The frame of this species' vectors and tensors: the `nuclear_positions`, the `hyperfine`, `efg` and `zfs` components
and the zfs `axis`. With the rotation off, or a species without these, it is not used (reported). A detuning has no
direction.

| Value | Meaning |
|---|---|
| `"file"` | the bath file's frame (e.g. the crystal frame): turned with the positions |
| `"qubit"` | the computational frame (z = `qubit_axis`): used as written |

```json
"frame": "file"
```

<a id="df_onsite"></a>
## 4. The on-site nuclear spins

The lists of `types`, `spins`, `gyros` and `quadrupole_moments` have one entry per on-site spin. The indexed lists of
`nuclear_positions`, `hyperfine` and `efg` have one entry per on-site spin and axis state,
`n_onsite_spins` × `n_axis_states` in all.

<a id="tag_df_types"></a>
### types

| Type | Default |
|---|---|
| list of strings | required with `n_onsite_spins` > 0 |

The names of the on-site nuclear spins.

```json
"types": ["14N"]
```

<a id="tag_df_spins"></a>
### spins

| Type | Default |
|---|---|
| list of numbers | required with `n_onsite_spins` > 0 |

The spins of the on-site nuclear spins.

```json
"spins": [1.0]
```

<a id="tag_df_gyros"></a>
### gyros

| Type | Default | Unit |
|---|---|---|
| `{"values": [...], "unit": ...}`, or a list | required with `n_onsite_spins` > 0 | a plain list is rad kHz/G (= rad/ms/G) |

The gyromagnetic ratios of the on-site nuclear spins. The units are `radkHz/G`, `rad/ms/G`, `kHz/G`, `MHz/G`, `MHz/T`,
`Hz/T`, `radMHz/T`, `rad/us/T` and `rad/s/T`. Those without `rad` are γ/2π.

```json
"gyros": {"values": [3.0777051853], "unit": "MHz/T"}
```

<a id="tag_df_quadrupole_moments"></a>
### quadrupole_moments

| Type | Default | Unit |
|---|---|---|
| `{"values": [...], "unit": ...}`, or a list | `0` for every spin | a plain list is fm² (10<sup>−30</sup> m² = 10 mbarn) |

The quadrupole moments eQ of the on-site nuclear spins. The units are `fm^2`, `m^2`, `cm^2`, `barn`, `mbarn` and other
spellings of these.

```json
"quadrupole_moments": {"values": [20.44], "unit": "mbarn"}
```

<a id="tag_df_nuclear_positions"></a>
### nuclear_positions

| Type | Default | Unit |
|---|---|---|
| indexed list `[[state, name, [x, y, z]], ...]` | `[0, 0, 0]` | Å |

The position of each on-site nuclear spin relative to its electron, per axis state.

```json
"nuclear_positions": [[1, "14N", [0.0, 0.0, 1.5]], [2, "14N", [1.4, 0.0, -0.5]]]
```

<a id="tag_df_hyperfine"></a>
### hyperfine

| Type | Default | Unit |
|---|---|---|
| `{"values": [...], "unit": ...}`, or an indexed list | zero for an entry not given | MHz; `Hz`, `kHz`, `MHz`, `GHz` |

The hyperfine tensor between each on-site nuclear spin and its electron, per axis state. A tensor is 9 numbers (row
by row) or a 3 × 3 list.

```json
"hyperfine": {"values": [[1, "14N", [[81.312, 0, 0], [0, 81.312, 0], [0, 0, 114.0264]]]], "unit": "MHz"}
```

<a id="tag_df_efg"></a>
### efg

| Type | Default | Unit |
|---|---|---|
| `{"values": [...], "unit": ...}`, or an indexed list | zero for an entry not given (no quadrupole term) | Hartree/Bohr² (`a.u.`); `V/m^2`, `V/A^2` |

The electric field gradient tensor at each on-site nuclear spin, per axis state, in the layout of
[hyperfine](#tag_df_hyperfine). With `quadrupole_moments` it gives the quadrupole term.

```json
"efg": {"values": [[1, "14N", [[0.521910, 0, 0], [0, 0.522771, -0.000086], [0, -0.000086, -1.044681]]]], "unit": "Hartree/Bohr^2"}
```

<a id="df_electron"></a>
## 5. The electron spin

<a id="tag_df_zfs"></a>
### zfs

| Type | Default | Unit |
|---|---|---|
| `{"D", "axis", "unit"}`, `{"tensor", "unit"}`, or an indexed list | zero | MHz; the object forms take `Hz`, `kHz`, `MHz`, `GHz` |

The zero-field splitting of the electron spin, in one of three forms:

| Form | Meaning |
|---|---|
| `{"D": D, "axis": [x, y, z], "unit": ...}` | the axial tensor D (n̂ n̂<sup>T</sup> − 1/3) about the normalized [axis](#tag_df_axis) n̂, for every axis state |
| `{"tensor": [[...], [...], [...]], "unit": ...}` | one tensor (3 × 3, or 9 numbers row by row) for every axis state |
| `[[state, "e", [9 numbers]], ...]` | a tensor per axis state, in MHz |

The axial form has no E, since an axis does not fix the transverse x and y; for E other than 0 write the tensor.

```json
"zfs": {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"}
"zfs": {"tensor": [[0.0, 956.6667, 956.6667], [956.6667, 0.0, 956.6667], [956.6667, 956.6667, 0.0]], "unit": "MHz"}
"zfs": [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
```

<a id="tag_df_axis"></a>
### axis (zfs)

| Type | Default |
|---|---|
| `[x, y, z]` | required with `D` |

The symmetry axis of the D form of [zfs](#tag_df_zfs), written inside it, in the species' [frame](#tag_df_frame). It
is normalized when read. A species-level `axis` is refused with this place.

```json
"zfs": {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"}
```

<a id="tag_df_detuning"></a>
### detuning (species)

| Type | Default | Unit |
|---|---|---|
| `{"values": [...], "unit": ...}`, or an indexed list | `0` for every axis state not given | **MHz**; `Hz`, `kHz`, `MHz`, `GHz` |

A static shift of the electron spin per axis state, detuning × S<sub>z</sub>. It stands in for couplings that are
taken as static, such as **S**·**A**·**I** ≈ A m<sub>I</sub> S<sub>z</sub> of an on-site nucleus. A qubit's
[detuning](InputFile_Qubit.md#tag_detuning) is kHz.

```json
"detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}
```
