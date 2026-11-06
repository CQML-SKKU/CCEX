# Defect tags
<a id="DefectTags"></a>

The `Defect` block describes paramagnetic defects in the bath, such as P1 centres or NV centres. It adds what the bath
file and the gyrofile cannot hold. Old key names are in [Upgrading old inputs](InputFile_OldNames.md).

1. [How a defect is described](#DefectModel): where its information comes from, and what an axis state is.
2. [The keys of a species](#df_block): `name`, `n_axis_states`, `frame` and `zfs`.
3. [The nucleus on the site](#df_nucleus): as a static detuning or as an on-site spin.
4. [Defects with several orientations](#df_orientations): P1 and NV.

<a id="DefectModel"></a>
## 1. How a defect is described

The information on a defect comes from four places, and the `Defect` block is one of them:

| Where | What | For a P1 |
|---|---|---|
| the bath file ([configuration.bathfile](InputFile_General.md#tag_bathfile)) | The position and the name of each spin. | `1317.8 159.4 1842.6 P1` |
| the gyrofile ([configuration.gyrofile](InputFile_General.md#tag_gyrofile)) | The spin and γ of each name. | `P1 0.5 -17608.597050` |
| the `Defect` block | What a defect of that name is, with its axis states, the tensors or the detuning of each, and the nuclear spins on its site. | Four Jahn-Teller axes and one ¹⁴N. |
| [defect_axis_file](InputFile_General.md#tag_defect_axis_file), [defect_state_file](InputFile_General.md#tag_defect_state_file), or the seed | The axis state of each defect, and the states of its on-site spins. | This P1 on axis 2, that one on axis 4. |

A bath spin whose name in the bath file is the `name` of a species is that defect. It is an electron spin, with its
spin and γ from the gyrofile. The defect sits in one of several *axis states*, for example the four Jahn-Teller axes
of a P1 centre, and each axis state has its own tensors or detuning. Each defect gets its axis state from
[defect_axis_file](InputFile_General.md#tag_defect_axis_file) or draws it from the seed.

<a id="df_block"></a>
<a id="tag_Defect"></a>
<a id="tag_df_species"></a>
<a id="df_species_keys"></a>
<a id="tag_df_name"></a>
<a id="tag_df_n_axis_states"></a>
<a id="tag_df_frame"></a>
## 2. The keys of a species

`Defect` holds one key, `species`, a list with one entry per defect species. Without the block no bath spin is a
defect. The keys of a species:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `name` | string | required | The name of the defect, as in the bath file and the gyrofile. |
| `n_axis_states` | integer | `0` | The number of axis states. |
| `frame` | `"file"` or `"qubit"` | required with [frame.rotation](InputFile_General.md#tag_frame_rotation) on, for a species with positions or tensors | The frame of the species' positions and tensors. `"file"` is the bath file's frame (turned), and `"qubit"` is the computational frame (used as written). |
| `zfs` | object or indexed list, in MHz or with a `unit` (`Hz`, `kHz`, `MHz`, `GHz`) | zero | The zero-field splitting of the electron, for a defect of spin 1 or more, in one of the three forms below. For E other than 0, write the tensor. |

<a id="tag_df_forms"></a>
Values that depend on the axis state are written as an indexed list, `[[axis state, name, value], ...]`. The axis
state is one of 1 .. `n_axis_states`, and the name is `"e"` for the electron or the name of an on-site spin. The same
list can be wrapped with a unit, `{"values": [...], "unit": ...}`. A plain frequency in this block is MHz, while in the
Qubit block it is kHz.

<a id="tag_df_zfs"></a>
<a id="tag_df_axis"></a>
The three forms of `zfs`. The two object forms need `n_axis_states` of 1 or more, and `axis` is in the species' frame.

| Form | Meaning |
|---|---|
| `{"D": D, "axis": [x, y, z], "unit": ...}` | The axial tensor D (n̂ n̂<sup>T</sup> − 1/3) about the normalised `axis` n̂, for every axis state. |
| `{"tensor": [[...], [...], [...]], "unit": ...}` | One tensor, 3 × 3 or 9 numbers row by row, for every axis state. |
| `[[state, "e", [9 numbers]], ...]` | A tensor per axis state, in MHz. |

```json
"zfs": {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"}
"zfs": [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
```

<a id="df_nucleus"></a>
## 3. The nucleus on the site

The nucleus on the defect site, such as the ¹⁴N of a P1 or of an NV, is included in one of two ways. Both are compared
for P1 baths in H. Park et al., Adv. Opt. Mater. (2026), [doi:10.1002/adom.71146](https://doi.org/10.1002/adom.71146) (arXiv:2503.05404).

| | 3.1 Static detuning | 3.2 On-site spin |
|---|---|---|
| `n_onsite_spins` | 0 (the default) | 1 or more |
| The nucleus | a fixed shift of the electron in each axis state | a spin of its own, in every cluster with its electron |
| Its m<sub>I</sub> | part of the axis state | the state of the on-site spin ([defect_state_file](InputFile_General.md#tag_defect_state_file) or the seed) |
| Its flips | left out | included |
| Cost | as a lone electron spin | larger clusters |
| Name in Park et al. | effective P1 model, CCE-2 | full P1 model, CCE-2e2n |

In that work the Hahn-echo T2 of the two models differs by 17.5 %, 16.2 % and 6.4 % at 1, 10 and 100 ppm P1.

<a id="df_electron"></a>
<a id="tag_df_detuning"></a>
### 3.1 Static detuning

The nucleus is not a spin of the calculation. Each axis state includes a value of m<sub>I</sub>, and in each the
hyperfine field of the nucleus shifts the electron by a fixed detuning.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `detuning` | indexed list in **MHz**, or `{"values", "unit"}` with `Hz`, `kHz`, `MHz`, `GHz` | `0` | The static shift detuning × S<sub>z</sub> of the electron in each axis state, such as the A m<sub>I</sub> of the nucleus. |

The P1 bath of [Basic 06](../../example/tutorial/basic/06_defect_bath/README.md) has 12 axis states, the 4 Jahn-Teller
axes times the 3 values of m<sub>I</sub> of the ¹⁴N:

```json
"Defect": {"species": [{
    "name": "P1", "n_axis_states": 12,
    "detuning": [[1, "e",  84.946496], [2, "e",  84.946496], [3,  "e",  84.946496],
                 [4, "e",   0.0],      [5, "e",   0.0],      [6,  "e",   0.0],
                 [7, "e", -84.946496], [8, "e", -84.946496], [9,  "e", -84.946496],
                 [10, "e", 114.0264],  [11, "e",  0.0],      [12, "e", -114.0264]]}]}
```

<a id="df_onsite"></a>
<a id="tag_df_n_onsite_spins"></a>
### 3.2 On-site spin

The nucleus is a spin of its own, with its hyperfine and quadrupole tensors, and joins its electron in every cluster.
Its state is drawn from the seed or read from [defect_state_file](InputFile_General.md#tag_defect_state_file). This is
also the approach of H. Park et al., npj Quantum Inf. 8, 95 (2022).

The lists of `types`, `spins`, `gyros` and `quadrupole_moments` have one entry per on-site spin. The indexed lists of
`nuclear_positions`, `hyperfine` and `efg` have one entry per on-site spin and axis state, `n_onsite_spins` ×
`n_axis_states` in all.

<a id="tag_df_types"></a>
<a id="tag_df_spins"></a>
<a id="tag_df_gyros"></a>
<a id="tag_df_quadrupole_moments"></a>
<a id="tag_df_nuclear_positions"></a>
<a id="tag_df_hyperfine"></a>
<a id="tag_df_efg"></a>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `n_onsite_spins` | integer | `0` | The number of on-site nuclear spins. |
| `types` | list of strings | required | The names of the on-site spins. |
| `spins` | list of numbers | required | Their spins. |
| `gyros` | list in rad kHz/G, or `{"values", "unit"}` | required | Their gyromagnetic ratios. |
| `quadrupole_moments` | list in fm² (10<sup>−30</sup> m² = 10 mbarn), or `{"values", "unit"}` | `0` | Their quadrupole moments eQ. |
| `nuclear_positions` | indexed list `[[state, name, [x, y, z]], ...]` in Å | `[0, 0, 0]` | Each spin's position relative to its electron, per axis state. |
| `hyperfine` | indexed list in MHz, or `{"values", "unit"}` | zero for an entry not given | The hyperfine tensor with its electron, per axis state, as 9 numbers row by row. The object form also takes a 3 × 3 list. |
| `efg` | indexed list in Hartree/Bohr², or `{"values", "unit"}` | zero for an entry not given | The electric field gradient, per axis state, written as `hyperfine`. With `quadrupole_moments` it gives the quadrupole term. |

The units of the object form:

- `gyros`: `radkHz/G`, `rad/ms/G`, `kHz/G`, `MHz/G`, `MHz/T`, `Hz/T`, `radMHz/T`, `rad/us/T`, `rad/s/T` (those
  without `rad` are γ/2π).
- `quadrupole_moments`: `fm^2`, `m^2`, `cm^2`, `barn`, `mbarn`.
- `hyperfine`: `Hz`, `kHz`, `MHz`, `GHz`.
- `efg`: `Hartree/Bohr^2` (`a.u.`), `V/m^2`, `V/A^2`.

[Advanced 05](../../example/tutorial/advanced/05_P1/README.md) runs the same P1 bath this way, with 4 axis states
(the m<sub>I</sub> are now states of the ¹⁴N spin):

```json
"Defect": {"species": [{
    "name": "P1", "n_onsite_spins": 1, "n_axis_states": 4,
    "types": ["14N"], "spins": [1], "gyros": [1.9337792], "quadrupole_moments": [2.044],
    "hyperfine": [[1, "14N", [  81.312000,    0.000000,    0.000000,
                                 0.000000,   81.312000,    0.000000,
                                 0.000000,    0.000000,  114.026400]],
                 [2, "14N", [ 103.121928,   12.591967,    8.903264,
                               12.591967,   88.581975,    5.140302,
                                8.903264,    5.140302,   84.946496]],
                 [3, "14N", [ 103.121923,  -12.591970,   -8.903263,
                              -12.591970,   88.581981,    5.140303,
                               -8.903263,    5.140303,   84.946496]],
                 [4, "14N", [  81.312000,   -0.000001,    0.000000,
                               -0.000001,  110.391904,  -10.280603,
                                0.000000,  -10.280603,   84.946496]]],
    "efg": [[1, "14N", [  0.521910,   0.000000,   0.000000,
                          0.000000,   0.522771,  -0.000086,
                          0.000000,  -0.000086,  -1.044681]],
           [2, "14N", [ -0.522214,  -0.603404,  -0.426335,
                        -0.603404,   0.174349,  -0.246178,
                        -0.426335,  -0.246178,   0.347865]],
           [3, "14N", [ -0.522354,   0.603323,   0.426364,
                         0.603323,   0.174489,  -0.246128,
                         0.426364,  -0.246128,   0.347865]],
           [4, "14N", [  0.522771,   0.000081,  -0.000029,
                         0.000081,  -0.870636,   0.492306,
                        -0.000029,   0.492306,   0.347865]]]}]}
```

<a id="df_orientations"></a>
## 4. Defects with several orientations

A defect whose axis can point in several directions can be written in two ways.

- **One species with axis states.** The directions are axis states of one species, and each defect draws its own. The
  P1 examples of section 3 treat the four Jahn-Teller axes this way.
- **One species per direction.** The bath file names each defect by its direction, and each species has its own
  axis. This suits an NV⁻ in the bath (spin 1, D = 2870 MHz), whose zero-field splitting follows its ⟨111⟩ direction.
  The example below has four species `NV_0` … `NV_3`, with the D form of `zfs` about each axis, in a bath file written
  in the cubic axes of diamond. The axis states are then left for the m<sub>I</sub> of the ¹⁴N, as in section 3.1.

```json
"frame": {"rotation": {"qubit_axis": [1, 1, 1]}},
"Defect": {"species": [
    {"name": "NV_0", "frame": "file", "n_axis_states": 3, "zfs": {"D": 2870.0, "axis": [ 1,  1,  1], "unit": "MHz"},
     "detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}},
    {"name": "NV_1", "frame": "file", "n_axis_states": 3, "zfs": {"D": 2870.0, "axis": [ 1, -1, -1], "unit": "MHz"},
     "detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}},
    {"name": "NV_2", "frame": "file", "n_axis_states": 3, "zfs": {"D": 2870.0, "axis": [-1,  1, -1], "unit": "MHz"},
     "detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}},
    {"name": "NV_3", "frame": "file", "n_axis_states": 3, "zfs": {"D": 2870.0, "axis": [-1, -1,  1], "unit": "MHz"},
     "detuning": {"values": [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]], "unit": "MHz"}}]}
```

The gyrofile then gives each of the four names spin 1 and the electron's γ (`NV_0 1.0 -17608.597050`, and so on).
