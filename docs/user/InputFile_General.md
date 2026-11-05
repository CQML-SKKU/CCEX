# General tags
<a id="GeneralTags"></a>

The method, the numbers of the calculation, the files of the bath, the frame, and the `hamiltonian` block. Keys that need
explaining have a section of their own, and the simple keys of each group are in one table. Old key names are in
[Upgrading old inputs](InputFile_OldNames.md).

1. [Calculation](#gen_calc): `method`, `quantity`, `order`, `bfield`, `rbath`, `rdip`, `deltat`, `nstep`.
2. [configuration](#tag_configuration): the bath files, the states and the seed.
3. [frame](#tag_frame): the rotation of the bath file's frame.
4. [hamiltonian](#tag_hamiltonian): the block of the Hamiltonian's terms, one section per method.
    1. [common: mean fields](#gen_mf): the static fields of the sampled bath states.
    2. [common: hyperfine](#gen_hf): the hyperfine tensors, point dipole or DFT.
    3. [common: quadrupole](#gen_qd): the quadrupole tensors from DFT files.
    4. [cCCE](#gen_ham_method): the terms `hfmedi`, `knight` and `delta_A` (read by cCCE and meCCE).
    5. [meCCE and megCCE](#gen_ham_me): the jump operators.
    6. [noise](#gen_ham_noise): the spectrum and what method noise computes.
    7. [pseudospin](#gen_ham_pseudospin): how the channels of a pair are evolved.
5. [Reference](#gen_reference): keys read and not used, and how to turn a part off.

<a id="gen_calc"></a>
## 1. Calculation

<a id="tag_method"></a>
### method

| Type | Default |
|---|---|
| string (case is not significant) | `"cCCE"` |

The CCE method. How the bath is cut into clusters is a separate key,
[cluster.mode](InputFile_Cluster.md#tag_cluster_mode). The mode `pcce` runs with cCCE, gCCE, meCCE and megCCE, `bundle`
with cCCE and gCCE, and `dsj`, `itb` and `dsjitb` with cCCE only.

| Value | Method |
|---|---|
| `cCCE` | Conventional CCE. The qubit is projected out into its two branches. |
| `gCCE` | Generalized CCE. The qubit is in the cluster Hamiltonian. |
| `meCCE` | cCCE with Lindblad jump operators on the bath spins, on the CPU and for one qubit ([master-equation CCE](../methods/master-equation.md)). |
| `megCCE` | gCCE with jump operators on the bath spins and on the qubit, on the CPU. |
| `noise` | The Gaussian coherence exp(−χ) of the qubit's phase noise, for one qubit and `cluster.mode` standard ([spectrum](#tag_noise_spectrum)). |
| `pseudospin` | The pairs of cCCE in the pseudospin model of nuclear flip-flops. It needs one qubit, the ensemble, order 2, `bfield` `[0, 0, B]` and `execution.propagator` eigen, and runs on the CPU ([channels](#tag_pseudospin_channels)). |

```json
"method": "gCCE"
"method": "cCCE", "cluster": {"mode": "pcce"}
```

<a id="tag_quantity"></a>
### quantity

| Type | Default |
|---|---|
| string | `"coherence"` |

| Value | Computes |
|---|---|
| `coherence` | The decay curve L(t) between \|α⟩ and \|β⟩, the usual choice. |
| `dm` | The qubit density matrix, for populations or several qubit levels. It needs gCCE with `cluster.mode` standard or pcce ([The density matrix](OutputFiles.md#out_dm)). |

```json
"quantity": "dm"
```

<a id="tag_bfield"></a>
### bfield

| Type | Default | Unit |
|---|---|---|
| number (Bz) or `[Bx, By, Bz]` | required (`0` for no field) | G, or `{"value", "unit"}` with `G`, `mT`, `T` ([Units](InputFile.md#Units)) |

The magnetic field, in the computational frame. With [frame.rotation](#tag_frame_rotation) on, only `[0, 0, B]` is
taken, a field along `qubit_axis`.

```json
"bfield": 500
"bfield": {"value": 50, "unit": "mT"}
```

<a id="tag_order"></a>
<a id="tag_rbath"></a>
<a id="tag_rdip"></a>
<a id="tag_deltat"></a>
<a id="tag_nstep"></a>
<a id="tag_rdip_ladder"></a>
### The other numbers

| Key | Type | Default | Meaning |
|---|---|---|---|
| `order` | integer | required (2 for method pseudospin, 1 for noise with a model spectrum) | The largest cluster, in bath spins (standard), partitions (pcce) or bundles (bundle). |
| `rbath` | number or `[inner, outer]`, Å | required (inner 0) | The bath spins within this distance of the qubit are kept. `inner` leaves out the closer ones. |
| `rdip` | number or `[inner, outer]`, Å | required, except with `cluster.mode` bundle whose `connections` is not `"distance"` | Two bath spins closer than this can be in one cluster. `inner` leaves out the closer pairs. |
| `deltat` | number, ms | required | The time step. Each time is the total time of the sequence, so a Hahn echo of t has its π pulse at t/2. |
| `nstep` | integer | required | The number of time points, 0, `deltat`, ..., (`nstep` − 1) × `deltat`. |
| `rdip_ladder` | (internal) | none | Written by the [automatic convergence](InputFile_AutoConvergence.md) into its runs, not for a hand-written input. |

`rbath`, `rdip` and `deltat` also take `{"value", "unit"}` ([Units](InputFile.md#Units)).

```json
"order": 2, "rbath": 50, "rdip": 10, "deltat": 0.001, "nstep": 2000
```

<a id="gen_config"></a>
## 2. configuration

<a id="tag_configuration"></a>
### configuration

| Type | Default |
|---|---|
| object | none (without it, or without `bathfile`, the bath is empty and L(t) = 1) |

The bath files, their states and the random seed in one block. The defect model itself is
[Defect.species](InputFile_Defect.md#tag_df_species).

```json
"configuration": {
    "path": "data",                           // optional: the folder of the file names below
    "k_range": [1, 10],                       // optional: one configuration per k
    "gyrofile": "gyro_13C",
    "bathfile": ["config_{k}"],
    "bath_state_file": ["cfg{k}_state{s}"],
    "s_range": [1, 5],                        // optional: sampled states; absent: the ensemble
    "seed": 1
}
```

<a id="tag_configuration_k_range"></a>
### configuration.k_range

| Type | Default |
|---|---|
| `[first, last]` | one configuration |

Runs the bath configurations first .. last one after another in the same run. Write [{k}](#tag_configuration_k) in
the file names, and it is replaced by each k. Each configuration writes its own outputs (`out_k3_wiDiv`, ...), and
configuration k uses the [seed](#tag_seed) + (k − first). It does not run with a `cluster.adaptive` block.

```json
"k_range": [1, 10], "bathfile": ["config_{k}"]      // config_1 ... config_10
```

<a id="tag_configuration_k"></a>
### {k} and {s} in file names

With `k_range`, `{k}` in the entries of `bathfile`, `bath_state_file`, `defect_axis_file`, `defect_state_file` and
`pcce.cache.file` (and in `gyrofile`) is replaced by k. An entry without `{k}` is the same for every k, so a common bath file is merged
into every configuration. `{s}` in the state file names is the state s of `s_range`, and a name without `{s}` has s
appended.
`{{k}}` and `{{s}}` give the literal text.

```json
"bathfile": ["config_{k}"], "bath_state_file": ["cfg{k}_state{s}"]
// with k_range [1, 2] and s_range [1, 3]: reads cfg1_state1 ... cfg2_state3
```

<a id="tag_gyrofile"></a>
### configuration.gyrofile

| Type | Default |
|---|---|
| string (a file, in [configuration.path](#tag_configuration_path)) | none (the built-in isotope table) |

The file of the bath spin species, with one line per species that holds its name, spin and gyromagnetic ratio
(rad/ms/G). A species it
does not name is taken from the built-in isotope table, whose names are mass number + symbol (`13C`, `14N`, `29Si`,
... in [IsotopeTable](Tools.md#tag_tool_isotopetable)). Species that are not nuclear isotopes (`P1`, an electron spin)
must be in the gyrofile.

```
13C 0.5 6.728284
P1  0.5 -17608.597050
```

```json
"configuration": {"gyrofile": "gyro_13C", ...}
```

<a id="tag_bathfile"></a>
### configuration.bathfile

| Type | Default |
|---|---|
| list of strings (files), also for one file | required |

The files of the bath's spins. Several files are merged into one bath. Each file starts with a line of four fields
(not read), then one line per spin with x, y, z in Å and the species name.

```
2.0 0.0 0.0 0.0
1.0 2.0 3.0 13C
4.0 5.0 6.0 P1
```

```json
"bathfile": ["bath_13C", "bath_P1"]
```

<a id="tag_bath_state_file"></a>
### configuration.bath_state_file

| Type | Default |
|---|---|
| list of strings, one per bath file | none (the states are drawn from [seed](#tag_seed)) |

The files of the sampled bath states, read with [s_range](#tag_s_range) only. State s of bath file j is read from name
j with `{s}` replaced by s (or s appended). An entry `""` draws its bath file's states from the seed.

<a id="tag_configuration_files"></a>
**File format** (also of the two defect files below). The first line holds one entry, by convention the bath file's
spin count, and is skipped. Do not copy the bath file's four-field first line. Then comes one line per spin line of the
bath file:

| Key | Each line |
|---|---|
| `bath_state_file` | The m of the spin (S, S−1, ..., −S). |
| `defect_axis_file` | The axis state, 1 .. n_axis_states, of a Defect spin, and 0 for a spin of another species. |
| `defect_state_file` | The m of each on-site spin of a Defect spin, space-separated. |

[StateFiles](Tools.md#tag_tool_statefiles) writes all three.

```json
"bath_state_file": ["state{s}"]                  // state1, state2, ...
```

<a id="tag_s_range"></a>
### configuration.s_range

| Type | Default |
|---|---|
| `[first, last]`, 1 ≤ first ≤ last | none (the ensemble) |

Without it the run computes the unpolarised bath (the ensemble), one output. With it, each state s from first to last
is a pure bath state, from [bath_state_file](#tag_bath_state_file) or drawn from [seed](#tag_seed), with its own
`_state<s>` output. Average them afterwards. Use states for an electron-spin bath (P1). State s is the same in
`[1, 20]` as in `[11, 20]`, so a run can be split.

```json
"configuration": {"s_range": [1, 10], ...}
```

<a id="tag_seed"></a>
### configuration.seed

| Type | Default |
|---|---|
| integer, or a list (one per configuration) | generated and printed in the log when the run draws random numbers |

The seed of the random numbers. The same seed draws the same states. With [k_range](#tag_configuration_k_range) a
number is the first configuration's seed, and a list gives each configuration its own. Before 2.0.0 a number was every
configuration's seed ([Upgrading old inputs](InputFile_OldNames.md)).

```json
"configuration": {"seed": 1241798, ...}
```

<a id="tag_configuration_path"></a>
<a id="tag_bath_shift"></a>
<a id="tag_defect_axis_file"></a>
<a id="tag_defect_state_file"></a>
### The other keys of the block

| Key | Type | Default | Meaning |
|---|---|---|---|
| `path` | string (a folder) | `"."` | Put in front of `gyrofile` and the per-bath-file names, not of the files of other blocks. |
| `bath_shift` | list of `[x, y, z]` (Å), one per bath file | `[0, 0, 0]` each | Shifts each bath file's sites, for example around several qubits. |
| `defect_axis_file` | list of strings, one per bath file | drawn from `seed` | The axis state of each Defect spin, in the [file format](#tag_configuration_files) above. |
| `defect_state_file` | list of strings, one per bath file | drawn from `seed` | The states of the defects' on-site nuclear spins, read with `s_range` as `bath_state_file` is. |

```json
"path": "../bath/14N_1ppm", "bathfile": ["bath_DiaP1_1ppm_1"]     // reads ../bath/14N_1ppm/bath_DiaP1_1ppm_1
```

<a id="tag_configuration_missing"></a>
A file that the input names and the run reads has to open. Otherwise the run stops before anything is calculated and
names the key and the file. A file the run does not read (a state file in the ensemble) is not opened.

<a id="tag_frame"></a>
## 3. frame

<a id="tag_frame_rotation"></a>
### frame.rotation

| Type | Default |
|---|---|
| object (the keys below) | no rotation without the block |

CCEX calculates with the qubit's quantization axis as z. When the bath file is written in another frame (the lab or
crystal frame), for example the cubic axes of diamond where an NV points along [111], give that direction as
`qubit_axis`. CCEX then turns every position about the first qubit, so that `qubit_axis` becomes z. `bfield` is given in the
turned frame, so `[0, 0, B]` is a field B along the NV axis.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `qubit_axis` | three numbers | `[0, 0, 1]` | The qubit's axis, in the bath file's frame. |
| `enabled` | boolean | `true` | `false` turns a written block off. |

A tensor written in the input or in a tensor file has a direction too, so with the rotation on it has to say which
frame it is written in, with a key `frame`:

- `"file"`: in the bath file's frame (the lab or crystal frame). CCEX turns it with the positions.
- `"qubit"`: already in the turned frame (z = `qubit_axis`). CCEX uses it as written.

This applies to [hyperfine.frame](#tag_hyperfine_frame), [quadrupole.frame](#tag_quadrupole_frame),
[Qubit.intmap[].frame](InputFile_Qubit.md#tag_intmap_frame) and [Defect.species[].frame](InputFile_Defect.md#tag_df_frame).
Point-dipole tensors are computed from the turned positions and need no `frame`. The details are in
[The coordinate frame rotation](../methods/coordinate-frame-rotation.md).

`qubit_axis` defines the rotation once for the whole input. Each `frame` only says how one tensor is written. For
example, an NV along [111] in a bath file written in the cubic axes of diamond, with its zero-field splitting
D = 2880 MHz as a `Qubit.intmap` entry, can give the splitting in either frame:

```json
"frame": {"rotation": {"qubit_axis": [1, 1, 1]}},
"Qubit": {"qubit": [{"file": "defect"}],
          "intmap": [{"between": ["q0", "q0"], "unit": "MHz", "frame": "qubit",
                      "tensor": [[-960, 0, 0], [0, -960, 0], [0, 0, 1920]]}]}
```

Here the tensor is the usual diagonal form, with z along the NV axis, so it is already in the turned frame. The same
splitting written in the cubic axes, with its axis along [111], is

```json
"intmap": [{"between": ["q0", "q0"], "unit": "MHz", "frame": "file",
            "tensor": [[0, 960, 960], [960, 0, 960], [960, 960, 0]]}]
```

and CCEX turns it into the first form. A DFT hyperfine file written in the cubic axes is given the same way:

```json
"hamiltonian": {"common": {"hyperfine": {"dipole": "dft", "file": "Afile", "frame": "file"}}}
```

<a id="tag_hamiltonian"></a>
## 4. hamiltonian

The block of the terms of the Hamiltonian, in one section per method.

| Type | Default |
|---|---|
| object of sections | every key at its default without the block |

The terms of the Hamiltonian that are chosen, not given by the files. A setting that is on in the section of another
method is refused (a `hamiltonian.cCCE` section in a gCCE or megCCE run, for example).

| Section | Read by | Holds |
|---|---|---|
| `common` | every method | The [mean fields](#gen_mf), the [hyperfine](#gen_hf) and the [quadrupole](#gen_qd) tensors. |
| `cCCE` | cCCE and meCCE | [knight](#tag_knight), [hfmedi](#tag_hfmedi), [delta_A](#tag_delta_A) |
| `gCCE` | gCCE | nothing yet |
| `meCCE`, `megCCE` | that method | [jump_operators](#tag_jump_operators) |
| `noise` | method noise | [spectrum](#tag_noise_spectrum) and the other noise keys |
| `pseudospin` | method pseudospin | [channels](#tag_pseudospin_channels) |

```json
"hamiltonian": {
    "common": {"bath_bath_mean_field": "z", "hyperfine": {"dipole": "dft", "contact": true, "file": "Afile"}},
    "cCCE": {"knight": false}
}
```

<a id="gen_mf"></a>
### 4.1 common: mean fields

With sampled bath states ([s_range](#tag_s_range)), every bath spin outside a cluster sits in a fixed sublevel, and
its static field still acts on the cluster. A cluster's calculation cannot include those spins, so CCEX adds their
field as a mean field. Two independent keys of `hamiltonian.common` say where it is added, and they can be set together.
In the ensemble (no `s_range`) there is no such field, and neither key is used.

| Key | The field of the spins outside the cluster on | Read by | Values | Default |
|---|---|---|---|---|
| [qubit_bath_mean_field](#tag_qubit_bath_mean_field) | the qubit (the Overhauser field), through the hyperfine coupling | gCCE, megCCE (the qubit is in their Hamiltonian) | `true`, `false` | `false` |
| [bath_bath_mean_field](#tag_bath_bath_mean_field) | the bath spins of the cluster, through the bath-bath coupling | every method | `"z"`, `"xyz"`, `"off"` | `"z"` |

```json
"hamiltonian": {"common": {"qubit_bath_mean_field": true, "bath_bath_mean_field": "xyz"}}
```

<a id="tag_qubit_bath_mean_field"></a>
#### qubit_bath_mean_field

`true` adds to the qubit the hyperfine field of the bath spins outside the cluster. With `false` the qubit feels the
spins of the cluster only. cCCE and meCCE project the qubit out and do not use the key.

<a id="tag_bath_bath_mean_field"></a>
#### bath_bath_mean_field

How much of the field of the spins outside the cluster each spin of the cluster feels.

| Value | Field |
|---|---|
| `"z"` | The z part only. |
| `"xyz"` | The whole vector. It can help where the z-only field makes the expansion diverge (clustered P1 baths). |
| `"off"` | None. The cluster is alone in its sampled state. |

[Basic 06](../../example/tutorial/basic/06_defect_bath/README.md) runs the three on a P1 bath.

<a id="gen_hf"></a>
### 4.2 common: hyperfine

The hyperfine tensors between the qubit and the bath spins.

<a id="tag_hyperfine"></a>
#### hyperfine

| Type | Default |
|---|---|
| object, in `hamiltonian.common` | `{"dipole": "point", "contact": false}` |

Where the hyperfine tensors between the qubit and the bath spins come from. [file](#tag_hyperfine_file) is required
when `dipole` is `"dft"` or `contact` is `true`.

| dipole | contact | Tensor |
|---|---|---|
| `"point"` | `false` | Point dipole from the positions (the default). |
| `"point"` | `true` | Fermi contact from the file + point dipole. |
| `"dft"` | `false` | The DFT dipolar tensor of the file. |
| `"dft"` | `true` | Fermi contact + DFT dipolar tensor. |

```json
"hamiltonian": {"common": {"hyperfine": {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}}}
```

<a id="tag_hyperfine_file"></a>
#### hyperfine.file

| Type | Default |
|---|---|
| string (a file) | required when `dipole` is `"dft"` or `contact` is `true` |

The hyperfine tensor file written by [SpinTensorReader](Tools.md#tag_tool_spintensorreader) (format v3). It holds one
row per atom with its position relative to the defect, the isotropic part and the nine components (MHz, per unit
nuclear g-factor). A bath spin outside the range the file covers gets a point-dipole tensor. One qubit only.

```json
"hyperfine": {"dipole": "dft", "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}
```

<a id="tag_hyperfine_missing_row"></a>
<a id="tag_hyperfine_frame"></a>
The other keys of `hyperfine`:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `missing_row` | `"stop"` or `"point_dipole"` | `"stop"` | What happens to a spin inside the file's range with no row. `"stop"` stops and prints its position, and `"point_dipole"` gives it a point-dipole tensor. |
| `frame` | `"file"` or `"qubit"` | required with `frame.rotation` and `dipole` `"dft"` | The file's components are in the bath file's frame (turned) or already computational (used as written). |

```json
"hyperfine": {"dipole": "dft", "file": "Afile", "missing_row": "point_dipole", "frame": "file"}
```

<a id="gen_qd"></a>
### 4.3 common: quadrupole

The quadrupole (electric field gradient) tensors of the bath spins.

<a id="tag_quadrupole"></a>
#### quadrupole

| Type | Default |
|---|---|
| `"off"`, or an object, in `hamiltonian.common` | `"off"` |

Where the quadrupole (electric field gradient) tensors of the bath spins come from. A bath spin inside the range of
[file](#tag_quadrupole_file) takes its row. Outside it:

- `"plain"` takes the species default of `file`.
- `"strained"` folds the position by [cell](#tag_quadrupole_cell) into the range of
  [file_without_defect](#tag_quadrupole_file_without_defect) and takes that file's row.
- `"bubble"` first looks the position up in `file_without_defect` as it is, then folded. It also uses
  `file_without_defect` for a spin inside `file`'s range that has no row there.

A position that still falls outside every range gets the species default of `file`. One qubit only.

| structure | Files |
|---|---|
| `"plain"` | `file`, computed with the defect. |
| `"strained"`, `"bubble"` | `file`, `file_without_defect` and `cell`. |

```json
"hamiltonian": {"common": {"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}}}
```

<a id="tag_quadrupole_file"></a>
#### quadrupole.file

| Type | Default |
|---|---|
| string (a file) | required with every `structure` |

The quadrupole tensor file computed with the defect, written by [SpinTensorReader](Tools.md#tag_tool_spintensorreader).
It has the layout of [hyperfine.file](#tag_hyperfine_file) with the nine EFG components (Hartree/Bohr²) and no
isotropic column, one `eQ___` line per species in place of the g-factor, and no code tag or provenance line.

```json
"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}
```

<a id="tag_quadrupole_file_without_defect"></a>
<a id="tag_quadrupole_cell"></a>
<a id="tag_quadrupole_frame"></a>
The other keys of `quadrupole`:

| Key | Type | Default | Meaning |
|---|---|---|---|
| `file_without_defect` | string (a file) | required with `"strained"`, `"bubble"` | The same structure without the defect. It gives the tensors of the spins outside `file`'s range. |
| `cell` | three numbers, Å | required with `"strained"`, `"bubble"` | The lattice parameters. Positions are folded into the range of `file_without_defect` along every axis whose value is not 0, and a range wider than the value is refused. |
| `frame` | `"file"` or `"qubit"` | required with `frame.rotation` and a quadrupole structure | As `hyperfine.frame`, for the quadrupole files. |

```json
"quadrupole": {"structure": "strained", "file": "Q_with_defect", "file_without_defect": "Q_without_defect", "cell": [2.5, 2.5, 0]}
```

<a id="gen_ham_method"></a>
### 4.4 cCCE

The section `hamiltonian.cCCE`, read by cCCE and meCCE.

<a id="tag_hfmedi"></a>
<a id="tag_knight"></a>
#### hfmedi and knight

| Key | Type | Default | Meaning |
|---|---|---|---|
| `hfmedi` | boolean, in `hamiltonian.cCCE` | `false` | The hyperfine-mediated coupling between bath spins through the qubit, in full for each branch (N. Zhao et al., Nature Nanotechnology 6, 242 (2011)). |
| `knight` | boolean, in `hamiltonian.cCCE` | `false` | The Knight-field shift of the bath spins under a field across the qubit axis. It needs a qubit with a zero-field splitting. |

The terms are in [cCCE and gCCE](../methods/cce.md), and
[Advanced 09](../../example/tutorial/advanced/09_hamiltonian_terms/README.md) runs them.

```json
"hamiltonian": {"cCCE": {"hfmedi": true}}
```

<a id="tag_delta_A"></a>
#### delta_A

| Type | Default |
|---|---|
| boolean, in `hamiltonian.cCCE` | `false` |

Adds the second-order hyperfine correction of each bath spin, δA<sub>n</sub> = −m<sub>s</sub>(A<sub>x</sub> × A<sub>y</sub>)/(2Δ)
(J. R. Maze et al., PRB 78, 094303 (2008)). It matters for nuclei close to the defect at low field, needs a single
spin-1 qubit with a zero-field splitting, and does not go with [hfmedi](#tag_hfmedi), which already contains it. Its
conditions are in [cCCE and gCCE](../methods/cce.md#cce_delta_A).

```json
"hamiltonian": {"cCCE": {"delta_A": true}}
```

<a id="gen_ham_me"></a>
### 4.5 meCCE and megCCE

The section `hamiltonian.meCCE` or `hamiltonian.megCCE`, the one of the run's method.

<a id="tag_jump_operators"></a>
#### jump_operators

| Type | Default | Unit |
|---|---|---|
| list of `{"bath_name", "operator", "rate"}`, in `hamiltonian.meCCE` or `hamiltonian.megCCE` | none (no dissipation) | `rate` in MHz (≥ 0), or `{"value", "unit"}` |

The Lindblad jump operators of the master-equation methods. Each acts as C = √γ S<sub>op</sub> on every spin named
`bath_name`, with γ = `rate` × 2π·10³ rad/ms. `operator` is `-`, `+`, `z`, `x` or `y`. `bath_name` is a bath spin name
(case not significant), or `"qubit"` or the qubit's name for a jump of the qubit (megCCE only). See [master-equation CCE](../methods/master-equation.md).

```json
"method": "megCCE",
"hamiltonian": {"megCCE": {"jump_operators": [
    {"bath_name": "e",     "operator": "-", "rate": 1.5915e-4},
    {"bath_name": "qubit", "operator": "z", "rate": 1.5915e-4}]}}
```

<a id="gen_ham_noise"></a>
### 4.6 noise

The section `hamiltonian.noise`, read by method noise.

<a id="tag_noise_spectrum"></a>
#### spectrum

| Type | Default |
|---|---|
| `"bath"`, or a model object, in `hamiltonian.noise` | `"bath"` |

Where the noise of the qubit frequency comes from, for method noise. The coherence is exp(−χ) in the Gaussian
approximation ([Semiclassical noise](../methods/noise.md)), which needs many weakly coupled spins.

| Value | Noise |
|---|---|
| `"bath"` | The bath of the run, from the modes of its cCCE clusters. |
| `{"model": "OU", "sigma": σ, "tau_c": τ_c}` | Ornstein–Uhlenbeck noise (σ in kHz, τ_c in ms). |
| `{"model": "power", "amplitude": A, "beta": β, "f_min": f_min, "f_max": f_max}` | S(f) = A (f/1 kHz)<sup>−β</sup> on [f_min, f_max] (kHz). |

A model run reads no bath, so `order`, `bfield`, `rbath`, `rdip`, `Qubit`, `configuration` and `Defect` may be left out.

```json
{"method": "noise", "hamiltonian": {"noise": {"spectrum": {"model": "OU", "sigma": 10, "tau_c": {"value": 50, "unit": "us"}}}},
 "deltat": 0.002, "nstep": 200, "pulse": {"count": 1}}
```

<a id="tag_noise_pseudo_secular"></a>
#### pseudo_secular

| Type | Default |
|---|---|
| boolean, in `hamiltonian.noise` (with `spectrum` `"bath"`) | `false` |

With `false` only the secular hyperfine term A<sub>zz</sub>I<sub>z</sub> counts as frequency noise. With `true`,
A<sub>zx</sub>I<sub>x</sub> + A<sub>zy</sub>I<sub>y</sub> also counts, the coupling cCCE uses. Which is the better
approximation has not been settled.

```json
"hamiltonian": {"noise": {"pseudo_secular": true}}
```

<a id="tag_noise_compute"></a>
#### compute

| Type | Default |
|---|---|
| string, in `hamiltonian.noise` (with `spectrum` `"bath"`) | `"coherence"` |

| Value | Computes | Writes |
|---|---|---|
| `"coherence"` | The coherence exp(−χ) from each cluster's eigenmodes. | `_wiDiv` (and `_noDiv`) |
| `"spectrum"` | The correlation C(t), its spectrum S(f) by FFT ([fft](#tag_noise_fft)), and the coherence from S(f). | `_noise_correlation`, `_noise_spectrum`, `_wiDiv` |
| `"correlation"` | C(t) at the run's times, as PyCCE's `CCENoise` does. | `_noise_correlation` |

C(t) is in kHz² (t in ms), S(f) two-sided in kHz²/kHz (f in kHz). Without `s_range`, `"correlation"` and `"spectrum"` give the same
C(t), and the coherence of `"spectrum"` approaches that of `"coherence"` as the `fft` grid is refined.

```json
"hamiltonian": {"noise": {"compute": "correlation"}}
```

<a id="tag_noise_fft"></a>
#### fft

| Type | Default |
|---|---|
| object, in `hamiltonian.noise` (with `compute` `"spectrum"`) | every key has a default |

The time grid of C(t) and the FFT.

| Key | Meaning | Default |
|---|---|---|
| `time_max` | The end of the grid (ms). | the run's last time |
| `nstep` | The number of intervals. | chosen from the bath's fastest modes |
| `window` | `"none"`, `"hann"`, `"tukey"` or `"exp"` on C(t). | `"none"` |
| `pad_factor` | The FFT length in units of the grid. | 1 |
| `subtract_mean` | Subtract the mean of the tail of C(t) first. | `false` |
| `clamp_negative` | Set S(f) < 0 to 0. | `false` |

For a finer S(f), give a longer `time_max`.

```json
"hamiltonian": {"noise": {"compute": "spectrum", "fft": {"time_max": 3.0, "nstep": 15000, "window": "hann"}}}
```

<a id="gen_ham_pseudospin"></a>
### 4.7 pseudospin

The section `hamiltonian.pseudospin`, read by method pseudospin.

<a id="tag_pseudospin_model"></a>
<a id="tag_pseudospin_channels"></a>
#### channels

| Type | Default |
|---|---|
| string, in `hamiltonian.pseudospin` | `"grouped"` |

How [method](#tag_method) pseudospin evolves the channels of a pair. A channel is an element of the pair's
Hamiltonian that joins two product states differing in both spins. It is a flip-flop, a flip-flip, or with a mixing
quadrupole a double-quantum element. Each channel has, in each qubit state (α, β):

- Ω<sub>α</sub>, Ω<sub>β</sub>: the coupling of its two states, 2\|⟨m<sub>1</sub>, m<sub>2</sub>\|H\|m<sub>1</sub>′, m<sub>2</sub>′⟩\| (kHz).
- Δ<sub>α</sub>, Δ<sub>β</sub>: the signed energy difference E(m<sub>1</sub>, m<sub>2</sub>) − E(m<sub>1</sub>′, m<sub>2</sub>′) (kHz).

The run writes every channel with its Ω and Δ to `<output.file>_pseudospin_channels`. The model is explained in
[the pseudospin model](../methods/pseudospin.md).

| Value | Evolves |
|---|---|
| `"grouped"` | The blocks the channels join the product states into, each evolved whole. It is exact within the cut Hamiltonian. |
| `"pairwise"` | Every channel alone, as the 2 × 2 block of its two states (coupling Ω, gap Δ). It misses what couples channels that share a state. |

```json
"method": "pseudospin", "hamiltonian": {"pseudospin": {"channels": "pairwise"}}
```

<a id="gen_reference"></a>
## 5. Reference

<a id="tag_not_used"></a>
### Keys read and not used

A key this run reads but does not use is reported once (`Warning: input key "<key>" is not used by this run: <why>`),
and the run goes on without it. The common cases:

| Key | Not used when |
|---|---|
| `Qubit.psi0` | cCCE, meCCE |
| `hamiltonian.common.qubit_bath_mean_field` | cCCE, meCCE, and the ensemble |
| `hamiltonian.common.bath_bath_mean_field` | the ensemble |
| `configuration.bath_state_file`, `defect_state_file` | the ensemble (no `s_range`) |
| `configuration.defect_axis_file`, `defect_state_file` | no Defect species |
| the `frame` key of a tensor (hyperfine, quadrupole, intmap entry, Defect species) | `frame.rotation` off |
| `execution.gpu_options` | `main.out` (built without CUDA) |
| the bath keys (`order`, `bfield`, `rbath`, `rdip`, `Qubit`, `configuration`, `Defect`) | method noise with a model spectrum |

A key that no part of the run reads is reported as not read, and a key that looks like a typo of a key that is read
stops the run.

<a id="tag_turn_off"></a>
### Turning a part off

A block written without `"enabled"` is on.

| Key | Off | Without the key |
|---|---|---|
| `hamiltonian.common.quadrupole` | `"off"` | off |
| `hamiltonian.common.bath_bath_mean_field` | `"off"` | `"z"` |
| `cluster.pcce.cache` | `{"mode": "off"}` | off |
| `frame.rotation` | `{"enabled": false}` | off |
| `execution.gcce_refinement` | `{"enabled": false}` | off |
| `execution.column_evolution` | `{"enabled": false}` | on |
| `output.cluster_contributions` | `{"enabled": false}` | off |
| `auto_convergence` | `{"enabled": false}` | off |
