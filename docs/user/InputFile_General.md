# General tags
<a id="GeneralTags"></a>

The method, the numbers of the calculation, the files of the bath, the frame, and the `hamiltonian` block. Each key
has a small table (its type, default and unit), what it does, and an example. Old key names are in
[Old names](InputFile_OldNames.md).

1. [Calculation](#gen_calc): `method`, `quantity`, `order`, `bfield`, `rbath`, `rdip`, `deltat`, `nstep`.
2. [configuration](#tag_configuration): the bath files, the states and the seed.
3. [frame](#tag_frame): the rotation of the bath file's frame.
4. [hamiltonian: common](#tag_hamiltonian): the mean fields, the hyperfine and the quadrupole tensors.
5. [hamiltonian: per method](#gen_ham_method): the cCCE terms, pseudospin, the jump operators and the noise.
6. [Reference](#gen_reference): keys read and not used, and how to turn a part off.

<a id="gen_calc"></a>
## 1. Calculation

<a id="tag_method"></a>
### method

| Type | Default |
|---|---|
| string (case is not significant) | `"cCCE"` |

The CCE method. How the bath is cut into clusters is a separate key,
[cluster.mode](InputFile_Cluster.md#tag_cluster_mode).

| Value | Method |
|---|---|
| `cCCE` | conventional CCE: the qubit is projected out into its two branches |
| `gCCE` | generalized CCE: the qubit is in the cluster Hamiltonian |
| `meCCE` | cCCE with Lindblad jump operators on the bath spins ([master-equation CCE](../methods/master-equation.md)) |
| `megCCE` | gCCE with jump operators on the bath spins and on the qubit |
| `noise` | the Gaussian coherence exp(−χ) of the qubit's phase noise ([spectrum](#tag_noise_spectrum)) |
| `pseudospin` | the pairs of cCCE in the pseudospin model of nuclear flip-flops ([channels](#tag_pseudospin_channels)) |

Which `cluster.mode` each method takes:

- `pcce` runs with cCCE, gCCE, meCCE and megCCE;
- `bundle` with cCCE and gCCE;
- `dsj`, `itb` and `dsjitb` with cCCE only.

Limits:

- `meCCE`: CPU, one qubit, ideal π pulses from `pulse.count`; qubit jumps need megCCE.
- `megCCE`: CPU. The qubit sequence is used as given (Manual and any rotation), with `pulse.per_qubit` false and
  `pulse.detuning_factor` 1. DEER bath pulses run with it.
- `noise`: cCCE clusters with `cluster.mode` standard, one qubit, ideal π pulses, `pulse.per_qubit` false,
  `execution.propagator` eigen, `quantity` coherence; no `hamiltonian.cCCE` terms and no
  `execution.gpu_options.ranks_per_gpu`. With `execution.backend` `cuda` the sums over the bath's modes run on the GPU.
- `pseudospin`: one qubit, the ensemble, order 2, `cluster.mode` standard, the field along the qubit axis, ideal π
  pulses, `pulse.per_qubit` false, `execution.propagator` eigen, CPU; no `hamiltonian.cCCE` terms, no Defect species with
  on-site spins and no `ranks_per_gpu`. It writes the pairs' product alone, without the single spins, and a table of every channel
  (`<output.file>_pseudospin_channels`). It is for reading a decay rather than for the best curve
  ([the pseudospin model](../methods/pseudospin.md)).

```json
"method": "gCCE"
"method": "cCCE", "cluster": {"mode": "pcce"}
```

<a id="tag_quantity"></a>
### quantity

| Type | Default |
|---|---|
| string | `"coherence"` |

What the run computes.

| Value | Computes |
|---|---|
| `coherence` | the decay curve L(t) between \|α⟩ and \|β⟩, the usual choice |
| `dm` | the qubit density matrix, for populations or several qubit levels: gCCE with `cluster.mode` standard or pcce ([The density matrix](OutputFiles.md#out_dm)) |

The other methods compute `coherence` only. The noise calculation is [method](#tag_method) `noise`, not a quantity.

```json
"quantity": "dm"
```

<a id="tag_order"></a>
### order

| Type | Default |
|---|---|
| integer | required (method pseudospin: 2; noise with a model spectrum: 1) |

The size of the largest cluster: counted in bath spins (`cluster.mode` standard), in partitions (pcce) or in bundles
(bundle).

```json
"order": 2
```

<a id="tag_bfield"></a>
### bfield

| Type | Default | Unit |
|---|---|---|
| number (Bz) or `[Bx, By, Bz]` | required (`0` for no field) | G; or `{"value", "unit"}` with `G`, `mT`, `T` ([Units](InputFile.md#Units)) |

The magnetic field, in the computational frame. With [frame.rotation](#tag_frame_rotation) on, the computational z is
`qubit_axis`, so `[0, 0, B]` (or B) is a field along `qubit_axis`, and a transverse part is refused.

```json
"bfield": 500
"bfield": [0, 0, 500]
"bfield": {"value": 50, "unit": "mT"}
```

<a id="tag_rbath"></a>
### rbath

| Type | Default | Unit |
|---|---|---|
| number, or `[inner, outer]` | required (inner 0) | Å; or `{"value", "unit"}` with `A`, `nm` |

The bath spins within this distance of the qubit are kept. With `[inner, outer]`, the spins closer than `inner` are
left out too.

```json
"rbath": 50
"rbath": [2, 50]
```

<a id="tag_rdip"></a>
### rdip

| Type | Default | Unit |
|---|---|---|
| number, or `[inner, outer]` | required, except with `cluster.mode` bundle (inner 0) | Å; or `{"value", "unit"}` with `A`, `nm` |

Two bath spins closer than this are connected and can be in one cluster. With `[inner, outer]`, the pairs closer than
`inner` are left out too. With `cluster.mode` bundle the key may be left out; the connections are then chosen by a
coherence check ([cluster.mode](InputFile_Cluster.md#tag_cluster_mode)).

```json
"rdip": 10
"rdip": {"value": 1.0, "unit": "nm"}
```

<a id="tag_deltat"></a>
### deltat

| Type | Default | Unit |
|---|---|---|
| number | required | ms; or `{"value", "unit"}` with `s`, `ms`, `us`, `ns` |

The time step. The outputs hold the times 0, `deltat`, ..., (`nstep` − 1) × `deltat`. Each time is the total time of
the sequence: a Hahn echo of time t has its π pulse at t/2.

```json
"deltat": 0.001
"deltat": {"value": 1, "unit": "us"}
```

<a id="tag_nstep"></a>
### nstep

| Type | Default |
|---|---|
| integer | required |

The number of time points, from t = 0.

```json
"nstep": 2000
```

<a id="tag_rdip_ladder"></a>
### rdip_ladder (internal)

The [automatic convergence](InputFile_AutoConvergence.md) writes this key into the inputs of its runs; a hand-written
input does not need it. It gives extra rdip values below `rdip`, computed from the same order-2 run, so that one run
gives several rdip values (`<output.file>_rdip<r>_...`). It needs order 2, cCCE with `cluster.mode` standard or pcce or
gCCE with standard, `quantity` and `output.savemode` coherence, and no `cluster.selection.max_per_order` or
`output.cluster_contributions`.

<a id="gen_config"></a>
## 2. configuration

<a id="tag_configuration"></a>
### configuration

| Type | Default |
|---|---|
| object | required: [bathfile](#tag_bathfile) has no default (the other keys have their own) |

The bath configuration in one block: the files of the bath, its states and the random seed. Every file the block
names and the run reads has to open ([files that do not open](#tag_configuration_missing)). The keys in their old
top-level places are refused with their place in the block ([Old names](InputFile_OldNames.md)). The defect model
itself is [Defect.species](InputFile_Defect.md#tag_df_species).

```json
"configuration": {
    "path": "data",                           // optional: the folder of the file names below
    "k_range": [1, 10],                       // optional: one configuration per k, run one after another
    "gyrofile": "gyro_13C",
    "bathfile": ["config_{k}"],
    "bath_shift": [[0, 0, 0]],
    "bath_state_file": ["cfg{k}_state{s}"],
    "defect_axis_file": ["cfg{k}_axis"],      // with Defect species in the bath
    "defect_state_file": ["cfg{k}_onsite{s}"],
    "s_range": [1, 5],                        // optional: sampled states; absent: the ensemble
    "seed": 1
}
```

<a id="tag_configuration_path"></a>
### configuration.path

| Type | Default |
|---|---|
| string (a folder) | `"."` (the working directory) |

The folder of the block's file names, put in front of `gyrofile` and of every entry of `bathfile`, `bath_state_file`,
`defect_axis_file` and `defect_state_file`. It is relative to the working directory, or absolute. A file name that is
itself absolute is used as it is. It does not apply to the files of other blocks (`Qubit.qubit[].file`, the tensor
files of `hamiltonian`), which stay relative to the working directory.

```json
"path": "../bath/14N_1ppm", "bathfile": ["bath_DiaP1_1ppm_1"]     // reads ../bath/14N_1ppm/bath_DiaP1_1ppm_1
```

<a id="tag_configuration_k_range"></a>
### configuration.k_range

| Type | Default |
|---|---|
| `[first, last]`, two nonnegative integers | none: one configuration, `{k}` is not replaced |

Several bath configurations in one run. Every k from `first` to `last` (both included) is a configuration of its own,
calculated one after another in the same execution, with [{k}](#tag_configuration_k) in the file names replaced by k.
Each configuration writes its own outputs, `<output.file>_k<k>_...` (e.g. `out_k3_wiDiv`, `out_k3_state1_noDiv`); no
average is written. The log has a "Configuration k = ..." header before each one.

Seeds: configuration k uses seed + (k − first), and a [seed](#tag_seed) list gives each configuration its own. Each
configuration draws its own states, so a run of its files alone with its seed gives it byte for byte. Without `seed`,
the log prints each one: "[seed] configuration k=...: seed S".

Refused, before anything is calculated:

- `k_range` with no `{k}` in any entry of `bathfile`, `bath_state_file`, `defect_axis_file` or `defect_state_file`;
- a file of any k that cannot be read (every bath file and axis file with `{k}`, and the state files of every state
  of [s_range](#tag_s_range));
- a `k_range` that is not two nonnegative integers with first ≤ last;
- an enabled pCCE partition cache ([pcce.cache](InputFile_Cluster.md#tag_pcce_cache) `mode` read or readwrite) whose
  `file` has no `{k}` (one cache file per configuration, e.g. `"partition_{k}.cache"`);
- `k_range` together with a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block (refine one configuration at a
  time).

```json
"k_range": [1, 10], "bathfile": ["config_{k}"]      // runs config_1 ... config_10
```

<a id="tag_configuration_k"></a>
### {k} and {s} in file names

With [k_range](#tag_configuration_k_range), `{k}` is replaced by the configuration's k in the entries of `bathfile`,
`bath_state_file`, `defect_axis_file`, `defect_state_file` and in `pcce.cache.file`.

- An entry without `{k}` is the same for every k. A bath file without `{k}` is merged into every configuration:
  `["nv1_{k}", "common"]` is, for each k, the bath of `nv1_<k>` and `common` together.
- `{{k}}` gives a literal `{k}`.
- Without `k_range` nothing is replaced: a file actually named `config_{k}` is read as it is.

`{s}` in the entries of [bath_state_file](#tag_bath_state_file) and [defect_state_file](#tag_defect_state_file) is the
state s of [s_range](#tag_s_range), by the same rule (`{{s}}` gives a literal `{s}`). A name without `{s}` has the
state appended: `"state_"` reads `state_1`, `state_2`, ....

```json
"bathfile": ["config_{k}"], "bath_state_file": ["cfg{k}_state{s}"]
// with k_range [1, 2] and s_range [1, 3]: reads cfg1_state1 ... cfg2_state3
```

<a id="tag_gyrofile"></a>
### configuration.gyrofile

| Type | Default |
|---|---|
| string (a file, in [configuration.path](#tag_configuration_path)) | none: the built-in isotope table |

The file of the bath spin species, one line per species: name, spin and gyromagnetic ratio (rad/ms/G).

```
13C 0.5 6.728284
P1  0.5 -17608.597050
```

A bath species that the gyrofile does not name, or every species without a gyrofile, is taken from the built-in
isotope table. The table holds every nuclear isotope (mass number > 0) with a positive spin of the EasySpin isotope
database (magnetic moments from N. Stone, IAEA INDC(NDS)-0658, 2014), with γ = g<sub>n</sub> μ<sub>N</sub>/ħ (CODATA
2018). Its names are mass number + symbol (`1H`, `11B`, `13C`, `14N`, `29Si`, ...); it is made by
[IsotopeTable](Tools.md#tag_tool_isotopetable).

- A gyrofile entry wins over the table. A value that differs from the table by more than 0.1 % is noted in the log.
- Species that are not nuclear isotopes (`P1`, `Vp`, an electron spin) must be in the gyrofile.
- A name found in neither stops the run. The log says where each species came from.

```json
"configuration": {"gyrofile": "gyro_13C", ...}
```

<a id="tag_bathfile"></a>
### configuration.bathfile

| Type | Default |
|---|---|
| list of strings (files, in [configuration.path](#tag_configuration_path); `{k}` with [k_range](#tag_configuration_k_range)); a bare string is refused | required |

The files of the bath's spins. Several files are merged into one bath, in their order; the log names the merge
("Bath = fileA (n spins) + fileB (m spins)"). Each file is a count line, then one line per spin:

```
2.0 0.0 0.0 0.0     four fields, not read (by convention the number of spins, then 0 0 0)
1.0 2.0 3.0 13C     x, y, z in Å, and the species name (as in the gyrofile)
4.0 5.0 6.0 P1
```

Two spins of different files close to each other (after [bath_shift](#tag_bath_shift)):

- within 0.5 Å: the log warns, naming the files and positions, and the calculation goes on. These are two spins at
  one site, e.g. species files made independently on one lattice, or two configurations merged by mistake (give each
  its own run, [k_range](#tag_configuration_k_range));
- at the same position (closer than 10<sup>-6</sup> Å): refused before anything is calculated, since their coupling is
  infinite.

```json
"bathfile": ["bath_13C"]
"bathfile": ["bath_13C", "bath_P1"]
```

<a id="tag_bath_shift"></a>
### configuration.bath_shift

| Type | Default | Unit |
|---|---|---|
| list of `[x, y, z]`, one per bath file (a flat `[1, 2, 3]` is refused) | `[0, 0, 0]` for every bath file | Å; or `{"value", "unit"}` with `A`, `nm` |

Shifts each bath file's sites by its `[x, y, z]`, e.g. when several bath files describe the surroundings of different
qubits. The list has exactly as many rows as `bathfile` has entries.

```json
"bath_shift": [[1.0, 2.0, 3.0]]
"bath_shift": [[0, 0, 0], [4.0, 5.0, 6.0]]
```

<a id="tag_bath_state_file"></a>
### configuration.bath_state_file

| Type | Default |
|---|---|
| list of strings, one per bath file (`{k}`, `{s}`); another length or a bare string is refused | none: the states are drawn from [seed](#tag_seed) |

The files of the bath states, one name per bath file, in the order of `bathfile`. They are read with
[s_range](#tag_s_range) only; in the ensemble the key is not used (the run says so).

- The state s of the spins of bath file j is read from name j with `{s}` replaced by s, or s appended
  ([{s}](#tag_configuration_k)).
- An entry `""` is a bath file without a state file. The whole bath is drawn from the seed first, and the files given
  then set the states of their own bath files. `""` for every bath file is the same as leaving the key out. The same
  rule holds for [defect_axis_file](#tag_defect_axis_file) and [defect_state_file](#tag_defect_state_file).
- A file that does not open stops the run, naming the key and the file.

<a id="tag_configuration_files"></a>
**File formats of the three keys.** The first line holds one entry, by convention the bath file's spin count, and is
skipped. Do not copy the whole four-column first line of the bath file: its extra numbers would be read as the first
spins' values. Then comes one line per spin line of the bath file, in its order:

| Key | Each line |
|---|---|
| `bath_state_file` | the m of the spin (one of S, S−1, ..., −S) |
| `defect_axis_file` | the axis state, 1 .. n_axis_states, of a Defect spin; 0 for a spin of another species (another value is refused) |
| `defect_state_file` | the m of each on-site spin of a Defect spin, space-separated, in the order of the species' types (any number for a line of another species) |

```
 5.0     the first line (not read)
-0.5     m of the spin of the bath file's second line
 0.5
-0.5
-0.5
```

[StateFiles](Tools.md#tag_tool_statefiles) writes the files of all three keys that an input names, by the same rule as
the automatic convergence.

```json
"bath_state_file": ["state{s}"]                  // state1, state2, ... ; ["state"] reads the same files
"bath_state_file": ["stateA_{s}", "stateB_{s}"]  // two bath files
```

<a id="tag_defect_axis_file"></a>
### configuration.defect_axis_file

| Type | Default |
|---|---|
| list of strings, one per bath file (`{k}`); another length or a bare string is refused | none: drawn from [seed](#tag_seed) |

The axis state (the principal axis) of each bath spin that is a [Defect species](InputFile_Defect.md#tag_df_species).
File j holds the lines of bath file j ([file formats](#tag_configuration_files)), and the files are read one after
another in the merged bath's order. `""` draws the axes of its bath file from the seed. Without Defect species the key
is not used.

```json
"bathfile": ["bath_1"], "defect_axis_file": ["axis_1"]
```

<a id="tag_defect_state_file"></a>
### configuration.defect_state_file

| Type | Default |
|---|---|
| list of strings, one per bath file (`{k}`, `{s}`); another length or a bare string is refused | none: drawn from [seed](#tag_seed) |

The states of the defects' on-site nuclear spins, read as [bath_state_file](#tag_bath_state_file) is: one name per
bath file, state s from name j with `{s}` replaced (or s appended). `""` draws its bath file's from the seed. In the
ensemble, or without Defect species, the key is not used.

```
 5.0     the first line (not read)
-1.0     m of the 14N of the P1 on the bath file's second line
 0.0
 1.0
-1.0
```

```json
"bathfile": ["bath_1"], "s_range": [1, 4], "defect_state_file": ["onsite{s}"]
```

<a id="tag_s_range"></a>
### configuration.s_range

| Type | Default |
|---|---|
| `[first, last]`, 1 ≤ first ≤ last (a list also for one state: `[3, 3]`) | none: the ensemble |

The sampled bath states.

| Value | Runs |
|---|---|
| absent | the ensemble: the average over the unpolarised bath, one output |
| `[first, last]` | the states first .. last, each a pure bath state (every spin in one sublevel, from [bath_state_file](#tag_bath_state_file) or drawn from [seed](#tag_seed)), each with its own `_state<s>` output; the static field of the other bath spins is added ([bath_bath_mean_field](#tag_bath_bath_mean_field)) |

Use states for an electron-spin bath (P1), where one bath state differs much from another, or for the spread between
single centres; their average is left to the analysis. States drawn from the seed are drawn for 1 .. last and
calculated for first .. last, so `[11, 20]` and `[1, 10]` together are `[1, 20]` split over two runs. With state files,
only the files of first .. last are read.

```json
"configuration": {"s_range": [1, 10], ...}
```

<a id="tag_seed"></a>
### configuration.seed

| Type | Default |
|---|---|
| integer, or a list of integers (one per configuration) | generated on rank 0 and printed (`[seed] ...`) when the run draws random numbers ([s_range](#tag_s_range), Defect species in the input, a pCCE partition); otherwise not used |

The seed of the random numbers; the same seed draws the same states again. A number is the seed of the first
configuration of [k_range](#tag_configuration_k_range), and configuration k uses seed + (k − first); without
`k_range` it is the run's seed. A list gives one seed per configuration: as many as `k_range` has, or one without it.
Before 2.0.0 a number was every configuration's seed; the converter writes such an input as a list of the same number
([Old names](InputFile_OldNames.md)).

```json
"configuration": {"seed": 1241798, ...}
"configuration": {"k_range": [1, 3], "seed": [7, 7, 9], ...}
```

<a id="tag_configuration_missing"></a>
### Files that do not open

A file that the input names and the run reads has to open; otherwise the run stops before anything is calculated,
naming the key and the file. The keys:

- `configuration.bathfile`, `gyrofile`;
- `bath_state_file` and `defect_state_file` (every state of `s_range`);
- `defect_axis_file` (with Defect species);
- `Qubit.qubit[].file`;
- `hamiltonian.common.hyperfine.file`, `hamiltonian.common.quadrupole.file` and `file_without_defect`;
- `cluster.pcce.cache.file` in mode `read`.

`""` in a list of per-bath-file names is no file: that bath file's states or axes are drawn from the seed. A file the
run does not read (a state file in the ensemble) is not opened; the key is reported as not used.

<a id="tag_frame"></a>
## 3. frame

The `frame` block has one key, [rotation](#tag_frame_rotation).

<a id="tag_frame_rotation"></a>
### frame.rotation

| Type | Default |
|---|---|
| object `{"enabled", "qubit_axis"}` | no block: no rotation |

| Key | Type | Values | Default |
|---|---|---|---|
| `enabled` | boolean | `true`, `false` | `true` in a given block |
| `qubit_axis` | three numbers | the qubit's quantization axis in the bath file's frame, any nonzero vector except along [0, 0, −1] | [0, 0, 1] |

Turns the frame of the bath file (the lab or crystal frame) so that `qubit_axis` becomes the computational z axis. Use
it when the bath file is written in a crystal frame (e.g. the cubic axes of diamond) and the qubit points elsewhere
(e.g. an NV along [111]). A given block is on; `"enabled": false` turns it off. The method is explained in
[The coordinate frame rotation](../methods/coordinate-frame-rotation.md).

- Every position (bath spins and qubits) is turned about the first qubit, and every tensor written in the file frame
  with it.
- [bfield](#tag_bfield) is read in the computational frame: `[0, 0, B]` is a field along `qubit_axis`.
- The qubit states (\|α⟩, \|β⟩, psi0) and the pulse axes stay in the computational frame. Its x is
  `qubit_axis` × [0, 0, 1], which sets the preparation phase; with gCCE and a qubit ZFS with E ≠ 0 the coherence
  depends on that choice.

Each set of directional numbers says which frame it is written in, with a key `frame` beside it: `"file"` (the bath
file's frame, turned with the positions) or `"qubit"` (the computational frame, used as written). With the rotation
on, `frame` is required at [hyperfine.frame](#tag_hyperfine_frame) (with `"dipole": "dft"`),
[quadrupole.frame](#tag_quadrupole_frame), [Qubit.intmap[].frame](InputFile_Qubit.md#tag_intmap_frame) (a qubit's
ZFS is a self entry there) and [Defect.species[].frame](InputFile_Defect.md#tag_df_frame) (a species with
`nuclear_positions`, `hyperfine`, `efg` or `zfs`). With the rotation off, a written `frame` is not used (reported).

The old keys `bath_axis`, `reference_qubit`, `qubit_position_frame`, `hf_tensor_frame` and `qd_tensor_frame` are
refused ([Old names](InputFile_OldNames.md)).

```json
"frame": {"rotation": {"qubit_axis": [1, 1, 1]}},
"bfield": 500,
"hamiltonian": {"common": {"hyperfine": {"dipole": "dft", "file": "Afile", "frame": "file"}}}
```

<a id="tag_hamiltonian"></a>
## 4. hamiltonian: common

### hamiltonian

| Type | Default |
|---|---|
| object of sections | no block: every key at its default |

The terms of the Hamiltonian that are chosen, not given by the files. Every key has a default, so the block may be left
out.

| Section | Read by | Holds |
|---|---|---|
| `common` | every method | the mean fields, `hyperfine`, `quadrupole` (this group) |
| `cCCE` | cCCE, and meCCE (it builds the cCCE's Hamiltonian) | [knight](#tag_knight), [hfmedi](#tag_hfmedi), [delta_A](#tag_delta_A) |
| `gCCE` | gCCE | nothing: an empty section is allowed, a key in it is refused |
| `meCCE`, `megCCE` | that method | [jump_operators](#tag_jump_operators) |
| `noise` | method noise | [spectrum](#tag_noise_spectrum) and the other noise keys |
| `pseudospin` | method pseudospin | [channels](#tag_pseudospin_channels) |

A setting in the section of another method is refused: a `cCCE` setting in a gCCE or megCCE run, a `meCCE` or `megCCE`
section in a run of another method, and a `noise` or `pseudospin` section in a run of another method. Without the block: no mean field on
the qubit, the z mean field between bath spins, point-dipole hyperfine, no quadrupole, no knight, hfmedi or delta_A, no
jump operators.

```json
"hamiltonian": {
    "common": {
        "qubit_bath_mean_field": false,
        "bath_bath_mean_field": "z",
        "hyperfine": {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3", "missing_row": "stop"},
        "quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}
    },
    "cCCE": {"knight": false, "hfmedi": false, "delta_A": false}
}
```

<a id="tag_qubit_bath_mean_field"></a>
### qubit_bath_mean_field

| Type | Default |
|---|---|
| boolean, in `hamiltonian.common` | `false` |

With [s_range](#tag_s_range), adds the Overhauser field of the bath on the qubit: the static hyperfine field of the
sampled bath states. In each cluster this is the field of every bath spin outside it, since the cluster's own spins act
on the qubit through their coupling. It enters the qubit Hamiltonian, which gCCE and megCCE build; with cCCE and meCCE,
and in the ensemble, it is not used (reported).

```json
"hamiltonian": {"common": {"qubit_bath_mean_field": true}}
```

<a id="tag_bath_bath_mean_field"></a>
### bath_bath_mean_field

| Type | Default |
|---|---|
| string, in `hamiltonian.common` | `"z"` |

With [s_range](#tag_s_range), each cluster feels the static field of every bath spin outside it, from that spin's
sampled state m<sub>j</sub>. This key says how; in the ensemble it is not used.

| Value | Field of the spins outside the cluster |
|---|---|
| `"z"` | the z part only, J<sub>zz</sub> m<sub>j</sub> I<sub>z</sub> |
| `"xyz"` | the whole vector J·(0, 0, m<sub>j</sub>), so its x and y parts too; a cluster removes its own members' field and keeps their full coupling instead |
| `"off"` | none: a cluster is the cluster alone in its sampled state |

The x and y parts matter for a bath spin whose own z field nearly cancels (a close neighbour's J<sub>zz</sub> m<sub>j</sub>
against its Zeeman term, as in clustered P1 baths); there the z-only field can make the expansion diverge. `"xyz"`
removed such a divergence but was not more accurate elsewhere, so it is not the default. The transverse field is used
only between main bath spins; a term with a defect sub-spin (e.g. a P1's <sup>14</sup>N) stays z only.
[Basic 06](../../example/tutorial/basic/06_defect_bath/README.md) runs the three on a P1 bath.

```json
"hamiltonian": {"common": {"bath_bath_mean_field": "xyz"}}
```

<a id="tag_hyperfine"></a>
### hyperfine

| Type | Default |
|---|---|
| object `{"dipole", "contact", "file", "missing_row", "frame"}`, in `hamiltonian.common` | `{"dipole": "point", "contact": false}` |

Where the hyperfine tensors between the qubit and the bath spins come from.

| dipole | contact | Tensor |
|---|---|---|
| `"point"` | `false` | point dipole from the positions (the default) |
| `"point"` | `true` | Fermi contact from the file + point dipole |
| `"dft"` | `false` | the DFT dipolar tensor of the file |
| `"dft"` | `true` | Fermi contact + DFT dipolar tensor |

[file](#tag_hyperfine_file) is required when `dipole` is `"dft"` or `contact` is `true`, and
[missing_row](#tag_hyperfine_missing_row) is then read too. With [frame.rotation](#tag_frame_rotation) on,
[frame](#tag_hyperfine_frame) gives the frame of the file's tensor components.

```json
"hamiltonian": {"common": {"hyperfine": {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}}}
```

<a id="tag_hyperfine_file"></a>
### hyperfine.file

| Type | Default |
|---|---|
| string (a file) | required when `dipole` is `"dft"` or `contact` is `true` |

The hyperfine tensor file, written from a DFT run by [SpinTensorReader](Tools.md#tag_tool_spintensorreader). After a
header it holds one line per species with its g-factor and an `etc____` line (read, but not used: a spin outside the
range gets a point-dipole tensor), the range of relative positions it covers (`MinDif[A]`, `MaxDif[A]`), and one row per atom: species, position
relative to the defect (Å), the isotropic part and the nine components (MHz, per unit nuclear g-factor; CCEX
multiplies them by the species' g-factor).

```
VASP                                         the DFT code (format v3)
# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = VASP OUTCAR, S = 1
AtomName   r_1i_X  r_1i_Y  r_1i_Z  iso[MHz]  Axx Axy Axz Ayx Ayy Ayz Azx Azy Azz
10B  g-factor___:  0.600215 [No-unit]       one line per species
10B  etc____:      ---  ---  ---  0.00000    read, not used
MinDif[A]:  -11.31  -6.53  0.00              the range of relative positions the file covers
MaxDif[A]:   10.06   6.53  0.00
B   -10.05626  -4.36351  0.00000  0.000000  0.030490 ... -0.022870   one row per atom
```

Format v3: line 1 is `VASP` or `QUANTUMESPRESSO`, and line 2 the provenance line shown above. The numbers are the
hyperfine of the spin-S defect, read as written (SpinTensorReader divides Quantum ESPRESSO's S = 1/2 coupling by 2S;
VASP's is already per S). Refused: a first line that is not the tag of the source, the provenance line elsewhere than
line 2 or twice, an unknown format or source, and an S other than the qubit spin. A file without the provenance line
(older than v3) is read by the earlier rule with a warning (`VASP` × 1, otherwise × 0.5/S); remake such files with
SpinTensorReader.

A bath spin is matched to a row by its position relative to the qubit. A spin outside the range gets a point-dipole
tensor; one inside the range with no row is decided by [missing_row](#tag_hyperfine_missing_row). The positions are
in the bath file's frame; the frame of the components is [hyperfine.frame](#tag_hyperfine_frame). One qubit only.

```json
"hyperfine": {"dipole": "dft", "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}
```

<a id="tag_hyperfine_missing_row"></a>
### hyperfine.missing_row

| Type | Default |
|---|---|
| string | `"stop"` |

What to do with a bath spin inside the range of the hyperfine file that has no row in it. Spins outside the range
always get a point-dipole tensor.

| Value | Does |
|---|---|
| `"stop"` | stops the run and prints the position that was not found |
| `"point_dipole"` | gives the spin a point-dipole tensor and goes on |

```json
"hyperfine": {"dipole": "dft", "file": "Afile", "missing_row": "point_dipole"}
```

<a id="tag_hyperfine_frame"></a>
### hyperfine.frame

| Type | Default |
|---|---|
| string | required with [frame.rotation](#tag_frame_rotation) on and `dipole` `"dft"` |

The frame the tensor components of [hyperfine.file](#tag_hyperfine_file) are written in. The position columns of the
file are always in the bath file's frame.

| Value | The components are |
|---|---|
| `"file"` | in the bath file's frame: turned with the positions, R T R<sup>T</sup> |
| `"qubit"` | already in the computational frame (z = `qubit_axis`): used as written; a spin the file does not cover still gets a turned point-dipole tensor |

With `dipole` `"point"` the tensors come from the positions and are always turned (the file then gives only the
isotropic term), so a written frame is not used. With the rotation off, a written frame is not used either (a
warning).

```json
"hyperfine": {"dipole": "dft", "file": "Afile", "frame": "file"}
```

<a id="tag_quadrupole"></a>
### quadrupole

| Type | Default |
|---|---|
| `"off"`, or object `{"structure", "file", "file_without_defect", "cell", "frame"}`, in `hamiltonian.common` | `"off"` |

Where the quadrupole (electric field gradient) tensors of the bath spins come from. An object reads a DFT
[file](#tag_quadrupole_file).

| structure | Files |
|---|---|
| (`"off"`) | no quadrupole interaction |
| `"plain"` | one DFT file computed with the defect, `file` |
| `"strained"` | `file` with the defect and `file_without_defect` without it; `cell` required |
| `"bubble"` | the same two files and `cell` |

- A bath spin inside the range of `file` takes the row at its position relative to the defect; a spin inside the
  range with no row stops the run.
- Outside the range, `"plain"` takes the species default of the file (`etc`). `"strained"` and `"bubble"` take the
  tensor of [file_without_defect](#tag_quadrupole_file_without_defect), with the relative position folded by
  [cell](#tag_quadrupole_cell) along every axis whose cell parameter is not zero.

One qubit only. With [frame.rotation](#tag_frame_rotation) on, [frame](#tag_quadrupole_frame) gives the frame of the
files' components.

```json
"hamiltonian": {"common": {"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}}}
```

<a id="tag_quadrupole_file"></a>
### quadrupole.file

| Type | Default |
|---|---|
| string (a file) | required with every `structure` |

The quadrupole tensor file computed with the defect, written by [SpinTensorReader](Tools.md#tag_tool_spintensorreader).
It has the layout of [hyperfine.file](#tag_hyperfine_file), with these differences: the nine electric field gradient
components (Hartree/Bohr²) in place of the hyperfine ones and no isotropic column; one `eQ___` line per species (its
quadrupole moment) in place of the g-factor; an `etc____` line with the nine components, used outside the range; and no
code tag or provenance line at the top.

```json
"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}
```

<a id="tag_quadrupole_file_without_defect"></a>
### quadrupole.file_without_defect

| Type | Default |
|---|---|
| string (a file) | required with `"strained"` and `"bubble"` |

The quadrupole tensor file of the same structure computed without the defect, in the layout of
[quadrupole.file](#tag_quadrupole_file). It supplies the tensors of the spins outside the range of the file with the
defect.

```json
"quadrupole": {"structure": "strained", "file": "Q_defect", "file_without_defect": "Q_pristine", "cell": [2.5, 2.5, 0]}
```

<a id="tag_quadrupole_cell"></a>
### quadrupole.cell

| Type | Default | Unit |
|---|---|---|
| three numbers | required with `"strained"` and `"bubble"` | Å; or `{"value", "unit"}` with `A`, `nm` |

The lattice parameters of the strained or bubble structure.

```json
"cell": [1.0, 2.0, 3.0]
```

<a id="tag_quadrupole_frame"></a>
### quadrupole.frame

| Type | Default |
|---|---|
| string | required with [frame.rotation](#tag_frame_rotation) on and a quadrupole structure |

The frame the components of [quadrupole.file](#tag_quadrupole_file) and
[file_without_defect](#tag_quadrupole_file_without_defect) are written in. With the rotation off, a written frame is not
used (a warning).

| Value | The components are |
|---|---|
| `"file"` | in the bath file's frame: turned with the positions |
| `"qubit"` | already in the computational frame: used as written |

```json
"quadrupole": {"structure": "plain", "file": "Qfile", "frame": "file"}
```

<a id="gen_ham_method"></a>
## 5. hamiltonian: per method

<a id="tag_hfmedi"></a>
### hfmedi

| Type | Default |
|---|---|
| boolean, in `hamiltonian.cCCE` | `false` |

Adds the hyperfine-mediated interaction between bath spins, the second-order coupling through the qubit (N. Zhao et
al., Nature Nanotechnology 6, 242 (2011), Supplementary, Eq. 9). CCEX adds the full second-order correction, for each
qubit branch α on its own:

$$H_\text{med}^{(\alpha)} = \sum_{\beta\neq\alpha} \frac{E_{\alpha\beta} E_{\alpha\beta}^\dagger}{\omega_\alpha - \omega_\beta}, \qquad E_{\alpha\beta} = \sum_i \langle\alpha|\mathbf{S}|\beta\rangle\cdot\mathsf{A}_i\cdot\mathbf{I}_i$$

Here α is the eigenstate of the qubit Hamiltonian closest to the branch state, and β runs over every other eigenstate.
Expanded, the sum holds both orders of every pair (i, j) and the i = j terms of each bath spin; it is Hermitian.
[Advanced 06](../../example/tutorial/advanced/06_terms_and_models/README.md) runs with and without it.

```json
"hamiltonian": {"cCCE": {"hfmedi": true}}
```

<a id="tag_knight"></a>
### knight

| Type | Default |
|---|---|
| boolean, in `hamiltonian.cCCE` | `false` |

Adds the Knight-field shift of the bath spins that the qubit induces under a transverse magnetic field (same reference
as [hfmedi](#tag_hfmedi), SI Eq. 14). For each bath spin:

$$H_K^{(\alpha)} = -\gamma_n \mathbf{B}\cdot(\mathsf{g} - 1)\cdot\mathbf{I} = -(3|m_s| - 2)\,(\gamma_e/\Delta)\,\mathbf{B}\cdot\tilde{\mathsf{A}}\cdot\mathbf{I}$$

with Δ = 3/2 of the zz entry of the qubit's zero-field splitting (its self term in `Qubit.intmap`), Ã the hyperfine
tensor with its z row set to zero, and \|m<sub>s</sub>\| = ⟨S<sub>z</sub><sup>2</sup>⟩ of the branch's eigenstate.
The term is zero for a field along z. A qubit without a zero-field splitting is refused, since the term divides by Δ
(`Error: hamiltonian.cCCE.knight divides by the zero-field splitting of the qubit ...`).
[Advanced 06](../../example/tutorial/advanced/06_terms_and_models/README.md) runs with and without it in a tilted field.

```json
"hamiltonian": {"cCCE": {"knight": true}}
```

<a id="tag_delta_A"></a>
### delta_A

| Type | Default |
|---|---|
| boolean, in `hamiltonian.cCCE` | `false` |

Adds the second-order hyperfine correction δA of each bath spin (J. R. Maze, J. M. Taylor and M. D. Lukin, PRB 78,
094303 (2008), Appendix). For each qubit branch and each bath spin n, CCEX adds δA<sub>n</sub>(m<sub>s</sub>)·I<sub>n</sub>
to the branch's bath Hamiltonian:

$$\delta\mathbf{A}_n(m_s) = -\frac{m_s}{2\Delta}\,(\mathbf{A}_x\times\mathbf{A}_y)$$

with A<sub>x</sub>, A<sub>y</sub> the x and y rows of the spin's hyperfine tensor, Δ = 3/2 of the zz entry of the
qubit's zero-field splitting, and m<sub>s</sub> the branch's m<sub>s</sub>. The term is of size \|A\|²/Δ, so it
matters for nuclei close to the defect at low field. The coefficient −m<sub>s</sub> differs from Maze's Eq. (A4); the
derivation and the comparison are in [cCCE and gCCE](../methods/cce.md#cce_delta_A).

- It needs a single spin-1 qubit with a zero-field splitting, and does not run with `execution.propagator` sparse.
- It is refused together with [hfmedi](#tag_hfmedi), which already contains it; use delta_A for the same-nucleus
  correction without the full hyperfine-mediated term.
- It assumes a field along the qubit axis with γ<sub>e</sub>B ≪ Δ. CCEX warns when the qubit eigenstate closest to
  \|α⟩ or \|β⟩ has less than 0.99 of its weight on that \|m<sub>s</sub>⟩.

```json
"hamiltonian": {"cCCE": {"delta_A": true}}
```

<a id="tag_pseudospin_model"></a>
<a id="tag_pseudospin_channels"></a>
### channels (pseudospin)

| Type | Default |
|---|---|
| string, in `hamiltonian.pseudospin` | `"grouped"` |

How [method](#tag_method) pseudospin gives a pair's term P, and how its channels are evolved. Each pair's Hamiltonian
in each qubit state is written in the product basis of the two spins' own eigenstates and cut to its diagonal and the
elements between product states that differ in both spins: the flip-flops, the flip-flips, and with a mixing
quadrupole the double-quantum ones. These elements are the channels. The α and β product states of one index are taken
as one state (their overlap, the pair's share of the ESEEM, is left out), and with the nuclear state I/d each of the d
product states weighs 1/d.

| Value | Evolves |
|---|---|
| `"grouped"` | the blocks the channels join the product states into (a chain of channels is one block), each block of k states whole: L<sub>pair</sub> = L<sub>diag</sub> + Σ<sub>B</sub> (k<sub>B</sub>/d)(l<sub>B</sub> − l<sub>B,diag</sub>). Exact within the cut Hamiltonian |
| `"pairwise"` | every channel alone as a two-level system with the fields (Ω<sub>α,β</sub>, 0, Δ<sub>α,β</sub>): L<sub>pair</sub> = L<sub>diag</sub> + (2/d) Σ<sub>c</sub> (l<sub>c</sub> − l<sub>c,diag</sub>). It misses what couples the channels of a chain, and is kept to compare with that picture |

The two agree when no two channels share a state. Then P = L<sub>pair</sub> / (L<sub>1</sub> L<sub>2</sub>), with the
single spins from their eigenvalues alone; for two spins 1/2 of one species with the flip-flop alone this is Seo's eq.
23, P = 1 − K sin²(√(ΔA² + D²) t/4) sin²(D t/4), K = ΔA²/(ΔA² + D²). Without a pulse the single spins can pass through
zero; such a time point is reported as `[numerical-invalid]`. Every channel's Ω and Δ (kHz) are written to
`<output.file>_pseudospin_channels`. The model, its sources (Seo et al., Nat. Commun. 7, 12935 (2016); Lee, Park, Seo,
npj 2D Mater. Appl. 6, 60 (2022)) and its accuracy are in [the pseudospin model](../methods/pseudospin.md).

```json
"method": "pseudospin", "hamiltonian": {"pseudospin": {"channels": "pairwise"}}
```

<a id="tag_jump_operators"></a>
### jump_operators

| Type | Default | Unit |
|---|---|---|
| list of `{"bath_name", "operator", "rate"}`, in `hamiltonian.meCCE` or `hamiltonian.megCCE` | none: no dissipation | `rate` in MHz (≥ 0); or `{"value", "unit"}` with `Hz`, `kHz`, `MHz`, `GHz` |

The Lindblad jump operators of the master-equation methods, one entry per operator: C = √γ S<sub>op</sub> on every
spin named `bath_name`, with γ = `rate` × 2π·10³ rad/ms.

| Key | Value |
|---|---|
| `operator` | `-` (S<sub>-</sub>), `+` (S<sub>+</sub>), `z`, `x` or `y`, of a spin of any S |
| `bath_name` | the name of bath spins (as in the bath file, or of the Defect sub-spins; case not significant), or `"qubit"` / the qubit's name for a jump of the qubit (megCCE only, one qubit) |
| `rate` | in MHz, ≥ 0 |

Refused: a `bath_name` that is neither (checked against the whole bath; the message lists the names it has), a
`bath_name` that is both the qubit's and a bath spin's, and the section in a run of another method. Without the key,
or with an empty list, the run is the rate-0 limit (cCCE for meCCE, gCCE for megCCE). The equations are in
[master-equation CCE](../methods/master-equation.md), and
[Advanced 06](../../example/tutorial/advanced/06_terms_and_models/README.md) runs S₋ and S₊ on ¹³C.

```json
"method": "megCCE",
"hamiltonian": {"megCCE": {"jump_operators": [
    {"bath_name": "e",     "operator": "-", "rate": 1.5915e-4},
    {"bath_name": "e",     "operator": "+", "rate": 1.5915e-4},
    {"bath_name": "qubit", "operator": "z", "rate": 1.5915e-4}]}}
```

<a id="tag_noise_spectrum"></a>
### spectrum (noise)

| Type | Default |
|---|---|
| `"bath"`, or a model object, in `hamiltonian.noise` (method noise only) | `"bath"` |

Where the noise of the qubit transition frequency comes from. The coherence is $W(T) = e^{-\chi(T)}$ with
$\chi = \langle\varphi^2\rangle/2$, where $\varphi = \int f(t)\,\delta\omega(t)\,dt$ is the qubit's phase and
$f = \pm 1$ flips at each ideal π pulse. This Gaussian approximation needs many weakly coupled noise sources; it has no
back-action of the qubit on the bath, so strongly coupled spins need cCCE or gCCE.

| Value | Noise |
|---|---|
| `"bath"` | the bath of the run: each cCCE cluster gives $\chi_C = \frac{1}{2D}\sum_{ab} \lvert\langle a\rvert\delta\omega\lvert b\rangle\rvert^2 \lvert Y(\omega_{ab})\rvert^2$, $Y(\omega) = \int f(t) e^{i\omega t} dt$, over the eigenstates of $H_\text{bath} + p \sum_i \mathbf{A}_{z,i}\cdot\mathbf{I}_i$, $p = (\lvert m_\alpha\rvert + \lvert m_\beta\rvert)/2$, at infinite temperature, with $\delta\omega = (m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$ |
| `{"model": "OU", "sigma": σ, "tau_c": τ_c}` | Ornstein–Uhlenbeck noise, $\langle\delta\nu(t)\delta\nu(0)\rangle = \sigma^2 e^{-\lvert t\rvert/\tau_c}$ (σ in kHz, τ_c in ms; both take a unit); χ exact for any ideal π sequence |
| `{"model": "power", "amplitude": A, "beta": β, "f_min": f_min, "f_max": f_max}` | $S(f) = A\,(f/1\,\text{kHz})^{-\beta}$ (kHz²/kHz, two-sided) on [f_min, f_max] (kHz, both take a unit), $\chi = 4\pi^2 \int S(f) \lvert Y(2\pi f)\rvert^2 df$ by Gauss–Legendre; the work grows as $f_\text{max} T$ |

With `"bath"`, the χ<sub>C</sub> add over the CCE terms, so no e<sup>−χ<sub>C</sub></sup> of one cluster can
underflow. With `s_range`, a bath state enters only through the mean field of the other bath spins. For a Hahn echo
this is `semi-cce`'s `mode_sum`; [pseudo_secular](#tag_noise_pseudo_secular) adds the pseudo-secular terms to δω. With
method noise, `output.savemode` `cluster_coherence` and `cluster_density_matrix`, `rdip_ladder` and
`output.cluster_contributions` are refused.

A model run reads no bath: `configuration`, the Defect block and the clusters are not used, and `order`, `bfield`,
`rbath`, `rdip` and `Qubit` may be left out (when given, the log warns that they are not used).

```json
{"method": "noise", "hamiltonian": {"noise": {"spectrum": {"model": "OU", "sigma": 10, "tau_c": {"value": 50, "unit": "us"}}}},
 "deltat": 0.002, "nstep": 200, "pulse": {"count": 1}}
```

<a id="tag_noise_pseudo_secular"></a>
### pseudo_secular (noise)

| Type | Default |
|---|---|
| boolean, in `hamiltonian.noise` (with [spectrum](#tag_noise_spectrum) `"bath"`) | `false` |

Which terms of the hyperfine coupling between the qubit and a bath spin, $m_s (A_{zx} I_x + A_{zy} I_y + A_{zz} I_z)$,
count as the qubit's frequency noise.

| Value | δω |
|---|---|
| `false` | the secular term alone, $(m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$ (`semi-cce`'s) |
| `true` | also the pseudo-secular terms, $(m_\alpha - m_\beta) \sum_i (A_{zx,i} I_{x,i} + A_{zy,i} I_{y,i} + A_{zz,i} I_{z,i})$, the coupling cCCE puts into its two Hamiltonians |

The bath Hamiltonian is the same either way. Which of the two is the better approximation has not been settled. With a
model spectrum the key is not used.

```json
"hamiltonian": {"noise": {"spectrum": "bath", "pseudo_secular": true}}
```

<a id="tag_noise_compute"></a>
### compute (noise)

| Type | Default |
|---|---|
| string (case is not significant), in `hamiltonian.noise` (with [spectrum](#tag_noise_spectrum) `"bath"`) | `"coherence"` |

What the run computes from the bath's noise $B = \sum_i A_{zz,i} I_{z,i}$ ($\mathbf{A}_{z,i}\cdot\mathbf{I}_i$ with
[pseudo_secular](#tag_noise_pseudo_secular)), and which files it writes.

| Value | Computes | Writes |
|---|---|---|
| `"coherence"` | χ<sub>C</sub> of each cluster from its eigenmodes, the coherence exp(−χ) at the run's times; exact for ideal π pulses | `_wiDiv` (and `_noDiv`) |
| `"spectrum"` | `semi-cce`'s `psd`: C(t) on the grid of [fft](#tag_noise_fft), S(f) by FFT, and $\chi = 2 (m_\alpha - m_\beta)^2 \int S(f)\, F(2\pi f T)/f^2\, df$ | `_noise_correlation`, `_noise_spectrum`, `_wiDiv` |
| `"correlation"` | PyCCE's noise (`CCENoise`): C(t) at the run's times | `_noise_correlation` |

$C(t) = \operatorname{Re}\operatorname{Tr}[\rho\,\{B(t), B\}]/2$ in kHz², summed over the CCE terms
(`<output.file>[_state<s>]_noise_correlation`: t in ms, C). S(f) is two-sided, in kHz²/kHz (`_noise_spectrum`: f in
kHz, S). The state ρ is 1/D (infinite temperature) for `"coherence"` and `"spectrum"`, where a sampled state enters
through the mean field alone; for `"correlation"` with `s_range` it is the cluster's pure state, as in PyCCE. Without
`s_range` the three agree, and `"spectrum"` gives the coherence of `"coherence"` to the accuracy of its grid. A model
spectrum computes the coherence alone.

```json
"hamiltonian": {"noise": {"compute": "correlation"}}
```

<a id="tag_noise_fft"></a>
### fft (noise)

| Type | Default |
|---|---|
| object, in `hamiltonian.noise` (with [compute](#tag_noise_compute) `"spectrum"`) | every key has a default |

The time grid of C(t) and the FFT, `semi-cce`'s `psd` keys.

| Key | Meaning | Default |
|---|---|---|
| `time_max` | the end of the grid (ms, takes a unit): C(t) at t = 0, Δt, ..., `time_max` | the run's last time, `deltat` × (`nstep` − 1) |
| `nstep` | its intervals (`nstep` + 1 points, Δt = `time_max`/`nstep`) | from the bath's modes (below) |
| `window` | `"none"`, `"hann"`, `"tukey"` (the last half tapered) or `"exp"` (exp(−3t/`time_max`)) on C(t) | `"none"` |
| `pad_factor` | the FFT length is `pad_factor` × (2 `nstep` + 1), made odd; the padding holds the mean of the last fifth of C(t) | 1 |
| `subtract_mean` | subtract that mean from C(t) first (and pad with 0) | `false` |
| `clamp_negative` | set S(f) < 0 to 0 before the coherence | `false` |

S(f) is Δt Re FFT of C(t) made even, at f = k/(M Δt) for an FFT of length M: the highest frequency is 1/(2Δt) and the
resolution 1/(M Δt). Δt has to resolve the fastest bath frequency (the ¹³C Larmor period is 1.9 µs at 500 G), and
`time_max` has to cover the times of the coherence. Without `nstep`, Δt gives 8 points per period of the fastest modes
left after dropping at most 10⁻³ of their weight, and the log prints the values taken. All the modes then stay in
memory to the end of the run (24 bytes each; a cluster has dimension² of them); give `nstep` if there are too many. The
run stops if `nstep` would pass 10⁸. The default `time_max` is what the coherence needs; for a finer S(f), or with a
`window`, give a longer one.

```json
"hamiltonian": {"noise": {"compute": "spectrum", "fft": {"time_max": 3.0, "nstep": 15000, "window": "hann"}}}
```

<a id="gen_reference"></a>
## 6. Reference

<a id="tag_not_used"></a>
### Keys read and not used

A key this run reads but does not use is reported once, as `Warning: input key "<key>" is not used by this run: <why>`.
The calculation goes on, and the result is the one without the key.

| Key | Not used when |
|---|---|
| `Qubit.psi0` | cCCE, meCCE (the qubit is projected out: alpha and beta only) |
| `hamiltonian.common.qubit_bath_mean_field` | cCCE, meCCE, and the ensemble (no sampled state) |
| `hamiltonian.common.bath_bath_mean_field` | the ensemble |
| `configuration.bath_state_file`, `defect_state_file` | the ensemble (no `s_range`) |
| `configuration.defect_axis_file`, `defect_state_file` | no Defect species |
| `pulse.duration` | gCCE with `pulse.detuning_factor` 1 (cCCE, meCCE and megCCE refuse it) |
| `frame` (any) | `frame.rotation` off; a Defect species' `frame` without directional numbers |
| `cluster.selection.max_per_order` | an entry 0-2 other than 0, or every entry with dsj, itb, dsjitb |
| `execution.gpu_options` | `main.out` (built without CUDA; `main_gpu.out` uses them) |
| `cluster.adaptive.policy` | always: `merge` is the only policy |
| `hamiltonian.noise.fft` | `compute` other than `"spectrum"` |
| `hamiltonian.noise.pseudo_secular` | a model spectrum |
| `execution.backend` `cuda` | method noise with a model spectrum |
| `order`, `bfield`, `rbath`, `rdip`, `Qubit`, `configuration`, `Defect` | method noise with a model spectrum |

A key that no part of the run reads is reported as `Warning: input key "<key>" is not read by this run and is ignored`,
and a key that looks like a typo of a key that is read stops the run.

<a id="tag_turn_off"></a>
### Turning a part off

How each part that can be switched off is switched off; leaving the key out gives the last column. A block written
without `"enabled"` is on.

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
