# Input file
<a id="InputFile"></a>

The input file is JSON and sets every option of a run. This page shows how it is written and lists every key by group;
each group has a page of its own, and every key below links to its entry there. `main.out -e` writes a template with
every key ([Getting started](GettingStarted.md#OwnInput)). An input of an earlier version can be converted with the
table of [Old names](InputFile_OldNames.md).

1. [How to write the input file](#HowToWrite): the format, comments, types, and what the run reports about keys.
2. [Units](#Units): the unit of every number, and the keys that take another unit.
3. [The groups of keys](#KeyGroups): General, Qubit, Cluster, Execution, Pulse, Defect, Output and Automatic
   convergence.

<a id="HowToWrite"></a>
## 1. How to write the input file

The input is a JSON object of keys and values. A `//`, `#` or `!` outside a string starts a comment to the end of the
line; inside a string these characters are ordinary (`"file": "run#1"`). The file may be named `.json` or `.jsonc`.

```
// a comment line
{
    "method" : "cCCE",      // a comment after a value
    # a comment line of another kind
    ! and of a third
    "output" : {"file": "run#1"}
}
```

A key that is left out takes its default. Keys are read without regard to case. A key that is present must have the
type its entry gives; a wrong type stops the run and names the key (with its path in a block, e.g. `pulse.count`).
Examples of wrong types: `"knight": "true"`, a quoted number, `"nstep": 10.5`. A switch is `true` or `false`, not 0 or
1.

What the run reports about keys:

- A key that no part of the run reads gets `Warning: input key ... is not read`. When it is one or two letters away
  from a key that is absent from the same block, the run stops instead: `"s_rang"` would otherwise leave `s_range` out
  and compute the ensemble.
- A key the run reads but does not use gets `Warning: input key "..." is not used by this run: <why>`
  ([the cases](InputFile_General.md#tag_not_used)).
- A file the input names that does not open stops the run ([files](InputFile_General.md#tag_configuration_missing)).
- A key of an old name is refused with its current name, e.g. `Error: input key "npulse" is now "pulse.count"`
  ([Old names](InputFile_OldNames.md)).

A whole input:

```
{
    "method"       : "cCCE",
    "cluster"      : {"mode": "standard"},
    "quantity"     : "coherence",
    "Qubit"        : {"qubit": [{"file": "../bath/bath_Diamond_NV_natab_defect"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
    "configuration": {"gyrofile": "../bath/gyro_13C", "bathfile": ["../bath/bath_Diamond_NV_natab_1"], "seed": 1},
    "order"        : 2,
    "bfield"       : 500,
    "rbath"        : 50,
    "rdip"         : 10,
    "deltat"       : 0.001,
    "nstep"        : 2000,
    "pulse"        : {"count": 1},
    "hamiltonian"  : {"common": {"bath_bath_mean_field": "z"}, "cCCE": {"knight": false}},
    "output"       : {"path": "results"}
}
```

<a id="Units"></a>
## 2. Units

A plain number is in the unit its key names (the table below). The keys of the table also take
`{"value": ..., "unit": "..."}`. The value is what the plain form takes: a number, or a list for `bfield`, `rbath`,
`rdip`, `bath_shift`, `xyz` and `energy_shift`.

| Key | Unit of a plain number | Units of `{"value", "unit"}` |
|---|---|---|
| `deltat`, `auto_convergence.time_window.value` | ms | `s`, `ms`, `us`, `ns` |
| `pulse.duration` | ns | `s`, `ms`, `us`, `ns` |
| `bfield` | G | `G`, `mT`, `T` |
| `rbath`, `rdip`, `configuration.bath_shift`, `Qubit.qubit[].xyz` | Å | `A` (Å), `nm` |
| `pulse.bath.energy_shift`, `pulse.bath.bandwidth`, `hamiltonian.<method>.jump_operators[].rate` | MHz | `Hz`, `kHz`, `MHz`, `GHz` |
| `Qubit.qubit[].detuning` | **kHz** (a Defect detuning is MHz) | `Hz`, `kHz`, `MHz`, `GHz` |
| `cluster.bundle.merge_cutoff` | rad/ms | `rad/s`, `rad/ms`, `rad/us` |
| `hamiltonian.common.quadrupole.cell` | Å | `A` (Å), `nm` |
| `Qubit.qubit[].gyro` | rad kHz/G | the units of Defect `gyros` (below); a unit other than rad kHz/G multiplies by its factor (2π for the plain-frequency units) |

```
"deltat": {"value": 5, "unit": "us"},              # = "deltat": 0.005
"bfield": {"value": [0, 0, 50], "unit": "mT"},     # = "bfield": [0, 0, 500]
"rdip":   {"value": [0.1, 1.0], "unit": "nm"}      # = "rdip": [1, 10]
```

The object is turned into the plain number before the input is read: multiplied by a power of ten for a smaller unit
(`nm` to Å: × 10), divided for a larger one (`us` to ms: / 1000), so `{"value": 5, "unit": "us"}` is exactly `0.005`.

Refused, with the units of the key: an object with other keys, a unit not in its row, and a value that is not a number
or a list of numbers.

Units are case-sensitive and written as in the table (`mT`, not `MT`), and so are the units of the Defect block and of
`Qubit.intmap` (`"mhz"` is refused).

Some keys keep a form of their own:

- the Defect `gyros`, `quadrupole_moments`, `hyperfine`, `efg` and `detuning`: `{"values": ..., "unit": ...}`.
  `"value"`, the word of the other keys, is taken there too; write one of the two.
- the Defect `zfs`: `{"D", "axis", "unit"}` or `{"tensor", "unit"}`;
- a `Qubit.intmap` entry: its own `unit`.

The units do not apply to the numbers inside the bath, state and tensor files (Å), to the gyromagnetic ratios, or to
the time column of the output files (ms).

In `auto_convergence.time_window`, `"value"` is a number of ms or `{"value": ..., "unit": ...}`. A `time_window.unit`
beside it was removed and is refused.

<a id="FixedUnits"></a>
### Keys without a unit object

Every other number has one fixed unit, or its own `unit` sub-key, as the table says.

The zero-field splitting and the detuning have different units in the two blocks:

- in the `Qubit` block, **kHz**: a `Qubit.intmap` tensor without the entry's `unit`, and `Qubit.qubit[].detuning`;
- in the `Defect` block, **MHz**: `zfs` and `Defect.species[].detuning`.

For example, a `Qubit.intmap` self entry with the tensor [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]] and a
Defect `zfs` {"D": 2880, "unit": "MHz"} describe the same D = 2.88 GHz.

| Key or file | Unit | Other units |
|---|---|---|
| the γ column of `configuration.gyrofile` | rad/ms/G (= rad kHz/G) | none |
| `Qubit.intmap[].tensor` | kHz | the entry's `unit`: `Hz`, `kHz`, `MHz`, `GHz` |
| Defect `zfs` (`{D, axis, unit}`, `{tensor, unit}` or the indexed list), `hyperfine`, `detuning` | MHz | `unit` (in `{"values", "unit"}` for `hyperfine` and `detuning`): `Hz`, `kHz`, `MHz`, `GHz` |
| Defect `gyros` | rad kHz/G (= rad/ms/G) | `unit`: `radkHz/G`, `rad/ms/G`, `kHz/G`, `MHz/G`, `MHz/T`, `Hz/T`, `radMHz/T`, `rad/us/T`, `rad/s/T` |
| Defect `quadrupole_moments` | fm² (10⁻³⁰ m²) | `unit`: `fm^2`, `m^2`, `cm^2`, `barn`, `mbarn` and spellings of these |
| Defect `efg` | Hartree/Bohr² (atomic units) | `unit`: `Hartree/Bohr^2` (`a.u.`), `V/m^2`, `V/A^2` |
| Defect `nuclear_positions` | Å | none |
| bath file positions, the position file of `Qubit.qubit[].file` | Å | none |
| hyperfine tensor file (`hyperfine.file`) | MHz, per unit nuclear g-factor | none |
| quadrupole tensor file (`quadrupole.file`) | Hartree/Bohr² | none |
| the time column of the output files | ms | none |

<a id="KeyGroups"></a>
## 3. The groups of keys

<a id="GeneralTag"></a>
### General tags

[General tags](InputFile_General.md): the method, the numbers of the calculation, the files of the bath, the frame, and the `hamiltonian` block.

- [method](InputFile_General.md#tag_method)
- [quantity](InputFile_General.md#tag_quantity)
- [order](InputFile_General.md#tag_order)
- [bfield](InputFile_General.md#tag_bfield)
- [rbath](InputFile_General.md#tag_rbath) (a number or `[inner, outer]`)
- [rdip](InputFile_General.md#tag_rdip) (a number or `[inner, outer]`)
- [deltat](InputFile_General.md#tag_deltat)
- [nstep](InputFile_General.md#tag_nstep)
- [configuration](InputFile_General.md#tag_configuration): the bath configuration
    - [path](InputFile_General.md#tag_configuration_path), [k_range](InputFile_General.md#tag_configuration_k_range)
      and [{k}](InputFile_General.md#tag_configuration_k): several configurations in one run, outputs `<output.file>_k<k>_...`
    - [gyrofile](InputFile_General.md#tag_gyrofile), [bathfile](InputFile_General.md#tag_bathfile),
      [bath_shift](InputFile_General.md#tag_bath_shift), [bath_state_file](InputFile_General.md#tag_bath_state_file)
      (one per bath file)
    - [defect_axis_file](InputFile_General.md#tag_defect_axis_file),
      [defect_state_file](InputFile_General.md#tag_defect_state_file) (one per bath file)
    - [s_range](InputFile_General.md#tag_s_range) (absent: the ensemble), [seed](InputFile_General.md#tag_seed)
- [frame](InputFile_General.md#tag_frame): [rotation](InputFile_General.md#tag_frame_rotation),
  the qubit-aligned frame, on when the block is given (`enabled`, `qubit_axis`). Each set of directional numbers names
  its frame with a key `frame` (`"file"` or `"qubit"`) beside it.
- [hamiltonian](InputFile_General.md#tag_hamiltonian)
    - common: [qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field),
      [bath_bath_mean_field](InputFile_General.md#tag_bath_bath_mean_field),
      [hyperfine](InputFile_General.md#tag_hyperfine) (`dipole`, `contact`, [file](InputFile_General.md#tag_hyperfine_file),
      [missing_row](InputFile_General.md#tag_hyperfine_missing_row), [frame](InputFile_General.md#tag_hyperfine_frame)),
      [quadrupole](InputFile_General.md#tag_quadrupole) (`structure`, [file](InputFile_General.md#tag_quadrupole_file),
      [file_without_defect](InputFile_General.md#tag_quadrupole_file_without_defect), [cell](InputFile_General.md#tag_quadrupole_cell),
      [frame](InputFile_General.md#tag_quadrupole_frame))
    - cCCE: [knight](InputFile_General.md#tag_knight), [hfmedi](InputFile_General.md#tag_hfmedi),
      [delta_A](InputFile_General.md#tag_delta_A) (meCCE reads them too)
    - meCCE, megCCE: [jump_operators](InputFile_General.md#tag_jump_operators)
      ([master-equation.md](../methods/master-equation.md))
- [keys read and not used](InputFile_General.md#tag_not_used), [turning a part off](InputFile_General.md#tag_turn_off)

<a id="QubitArrayTags"></a>
### Qubit tags

On the page [Qubit tags](InputFile_Qubit.md).

The `Qubit` block: [qubit](InputFile_Qubit.md#tag_qubit), one entry per qubit, with its position written out (`xyz`) or
read from a file (`file`). Beside the list: `intmap`, `alpha_state`, `beta_state` and `psi0`. For example:

```
"Qubit" : {"qubit": [{"file": "file-path"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}}
```

- [Qubit](InputFile_Qubit.md#tag_Qubit)
- [qubit](InputFile_Qubit.md#tag_qubit)
    - [name](InputFile_Qubit.md#tag_name), [spin](InputFile_Qubit.md#tag_spin), [gyro](InputFile_Qubit.md#tag_gyro),
      [xyz](InputFile_Qubit.md#tag_xyz) or [file](InputFile_Qubit.md#tag_Qubit_file), [detuning](InputFile_Qubit.md#tag_detuning)
- [intmap](InputFile_Qubit.md#tag_intmap): [between](InputFile_Qubit.md#tag_between), [tensor](InputFile_Qubit.md#tag_tensor),
  [unit](InputFile_Qubit.md#tag_intmap_unit), [frame](InputFile_Qubit.md#tag_intmap_frame)
- [alpha_state](InputFile_Qubit.md#tag_alpha_state), [beta_state](InputFile_Qubit.md#tag_beta_state) (`{"ms": [...]}`, `{"index": n}` or `{"vector": [...]}`)
- [psi0](InputFile_Qubit.md#tag_psi0)
- the Overhauser field: [hamiltonian.common.qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field)

<a id="ClusterTags"></a>
### Cluster tags

On the page [Cluster tags](InputFile_Cluster.md).

- [order](InputFile_General.md#tag_order) and [method](InputFile_General.md#tag_method) (top level)
- the `cluster` block:
    - [mode](InputFile_Cluster.md#tag_cluster_mode) (`standard`, `pcce`, `bundle`, `dsj`, `itb`, `dsjitb`)
    - `selection`: [add_subclusters](InputFile_Cluster.md#tag_add_subclusters),
      [max_per_order](InputFile_Cluster.md#tag_max_clusters_per_order)
    - [pcce](InputFile_Cluster.md#tag_pcce): [partition_size](InputFile_Cluster.md#tag_pcce_partition_size),
      [trials](InputFile_Cluster.md#tag_pcce_trials), [max_iterations](InputFile_Cluster.md#tag_pcce_max_iterations),
      [kmeans_pp](InputFile_Cluster.md#tag_pcce_kmeans_pp), [print_trials](InputFile_Cluster.md#tag_pcce_print_trials),
      [partition](InputFile_Cluster.md#tag_pcce_partition),
      [connections](InputFile_Cluster.md#tag_pcce_connections), [cache](InputFile_Cluster.md#tag_pcce_cache)
    - [bundle](InputFile_Cluster.md#tag_bundle): [max_dim](InputFile_Cluster.md#tag_bundle_max_dim),
      [merge_cutoff](InputFile_Cluster.md#tag_bundle_merge_cutoff), [connections](InputFile_Cluster.md#tag_bundle_connections),
      [connection_budget](InputFile_Cluster.md#tag_bundle_connection_budget), [max_clusters](InputFile_Cluster.md#tag_bundle_max_clusters),
      [max_work](InputFile_Cluster.md#tag_bundle_max_work)
    - [adaptive](InputFile_Cluster.md#tag_adaptive) ([probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb),
      [work_file](InputFile_Cluster.md#tag_adaptive_work_file), internal)

<a id="BackendTags"></a>
### Execution tags

On the page [Execution tags](InputFile_Backend.md).

Where and how the cluster dynamics is computed, and what the run reports about it: the `execution` block.

- [backend](InputFile_Backend.md#tag_backend)
- [gpu_options](InputFile_Backend.md#tag_gpu_options) (min_dim, packed_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb, warp_max_dim, ranks_per_gpu)
- [propagator](InputFile_Backend.md#tag_propagator)
- [evolution](InputFile_Backend.md#tag_evolution)
- [evolution_update](InputFile_Backend.md#tag_evolution_update)
- [gcce_refinement](InputFile_Backend.md#tag_gcce_refinement) (writing the block turns it on; `enabled`, [max_sweeps](InputFile_Backend.md#tag_gcce_refinement_max_sweeps))
- [column_evolution](InputFile_Backend.md#tag_column_evolution) (on by default; `enabled`, `min_dim`, `columns`)
- [profile](InputFile_Backend.md#tag_profile)

On this page too, in other blocks: [cluster.adaptive.probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb) and
[output.cluster_contributions](InputFile_Output.md#tag_cluster_contributions) (`enabled`, [max_mb](InputFile_Output.md#tag_cluster_contribution_max_mb)).

<a id="PulseTags"></a>
### Pulse tags

On the page [Pulse tags](InputFile_Pulse.md).

- [pulse](InputFile_Pulse.md#tag_pulse) (required)
    - [count](InputFile_Pulse.md#tag_pulse_count)
    - [name](InputFile_Pulse.md#tag_pulse_name)
    - [sequence](InputFile_Pulse.md#tag_sequence)
    - [per_qubit](InputFile_Pulse.md#tag_pulse_per_qubit)
    - [duration](InputFile_Pulse.md#tag_pulse_duration)
    - [detuning_factor](InputFile_Pulse.md#tag_detuning_factor)
    - [bath](InputFile_Pulse.md#tag_pulse_bath) (DEER bath pulses, gCCE and megCCE; [deer.md](../methods/deer.md))
        - [count](InputFile_Pulse.md#tag_pulse_bath_count)
        - [sequence](InputFile_Pulse.md#tag_pulse_bath_sequence)
        - [defect](InputFile_Pulse.md#tag_pulse_bath_defect)
        - [energy_shift](InputFile_Pulse.md#tag_pulse_bath_energy_shift)
        - [bandwidth](InputFile_Pulse.md#tag_pulse_bath_bandwidth)
        - [levels](InputFile_Pulse.md#tag_pulse_bath_levels)

<a id="DefectArrayTags"></a>
### Defect tags

On the page [Defect tags](InputFile_Defect.md).

The `Defect` block of a paramagnetic defect bath: `species`, the list of the defect species (the model). Their axis and
state files are keys of the configuration
([configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file),
[configuration.defect_state_file](InputFile_General.md#tag_defect_state_file)). The full example is on the
[Defect page](InputFile_Defect.md).

- [Defect](InputFile_Defect.md#tag_Defect)
- [species](InputFile_Defect.md#tag_df_species)
    - [name](InputFile_Defect.md#tag_df_name)
    - [n_onsite_spins](InputFile_Defect.md#tag_df_n_onsite_spins)
    - [n_axis_states](InputFile_Defect.md#tag_df_n_axis_states)
    - [frame](InputFile_Defect.md#tag_df_frame)
    - [types](InputFile_Defect.md#tag_df_types)
    - [spins](InputFile_Defect.md#tag_df_spins)
    - [gyros](InputFile_Defect.md#tag_df_gyros)
    - [quadrupole_moments](InputFile_Defect.md#tag_df_quadrupole_moments)
    - [nuclear_positions](InputFile_Defect.md#tag_df_nuclear_positions)
    - [hyperfine](InputFile_Defect.md#tag_df_hyperfine)
    - [efg](InputFile_Defect.md#tag_df_efg)
    - [zfs](InputFile_Defect.md#tag_df_zfs) (its [axis](InputFile_Defect.md#tag_df_axis))
    - [detuning](InputFile_Defect.md#tag_df_detuning)

<a id="OutputTags"></a>
### Output tags

On the page [Output tags](InputFile_Output.md).

The `output` block:

- [path](InputFile_Output.md#tag_output_path) (the output folder, made when missing; `-o DIR` takes its place; default: the working directory)
- [file](InputFile_Output.md#tag_outfile) (default: the input file's name without `.json` or `.jsonc`)
- [savemode](InputFile_Output.md#tag_savemode) (`mode`, `nodiv`; old names are refused)
- [cluster_contributions](InputFile_Output.md#tag_cluster_contributions) (`enabled`, `max_mb`)
- [cluster_only](InputFile_Output.md#tag_cluster_only) (`true`: the cluster plan only)

<a id="AutoConvergenceTags"></a>
### Automatic convergence tags

On the page [Automatic convergence tags](InputFile_AutoConvergence.md). The block's
sections are `time_window`, `rdip`, `rbath`, `t2_fit` and `sampling`; the bath configurations come from the
`configuration` block. `rdip_ladder` is an internal key it writes into its runs' inputs
([rdip_ladder](InputFile_General.md#tag_rdip_ladder)).
