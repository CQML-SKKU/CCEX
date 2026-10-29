# Input file
<a id="InputFile"></a>

Every key of the input file, in eight groups. Each group is a page of its own, and
every key in the lists below links to its entry there. Inputs written for an earlier version use other key names and
places: the table of [old names](#OldNames) is at the end of this page, and
`example/tools/InputKeys/convert_input_keys.py` rewrites an old input. `main.out -e` writes an input template
(`ccex_example.jsonc`) with every key ([Getting started](GettingStarted.md)).

<a id="HowToWrite"></a>
## How to write the input file

The input file consists of key-value pairs where each key represents a specific setting.
Each parameter should be written in the following format:

```
{
    parameter_name_1 : value_1,
    parameter_name_2 : value_2,
                :           
    parameter_name_N : value_N
}
```

A `//`, `#` or `!` outside a JSON string starts a comment that runs to the end of the line, so a whole line or the rest
of a line after a value can be a comment. Inside a string (`"#note"`, `"file": "run#1"`, `"http://x"`, `"a//b"`) they
are ordinary characters. An input written with `//` comments can be named `.jsonc` (JSON with comments), which vim and
VS Code colour as JSON with comments; the default output name drops `.jsonc` as it drops `.json`.

```
// a comment line
{
    "method" : "cCCE",      // a comment after a value
    # a comment line of another kind
    ! and of a third
    "output" : {"file": "run#1"}
}
```

A key that is left out takes its default. A key that is present has to have the type its
entry gives: `"knight": "true"`, a quoted number, `"nstep": 10.5` or a string inside a
number array stops the run with the key named (with its path in a block, e.g. `pulse.count`). A switch is `true` or `false` (the 0 and 1 the old
keys took are not read by the new ones). Keys are read without regard to case.

A key that no part of the run reads is reported (`Warning: input key ... is not read`), and
when it is one or two letters away from a key that was looked up and is absent from the same
block, the run stops: `"s_rang"` would leave `s_range` out (the ensemble), `"alpha_stat"` in the `Qubit` block its `alpha_state`.
A key the run reads but does not use is reported with one line for every such key,
`Warning: input key "..." is not used by this run: <why>` ([the cases](InputFile_General.md#tag_not_used)), and a file
the input names that does not open stops the run ([files](InputFile_General.md#tag_configuration_missing)).
A key of an old name is refused with its current name (`Error: input key "npulse" is now "pulse.count"`) and
a pointer to the converter ([old names](#OldNames)).

A whole input in the new keys:

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
### Units

A number is in the unit its key names (the table). The keys of the table also take `{"value": ..., "unit": "..."}`,
the value being what the plain form takes (a number, or a list for `bfield`, `rbath`, `rdip`, `bath_shift`, `xyz`,
`energy_shift`):

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

The object is turned into the plain number before any reader sees it, by one exact operation: a multiplication by a power
of ten going to a smaller unit (`nm` to Å: × 10), a division by one going to a larger (`us` to ms: / 1000; 0.001 has no
exact binary form, x / 1000 is the double nearest the decimal value). So `{"value": 5, "unit": "us"}` is the number
`0.005` bit for bit, and a plain number is read as it is. An object with other keys, a unit not in its row, or a value
that is not a number or a list of numbers is refused with the units of the key. The units are written as in the table
(case counts: `mT`, not `MT`), and so are the units of the Defect block and `Qubit.intmap` (they were compared
without regard to case: `"mhz"` is refused with `"MHz"`, and the converter writes the table's spelling). The keys that
take `{"values": ..., "unit": ...}` (the Defect `gyros`, `quadrupole_moments`, `hyperfine`, `efg`, `detuning`; `"value"`
is taken there too, the word of the other keys, one of the two) and the Defect `zfs` (`{"D", "axis", "unit"}` or
`{"tensor", "unit"}`) and a `Qubit.intmap` entry (its own `unit`) keep that form. The units do not apply to the numbers inside the bath, state
and tensor files (Å), the gyromagnetic ratios, and the time column of the output files (ms). In
`auto_convergence.time_window`, `"value"` is a number of ms or `{"value": ..., "unit": ...}` (`time_window.unit` beside it was
removed and is refused).

<a id="FixedUnits"></a>
#### Keys without a unit object

Every other number has one fixed unit (or its own `unit` sub-key, as the table says). The same zero-field splitting is
**kHz** in the `Qubit` block (without the intmap entry's `unit`) and **MHz** in the `Defect` block: a `Qubit.intmap` self
entry with the tensor [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]] and a Defect `zfs` {"D": 2880, "unit": "MHz"}
describe the same D = 2.88 GHz. The same holds for a detuning: `Qubit.qubit[].detuning` is kHz, `Defect.species[].detuning`
MHz. A plain number other than 0 in `Qubit.qubit[].detuning`, and a nonzero `Qubit.intmap` tensor without `unit`, are noted
in the log as read in kHz.

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

<a id="GeneralTag"></a>
## General tags
The keys of this group are described in [InputFile_General.md](InputFile_General.md).

The method, the numbers of the calculation, the files of the bath, the frame, and the `hamiltonian` block.

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
  the qubit-aligned frame, on when the block is given (`enabled`, `qubit_axis`); each set of directional numbers says
  its frame with a key `frame` (`"file"` or `"qubit"`) beside it
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
## Qubit tags
The keys of this group are described in [InputFile_Qubit.md](InputFile_Qubit.md).

The `Qubit` block: [qubit](InputFile_Qubit.md#tag_qubit), one entry per qubit, with its position written out (`xyz`) or
read from a file (`file`), and `intmap`, `alpha_state`, `beta_state` and `psi0` beside the list.
<b>Example:</b> ``` "Qubit" : {"qubit": [{"file": "file-path"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}} ```
The full example is on the [Qubit page](InputFile_Qubit.md).

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
## Cluster tags
The keys of this group are described in [InputFile_Cluster.md](InputFile_Cluster.md).

- [order](InputFile_General.md#tag_order) and [method](InputFile_General.md#tag_method) (top level)
- the `cluster` block:
- [mode](InputFile_Cluster.md#tag_cluster_mode) (`standard`, `pcce`, `bundle`, `dsj`, `itb`, `dsjitb`)
- `selection`: [add_subclusters](InputFile_Cluster.md#tag_add_subclusters),
  [max_per_order](InputFile_Cluster.md#tag_max_clusters_per_order)
- [pcce](InputFile_Cluster.md#tag_pcce):[partition_size](InputFile_Cluster.md#tag_pcce_partition_size),
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
## Execution tags
The keys of this group are described in [InputFile_Backend.md](InputFile_Backend.md).

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
## Pulse tags
The keys of this group are described in [InputFile_Pulse.md](InputFile_Pulse.md).

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
## Defect tags
The keys of this group are described in [InputFile_Defect.md](InputFile_Defect.md).

The `Defect` block of a paramagnetic defect bath: `species`, the list of the defect species (the model). Their axis and
state files are the configuration's
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
## Output tags
The keys of this group are described in [InputFile_Output.md](InputFile_Output.md).

The `output` block:

- [path](InputFile_Output.md#tag_output_path) (the output folder, made when missing; `-o DIR` takes its place; default: the working directory)
- [file](InputFile_Output.md#tag_outfile) (default: the input file's name without `.json` or `.jsonc`)
- [savemode](InputFile_Output.md#tag_savemode) (`mode`, `nodiv`; old names are refused)
- [cluster_contributions](InputFile_Output.md#tag_cluster_contributions) (`enabled`, `max_mb`)
- [cluster_only](InputFile_Output.md#tag_cluster_only) (`true`: the cluster plan only)

<a id="AutoConvergenceTags"></a>
## Automatic convergence tags
The `auto_convergence` block (its sections `time_window`, `rdip`, `rbath`, `t2_fit` and
`sampling`; the bath configurations are the `configuration` block's) is described in
[InputFile_AutoConvergence.md](InputFile_AutoConvergence.md). `rdip_ladder` is an internal key it writes into its runs'
inputs ([rdip_ladder](InputFile_General.md#tag_rdip_ladder)).

<a id="OldNames"></a>
## Old names

Inputs written for an earlier version (v1.2.2 and before, or the adding_DEER / adding_MECCE branches) use other names
and places for many keys. An old name, value, form or place is refused with its current one
(`Error: input key "npulse" is now "pulse.count"`, `Error: input key "bathfile" is now "configuration.bathfile" ...`),
except where the table says otherwise. The converter rewrites old inputs, keeping their comments, and lists what it
changed:

```
python3 example/tools/InputKeys/convert_input_keys.py FILE...           # rewrite in place
python3 example/tools/InputKeys/convert_input_keys.py --check FILE...   # list only; exit 1 when a file has old keys
python3 example/tools/InputKeys/convert_input_keys.py --stdout FILE     # print the new text
```

It drops keys the run did not read (e.g. `qd_tensorfile` with `qd_readmode` 0) and the defaults `hf_readmode` 0 and
`qd_readmode` 0, and refuses an input it cannot express (a `CorrTotSpin` other than 0). An old key whose new name is
already in the input (`npulse` beside `pulse.count`; case not read) is kept once when the two values are the same, and
stops the file, before it is written, when they differ. Its tests are `example/tools/InputKeys/tests` (`python3 -m
unittest discover -s example/tools/InputKeys/tests`). The results do not change: rewritten in the current keys, the
inputs of the standing checks (check_numerics, check_input_forms, check_rotation) give the same output bytes.
`order`, `rbath`, `rdip`, `quantity`, `bfield`, `deltat`, `nstep` and `rdip_ladder` (an internal key of the automatic
convergence) are at the top level, and so are the keys of `Qubit.qubit` and `Qubit.intmap` entries and the
Defect species' `types`, `spins`, `gyros`, `efg`, `zfs` and `detuning`.

| Old | Now |
|---|---|
| `"method": "CCE"` | `"method": "cCCE"` |
| `"method": "pCCE"` | `"method": "cCCE", "cluster": {"mode": "pcce"}` |
| `"method": "dsj"` / `"itb"` / `"dsjitb"` | `"method": "cCCE", "cluster": {"mode": "dsj"}` / `"itb"` / `"dsjitb"` |
| `"method": "pgCCE"` / `"pmeCCE"` / `"pmegCCE"` (adding_MECCE) | `"method": "gCCE"` / `"meCCE"` / `"megCCE"` with `"cluster": {"mode": "pcce"}` |
| `cluster_mode` | `cluster.mode` |
| `"cluster_mode": "adaptive"` | `"cluster": {"mode": "bundle"}` |
| `Adaptive` block, top-level `adaptive` block | `cluster.adaptive` block (keys unchanged; `cluster.mode` pcce or bundle; keys are read without regard to case, so `Adaptive` inside `cluster` is read) |
| `Refine` block | `cluster.adaptive` block |
| `adaptive_work_file` | `cluster.adaptive.work_file` (written by the refinement driver) |
| `probe_memory_gb` | `cluster.adaptive.probe_memory_gb` (the converter removes it from an input without an `adaptive` block: it was read with the refinement only) |
| `rbathcut` | `rbath: [inner, outer]` (the inner radius) |
| `rdipcut` | `rdip: [inner, outer]` (the inner radius) |
| `gyrofile`, `bathfile`, `nstate`, `seed` (top level) | `configuration.gyrofile`, `configuration.bathfile`, `configuration.s_range` (as `configuration.nstate` below), `configuration.seed` |
| `configuration.nstate` N | `configuration.s_range: [1, N]`; `nstate` 0 (the ensemble): no `s_range` |
| a `configuration.seed` number with `k_range` | the number is now the first configuration's seed (k takes seed + k - first); the old meaning, the same seed for every k, is a list of that number, one per configuration, which the converter writes where the run draws (sampled states, Defect species, pcce) |
| `bathadjust`, `bath_shift` (top level) | `configuration.bath_shift` |
| `statefile`, `bath_state_file: "prefix"` (top level) | `configuration.bath_state_file: ["prefix"]`, one name per bath file (`{s}` the state, or the state appended) (the converter removes a single prefix beside several bath files: such a prefix was not read, the states were drawn at random) |
| `avaaxfile`, `Defect.axis_file: "file"` | `configuration.defect_axis_file: ["file"]`, one per bath file |
| `exstatefile`, `Defect.state_file: "prefix"` | `configuration.defect_state_file: ["prefix"]`, one per bath file |
| `addsubclus` (0/1 or true/false), `add_subclusters` (top level) | `cluster.selection.add_subclusters` (true/false) |
| `nk`, `max_clusters_per_order` | `cluster.selection.max_per_order` |
| `qubitfile` | `Qubit.qubit: [{"file": ...}]` (as `Qubit.file` below) |
| `qspin`, `qzfs`, `qalphams`, `qbetams` | the entry's `spin`, then as `Qubit.zfs`, `Qubit.alphams`, `Qubit.betams` below |
| `Qubit.file`, `Qubit.spin`, `Qubit.name` (the one-qubit file form) | an entry of `Qubit.qubit`: `[{"file": ..., "spin": ..., "name": ...}]` (`file` in place of `xyz`) |
| `Qubit.qubit[].detuning: [x]` | `x` (a number, or `{"value", "unit"}`) |
| `Qubit.nqubit` | removed: the length of `Qubit.qubit` |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| `Qubit.alphaidx`, `Qubit.betaidx` | `Qubit.alpha_state: {"index": n}`, `Qubit.beta_state: {"index": n}` |
| `Qubit.psia`, `Qubit.psib` | `Qubit.alpha_state: {"vector": [...]}`, `Qubit.beta_state: {"vector": [...]}` (each side on its own; old inputs used `psia` only with `psib` and `alphaidx` only with `betaidx`, and a vector won over an index) |
| `Qubit.qubit[].alphams`, `Qubit.qubit[].betams` | `Qubit.alpha_state: {"ms": [m_1, m_2, ...]}`, `Qubit.beta_state: {"ms": [...]}`, one m<sub>s</sub> per qubit in the order of `Qubit.qubit` (the converter fills +1 / 0 for a qubit that wrote none; beside an `alpha_state` or `beta_state` already written, which CCEX used for the state, it removes them with a note) |
| `Qubit.alphams`, `Qubit.betams` (with `Qubit.file`) | `Qubit.alpha_state: {"ms": [m]}`, `Qubit.beta_state: {"ms": [m]}` |
| `Qubit.zfs` (with `Qubit.file`, kHz) | a `Qubit.intmap` self entry, `{"between": ["q0", "q0"], "tensor": ..., "unit": "kHz"}` (`q0` the qubit's name, its default) |
| `Qubit.frame` (with `Qubit.file`, development versions) | the `frame` of that self entry |
| `vector_mean_field` (false / true) | `hamiltonian.common.bath_bath_mean_field` (`"z"` / `"xyz"`; `"off"` is new) |
| `hf_readmode` 0 / 1 / 2 / 3 | `hamiltonian.common.hyperfine`: `dipole` `"point"` / `"point"` / `"dft"` / `"dft"`, `contact` false / true / false / true |
| `hf_tensorfile` | `hamiltonian.common.hyperfine.file` |
| `hf_ignore_oor` 0 / 1 | `hamiltonian.common.hyperfine.missing_row` `"stop"` / `"point_dipole"` |
| `qd_readmode` 0 | `hamiltonian.common.quadrupole: "off"` (the default) |
| `qd_readmode` 1 | not available (experimental constants) |
| `qd_readmode` 2 / 3 / 4 | `hamiltonian.common.quadrupole.structure` `"plain"` / `"strained"` / `"bubble"` |
| `qd_tensorfile` | `hamiltonian.common.quadrupole.file` |
| `qd_tensorfile_woqubit` | `hamiltonian.common.quadrupole.file_without_defect` |
| `qd_cellpara` | `hamiltonian.common.quadrupole.cell` |
| `knight` | `hamiltonian.cCCE.knight` |
| `hfmedi` | `hamiltonian.cCCE.hfmedi` |
| (none) | `hamiltonian.cCCE.delta_A` (new: Maze's second-order δA term, default `false` = the old calculation) |
| `DefectTotSpin` | removed: the qubit spin is used (also for the S of a v3 tensor file) |
| `CorrTotSpin` | removed: it served old-format VASP tensor files only; remake such a file with `example/tools/SpinTensorReader` (format v3) |
| `sK` | `cluster.pcce.partition_size` |
| `max_trial` | `cluster.pcce.trials` |
| `max_iter` | `cluster.pcce.max_iterations` |
| `kmeans_pp` | `cluster.pcce.kmeans_pp` |
| `iter_detail`, `pcce_detail_limit` | `cluster.pcce.print_trials` (as below) |
| `cluster.pcce.print_trials` true / false and `cluster.pcce.print_limit` L | `cluster.pcce.print_trials`: the number of trials printed, L (32 without `print_limit`) / 0 |
| `pcce_partition` | `cluster.pcce.partition` |
| `pcce_connections` | `cluster.pcce.connections` |
| `pcce_cache_mode` | `cluster.pcce.cache.mode` |
| `pcce_cache_file` | `cluster.pcce.cache.file` |
| `pcce`, `bundle` (top-level blocks) | `cluster.pcce`, `cluster.bundle` |
| `bundle_max_dim` | `cluster.bundle.max_dim` |
| `bundle_merge_cutoff` | `cluster.bundle.merge_cutoff` |
| `bundle_connections` | `cluster.bundle.connections` |
| `bundle_connection_budget` | `cluster.bundle.connection_budget` |
| `expansion_max_clusters` | `cluster.bundle.max_clusters` |
| `expansion_max_work` | `cluster.bundle.max_work` |
| `npulse` | `pulse.count` |
| `pulsename` | `pulse.name` |
| `sequence` | `pulse.sequence` |
| `pulseiter` | `pulse.per_qubit` |
| `pulse_time` | `pulse.duration` |
| `detuning_factor` | `pulse.detuning_factor` |
| `backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `column_evolution`, `profile` (top level) | `execution.backend`, `execution.gpu_options`, ... (the same names in `execution`) |
| `gpu_min_dim`, `gpu_batch_steps`, `gpu_cache_mb`, `gpu_batch_clusters`, `gpu_cluster_batch_mb` | `execution.gpu_options.min_dim`, `.evolution_batch_steps`, `.cache_mb`, `.evolution_batch_size`, `.cluster_batch_mb` |
| `gcce_refinement` (true/false) | `execution.gcce_refinement: {"enabled": true/false}` |
| `gcce_refinement_max_sweeps` | `execution.gcce_refinement.max_sweeps` |
| `outfile` | `output.file` |
| `savemode`, `save_nodiv` | `output.savemode` |
| `savemode` `"normal"` / `"info"` / `"all"` / `"allfull"` | `"coherence"` / `"coherence_tagged"` / `"cluster_coherence"` / `"cluster_density_matrix"` |
| `"save_nodiv": false` | `"savemode": {"mode": "coherence", "nodiv": false}` |
| `cluster_contributions` (true/false) | `output.cluster_contributions: {"enabled": true/false}` |
| `cluster_contribution_max_mb` | `output.cluster_contributions.max_mb` |
| `cluster_only` (top level) | `output.cluster_only` |
| `coordinate_frame_rotation`, `bath_coordinate_rotation` | `frame.rotation`, on when the block is given (`enabled` defaults to true; the old key defaulted to false, and the converter removes an old block that was off) |
| `defect_axis_reference` | removed (it only named the rotation's qubit axis in the log) |
| `frame.rotation.bath_axis` | removed: it only turned the computational x and y about `qubit_axis`; x is now `qubit_axis` x [0, 0, 1], the old default (the converter notes a value other than [0, 0, 1]) |
| `frame.rotation.reference_qubit` | removed: it only chose the rotation origin (a translation); the origin is the first qubit |
| `frame.rotation.qubit_position_frame` | removed: its one value was `"bath"`; qubit positions are always in the bath file's frame |
| `frame.rotation.hf_tensor_frame` / `qd_tensor_frame` | `hamiltonian.common.hyperfine.frame` / `hamiltonian.common.quadrupole.frame` |
| `Qubit.intmap[].tensor_frame` | `Qubit.intmap[].frame` |
| `Defect.species[].coordinate_frame` | `Defect.species[].frame` |
| `Defect.species[].axis` | `Defect.species[].zfs.axis`: `"zfs": {"D": ..., "axis": [...], "unit": ...}` (only the D form used it; the converter drops it elsewhere) |
| `Defect.species[].zfs.E` (it took 0 only) | removed: the D form has no E, `zfs.tensor` for an E other than 0 |
| a Defect or `Qubit.intmap` unit in another case, e.g. `"mhz"` | the table's spelling, `"MHz"` (units are case-sensitive) |
| the frame value `"bath"` | `"file"` |
| no frame, with the rotation on: a Defect species, a single qubit's `Qubit.intmap` entry, `Qubit.zfs` (read as written before) | refused; the converter writes the old meaning, `"frame": "qubit"` (for `Qubit.zfs` in its intmap self entry) |
| `magnetic_field_axis` | removed: write `bfield` as [0, 0, Bz] and turn the bath with `frame.rotation` |
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` |
| Defect entry `dfname` | `name` |
| Defect entry `apprx` | removed (it changed nothing) |
| Defect entry `naddspin`, `n_nuclear_spins` | `n_onsite_spins` |
| Defect entry `navaax`, `n_configurations` | `n_axis_states` |
| Defect entry `eqs` | `quadrupole_moments` |
| Defect entry `rxyzs` | `nuclear_positions` |
| Defect entry `hypf` | `hyperfine` |
| `auto_convergence.launcher` | removed: the runs use the `mpirun -n N` main.out was started with |
| `auto_convergence.configurations` | `auto_convergence.sampling`; its `bathfile`, `statefile`, `avaaxfile`, `exstatefile` are `configuration.bathfile`, `bath_state_file`, `defect_axis_file`, `defect_state_file` |
| `converge` of `auto_convergence.time_window`, `rdip`, `rbath`, `sampling` | removed: a section whose value (`value`; `sampling.configurations`) is given is fixed, one without it is found |
| `auto_convergence.time_window.unit` | `time_window.value: {"value": ..., "unit": ...}` |
| `auto_convergence.time_window.threshold` | `time_window.decay_level` |
| `auto_convergence.sampling.number`, `test_number` | `sampling.configurations`, `test_configurations` |
| `output.path`, `output.file` with `auto_convergence` (never used) | removed: the runs go to `auto_convergence.workdir` |
| `cluster.adaptive.output` | `cluster.adaptive.workdir` |
| `cluster.adaptive` keys spelled with `-`, e.g. `max-rounds` | `max_rounds` (`-` is the command line's) |
| `cluster.adaptive.policy: "merge"` | not needed (the only policy): removed by the converter, reported as not used by CCEX |
| `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift`, `bspin`, `balphams`, `bbetams` (adding_DEER) | `pulse.bath.count`, `pulse.bath.sequence`, `pulse.bath.defect`, `pulse.bath.energy_shift`, then as below (not rewritten by the converter) |
| `pulse.bath.spin` | removed: the driven spins' own S is used (it had to be the same) |
| `pulse.bath.alphams`, `pulse.bath.betams` | `pulse.bath.levels: [alphams, betams]` |
| `jump_operators` (top level, adding_MECCE) | `hamiltonian.meCCE.jump_operators` or `hamiltonian.megCCE.jump_operators`, the section of the run's method (not rewritten by the converter) |
| `hamiltonian.corrections`, `hamiltonian.dissipation` (development versions only) | `hamiltonian.cCCE`; `hamiltonian.meCCE` or `hamiltonian.megCCE`, the section of the run's method |
| jump `operator` `"+-"` / `"-+"` (adding_MECCE) | `"-"` / `"+"` |
| command-line options `-m`, `-q`, `-B`, `-I`, `-s`, `-a`, `-S`, `-N`, `-P` and the others that set a key | the input keys ([Getting started](GettingStarted.md)); the command line is `-f`, `-o DIR`, `-v`, `-e` and `-h` alone, and `-o` names the output folder, not the output file |
