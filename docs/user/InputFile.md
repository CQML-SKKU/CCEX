# Input File
<a id="InputFile"></a>

Every key of the input file, in eight groups. Each group is a page of its own, and
every key in the lists below links to its entry there. The keys were renamed and grouped into blocks on 2026-10-13
(A47); the [table of old and new names](#OldNames) is at the end of this page, and
`example/tools/InputKeys/convert_input_keys.py` rewrites an old input.

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

A `#` or `!` outside a JSON string starts a comment that runs to the end of the line, so a whole line or the rest of a
line after a value can be a comment. Inside a string (`"#note"`, `"outfile": "run#1"`) they are ordinary characters
(until 2026-10-08 the line was cut there too and such an input was refused).

```
# a comment line
{
    "method" : "cCCE",      # a comment after a value
    ! a comment line of the other kind
    "outfile" : "run#1"
}
```

A key that is left out takes its default. A key that is present has to have the type its
entry gives: `"knight": "true"`, a quoted number, `"nstep": 10.5` or a string inside a
number array stops the run with the key named (with its path in a block, e.g. `pulse.count`). Until 2026-10-02 such a
value was taken as missing and the default was used without a word. A switch is `true` or `false` (the 0 and 1 the old
keys took are not read by the new ones), and a qubit's `detuning` may still be written `[x]`. Keys are read without
regard to case.

A key that no part of the run reads is reported (`Warning: input key ... is not read`), and
when it is one or two letters away from a key that was looked up and is absent from the same
block, the run stops: `"nstates"` would leave `nstate` at its default, `"alphms"` in the `Qubit` block its `alphams`.
A key of the names before 2026-10-13 is refused with its new name (`Error: input key "npulse" is now "pulse.count"`) and
a pointer to the converter.

A whole input in the new keys:

```
{
    "method"       : "cCCE",
    "cluster_mode" : "standard",
    "quantity"     : "coherence",
    "Qubit"        : {"file": "../bath/bath_Diamond_NV_natab_defect", "alphams": 1, "betams": 0},
    "gyrofile"     : "../bath/gyro_13C",
    "bathfile"     : ["../bath/bath_Diamond_NV_natab_1"],
    "order"        : 2,
    "bfield"       : 500,
    "rbath"        : 50,
    "rdip"         : 10,
    "deltat"       : 0.001,
    "nstep"        : 2000,
    "pulse"        : {"count": 1},
    "hamiltonian"  : {"common": {"bath_bath_mean_field": "z"}, "cCCE": {"knight": false}},
    "outfile"      : "./results/CCE2_Diamond_NV_natab_1"
}
```

<a id="GeneralTag"></a>
## General tags
The keys of this group are described in [InputFile_General.md](InputFile_General.md).

The method, the numbers of the calculation, the files of the bath, and the `hamiltonian` block.

- [method](InputFile_General.md#tag_method)       
- [quantity](InputFile_General.md#tag_quantity)     
- [order](InputFile_General.md#tag_order)        
- [bfield](InputFile_General.md#tag_bfield)       
- [rbath](InputFile_General.md#tag_rbath) (a number or `[inner, outer]`)
- [rdip](InputFile_General.md#tag_rdip) (a number or `[inner, outer]`)
- [deltat](InputFile_General.md#tag_deltat)       
- [nstep](InputFile_General.md#tag_nstep)        
- [rdip_ladder](InputFile_General.md#tag_rdip_ladder)
- [nstate](InputFile_General.md#tag_nstate)        
- [seed](InputFile_General.md#tag_seed)       
- [gyrofile](InputFile_General.md#tag_gyrofile)      
- [bathfile](InputFile_General.md#tag_bathfile)      
- [bath_shift](InputFile_General.md#tag_bath_shift)
- [bath_state_file](InputFile_General.md#tag_bath_state_file)
- [coordinate_frame_rotation](InputFile_General.md#tag_coordinate_frame_rotation)
- [defect_axis_reference](InputFile_General.md#tag_defect_axis_reference)
- [hamiltonian](InputFile_General.md#tag_hamiltonian)
    - common: [qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field),
      [bath_bath_mean_field](InputFile_General.md#tag_bath_bath_mean_field),
      [hyperfine](InputFile_General.md#tag_hyperfine) (`dipole`, `contact`, [file](InputFile_General.md#tag_hyperfine_file),
      [missing_row](InputFile_General.md#tag_hyperfine_missing_row)),
      [quadrupole](InputFile_General.md#tag_quadrupole) (`structure`, [file](InputFile_General.md#tag_quadrupole_file),
      [file_without_defect](InputFile_General.md#tag_quadrupole_file_without_defect), [cell](InputFile_General.md#tag_quadrupole_cell))
    - cCCE: [knight](InputFile_General.md#tag_knight), [hfmedi](InputFile_General.md#tag_hfmedi),
      [delta_A](InputFile_General.md#tag_delta_A) (meCCE reads them too)
    - meCCE, megCCE: [jump_operators](InputFile_General.md#tag_jump_operators)
      ([master-equation.md](../methods/master-equation.md))
    - gCCE, noise: no settings yet (noise is a reserved method)

<a id="QubitArrayTags"></a>
## Qubit-related tags
The keys of this group are described in [InputFile_Qubit.md](InputFile_Qubit.md).

The `Qubit` block, in one of two forms: <br>
<b>Form 1.</b> [file](InputFile_Qubit.md#tag_Qubit_file): one qubit from a file, with `spin`, `zfs`, `alphams`,
`betams` beside it. <b>Example:</b> ``` "Qubit" : {"file": "file-path", "alphams": 1, "betams": 0} ``` <br>
<b>Form 2.</b> [qubit](InputFile_Qubit.md#tag_qubit): one or more qubits written out, with `intmap`, `alpha_state`,
`beta_state` and `psi0` beside the list. The two forms in one block are refused. The full example is on the
[Qubit page](InputFile_Qubit.md).

- [Qubit](InputFile_Qubit.md#tag_Qubit) 
- [file](InputFile_Qubit.md#tag_Qubit_file), [spin](InputFile_Qubit.md#tag_Qubit_spin), [zfs](InputFile_Qubit.md#tag_Qubit_zfs),
  [alphams](InputFile_Qubit.md#tag_Qubit_alphams), [betams](InputFile_Qubit.md#tag_Qubit_betams) (form 1)
- [qubit](InputFile_Qubit.md#tag_qubit) (form 2)
    - [name](InputFile_Qubit.md#tag_name), [spin](InputFile_Qubit.md#tag_spin), [gyro](InputFile_Qubit.md#tag_gyro),
      [xyz](InputFile_Qubit.md#tag_xyz), [detuning](InputFile_Qubit.md#tag_detuning), [alphams](InputFile_Qubit.md#tag_alphams),
      [betams](InputFile_Qubit.md#tag_betams)
- [intmap](InputFile_Qubit.md#tag_intmap): [between](InputFile_Qubit.md#tag_between), [tensor](InputFile_Qubit.md#tag_tensor),
  [unit](InputFile_Qubit.md#tag_intmap_unit), [tensor_frame](InputFile_Qubit.md#tag_tensor_frame)
- [alpha_state](InputFile_Qubit.md#tag_alpha_state), [beta_state](InputFile_Qubit.md#tag_beta_state)
- [psi0](InputFile_Qubit.md#tag_psi0)
- the Overhauser field: [hamiltonian.common.qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field)

<a id="ClusterTags"></a>
## Cluster-related tags
The keys of this group are described in [InputFile_Cluster.md](InputFile_Cluster.md).

- [order](InputFile_General.md#tag_order)        
- [method](InputFile_General.md#tag_method)       
- [cluster_mode](InputFile_Cluster.md#tag_cluster_mode) (`standard`, `pcce`, `bundle`, `dsj`, `itb`, `dsjitb`)
- [add_subclusters](InputFile_Cluster.md#tag_add_subclusters)   
- [max_clusters_per_order](InputFile_Cluster.md#tag_max_clusters_per_order)           
- [pcce](InputFile_Cluster.md#tag_pcce): [partition_size](InputFile_Cluster.md#tag_pcce_partition_size),
  [trials](InputFile_Cluster.md#tag_pcce_trials), [max_iterations](InputFile_Cluster.md#tag_pcce_max_iterations),
  [kmeans_pp](InputFile_Cluster.md#tag_pcce_kmeans_pp), [print_trials](InputFile_Cluster.md#tag_pcce_print_trials),
  [print_limit](InputFile_Cluster.md#tag_pcce_print_limit), [partition](InputFile_Cluster.md#tag_pcce_partition),
  [connections](InputFile_Cluster.md#tag_pcce_connections), [cache](InputFile_Cluster.md#tag_pcce_cache)
- [bundle](InputFile_Cluster.md#tag_bundle): [max_dim](InputFile_Cluster.md#tag_bundle_max_dim),
  [merge_cutoff](InputFile_Cluster.md#tag_bundle_merge_cutoff), [connections](InputFile_Cluster.md#tag_bundle_connections),
  [connection_budget](InputFile_Cluster.md#tag_bundle_connection_budget), [max_clusters](InputFile_Cluster.md#tag_bundle_max_clusters),
  [max_work](InputFile_Cluster.md#tag_bundle_max_work)
- [adaptive](InputFile_Cluster.md#tag_adaptive) ([work_file](InputFile_Cluster.md#tag_adaptive_work_file), internal)

<a id="BackendTags"></a>
## Backend tags
The keys of this group are described in [InputFile_Backend.md](InputFile_Backend.md).

Where and how the cluster dynamics is computed, and what the run reports about it.

- [backend](InputFile_Backend.md#tag_backend)
- [gpu_options](InputFile_Backend.md#tag_gpu_options) (min_dim, packed_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb, warp_max_dim)
- [probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb)
- [propagator](InputFile_Backend.md#tag_propagator)
- [evolution](InputFile_Backend.md#tag_evolution)
- [evolution_update](InputFile_Backend.md#tag_evolution_update)
- [gcce_refinement](InputFile_Backend.md#tag_gcce_refinement) (`enabled`, [max_sweeps](InputFile_Backend.md#tag_gcce_refinement_max_sweeps))
- [profile](InputFile_Backend.md#tag_profile)
- [cluster_contributions](InputFile_Backend.md#tag_cluster_contributions) (`enabled`, [max_mb](InputFile_Backend.md#tag_cluster_contribution_max_mb))

<a id="PulseTags"></a>
## Pulse-related tags
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
        - [spin](InputFile_Pulse.md#tag_pulse_bath_spin)
        - [alphams](InputFile_Pulse.md#tag_pulse_bath_alphams)
        - [betams](InputFile_Pulse.md#tag_pulse_bath_betams)

<a id="DefectArrayTags"></a>
## DefectArray-related tags
The keys of this group are described in [InputFile_Defect.md](InputFile_Defect.md).

The `Defect` block of a paramagnetic defect bath: `species`, the list of the defect species, and beside it the files
`axis_file` and `state_file`. The full example is on the [Defect page](InputFile_Defect.md).

- [Defect](InputFile_Defect.md#tag_Defect)       
- [axis_file](InputFile_Defect.md#tag_df_axis_file)
- [state_file](InputFile_Defect.md#tag_df_state_file)
- [species](InputFile_Defect.md#tag_df_species)
    - [name](InputFile_Defect.md#tag_df_name)       
    - [n_nuclear_spins](InputFile_Defect.md#tag_df_n_nuclear_spins)     
    - [n_configurations](InputFile_Defect.md#tag_df_n_configurations)       
    - [axis](InputFile_Defect.md#tag_df_axis)
    - [coordinate_frame](InputFile_Defect.md#tag_df_coordinate_frame)
    - [types](InputFile_Defect.md#tag_df_types)        
    - [spins](InputFile_Defect.md#tag_df_spins)        
    - [gyros](InputFile_Defect.md#tag_df_gyros)        
    - [quadrupole_moments](InputFile_Defect.md#tag_df_quadrupole_moments)          
    - [nuclear_positions](InputFile_Defect.md#tag_df_nuclear_positions)        
    - [hyperfine](InputFile_Defect.md#tag_df_hyperfine)         
    - [efg](InputFile_Defect.md#tag_df_efg)          
    - [zfs](InputFile_Defect.md#tag_df_zfs)          
    - [detuning](InputFile_Defect.md#tag_df_detuning)

<a id="OutputTags"></a>
## Output tags
The keys of this group are described in [InputFile_Output.md](InputFile_Output.md).

- [savemode](InputFile_Output.md#tag_savemode) (`mode`, `nodiv`; the old names and `save_nodiv` are refused)
- [outfile](InputFile_Output.md#tag_outfile)

<a id="AutoConvergenceTags"></a>
## Automatic convergence tags
The `auto_convergence` block (its sections `time_window`, `rdip`, `rbath`, `t2_fit` and
`configurations`) is described in [InputFile_AutoConvergence.md](InputFile_AutoConvergence.md).

<a id="OldNames"></a>
## Old and new names (2026-10-13)

The input keys were renamed and grouped on 2026-10-13 (A47; the decisions are in
(internal notes)). The results do not change: rewritten in the new keys, the inputs of the standing checks
(check_numerics, check_input_forms, check_rotation) give the same output bytes as before. An old name, value or form is refused with its new one; the converter rewrites old inputs,
keeping their comments, and lists what it changed:

```
python3 example/tools/InputKeys/convert_input_keys.py FILE...           # rewrite in place
python3 example/tools/InputKeys/convert_input_keys.py --check FILE...   # list only; exit 1 when a file has old keys
python3 example/tools/InputKeys/convert_input_keys.py --stdout FILE     # print the new text
```

It drops keys the run did not read (e.g. `qd_tensorfile` with `qd_readmode` 0) and the defaults `hf_readmode` 0 and
`qd_readmode` 0, and refuses an input it cannot express (a `CorrTotSpin` other than 0). Keys not in the table kept their
names (`quantity`, `order`, `bfield`, `deltat`, `nstep`, `nstate`, `seed`, `rdip_ladder`, `gyrofile`, `bathfile`,
`outfile`, `savemode`, `backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `probe_memory_gb`,
`profile`, `coordinate_frame_rotation`, `defect_axis_reference`, `auto_convergence`, the keys of `Qubit.qubit` and
`Qubit.intmap` entries, and the Defect species' `types`, `spins`, `gyros`, `efg`, `zfs`, `detuning`, `axis`,
`coordinate_frame`). The command-line options are the same; `-m` takes `cCCE` or `gCCE` and keeps the input's
`cluster_mode`.

| Before 2026-10-13 | Now |
|---|---|
| `"method": "CCE"` | `"method": "cCCE"` |
| `"method": "pCCE"` | `"method": "cCCE", "cluster_mode": "pcce"` |
| `"method": "dsj"` / `"itb"` / `"dsjitb"` | `"method": "cCCE", "cluster_mode": "dsj"` / `"itb"` / `"dsjitb"` |
| `"cluster_mode": "adaptive"` | `"cluster_mode": "bundle"` |
| `Adaptive` block | `adaptive` block (keys unchanged; `cluster_mode` pcce or bundle) |
| `adaptive_work_file` | `adaptive.work_file` (written by the refinement driver) |
| `rbathcut` | `rbath: [inner, outer]` (the inner radius) |
| `rdipcut` | `rdip: [inner, outer]` (the inner radius) |
| `bathadjust` | `bath_shift` |
| `statefile` | `bath_state_file` |
| `exstatefile` | `Defect.state_file` |
| `avaaxfile` | `Defect.axis_file` |
| `addsubclus` (0/1 or true/false) | `add_subclusters` (true/false) |
| `nk` | `max_clusters_per_order` |
| `qubitfile` | `Qubit.file` |
| `qspin`, `qzfs`, `qalphams`, `qbetams` | `Qubit.spin`, `Qubit.zfs`, `Qubit.alphams`, `Qubit.betams` (with `Qubit.file`) |
| `Qubit.nqubit` | removed: the length of `Qubit.qubit` |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| `Qubit.alphaidx`, `Qubit.betaidx` | `Qubit.alpha_state: {"index": n}`, `Qubit.beta_state: {"index": n}` |
| `Qubit.psia`, `Qubit.psib` | `Qubit.alpha_state: {"vector": [...]}`, `Qubit.beta_state: {"vector": [...]}` |
| `vector_mean_field` (false / true) | `hamiltonian.common.bath_bath_mean_field` (`"z"` / `"xyz"`; `"off"` is new) |
| `hf_readmode` 0 / 1 / 2 / 3 | `hamiltonian.common.hyperfine`: `dipole` `"point"` / `"point"` / `"dft"` / `"dft"`, `contact` false / true / false / true |
| `hf_tensorfile` | `hamiltonian.common.hyperfine.file` |
| `hf_ignore_oor` 0 / 1 | `hamiltonian.common.hyperfine.missing_row` `"stop"` / `"point_dipole"` |
| `qd_readmode` 0 | `hamiltonian.common.quadrupole: "off"` (the default) |
| `qd_readmode` 2 / 3 / 4 | `hamiltonian.common.quadrupole.structure` `"plain"` / `"strained"` / `"bubble"` |
| `qd_tensorfile` | `hamiltonian.common.quadrupole.file` |
| `qd_tensorfile_woqubit` | `hamiltonian.common.quadrupole.file_without_defect` |
| `qd_cellpara` | `hamiltonian.common.quadrupole.cell` |
| `knight` | `hamiltonian.cCCE.knight` |
| `hfmedi` | `hamiltonian.cCCE.hfmedi` |
| (none) | `hamiltonian.cCCE.delta_A` (new; `true` is the calculation before) |
| `DefectTotSpin` | removed: the qubit spin is used (also for the S of a v3 tensor file) |
| `CorrTotSpin` | removed: it served old-format VASP tensor files only; remake such a file with `example/tools/SpinTensorReader` (format v3) |
| `sK` | `pcce.partition_size` |
| `max_trial` | `pcce.trials` |
| `max_iter` | `pcce.max_iterations` |
| `kmeans_pp` | `pcce.kmeans_pp` |
| `iter_detail` | `pcce.print_trials` |
| `pcce_detail_limit` | `pcce.print_limit` |
| `pcce_partition` | `pcce.partition` |
| `pcce_connections` | `pcce.connections` |
| `pcce_cache_mode` | `pcce.cache.mode` |
| `pcce_cache_file` | `pcce.cache.file` |
| `bundle_max_dim` | `bundle.max_dim` |
| `bundle_merge_cutoff` | `bundle.merge_cutoff` |
| `bundle_connections` | `bundle.connections` |
| `bundle_connection_budget` | `bundle.connection_budget` |
| `expansion_max_clusters` | `bundle.max_clusters` |
| `expansion_max_work` | `bundle.max_work` |
| `npulse` | `pulse.count` |
| `pulsename` | `pulse.name` |
| `sequence` | `pulse.sequence` |
| `pulseiter` | `pulse.per_qubit` |
| `pulse_time` | `pulse.duration` |
| `detuning_factor` | `pulse.detuning_factor` |
| `gcce_refinement` (true/false) | `gcce_refinement: {"enabled": true/false}` |
| `gcce_refinement_max_sweeps` | `gcce_refinement.max_sweeps` |
| `cluster_contributions` (true/false) | `cluster_contributions: {"enabled": true/false}` |
| `cluster_contribution_max_mb` | `cluster_contributions.max_mb` |
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` (with `state_file` and `axis_file` beside it) |
| Defect entry `dfname` | `name` |
| Defect entry `apprx` | removed (it changed nothing) |
| Defect entry `naddspin` | `n_nuclear_spins` |
| Defect entry `navaax` | `n_configurations` |
| Defect entry `eqs` | `quadrupole_moments` |
| Defect entry `rxyzs` | `nuclear_positions` |
| Defect entry `hypf` | `hyperfine` |
| `-m CCE`, `-m pCCE`, ... | `-m cCCE`, `-m gCCE`, `-m meCCE` or `-m megCCE`; the cluster shape is the input's `cluster_mode` |

The DEER and master-equation keys of the adding_DEER / adding_MECCE branches were never keys of this code; they are
refused with their new names (`Error: input key "bnpulse" is now "pulse.bath.count"`), and the converter does not
rewrite them:

| adding_DEER / adding_MECCE | Now |
|---|---|
| `bnpulse` | `pulse.bath.count` |
| `bsequence` | `pulse.bath.sequence` |
| `bpulse_defect` | `pulse.bath.defect` |
| `bpulse_energy_shift` | `pulse.bath.energy_shift` |
| `bspin` | `pulse.bath.spin` |
| `balphams` | `pulse.bath.alphams` |
| `bbetams` | `pulse.bath.betams` |
| `jump_operators` (top level) | `hamiltonian.meCCE.jump_operators` or `hamiltonian.megCCE.jump_operators`, the section of the run's method |
| jump `operator` `"+-"` / `"-+"` | `"-"` / `"+"` |
| `"method": "pgCCE"` / `"pmeCCE"` / `"pmegCCE"` (also `-m`) | `"method": "gCCE"` / `"meCCE"` / `"megCCE"` with `"cluster_mode": "pcce"` |

The names of savemode changed the same day and are in [savemode](InputFile_Output.md#tag_savemode).
