# Upgrading old inputs
<a id="OldNames"></a>

Inputs written for an earlier version (v1.2.2 and before, or the development branches adding_DEER and adding_MECCE) use other
names and places for many keys. CCEX refuses an old name, value, form or place and names the current one, for example
`Error: input key "npulse" is now "pulse.count"`. The tables below list them all, by block.

1. [Converting an old input](#OldConvert): the converter.
2. [The old names](#OldTable): every old name and what it is now, by block.

<a id="OldConvert"></a>
## 1. Converting an old input

`example/tools/InputKeys/convert_input_keys.py` rewrites an old input, keeps its comments and lists what it changed:

```sh
python3 example/tools/InputKeys/convert_input_keys.py FILE...           # rewrite in place
python3 example/tools/InputKeys/convert_input_keys.py --check FILE...   # list only; exit 1 when a file has old keys
python3 example/tools/InputKeys/convert_input_keys.py --stdout FILE     # print the new text
```

It also drops keys the old run did not read (for example `qd_tensorfile` with `qd_readmode` 0) and the defaults
`hf_readmode` 0 and `qd_readmode` 0. It refuses an input it cannot express (a `CorrTotSpin` other than 0). An old key
whose new name is also in the input (`npulse` beside `pulse.count`) is kept once when the two values agree. When they
differ, the file is not written. A few rows of the table say "not rewritten by the converter". Change those by hand.

<a id="OldTable"></a>
## 2. The old names

These keys kept their names and places and are not in the table: `order`, `rbath`, `rdip`, `quantity`, `bfield`,
`deltat`, `nstep` and `rdip_ladder` at the top level, the keys of the `Qubit.qubit` and `Qubit.intmap` entries, and the
Defect species' `types`, `spins`, `gyros`, `efg`, `zfs` and `detuning`.

<a id="old_method"></a>
### method, cluster and the radii

| Old | Now |
|---|---|
| `"method": "CCE"` | `"method": "cCCE"` |
| `"method": "pCCE"` | `"method": "cCCE", "cluster": {"mode": "pcce"}` |
| `"method": "dsj"` / `"itb"` / `"dsjitb"` | `"method": "cCCE", "cluster": {"mode": "dsj"}` / `"itb"` / `"dsjitb"` |
| `"method": "pgCCE"` / `"pmeCCE"` / `"pmegCCE"` | `"method": "gCCE"` / `"meCCE"` / `"megCCE"` with `"cluster": {"mode": "pcce"}` |
| `cluster_mode` | `cluster.mode` |
| `"cluster_mode": "adaptive"` | `"cluster": {"mode": "bundle"}` |
| `Adaptive` block, top-level `adaptive` block | `cluster.adaptive` block (keys unchanged, `cluster.mode` bundle). Keys are read without regard to case, so `Adaptive` inside `cluster` is read. |
| `Refine` block | `cluster.adaptive` block |
| `adaptive_work_file` | `cluster.adaptive.work_file` (written by the refinement driver) |
| `probe_memory_gb` | `cluster.adaptive.probe_memory_gb` (the converter removes it from an input without an `adaptive` block, since it was read with the refinement only) |
| `rbathcut` | `rbath: [inner, outer]` (the inner radius) |
| `rdipcut` | `rdip: [inner, outer]` (the inner radius) |
| `addsubclus` (0/1 or true/false), `add_subclusters` (top level) | `cluster.selection.add_subclusters` (true/false) |
| `nk`, `max_clusters_per_order` | `cluster.selection.max_per_order` |
| `sK` | `cluster.pcce.partition_size` |
| `max_trial` | `cluster.pcce.trials` |
| `max_iter` | `cluster.pcce.max_iterations` |
| `kmeans_pp` | `cluster.pcce.kmeans_pp` |
| `iter_detail`, `pcce_detail_limit` | `cluster.pcce.print_trials` (as below) |
| `cluster.pcce.print_trials` true / false and `cluster.pcce.print_limit` L | `cluster.pcce.print_trials`, the number of trials printed: L (32 without `print_limit`) / 0 |
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

<a id="old_config"></a>
### configuration

| Old | Now |
|---|---|
| `gyrofile`, `bathfile`, `nstate`, `seed` (top level) | `configuration.gyrofile`, `configuration.bathfile`, `configuration.s_range` (as `configuration.nstate` below), `configuration.seed` |
| `configuration.nstate` N | `configuration.s_range: [1, N]`. For `nstate` 0 (the ensemble), no `s_range`. |
| a `configuration.seed` number with `k_range` | The number is now the seed of the first configuration, and configuration k takes seed + k − first. For the old meaning, the same seed for every k, give a list of that number. The converter writes it. |
| `bathadjust`, `bath_shift` (top level) | `configuration.bath_shift` |
| `statefile`, `bath_state_file: "prefix"` (top level) | `configuration.bath_state_file: ["prefix"]`, one name per bath file. The converter removes a single prefix beside several bath files, which the old run did not read. |
| `avaaxfile`, `Defect.axis_file: "file"` | `configuration.defect_axis_file: ["file"]`, one per bath file |
| `exstatefile`, `Defect.state_file: "prefix"` | `configuration.defect_state_file: ["prefix"]`, one per bath file |

<a id="old_qubit"></a>
### Qubit

| Old | Now |
|---|---|
| `qubitfile` | `Qubit.qubit: [{"file": ...}]` (as `Qubit.file` below) |
| `qspin`, `qzfs`, `qalphams`, `qbetams` | the entry's `spin`, then as `Qubit.zfs`, `Qubit.alphams`, `Qubit.betams` below |
| `Qubit.file`, `Qubit.spin`, `Qubit.name` (the one-qubit file form) | an entry of `Qubit.qubit`, `[{"file": ..., "spin": ..., "name": ...}]` (`file` in place of `xyz`) |
| `Qubit.qubit[].detuning: [x]` | `x` (a number, or `{"value", "unit"}`) |
| `Qubit.nqubit` | Removed. The number of qubits is the length of `Qubit.qubit`. |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| `Qubit.alphaidx`, `Qubit.betaidx` | `Qubit.alpha_state: {"index": n}`, `Qubit.beta_state: {"index": n}` |
| `Qubit.psia`, `Qubit.psib` | `Qubit.alpha_state: {"vector": [...]}`, `Qubit.beta_state: {"vector": [...]}`. |
| `Qubit.qubit[].alphams`, `Qubit.qubit[].betams` | `Qubit.alpha_state: {"ms": [m_1, m_2, ...]}`, `Qubit.beta_state: {"ms": [...]}`, one m<sub>s</sub> per qubit. The converter writes +1 / 0 for a qubit that gave none. |
| `Qubit.alphams`, `Qubit.betams` (with `Qubit.file`) | `Qubit.alpha_state: {"ms": [m]}`, `Qubit.beta_state: {"ms": [m]}` |
| `Qubit.zfs` (with `Qubit.file`, kHz) | a `Qubit.intmap` self entry, `{"between": ["q0", "q0"], "tensor": ..., "unit": "kHz"}` (`q0` the qubit's name, its default) |
| `Qubit.frame` (with `Qubit.file`, development versions) | the `frame` of that self entry |
| `Qubit.intmap[].tensor_frame` | `Qubit.intmap[].frame` |

<a id="old_ham"></a>
### hamiltonian

| Old | Now |
|---|---|
| `vector_mean_field` (false / true) | `hamiltonian.common.bath_bath_mean_field` (`"z"` / `"xyz"`, and `"off"` is new) |
| `hf_readmode` 0 / 1 / 2 / 3 | `hamiltonian.common.hyperfine` with `dipole` `"point"` / `"point"` / `"dft"` / `"dft"`, `contact` false / true / false / true |
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
| `DefectTotSpin` | Removed. The qubit spin is used (also for the S of a v3 tensor file). |
| `CorrTotSpin` | Removed. It served old-format VASP tensor files only. Remake such a file with `example/tools/SpinTensorReader` (format v3). |
| `jump_operators` (top level) | `hamiltonian.meCCE.jump_operators` or `hamiltonian.megCCE.jump_operators`, the section of the run's method (not rewritten by the converter) |
| `hamiltonian.corrections`, `hamiltonian.dissipation` (development versions only) | `hamiltonian.cCCE`, or `hamiltonian.meCCE` / `hamiltonian.megCCE` (the section of the run's method) |
| jump `operator` `"+-"` / `"-+"` | `"-"` / `"+"` |

<a id="old_pulse"></a>
### pulse

| Old | Now |
|---|---|
| `npulse` | `pulse.count` |
| `pulsename` | `pulse.name` |
| `sequence` | `pulse.sequence` |
| `pulseiter` | `pulse.per_qubit` |
| `pulse_time` | `pulse.duration` |
| `detuning_factor` | `pulse.detuning_factor` |
| `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift`, `bspin`, `balphams`, `bbetams` | `pulse.bath.count`, `pulse.bath.sequence`, `pulse.bath.defect`, `pulse.bath.energy_shift`, then as below (not rewritten by the converter) |
| `pulse.bath.spin` | Removed. The driven spins' own S is used (it had to be the same). |
| `pulse.bath.alphams`, `pulse.bath.betams` | `pulse.bath.levels: [alphams, betams]` |

<a id="old_exec"></a>
### execution

| Old | Now |
|---|---|
| `backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `column_evolution`, `profile` (top level) | `execution.backend`, `execution.gpu_options`, ... (the same names in `execution`) |
| `gpu_min_dim`, `gpu_batch_steps`, `gpu_cache_mb`, `gpu_batch_clusters`, `gpu_cluster_batch_mb` | `execution.gpu_options.min_dim`, `.evolution_batch_steps`, `.cache_mb`, `.evolution_batch_size`, `.cluster_batch_mb` |
| `gcce_refinement` (true/false) | `execution.gcce_refinement: {"enabled": true/false}` |
| `gcce_refinement_max_sweeps` | `execution.gcce_refinement.max_sweeps` |

<a id="old_output"></a>
### output

| Old | Now |
|---|---|
| `outfile` | `output.file` |
| `savemode`, `save_nodiv` | `output.savemode` |
| `savemode` `"normal"` / `"info"` / `"all"` / `"allfull"` | `"coherence"` / `"coherence_tagged"` / `"cluster_coherence"` / `"cluster_density_matrix"` |
| `"save_nodiv": false` | `"savemode": {"mode": "coherence", "nodiv": false}` |
| `cluster_contributions` (true/false) | `output.cluster_contributions: {"enabled": true/false}` |
| `cluster_contribution_max_mb` | `output.cluster_contributions.max_mb` |
| `cluster_only` (top level) | `output.cluster_only` |

<a id="old_frame"></a>
### frame and units

| Old | Now |
|---|---|
| `coordinate_frame_rotation`, `bath_coordinate_rotation` | `frame.rotation`, on when the block is given (`enabled` defaults to true). The old key defaulted to false, and the converter removes an old block that was off. |
| `defect_axis_reference` | Removed (it only named the rotation's qubit axis in the log). |
| `frame.rotation.bath_axis` | Removed. The computational x is now `qubit_axis` × [0, 0, 1], the old default. An input with another `bath_axis` gets this x axis, and with gCCE and a qubit ZFS that is not symmetric about `qubit_axis` its result can change. |
| `frame.rotation.reference_qubit` | Removed. It only chose the rotation origin (a translation). The origin is the first qubit. |
| `frame.rotation.qubit_position_frame` | Removed. Its one value was `"bath"`. Qubit positions are always in the bath file's frame. |
| `frame.rotation.hf_tensor_frame` / `qd_tensor_frame` | `hamiltonian.common.hyperfine.frame` / `hamiltonian.common.quadrupole.frame` |
| `Defect.species[].coordinate_frame` | `Defect.species[].frame` |
| `Defect.species[].axis` | `Defect.species[].zfs.axis`, as in `"zfs": {"D": ..., "axis": [...], "unit": ...}`. Only the D form used it, and the converter drops it elsewhere. |
| `Defect.species[].zfs.E` (it took 0 only) | Removed. The D form has no E. Use `zfs.tensor` for an E other than 0. |
| a Defect or `Qubit.intmap` unit in another case, e.g. `"mhz"` | the table's spelling, `"MHz"` (units are case-sensitive) |
| the frame value `"bath"` | `"file"` |
| no frame with the rotation on, for a Defect species, a single qubit's `Qubit.intmap` entry, `Qubit.zfs` (read as written before) | Refused. The converter writes the old meaning, `"frame": "qubit"` (for `Qubit.zfs` in its intmap self entry). |
| `magnetic_field_axis` | Removed. Write `bfield` as [0, 0, Bz] and turn the bath with `frame.rotation`. |

<a id="old_defect"></a>
### Defect

| Old | Now |
|---|---|
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` |
| Defect entry `dfname` | `name` |
| Defect entry `apprx` | Removed (it changed nothing). |
| Defect entry `naddspin`, `n_nuclear_spins` | `n_onsite_spins` |
| Defect entry `navaax`, `n_configurations` | `n_axis_states` |
| Defect entry `eqs` | `quadrupole_moments` |
| Defect entry `rxyzs` | `nuclear_positions` |
| Defect entry `hypf` | `hyperfine` |

<a id="old_auto"></a>
### auto_convergence and cluster.adaptive

| Old | Now |
|---|---|
| `auto_convergence.launcher` | Removed. The runs use the `mpirun -n N` main.out was started with. |
| `auto_convergence.configurations` | `auto_convergence.sampling`. Its file names are now `configuration.bathfile`, `bath_state_file`, `defect_axis_file` and `defect_state_file`. |
| `converge` of `auto_convergence.time_window`, `rdip`, `rbath`, `sampling` | Removed. A section whose value (`value`, or `sampling.configurations`) is given is fixed, one without it is found. |
| `auto_convergence.time_window.unit` | `time_window.value: {"value": ..., "unit": ...}` |
| `auto_convergence.time_window.threshold` | `time_window.decay_level` |
| `auto_convergence.sampling.number`, `test_number` | `sampling.configurations`, `test_configurations` |
| `output.path`, `output.file` with `auto_convergence` (never used) | Removed. The runs go to `auto_convergence.workdir`. |
| `cluster.adaptive.output` | `cluster.adaptive.workdir` |
| `cluster.adaptive` keys spelled with `-`, e.g. `max-rounds` | `max_rounds` (`-` is the command line's) |
| `cluster.adaptive.policy: "merge"` | Not needed, since it is the only policy. The converter removes it, and CCEX reports it as not used. |

<a id="old_cli"></a>
### command-line options

| Old | Now |
|---|---|
| command-line options `-m`, `-q`, `-B`, `-I`, `-s`, `-a`, `-S`, `-N`, `-P` and the others that set a key | The input keys. The command line is `-f`, `-o DIR`, `-v`, `-e` and `-h` alone, and `-o` names the output folder, not the output file. |
| `--auto_conv` | An `auto_convergence` block in the input. |
| `--auto_status` | `python3 -m autoconv status FILE` ([AutoConvergence](Tools.md#tag_tool_autoconvergence)). |
| `--adaptive` | A `cluster.adaptive` block in the input. |
