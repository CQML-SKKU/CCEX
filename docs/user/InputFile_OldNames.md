# Old names
<a id="OldNames"></a>

Inputs written for an earlier version (v1.2.2 and before, or the adding_DEER and adding_MECCE branches) use other
names and places for many keys. CCEX refuses an old name, value, form or place and names the current one, for example
`Error: input key "npulse" is now "pulse.count"`. The table below lists them all.

1. [Converting an old input](#OldConvert): the converter.
2. [The table](#OldTable): every old name and what it is now.

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
whose new name is also in the input (`npulse` beside `pulse.count`) is kept once when the two values agree; when they
differ, the file is not written. A few rows of the table say "not rewritten by the converter"; change those by hand.

<a id="OldTable"></a>
## 2. The table

These keys kept their names and places and are not in the table: `order`, `rbath`, `rdip`, `quantity`, `bfield`,
`deltat`, `nstep` and `rdip_ladder` at the top level; the keys of the `Qubit.qubit` and `Qubit.intmap` entries; and the
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
| command-line options `-m`, `-q`, `-B`, `-I`, `-s`, `-a`, `-S`, `-N`, `-P` and the others that set a key | the input keys; the command line is `-f`, `-o DIR`, `-v`, `-e` and `-h` alone, and `-o` names the output folder, not the output file |
