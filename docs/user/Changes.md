# Changes

What changed in each version, in the current key names. Every change, with its details and measurements, is in the
file `CHANGELOG.md` at the root of the CCEX repository.

## 2.1.0

Version 2.1.0, released 2026-11-03: what changed since 2.0.1 (2.0.1, 2026-10-22, fixed GPU gCCE and the builds of 2.0.0).

- **Two new methods.** `noise` gives the Gaussian coherence exp(−χ) of the qubit's phase noise, from the modes of the
  bath's cCCE clusters or from a model spectrum (Ornstein–Uhlenbeck or a power law), and with
  `hamiltonian.noise.compute` its correlation function C(t) and spectrum S(f); ported from semi-cce, its mode sums
  also run on the GPU ([method](InputFile_General.md#tag_method), tutorial
  [11_noise](../../example/tutorial/11_noise/README.md)). `pseudospin` gives the bath's pairs alone in the pseudospin
  model of nuclear flip-flops, for any nuclear spin, with a table of every channel's Ω and Δ
  ([channels](InputFile_General.md#tag_pseudospin_channels), tutorial
  [12_pseudospin](../../example/tutorial/12_pseudospin/README.md)).
- **GPU runs on a full node.** `execution.gpu_options.ranks_per_gpu` lets only some ranks of a node use the GPU while
  the others stay on the CPU, the clusters split between them by measured cost; a GPU batch of cCCE clusters is
  diagonalized in one call; for cCCE with `backend` `auto` or `cuda` the `gpu_options` an input leaves out are chosen
  at the start (a `[gpu-auto]` line in the log) ([Execution tags](InputFile_Backend.md#tag_backend)).
- **Builds.** An nvcc older than 13.2 is accepted for cards below compute capability 12.0 (CUDA 12.8 with MPS on a
  GH200).
- **Fixed.** The MKL builds write the same last digits from run to run.
- **Inputs.** No key of 2.0.1 was renamed or removed: a 2.0.1 input runs as it is. In `bin/main_gpu.out` with
  `backend` `auto` the automatic `gpu_options` change which clusters go to the card, so the outputs move within
  rounding.

## 2.0.0

Version 2.0.0, released 2026-10-17: what changed since v1.2.2, in the current key names.

### Summary

- **A new input format.** The keys are grouped in blocks (`configuration`, `Qubit`, `Defect`, `hamiltonian`, `pulse`,
  `cluster`, `execution`, `output`, `frame`, `auto_convergence`); inputs may be `.jsonc` files with `//` comments; many
  numbers take a unit (`{"value": 5, "unit": "us"}`); `main.out -e` writes a commented template with every key. The
  command line is only `-f`, `-o DIR`, `-v`, `-e`, `-h`. Old inputs are refused with the new name, and a converter
  rewrites them (see "Updating old inputs").
- **Faster.** The test inputs run up to about 35 times faster than v1.2.2 (P1 cCCE2 production: 442 s → 22 s on 64
  ranks), mostly with unchanged results; MKL is used for large clusters when the build finds it; a CUDA backend runs cCCE, gCCE
  and pCCE on a GPU.
- **New methods and options.** ME-CCE and ME-gCCE (Lindblad dissipation), DEER pulses on bath spins, pCCE with gCCE and
  the master-equation methods, bundle clustering with an adaptive refinement, the automatic convergence of rbath, rdip,
  the time window and the number of configurations, several bath configurations from one input.
- **Physics fixes that change results**: `hfmedi` (complete second-order term), `knight` (nuclear Zeeman vector, $|m_s|$
  of each branch), the static Overhauser field of gCCE with sampled states, the h-BN tutorial hyperfine tensors (read
  twice too large), and a memory corruption of pCCE with state files.
- **Safer inputs.** A misspelt key, a file that does not open, a combination a method cannot do or a run that does not
  fit in memory stops before the calculation, with one message on one rank; a key that is read but not used is
  reported.
- **Tutorials and documentation.** Eleven tutorials in learning order (`example/tutorial/00_quickstart` …
  `10_coordinate_rotation`), each with `run.sh` and `expected.png`; a rewritten manual (`docs/user/`), method pages
  (`docs/methods/`) and a doxygen site.

### Updating old inputs

`python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites inputs in place, keeping their comments and
listing what it changed (`--check` lists only, `--stdout` prints the new text). A rewritten input gives the same results
as before. The table below lists every renamed or removed key; [Old names](InputFile.md#OldNames) has the same table
with more detail.

| Old | New |
|---|---|
| **Method and cluster shape** | |
| `"method": "CCE"` | `"method": "cCCE"` |
| `"method": "pCCE"`, `"dsj"`, `"itb"`, `"dsjitb"` | `"method": "cCCE"` with `"cluster": {"mode": "pcce"}` (`"dsj"`, `"itb"`, `"dsjitb"`) |
| `"method": "pgCCE"`, `"pmeCCE"`, `"pmegCCE"` (adding_MECCE branch) | `gCCE`, `meCCE`, `megCCE` with `"cluster": {"mode": "pcce"}` |
| `cluster_mode`; its value `adaptive` | `cluster.mode`; `bundle` |
| `addsubclus`, `add_subclusters` | `cluster.selection.add_subclusters` (true/false) |
| `nk`, `max_clusters_per_order` | `cluster.selection.max_per_order` |
| `rbathcut`, `rdipcut` | `rbath: [inner, outer]`, `rdip: [inner, outer]` |
| **Bath configuration** | |
| `gyrofile`, `bathfile`, `seed` (top level) | `configuration.gyrofile`, `configuration.bathfile`, `configuration.seed` |
| `nstate` N (top level or in `configuration`) | `configuration.s_range: [1, N]`; `nstate` 0 (the ensemble): no `s_range` |
| a `seed` number with `k_range`, meaning the same seed for every k | a list of that number, one per configuration (a number is now the first configuration's seed) |
| `bathadjust`, `bath_shift` (top level) | `configuration.bath_shift` |
| `statefile`, `bath_state_file` (top level, a string) | `configuration.bath_state_file: [...]`, one per bath file (`{s}` the state) |
| `avaaxfile`, `Defect.axis_file` | `configuration.defect_axis_file: [...]`, one per bath file |
| `exstatefile`, `Defect.state_file` | `configuration.defect_state_file: [...]`, one per bath file |
| **Qubit** | |
| `qubitfile`, `Qubit.file`, `Qubit.spin`, `Qubit.name`, `qspin` | an entry of `Qubit.qubit`: `[{"file": ..., "spin": ..., "name": ...}]` (`file` in place of `xyz`) |
| `Qubit.nqubit` | removed: the length of `Qubit.qubit` |
| `qalphams`/`qbetams`, `Qubit.alphams`/`betams`, `Qubit.qubit[].alphams`/`betams` | `Qubit.alpha_state: {"ms": [m_1, m_2, ...]}`, `Qubit.beta_state: {"ms": [...]}` (one $m_s$ per qubit) |
| `Qubit.alphaidx`, `Qubit.betaidx` | `Qubit.alpha_state: {"index": n}`, `Qubit.beta_state: {"index": n}` |
| `Qubit.psia`, `Qubit.psib` | `Qubit.alpha_state: {"vector": [...]}`, `Qubit.beta_state: {"vector": [...]}` |
| `qzfs`, `Qubit.zfs` (kHz) and its `Qubit.frame` | a `Qubit.intmap` self entry `{"between": ["q0", "q0"], "tensor": ..., "unit": "kHz", "frame": ...}` |
| `Qubit.qubit[].detuning: [x]` | `x` (a number or `{"value", "unit"}`) |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| **Hamiltonian** | |
| `vector_mean_field` false / true | `hamiltonian.common.bath_bath_mean_field` `"z"` / `"xyz"` (`"off"` is new) |
| `hf_readmode` 0 / 1 / 2 / 3 | `hamiltonian.common.hyperfine`: `dipole` `"point"` / `"point"` / `"dft"` / `"dft"`, `contact` false / true / false / true |
| `hf_tensorfile`, `hf_ignore_oor` 0 / 1 | `hamiltonian.common.hyperfine.file`, `.missing_row` `"stop"` / `"point_dipole"` |
| `qd_readmode` 0; 1; 2 / 3 / 4 | `hamiltonian.common.quadrupole: "off"`; not available; `.structure` `"plain"` / `"strained"` / `"bubble"` |
| `qd_tensorfile`, `qd_tensorfile_woqubit`, `qd_cellpara` | `hamiltonian.common.quadrupole.file`, `.file_without_defect`, `.cell` |
| `knight`, `hfmedi` | `hamiltonian.cCCE.knight`, `hamiltonian.cCCE.hfmedi` |
| `jump_operators` (top level, adding_MECCE branch) | `hamiltonian.meCCE.jump_operators` or `hamiltonian.megCCE.jump_operators` (not rewritten by the converter) |
| jump `operator` `"+-"` / `"-+"` | `"-"` / `"+"` |
| `DefectTotSpin` | removed: the qubit spin is used |
| `CorrTotSpin` | removed: remake old-format VASP tensor files with `example/tools/SpinTensorReader` |
| `hf_cutoff` | removed |
| `magnetic_field_axis` | removed: write `bfield` as [0, 0, Bz] and turn the bath with `frame.rotation` |
| **Pulse** | |
| `npulse`, `pulsename`, `sequence`, `pulseiter`, `pulse_time`, `detuning_factor` | `pulse.count`, `pulse.name`, `pulse.sequence`, `pulse.per_qubit`, `pulse.duration`, `pulse.detuning_factor` |
| `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift` (adding_DEER branch) | `pulse.bath.count`, `.sequence`, `.defect`, `.energy_shift` (not rewritten by the converter) |
| `bspin`, `pulse.bath.spin` | removed: the driven spins' own S |
| `balphams`/`bbetams`, `pulse.bath.alphams`/`betams` | `pulse.bath.levels: [m_a, m_b]` |
| **Defect** | |
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` |
| `dfname`, `eqs`, `rxyzs`, `hypf` | `name`, `quadrupole_moments`, `nuclear_positions`, `hyperfine` |
| `naddspin`, `n_nuclear_spins` | `n_onsite_spins` |
| `navaax`, `n_configurations` | `n_axis_states` |
| `apprx` | removed (it changed nothing) |
| `coordinate_frame` | `frame` |
| species `axis` | `zfs.axis`: `"zfs": {"D": ..., "axis": [...], "unit": ...}` |
| `zfs.E` (0 only) | removed: `zfs.tensor` for E ≠ 0 |
| a unit in another case (`"mhz"`) | the exact spelling (`"MHz"`) |
| **pCCE and bundles** | |
| `sK`, `max_trial`, `max_iter`, `kmeans_pp` | `cluster.pcce.partition_size`, `.trials`, `.max_iterations`, `.kmeans_pp` |
| `iter_detail`, `pcce_detail_limit`, `print_trials` true/false + `print_limit` | `cluster.pcce.print_trials`: the number of trials printed |
| `pcce_partition`, `pcce_connections`, `pcce_cache_mode`, `pcce_cache_file` | `cluster.pcce.partition`, `.connections`, `.cache.mode`, `.cache.file` |
| `bundle_max_dim`, `bundle_merge_cutoff`, `bundle_connections`, `bundle_connection_budget` | `cluster.bundle.max_dim`, `.merge_cutoff`, `.connections`, `.connection_budget` |
| `expansion_max_clusters`, `expansion_max_work` | `cluster.bundle.max_clusters`, `.max_work` |
| top-level `pcce`, `bundle` blocks | `cluster.pcce`, `cluster.bundle` |
| `Adaptive`, `Refine` or top-level `adaptive` block | `cluster.adaptive` |
| `adaptive_work_file`, `probe_memory_gb` | `cluster.adaptive.work_file`, `cluster.adaptive.probe_memory_gb` |
| `cluster.adaptive.output` | `cluster.adaptive.workdir` |
| adaptive keys with `-` (`max-rounds`) | `_` (`max_rounds`) |
| `cluster.adaptive.policy` | removed (one policy) |
| **Execution** | |
| `backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `column_evolution`, `profile` (top level) | the same names in `execution` |
| `gpu_min_dim`, `gpu_batch_clusters`, `gpu_batch_steps`, `gpu_cluster_batch_mb`, `gpu_cache_mb` | `execution.gpu_options.min_dim`, `.evolution_batch_size`, `.evolution_batch_steps`, `.cluster_batch_mb`, `.cache_mb` |
| `gcce_refinement` true/false, `gcce_refinement_max_sweeps` | `execution.gcce_refinement: {"enabled": ..., "max_sweeps": ...}` |
| **Output** | |
| `outfile` | `output.file` |
| `savemode` `normal` / `info` / `all` / `allfull` | `output.savemode` `coherence` / `coherence_tagged` / `cluster_coherence` / `cluster_density_matrix` |
| `save_nodiv` | `output.savemode: {"mode": ..., "nodiv": ...}` |
| `cluster_contributions` true/false, `cluster_contribution_max_mb` | `output.cluster_contributions: {"enabled": ..., "max_mb": ...}` |
| `cluster_only` (top level) | `output.cluster_only` |
| **Frame** | |
| `coordinate_frame_rotation`, `bath_coordinate_rotation` | `frame.rotation`, on when the block is given (a block that was off is removed) |
| `hf_tensor_frame`, `qd_tensor_frame` | `hamiltonian.common.hyperfine.frame`, `hamiltonian.common.quadrupole.frame` |
| `Qubit.intmap[].tensor_frame` | `Qubit.intmap[].frame` |
| frame value `"bath"` | `"file"` |
| no frame with the rotation on (a Defect species, a single qubit's intmap entry, the qubit ZFS) | `"frame": "qubit"` (the old hidden default; now required) |
| `bath_axis`, `reference_qubit`, `qubit_position_frame`, `defect_axis_reference` | removed (x is `qubit_axis` × [0, 0, 1], the origin the first qubit) |
| **Automatic convergence** | |
| `auto_convergence.configurations` | `auto_convergence.sampling`; its file names are the `configuration` block's |
| `converge` (each section) | removed: a section whose value is given is fixed |
| `time_window.unit`, `time_window.threshold` | `time_window.value: {"value", "unit"}`, `time_window.decay_level` |
| `sampling.number`, `sampling.test_number` | `sampling.configurations`, `sampling.test_configurations` |
| `launcher`; `output.path`, `output.file` with `auto_convergence` | removed (the runs use `mpirun -n N` of `main.out` and go to `auto_convergence.workdir`) |
| **Command line** | |
| `-I`, `-N`, `-m`, `-q`, `-B` | `configuration.bathfile`, `configuration.s_range`, `method`, `quantity`, `bfield` |
| `-s`, `-S`, `-a` | `configuration.bath_state_file`, `defect_state_file`, `defect_axis_file` |
| `-P` | `output.cluster_only` |
| `--auto_conv`, `--adaptive` | the `auto_convergence` and `cluster.adaptive` blocks |
| `--auto_status` | `python -m autoconv status FILE` in `example/tools/AutoConvergence` |
| `-o FILE` (an output name) | `-o DIR` is the output folder; the name is `output.file` |

### The full list

Every change, with its details, the speed measurements and the comparison with v1.2.2, is in the file `CHANGELOG.md` at
the root of the CCEX repository.
