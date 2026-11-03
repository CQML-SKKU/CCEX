# check_refusals/: the inputs that must be refused are refused

    bash test/code_verification/check_refusals/run.sh     # needs bin/main.out

267 inputs that must be refused (4 of them outputs that cannot be written, 18 removed or wrong command-line uses),
five that must run, and `--help` (the same as `-h`). Each is run once to see that CCEX stops with the expected
message, most of them before it reads any tensor file. They are edits of the legacy inputs
`../check_input_forms/legacy.sh` writes (so it is run first), on the real systems: the h-BN bath with DFT tensor
files (S4), two NV centres with explicit `intmap` tensors (S5), the P1 defect bath (S3).

`../check_rotation/run.sh` refuses the same combinations on toy systems; this folder shows the guards fire on real
inputs and before any file is interpreted. What it writes goes to `../check_input_forms/_work/`, which is ignored.

The inputs of these cases are written with the old flat key names and put into the current blocks by
`../input_layout.py` (`flat` when read, `nest` when written); a case of an old place lists it in `__keep_old__`.

## Frames, removed keys and old places

| refused | because |
|---|---|
| with the rotation on, a DFT `hyperfine` file / a `quadrupole` without `hamiltonian.common.hyperfine.frame` / `quadrupole.frame`, a `Qubit.intmap` entry (one qubit or more, the `Qubit.file` ZFS included) without `frame`, a `Defect` species without `frame` | numbers with a direction have a frame, and the run must not guess it |
| a tilted `bfield` with the rotation on | the field is applied along z of the computational frame |
| `frame.rotation.bath_axis`, `reference_qubit`, `qubit_position_frame` (removed) | they changed nothing a result depends on |
| `frame.rotation.hf_tensor_frame` / `qd_tensor_frame`, `Qubit.intmap[].tensor_frame`, `Defect.species[].coordinate_frame`, the value `"bath"` | each is refused with its new name (`hamiltonian.common.hyperfine.frame`, `quadrupole.frame`, `frame`) and value (`"file"`) |
| `magnetic_field_axis`; `defect_axis_reference` (removed) | a key that changed nothing would look read |
| `zfs` given as D and E without `axis`; E non-zero with one axis | a zero-field splitting needs its axes |
| a Jahn-Teller axis index out of `1..n_axis_states` | it indexes the per-axis tensors |
| an unknown `unit` | the conversion would be silent |
| `bath_coordinate_rotation` (and `coordinate_frame_rotation`) | the rotation is `frame.rotation` |
| the old places (top-level `backend`, `propagator`, `cluster_mode`, `add_subclusters`, `max_clusters_per_order`, `probe_memory_gb`, `outfile`, `savemode`, the Defect species' `n_nuclear_spins`, `n_configurations`) | each is refused with its new place (`execution.`, `cluster.`, `output.`, `n_onsite_spins`, `n_axis_states`) |
| `output.path` with an absolute `output.file`; an empty `output.path` | two places for the outputs; no folder |

Without a row of their own (`run.sh` has a comment on each):
- the removed command-line options `-I -N -m -q -B -s -S -a -P --auto_conv --adaptive --auto_status`, each refused
  naming its input key; an unknown option; an absolute `output.file` with `-o DIR`; a file where `-o` needs a folder;
  `-e` with `-f`, and `-e` over an existing `ccex_example.jsonc` (the accepted `-o`/`-e` forms are
  `../check_input_forms/cli.sh`);
- a missing `Qubit.file`; `nstep` 0; a Manual sequence that does not sum to 1 or is malformed (a number fraction and
  an angle of 90.5 must run); a number outside the float it is stored in;
- the old per-qubit `alphams`/`betams` and `Qubit.alphams`/`betams`, `Qubit.zfs`, `Qubit.frame` (with the new form in
  the message); an `alpha_state`/`beta_state` with more than one or none of `ms`/`index`/`vector`; an `ms` list of the
  wrong length or with a non-number; a bare number or an `index` with `Qubit.file`; a `Qubit.file` intmap entry naming
  another qubit than `Qubit.name`;
- the phase of an `{"index"}` eigenvector (the largest component real and positive, on a tie of two equal components
  and a mixed NV state; these run, `a93_phase_*`);
- `intmap` tensors shorter than 3x3 or not numbers; no qubit in `Qubit.qubit`; malformed entries of the indexed form;
  arrays shorter than the input needs; a `configuration.defect_axis_file` axis past `n_axis_states`; a key of the
  wrong type; a misspelt key;
- `quantity` `noise`, in the file and as `-q noise` (nothing computes it); `Qubit.file` together with `Qubit.qubit`.

## Failures that must not end with exit status 0

| refused | because |
|---|---|
| an h-BN tensor file without its `g-factor___` lines, or without its `etc____` lines | the reader (`READ_Tensor_const`, `READ_Tensor_etc`) must not report that the condition cannot be read and then exit 0. These come after a tensor reader has started, so the "before any file is read" condition is not applied to them |
| a result file that cannot be written in full, in `normal` and in `info` mode | `Output_save` and `Output_save_info` must check every write, not only that the files opened. The result files of `check_numerics/cases/natab_cce2_hahn` are linked to `/dev/full` (every write fails with ENOSPC, as on a full disk). A file-size limit (`ulimit -f 4`) is not used because it also breaks Open MPI's own start-up |
| a `_noDiv` file that cannot be opened, in both modes | a directory in its place; `Output_save_info` must not write to a null file |

## The `factors` and `except` work modes of method D's regroup

Written by hand on the P1 input of the regroup regression
(`test/cluster_verification/check_adaptive/regression/inputs/regroup-p1-0.json`). Each has to be refused before any
state is calculated (no response file of a state is written).

| refused | because |
|---|---|
| a factors state that is not a state of the run (0 with `nstate` 5, or 6), a step outside `0..nstep-1`, a step that is not an integer | the factors are one state's point of the run |
| factors or except with `method` gCCE, with `order` 1, with `probes`, with `cluster_contributions`, with `cluster_only` | the modes are defined for cCCE bundle order 2 over an explicit partition, with nothing else in the product |
| `state`/`step` outside factors, `except` outside except | a key a mode does not read is not ignored |
| an empty except list, a term listed twice, a site set that is not a term of the partition's expansion | the driver works X out of the same exponent table; one that differs is a different table |

## `evolution_update`

(`docs/user/InputFile_Backend.md`), on inputs of `../check_numerics/cases/`, refused before the calculation starts.

| refused | because |
|---|---|
| a value other than `exact` or `recursive` (a string, or a number) | only the two are defined |
| `recursive` given with gCCE `propagator` = `expm` (`p1_manual_hahn`, `evolution` = `matrix`), or with `backend` = `cuda` (`p1_cce2_single`) | the recurrence exists on the CPU and in the packed kernel only (the GPU large-cluster path computes as `exact`), and `expm` forms no phases to carry over. Without the key these paths run `exact` |

The ensemble cCCE and gCCE `matrix` run with `recursive`, and `backend` = `auto` in this CPU-only binary computes every
cluster on the CPU; both are checked in `../check_input_forms/`. With `backend` `auto` or `cuda` in `bin/main_gpu.out`
a given `recursive` is refused: `test/gpu_verification/check_gpu/check_options.py`.

## `gpu_options`

(`docs/user/InputFile_Backend.md#tag_gpu_options`), on `../check_numerics/cases/`. Each of the five flat keys it
replaced (`gpu_min_dim` ...) is refused with the new form named. So are a block that is not an object, an unknown key
in it, `packed_max_dim` 3, 8 or a string, `min_dim` 0, `packed_max_dim` in `bin/main.out` (it needs the GPU binary)
and with gCCE (it is cCCE only). The accepted forms are `../check_input_forms/gpu_options.sh`.

## `propagator`

(`docs/user/InputFile_Backend.md`), on inputs of `../check_numerics/cases/`, refused before the calculation starts.
The malformed values are in `../check_input_forms/propagator.sh`.

| refused | because |
|---|---|
| `expm` with cCCE (`p1_cce2_single`) | the matrix exponential path exists in gCCE only |
| `sparse` with the ensemble (`p1_cce2_ensemble`), gCCE `evolution` = `matrix` (`p1_manual_hahn`), `backend` = `cuda` or `auto`, `evolution_update` = `recursive` (cCCE and gCCE), `knight` (`natab_knight_shift`), `hfmedi` (`hbn_hfmedi`), two qubits (`natab_gcce2_two_qubits`) | the Krylov path propagates one state vector on the CPU with the phases of `exact`, from the Hamiltonian terms it assembles itself |

## Hyperfine tensor files of format v3

(`docs/user/InputFile_General.md#tag_hyperfine_file`), on `hbn_cce2_hf3_qd0` with the h-BN v3 file edited once each:
a first line that is not the tag of the provenance line's source (both ways), the provenance line not on line 2 or
twice, an unknown format version or source, no S or an S that is not a half-integer, and S other than the qubit spin
(both ways). What such files are read as is `../check_input_forms/tensor_format.sh`.

## `savemode`

(`docs/user/InputFile_Output.md#tag_savemode`), on `p1_cce2_single`. The old names `normal`, `info`, `all`, `allfull`
(also as `{"mode": "all"}`) and the top-level `save_nodiv` are refused with the new form named
(`savemode "all" is now "cluster_coherence"`, `save_nodiv is now savemode.nodiv`). So are an unknown key of the block,
a `savemode`, `mode` or `nodiv` of the wrong type and an unknown mode. The accepted block is
`../check_input_forms/nodiv.sh`.

## Renamed keys and invalid combinations

(`docs/user/InputFile.md#OldNames`) Each case puts one old key, old value or old form back into a converted
`check_numerics` input, and the run must stop naming the new key and the converter. Also refused: gCCE with
`cluster_mode` dsj, itb or dsjitb; a `hamiltonian.cCCE` setting in a gCCE run; the `adaptive` block outside
`cluster_mode` pcce or bundle and the `pcce` block outside pcce; `delta_A` (the Maze δA term) without a spin-1 qubit's
zero-field splitting, or with `hfmedi`; an unknown `bath_bath_mean_field` value; and both forms of a qubit state at
once. What the two new settings compute is `../check_analytic/`.

## Units and memory

A key's `{"value", "unit"}` with an unknown unit (also in the wrong case, `MT`), with `values` for `value`, a string
value or a third key is refused with the units of the key. The memory check refuses with `SLURM_MEM_PER_NODE=1`
(1 MB per node) after the clusters and before the main calculation.

## A run that has to go through

`checked_merged`: bundle clustering with `rdip` omitted (`connections` "checked") and bundles of several sites
(`max_dim` 8, `merge_cutoff` 1) on `../../cluster_verification/check_bundle/checked_ensemble.json`. It once stopped
with a segmentation fault; it has to end with its `[adaptive-check]` line and a coherence file.
