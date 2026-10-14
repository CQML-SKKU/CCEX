# check_refusals/ — the inputs that must be refused are refused

    bash test/code_verification/check_refusals/run.sh     # needs bin/main.out

246 inputs that must be refused (4 of them outputs that cannot be written), eight that must run, and `--help` (the same as `-h`), each run once to see
that CCEX stops with the expected message, most of them before it reads any tensor
file. They are edits of the legacy inputs `../check_input_forms/legacy.sh`
writes (so it is run first), on the real systems: the h-BN bath with DFT
tensor files (S4), two NV centres with explicit `intmap` tensors (S5), the P1
defect bath (S3).

| refused | because |
|---|---|
| a DFT `hyperfine` file / a `quadrupole` with the rotation on and no `hf_tensor_frame` / `qd_tensor_frame` | a tensor file has a frame, and the run must not guess it |
| a tilted `bfield` with the rotation on | the field is applied along z of the computational frame |
| `qubit_position_frame` other than `bath` | only one meaning is defined |
| two qubits without `reference_qubit`, or with a name that is not a qubit | the rotation is about one qubit |
| an `intmap` tensor without its frame | same as the tensor files |
| `magnetic_field_axis`, `defect_axis_reference` without the rotation | keys that belong to the rotation block alone |
| `zfs` given as D and E without `axis`; `axis` without `coordinate_frame`; object-form `zfs` without `coordinate_frame`; E non-zero with one axis | a zero-field splitting needs its axes and their frame |
| a Jahn-Teller axis index out of `1..n_configurations` | it indexes the per-axis tensors |
| an unknown `unit` | the conversion would be silent |
| `coordinate_frame_rotation` and `bath_coordinate_rotation` both present | the same option under two names |

Without a row of their own (`run.sh` has a comment on each): `-N`/`-B` that are not numbers, `-P` with no output name, and the command line
winning over the file in either order; a missing `Qubit.file`; `nstep` 0; a Manual sequence that does not sum to 1 or
is malformed (a number fraction and an angle of 90.5 must run); a number outside the float it is stored in;
`alpha_state`/`beta_state` `{"index"}` naming other states than the written `alphams`/`betams` (D43; one that agrees must run);
`intmap` tensors shorter than 3x3 or not numbers (D28); no qubit in `Qubit.qubit`; malformed entries of
the indexed form (D28); arrays shorter than the input needs (D42); a `Defect.axis_file` axis past `n_configurations`; a key of the
wrong type (D18); a misspelt key (D24); `quantity` `noise`, in the file and as `-q noise` (nothing computes it; review 20261013-doxygen D22); `Qubit.file` together with `Qubit.qubit` (until 2026-10-13 `qubitfile` with a `Qubit` block, and the block was ignored; review 20261013-readme R01).

Since 2026-10-01 it also checks failures that used to end with exit status 0 (or, for the second file in info mode, crash):

| refused | because |
|---|---|
| an h-BN tensor file without its `g-factor___` lines, or without its `etc____` lines | the reader said the condition could not be read and then exited 0 (`READ_Tensor_const`, `READ_Tensor_etc`). These come after a tensor reader has started, so the "before any file is read" condition is not applied to them |
| a result file that cannot be written in full, in `normal` and in `info` mode | `Output_save` and `Output_save_info` checked only that the files opened. With its result files linked to `/dev/full` (every write fails with ENOSPC, as on a full disk), `check_numerics/cases/natab_cce2_hahn` ended with exit 0. Until 2026-10-02 the case used a 4 KB file-size limit (`ulimit -f 4`), which also broke Open MPI's own start-up |
| a `_noDiv` file that cannot be opened, in both modes | a directory in its place; `Output_save_info` checked only the `_wiDiv` file and would have written to a null file |

Since 2026-10-11 it also checks the requests of the `factors` and `except` work modes of method D's regroup
((internal notes) 4절, `-stage3b.md` 1절), written by hand on the P1 input of the regroup
regression (`test/cluster_verification/check_adaptive/regression/inputs/regroup-p1-0.json`). Each has to be refused
before any state is calculated (no response file of a state is written):

| refused | because |
|---|---|
| a factors state that is not a state of the run (0 with `nstate` 5, or 6), a step outside `0..nstep-1`, a step that is not an integer | the factors are one state's point of the run |
| factors or except with `method` gCCE, with `order` 1, with `probes`, with `cluster_contributions`, with `-P` | the modes are defined for cCCE bundle order 2 over an explicit partition, with nothing else in the product |
| `state`/`step` outside factors, `except` outside except | a key a mode does not read is not ignored |
| an empty except list, a term listed twice, a site set that is not a term of the partition's expansion | the driver works X out of the same exponent table; one that differs is a different table |

Since 2026-10-11 it also checks `evolution_update` (`docs/user/InputFile_Backend.md`), on inputs of
`../check_numerics/cases/`, refused before the calculation starts:

| refused | because |
|---|---|
| a value other than `exact` or `recursive` (a string, or a number) | only the two are defined |
| `recursive` given with gCCE `propagator` = `expm` (`p1_manual_hahn`, `evolution` = `matrix`), or with `backend` = `cuda` (`p1_cce2_single`) | the recurrence exists on the CPU and in the packed kernel only (the GPU large-cluster path computes as `exact`), and `expm` forms no phases to carry over. Without the key these paths run `exact` (A29, 2026-10-12). The ensemble cCCE and gCCE `matrix` were refused until A18 (2026-10-11) and run now (`../check_input_forms/`); `backend` = `auto` in this CPU-only binary computes every cluster on the CPU and was refused until A36 (2026-10-12), it runs now (`../check_input_forms/` E5). With `backend` `auto` or `cuda` in `bin/main_gpu.out` a given `recursive` is refused: `test/gpu_verification/check_gpu/check_options.py` |

Since 2026-10-12 it also checks the `gpu_options` block (`docs/user/InputFile_Backend.md#tag_gpu_options`, A26), on
`../check_numerics/cases/`: each of the five flat keys it replaced (`gpu_min_dim` ...) is refused with the new form
named, and so are a block that is not an object, an unknown key in it, `packed_max_dim` 3, 8 or a string, `min_dim` 0,
`packed_max_dim` in `bin/main.out` (it needs the GPU binary) and with gCCE (it is cCCE only). The accepted forms are
`../check_input_forms/gpu_options.sh`.

Since 2026-10-12 it also checks `propagator` (A2, `docs/user/InputFile_Backend.md`), on inputs of `../check_numerics/cases/`,
refused before the calculation starts (the malformed values are in `../check_input_forms/propagator.sh`):

| refused | because |
|---|---|
| `expm` with cCCE (`p1_cce2_single`) | the matrix exponential path exists in gCCE only |
| `sparse` with the ensemble (`p1_cce2_ensemble`), gCCE `evolution` = `matrix` (`p1_manual_hahn`), `backend` = `cuda` or `auto`, `evolution_update` = `recursive` (cCCE and gCCE), `knight` (`natab_knight_shift`), `hfmedi` (`hbn_hfmedi`), two qubits (`natab_gcce2_two_qubits`) | the Krylov path propagates one state vector on the CPU with the phases of `exact`, from the Hamiltonian terms it assembles itself |

Since 2026-10-13 it also checks hyperfine tensor files of format v3 (`docs/user/InputFile_General.md#tag_hyperfine_file`), on
`hbn_cce2_hf3_qd0` with the h-BN v3 file edited once each: a first line that is not the tag of the provenance line's source
(both ways), the provenance line not on line 2 or twice, an unknown format version or source, no S or an S that is not a
half-integer, and S other than the qubit spin (both ways; `DefectTotSpin` until 2026-10-13). What such files are read as is `../check_input_forms/tensor_format.sh`.

Since 2026-10-13 it also checks the `savemode` block (`docs/user/InputFile_Output.md#tag_savemode`), on
`p1_cce2_single`: the old names `normal`, `info`, `all`, `allfull` (also as `{"mode": "all"}`) and the top-level
`save_nodiv` are refused with the new form named (`savemode "all" is now "cluster_coherence"`, `save_nodiv is now
savemode.nodiv`), and so are an unknown key of the block, a `savemode`, `mode` or `nodiv` of the wrong type and an unknown
mode. The accepted block is `../check_input_forms/nodiv.sh`.

Since 2026-10-13 (A47) it also checks the renamed input keys (`docs/user/InputFile.md#OldNames`): each case puts one old key,
old value or old form back into a converted `check_numerics` input, and the run must stop naming the new key and the
converter; gCCE with `cluster_mode` pcce, dsj, itb or dsjitb (from the input and from `-m`); a `hamiltonian.cCCE` setting in a
gCCE run; the `adaptive` block outside `cluster_mode` pcce or bundle and the `pcce` block outside pcce; `delta_A` false with a
Defect species or with `qubit_bath_mean_field` true; an unknown `bath_bath_mean_field` value; and both forms of a qubit
state at once. What the two new settings compute is `../check_analytic/`.

`../check_rotation/run.sh` refuses the same combinations on toy systems; this
folder shows the guards fire on real inputs and before any file is interpreted.
What it writes goes to `../check_input_forms/_work/`, which is ignored.
