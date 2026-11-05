# check_refusals/

Checks that the inputs which must be refused are refused, with the expected message, and most of them before any
tensor file is read. The cases are edits of real inputs: the legacy inputs of `../check_input_forms/legacy.sh` (S3,
S4, S5), the cases of `../check_numerics/cases/` and the P1 regroup input of `check_adaptive`.

```
bash test/code_verification/check_refusals/run.sh     # bin/main.out or CCEX_BIN; runs ../check_input_forms/legacy.sh first
```

It writes to `../check_input_forms/_work/`, which git ignores.

## What it checks

| Group | Refused | Expected |
|---|---|---|
| Frames | with the rotation on: directional input (DFT `hyperfine`, `quadrupole`, `intmap`, `Defect`) without `frame`; a tilted `bfield` | message, before any tensor file is read |
| Old keys | removed or renamed keys, values and places (`frame.rotation` keys, top-level `backend`, `outfile`, old `alphams`, old `savemode` names, flat `gpu_*` keys, ...) | message naming the new form |
| Command line | removed options `-I -N -m -q -B -s -S -a -P --auto_conv --adaptive --auto_status`, unknown options, `-o` and `-e` misuse | message; `--help` equals `-h` |
| Malformed values | missing axes, out-of-range indices, unknown `unit`, wrong types, misspelt keys, short arrays, bad Manual sequences, `nstep` 0, bad qubit states | message |
| Combinations | gCCE with dsj or itb; `hamiltonian.cCCE` in gCCE; `adaptive` outside bundle; `delta_A` without a spin-1 ZFS; `quantity` `noise`; bad `evolution_update`, `propagator` and `gpu_options` pairings | message, before the calculation |
| v3 tensor files | wrong tag or provenance line, unknown version or source, bad S | message |
| Regroup `factors`, `except` | states or steps outside the run, gCCE, order 1, keys outside their mode | before any state is calculated |
| Must not exit 0 | h-BN tensor file without `g-factor___` or `etc____`; result files on `/dev/full`; a directory in place of `_noDiv` | non-zero exit |
| Memory | `SLURM_MEM_PER_NODE=1` | after the clusters, before the main calculation |
| Must run | `{"index"}` eigenvector phases, a 90.5 angle in Manual, `checked_merged` (bundle, `rdip` omitted, `max_dim` 8, `merge_cutoff` 1) | finishes; `checked_merged` prints `[adaptive-check]` and a coherence file |

`run.sh` has a comment on every case. The inputs are written with the old flat key names and put into the current
blocks by `../input_layout.py`; a case that needs an old place lists it in `__keep_old__`. The accepted forms are in
`../check_input_forms/`; the GPU binary's refusals in `../../gpu_verification/check_gpu/check_options.py`.

## Files

| File | Contents |
|---|---|
| `run.sh` | writes every case, runs it, and checks the exit status and message |
