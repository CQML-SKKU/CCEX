# check_input_forms/ — the new input forms give the physics of the legacy ones

    bash test/code_verification/check_input_forms/run.sh        # needs bin/main.out (bash do_compile.sh wsl-gnu), a python3 with numpy, scipy and matplotlib (test/requirements.txt)

Stages, in order: `legacy` (writes the legacy-form input of every system to
`_work/`, from `inputs/`, `../bath_data/` and tutorial 03), `axis2`
(legacy form beside new form, one binary), `errcorr` (the analyzer reads the
outputs), `gyro` (the built-in isotope table: no gyrofile gives the
bytes of a gyrofile holding the table's value, a species missing from the gyrofile comes from the
table, a species in neither is refused), `comments` ('#' and '!' are comments only outside JSON
strings: a "#note" key and an outfile with '#' give the plain input's curves), `nodiv` (`savemode` `{"mode": "coherence", "nodiv": false}` writes the
default run's `_wiDiv` and no `_noDiv`, as `save_nodiv` false did; `nodiv` true writes the default's files byte for byte; a string nodiv is refused), `allfull` (`savemode` cluster_density_matrix writes every file cluster_coherence writes, byte for byte, and with gCCE one
full density-matrix file per cluster named after it, no two alike; cCCE and gCCE, 1 and 3 ranks the same files), `evolution_update` (an explicit `exact` gives the
`../check_numerics` references, which were recorded with `exact`, within that suite's 1e-9 — skipped on a CI runner, since the references are one machine's digits; `recursive` runs on single-state cCCE2 and gCCE2 vector, P1 10 ppm, on the ensemble cCCE2 (Hahn, CPMG-8,
Manual) and gCCE2 `evolution` = `matrix` (ensemble Hahn, one-state CPMG-8), and stays within 1e-8 of `exact` in the raw
`_noDiv`; without the key (A29) single-state and ensemble cCCE2 and gCCE2 vector and matrix give the bytes of an explicit
`recursive` and the log says `recursive`, chosen by default - the path is judged by the log, the value the time evolutions
dispatch on; that `exact` writes other bytes is only noted since A56, as on another machine it need not -, and gCCE `propagator` = `expm` those of an explicit `exact`, not a refusal; `backend` = `auto`, all on the CPU in
this binary, gives the bytes of `backend` = `cpu`, with or without a given `recursive`, A36), `gpu_options` (no
block, an empty one, every default written out and every key changed give the same bytes on a CPU run of single-state
cCCE2 and gCCE2, and the log shows the values read; A26, 2026-10-12), `propagator` (A2: `"eigen"` and `{"method": "eigen"}` give the default run's bytes, `"sparse"` and `{"method": "sparse",
"min_dim": 4096}` those of an explicit `exact` on cCCE2 and gCCE2 inputs whose clusters are all below 4096 (sparse runs
`exact` without the key); `{"method": "sparse", "min_dim": 1}` runs every cluster through the Krylov expansion and stays
within 1e-8 of eigen `exact` in the raw `_noDiv` and in `_wiDiv` at |L| <= 1 (P1 and NV-13C cCCE2, P1 gCCE2 with zero-field
splitting); malformed values are refused), `tensor_format` (2026-10-13: on the five h-BN cases of `../check_numerics`, the hyperfine
file of format v3 - tag line and provenance line, read x 1 - gives the bytes of the old untagged v1 file read x 0.5; the
VASP-tagged v1 file of the tutorial until then is read as before, x 1, and gives the bytes of `head-20261012-aa883d0d4`, recorded
with it - skipped on a CI runner, since that reference is one machine's digits (A56); old files warn in the log; values of 10 characters or more, written by SpinTensorReader's VASP and QE writers in
their format and in the tab-glued one before 2026-10-13, are read as written: `example/tools/SpinTensorReader/tests/check_wide_fields.py`;
and a VASP Q file alone from an OUTCAR with magnetization 0 or none is written, while an A file from it is refused:
`example/tools/SpinTensorReader/tests/check_qfile_spin.py`), `autoconv` (`auto_convergence.enabled`: `false` gives
the bytes of the input without the block, an enabled block makes `main.out -f` start the runner, once under `mpirun -n 2`,
a value that is not a boolean and another option beside it are refused, `--auto_status` reads the workdir, and a tool
file or folder writable by others is refused). The refused key combinations are `../check_refusals/`, the
rotation `../check_rotation/`; the comparison against v1.1.0 and
v1.2.2 is `../check_numerics/`. Nothing here needs an old binary.

Exit status 0 only if every stage passed. Everything written goes to `_work/`,
which is gitignored.

## What this answers

`test/code_verification/check_rotation/` shows that the rotation does what it says. It cannot show that
the branch left v1.1.0 alone, because every comparison it makes is between two
runs of the *same* binary. That is what this suite is for.

## Why it takes three axes

The obvious check -- run the new input under v1.1.0 and compare -- is not
available: v1.1.0's parser rejects every key this branch adds. So the answer is
assembled from comparisons that can each actually be made.

| axis | compares | script | verdict |
|---|---|---|---|
| 1 | ~~the same **legacy** input, v1.1.0 binary vs this one~~ — retired 2026-09-24: since the fixes of 2026-09-19 (deltat as a double, the even-pulse conjugate) the numbers are not byte-identical to v1.1.0 on purpose; `../check_numerics/` compares against v1.1.0 and v1.2.2 with a tolerance and names those differences | (`axis1.sh` in git history up to 0ef682b54) | — |
| 2 | **legacy form vs new form** of the same physics, one binary | `axis2.sh` | byte-identical, or a reported deviation where a unit conversion sits between them |
| 3 | a **hand-transformed** input with the rotation off vs the original with it on | `../check_rotation/run.sh` `[A/B]` | numeric agreement |

Axes 1 and 2 composed: if legacy(v1.1.0) == legacy(dev) and legacy(dev) ==
new-form(dev), then new-form(dev) == legacy(v1.1.0). With axis 1 retired, the
legacy side is held to v1.1.0 and v1.2.2 by `../check_numerics/` instead, within
its tolerance and named differences.

## The systems

Existing fixtures, reduced. Equivalence needs reproducibility, not convergence,
so `rbath` and `nstep` are cut until a case runs in seconds on one core. What
each was cut to, and from what, is in the generator next to the numbers.

| | system | reaches | wall (1 core) |
|---|---|---|---|
| S1a/b | NV + natural-abundance 13C, CCE-2, ensemble and sampled bath | the ordinary CCE path, point-dipole hyperfine | 0.3 s |
| S2a | the same, gCCE-2 | gCCE | 0.35 s |
| S2b/c | NV + P1 10 ppm, CCE-2 and gCCE-2 | `Defect` with `n_nuclear_spins = 0`, `n_configurations = 12`, `Defect.axis_file` | 0.3 / 1.0 s |
| S3 | P1 bath, `Defect` `n_nuclear_spins = 1` (14N), `n_configurations = 4` | `gyros` `quadrupole_moments` `hyperfine` `efg` `zfs` `detuning`, state and extra-state files | 1.4 s |
| S4 | hBN V_B, `hyperfine` `{"dipole": "dft", "contact": true}` + `quadrupole` `"plain"` | both DFT tensor readers, the vertex boundary, `hf/qd_tensor_frame` | 0.44 s |
| S4b | the same over `Afile` instead of `Afile_vertex` | `CheckBD_Range`, the axis-aligned MinDif/MaxDif box | 0.43 s |
| S5 | two NV, gCCE, explicit `intmap` self/ZFS tensors | `intmap` `tensor_frame`, `reference_qubit`, multi-qubit rotation | 0.29 s |

S3's fixture is written for 24 ranks at `rbath` 850 (469 spins, 80 s). At
`rbath` 360 it keeps 45 spins and 518 pair clusters, which still reaches every
Defect code path. S4's is cut from `rbath` 20 to 6: 3897 spins and 121 s become
101 spins and 0.44 s, well inside the region the A and Q files cover, so both
readers stay on their file path rather than falling back.

S4 and S4b differ only in which boundary test the tensor file selects, and a file
picks exactly one: `Afile_vertex` carries `v1..v8` vertices and reaches
`CheckBD_vertex`, `Afile` carries `MinDif[A]` / `MaxDif[A]` and reaches
`CheckBD_Range`. Without S4b the axis-aligned box -- the one thing that makes a
rotation of the tensor file's own coordinates impossible -- would never run.

S4b needs `hyperfine.missing_row = "point_dipole"` and no quadrupole: both boxes are wider than the set
of atoms their files list, so a spin can sit inside the box with no row of its own.
The HF reader has a flag for that and the QD reader does not, which leaves the QD
range boundary uncovered. So are the quadrupole structures `strained` and `bubble`, which upstream also
records as untested rather than verified.

## What is compared

The **raw** `_noDiv` / `_wiDiv`. Not the error-corrected file and not the
ensemble average: both smooth, and a real difference could hide under them.

`errcorr.sh` is the other half of that decision. It runs the analyzer on each
baseline and reports `max|L|` before and after correction, so the suite says out
loud whether the systems it compared are numerically sound. They mostly are;
`S2a` is not -- gCCE-2 on that 13C bath produces `max|L| = 4.2` at two points,
which the correction takes back to 1. That is the ordinary CCE division
divergence, it is present in v1.1.0 exactly as it is here, and it does not
disqualify the system as an equivalence fixture. It does mean the raw file is
not physics.

## Unit conversions

`gyros`, `eqs`, `hypf`, `efg`, `zfs` and `detuning` are each given twice: once in
the unit whose scale factor is exactly 1 -- where byte identity with the legacy
array is required -- and once in a second unit, where it is not. Three of the
four second units come back byte-equal anyway. `efg` in `V/angstrom^2` lands
1e-10 away, which is one unit in the last printed place of a file written to ten
decimals.

## Refusals

`../check_refusals/run.sh` re-checks the refused combinations on the physical fixtures rather
than the toy ones, and requires two things of each: the run stops with the
expected message, and neither DFT tensor reader announced itself first. The
second half is what would regress quietly if a guard ever drifted later into the
run -- the refusal would still happen, but only after a file had been read in a
frame it was not written in.
