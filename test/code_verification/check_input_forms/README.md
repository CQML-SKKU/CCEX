# check_input_forms/: the new input forms give the physics of the legacy ones

    bash test/code_verification/check_input_forms/run.sh        # needs bin/main.out (bash do_compile.sh wsl-gnu), a python3 with numpy, scipy and matplotlib (test/requirements.txt)

Exit status 0 only if every stage passed. Everything written goes to `_work/`, which is gitignored. Nothing here
needs an old binary. The refused key combinations are `../check_refusals/`, the rotation `../check_rotation/`, the
comparison against v1.1.0 and v1.2.2 `../check_numerics/`.

## Stages

They run in this order.

| stage | what it checks |
|---|---|
| `legacy` | writes the legacy-form input of every system to `_work/`, from `inputs/`, `../bath_data/` and tutorial 03 |
| `axis2` | legacy form beside new form, one binary |
| `errcorr` | the analyzer reads the outputs |
| `gyro` | the built-in isotope table: no gyrofile gives the bytes of a gyrofile holding the table's value; a species missing from the gyrofile comes from the table; a species in neither is refused |
| `comments` | '#' and '!' are comments only outside JSON strings: a "#note" key and an outfile with '#' give the plain input's curves |
| `nodiv` | `savemode` `{"mode": "coherence", "nodiv": false}` writes the default run's `_wiDiv` and no `_noDiv`; `nodiv` true writes the default's files byte for byte; a string nodiv is refused |
| `allfull` | `savemode` cluster_density_matrix writes every file cluster_coherence writes, byte for byte, and with gCCE one full density-matrix file per cluster named after it, no two alike; cCCE and gCCE, 1 and 3 ranks write the same files |
| `evolution_update` | see below |
| `column_evolution` | see below |
| `gpu_options` | no block, an empty one, every default written out and every key changed give the same bytes on a CPU run of single-state cCCE2 and gCCE2; the log shows the values read |
| `propagator` | see below |
| `tensor_format` | see below |
| `cli` | see below |
| `configuration` | see below |
| `autoconv` | `auto_convergence.enabled` `false` gives the bytes of the input without the block; an enabled block makes `main.out -f` start the runner, once under `mpirun -n 2`; a value that is not a boolean and another option beside it are refused; `python -m autoconv status` reads the workdir; a tool file or folder writable by others is refused |
| `units` | a key written as `{"value": ..., "unit": ...}` in another unit gives the bytes of the plain number (deltat in us, bfield in mT, rbath and rdip in nm, the qubit's xyz in nm and detuning in Hz, bath_shift in nm, a meCCE jump rate in kHz, pulse.duration in us); every key in its own unit written out gives the plain input's bytes. DEER energy_shift and bandwidth are in `../check_deer` |
| `pcce_rest` | pCCE with spins left over by `partition_size` and states read from `bath_state_file` finishes, byte for byte the run whose files lack the left-over spins' lines. With `CCEX_BIN` a `CCEX_SANITIZE=1` build, a write into the deleted spins stops the run |

**`evolution_update`.** An explicit `exact` gives the `../check_numerics` references, which were recorded with
`exact`, within that suite's 1e-9. This part is skipped on a CI runner, since the references are one machine's
digits. `recursive` runs on single-state cCCE2 and gCCE2 vector (P1 10 ppm), on the ensemble cCCE2 (Hahn, CPMG-8,
Manual) and on gCCE2 `evolution` = `matrix` (ensemble Hahn, one-state CPMG-8). It stays within 1e-8 of `exact` in
the raw `_noDiv`. Without the key, single-state and ensemble cCCE2 and gCCE2 vector and matrix give the bytes of an
explicit `recursive`, and the log says `recursive`, chosen by default. The path is judged by the log, the value the
time evolutions dispatch on. That `exact` writes other bytes is only noted, since on another machine it need not.
gCCE `propagator` = `expm` without the key gives the bytes of an explicit `exact`, not a refusal. `backend` = `auto`,
all on the CPU in this binary, gives the bytes of `backend` = `cpu`, with or without a given `recursive`.

**`column_evolution`.** On by default. Without the block the run gives the bytes of an explicit
`{"enabled": true, "min_dim": 256, "columns": 64}`, and the log prints both values. `{"enabled": false}` gives the bytes
of a run whose `min_dim` no cluster reaches. The default stays within 1e-10 of that in the raw `_noDiv` (on WSL the
same written digits; the path is judged by the log). The case is p1_cce2e2n at order 4 with one cluster of
dimension 1296.

**`propagator`.** `"eigen"` and `{"method": "eigen"}` give the default run's bytes. `"sparse"` and
`{"method": "sparse", "min_dim": 4096}` give those of an explicit `exact` on cCCE2 and gCCE2 inputs whose clusters are
all below 4096 (sparse runs `exact` without the key). `{"method": "sparse", "min_dim": 1}` runs every cluster through
the Krylov expansion and stays within 1e-8 of eigen `exact` in the raw `_noDiv` and in `_wiDiv` at |L| <= 1 (P1 and
NV-13C cCCE2, P1 gCCE2 with zero-field splitting). Malformed values are refused.

**`tensor_format`.** On the five h-BN cases of `../check_numerics`, the hyperfine file of format v3 (tag line and
provenance line, read x 1) gives the bytes of the untagged v1 file read x 0.5. The VASP-tagged v1 file of the tutorial
is read x 1 and gives `head-20261012-aa883d0d4`, recorded with it (_noDiv/dm within 1e-12 relative). That part is
skipped on a CI runner, since the reference is one machine's digits. Old files warn in the log. Values of 10
characters or more, written by SpinTensorReader's VASP and QE writers in their format and in the older tab-glued one,
are read as written: `example/tools/SpinTensorReader/tests/check_wide_fields.py`. A VASP Q file alone from an OUTCAR
with magnetization 0 or none is written, while an A file from it is refused:
`example/tools/SpinTensorReader/tests/check_qfile_spin.py`.

**`cli`.**
- `-e` writes the template `ccex_example.jsonc` alone, with and without `-o DIR`.
- A sed copy with its `bathfile` pointed at 00_quickstart's bath (the path written with a `//` inside the string)
  runs to the end with no unread key and writes `run_*`. The `//` comments are read and `.jsonc` leaves the output
  name.
- Every trailing comment of the template starts at column 41.
- `example_keys.py` holds the example to the input readers' source. Every key `option.cpp` and `main.cpp` look up with
  `cJSON_GetObjectItem`/`readBlock` is written as `"key":`. Every other identifier-like literal of `option.cpp` appears
  in it, except the refused old names and a short list of internal names. Every option of the refinement and every
  key of the automatic convergence is written as `"key":`, and every `"key":` of the example is one of those names.
- `example_autoconv.py` uncomments the example's `auto_convergence` lines and has the runner read every choice of the
  alternative lines (a key written twice: found or fixed). All are accepted but the one with every section fixed.
- The `-o` rules: no `output.file` gives the input name; `-o DIR` (made, two levels) holds it; a relative
  `output.file` goes inside `DIR` with `-o` before or after `-f`; an absolute one without `-o` is that path;
  `output.cluster_only` true writes `DIR/<name>.clusters.json` and no curve; `output.path` (made, two levels) holds
  the outputs, and `-o DIR` wins over it.

**`configuration`.** `configuration.k_range` runs every configuration k, each `<outfile>_k<k>` byte for byte the run of
its files alone, with random axes and states and with state, axis and defect state files per k. Two bath files with
their per-bath-file lists give the bytes of the files merged. `{{k}}` and `{k}` without `k_range` are literal. There
is a pCCE cache file per k. Refused before any calculation: `k_range` without `{k}`, a missing file of one k, a bad
`k_range`, a list of the wrong length or a string for a list, an enabled cache without `{k}`, the old top-level
`bathfile` and `nstate`.

## What this answers

`test/code_verification/check_rotation/` shows that the rotation does what it says. It cannot show that the code
left v1.1.0 alone, because every comparison it makes is between two runs of the *same* binary. This suite does that.

## Why it takes three axes

The obvious check, running the new input under v1.1.0, is not available: v1.1.0's parser rejects every new key. So
the answer is assembled from comparisons that can each be made.

| axis | compares | script | verdict |
|---|---|---|---|
| 1 | the same **legacy** input, v1.1.0 binary vs this one | not run | the numbers differ from v1.1.0 on purpose (deltat as a double, the even-pulse conjugate); `../check_numerics/` compares against v1.1.0 and v1.2.2 with a tolerance and names those differences |
| 2 | **legacy form vs new form** of the same physics, one binary | `axis2.sh` | byte-identical, or a reported deviation where a unit conversion sits between them |
| 3 | a **hand-transformed** input with the rotation off vs the original with it on | `../check_rotation/run.sh` `[A/B]` | numeric agreement |

If legacy(v1.1.0) == legacy(dev) and legacy(dev) == new-form(dev), then new-form(dev) == legacy(v1.1.0). The legacy
side is held to v1.1.0 and v1.2.2 by `../check_numerics/`, within its tolerance and named differences.

## The systems

Existing fixtures, reduced. Equivalence needs reproducibility, not convergence, so `rbath` and `nstep` are cut until
a case runs in seconds on one core. What each was cut to, and from what, is in the generator next to the numbers.

| | system | reaches | wall (1 core) |
|---|---|---|---|
| S1a/b | NV + natural-abundance 13C, CCE-2, ensemble and sampled bath | the ordinary CCE path, point-dipole hyperfine | 0.3 s |
| S2a | the same, gCCE-2 | gCCE | 0.35 s |
| S2b/c | NV + P1 10 ppm, CCE-2 and gCCE-2 | `Defect` with `n_nuclear_spins = 0`, `n_configurations = 12`, `Defect.axis_file` | 0.3 / 1.0 s |
| S3 | P1 bath, `Defect` `n_nuclear_spins = 1` (14N), `n_configurations = 4` | `gyros` `quadrupole_moments` `hyperfine` `efg` `zfs` `detuning`, state and extra-state files | 1.4 s |
| S4 | hBN V_B, `hyperfine` `{"dipole": "dft", "contact": true}` + `quadrupole` `"plain"` | both DFT tensor readers, the vertex boundary, `hyperfine.frame` / `quadrupole.frame` | 0.44 s |
| S4b | the same over `Afile` instead of `Afile_vertex` | `CheckBD_Range`, the axis-aligned MinDif/MaxDif box | 0.43 s |
| S5 | two NV, gCCE, explicit `intmap` self/ZFS tensors | `intmap` `frame`, multi-qubit rotation (about the first qubit) | 0.29 s |

S3's fixture is written for 24 ranks at `rbath` 850 (469 spins, 80 s). At `rbath` 360 it keeps 45 spins and 518 pair
clusters, which still reaches every Defect code path. S4's is cut from `rbath` 20 to 6: 3897 spins and 121 s become
101 spins and 0.44 s, well inside the region the A and Q files cover, so both readers stay on their file path rather
than falling back.

S4 and S4b differ only in which boundary test the tensor file selects, and a file picks exactly one: `Afile_vertex`
carries `v1..v8` vertices and reaches `CheckBD_vertex`; `Afile` carries `MinDif[A]` / `MaxDif[A]` and reaches
`CheckBD_Range`. Without S4b the axis-aligned box would never run. That box is the one thing that makes a rotation of
the tensor file's own coordinates impossible.

S4b needs `hyperfine.missing_row = "point_dipole"` and no quadrupole: both boxes are wider than the set of atoms their
files list, so a spin can sit inside the box with no row of its own. The HF reader has a flag for that and the QD
reader does not, which leaves the QD range boundary uncovered. So are the quadrupole structures `strained` and
`bubble`, which upstream also records as untested.

## What is compared

The **raw** `_noDiv` / `_wiDiv`. Not the error-corrected file and not the ensemble average: both smooth, and a real
difference could hide under them.

`errcorr.sh` is the other half of that decision. It runs the analyzer on each baseline and reports `max|L|` before and
after correction, so the suite says whether the systems it compared are numerically sound. Most are. `S2a` is not:
gCCE-2 on that 13C bath produces `max|L| = 4.2` at two points, which the correction takes back to 1. That is the
ordinary CCE division divergence, present in v1.1.0 exactly as here, and it does not disqualify the system as an
equivalence fixture. It does mean the raw file is not physics.

## Unit conversions

`gyros`, `eqs`, `hypf`, `efg`, `zfs` and `detuning` are each given twice: once in the unit whose scale factor is
exactly 1, where byte identity with the legacy array is required, and once in a second unit, where it is not. Three
of the four second units come back byte-equal anyway. `efg` in `V/angstrom^2` lands 1e-10 away, one unit in the last
printed place of a file written to ten decimals.

## Refusals

`../check_refusals/run.sh` re-checks the refused combinations on the physical fixtures rather than the toy ones. It
requires two things of each: the run stops with the expected message, and neither DFT tensor reader announced itself
first. The second half would regress quietly if a guard drifted later into the run: the refusal would still happen,
but only after a file had been read in a frame it was not written in.
