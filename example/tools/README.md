# example/tools/

Four tools that make what CCEX reads (`IsotopeTable/` its built-in isotope table), the runner of the automatic convergence, which runs CCEX itself, and the converter of inputs written with the keys before 2026-10-13.

| Tool | Makes | Run | Needs |
|---|---|---|---|
| `BathGenerator/` | a bath file (site positions and isotopes) from a POSCAR unit cell, drawing isotopes with the probabilities in a `prob` file; options for 2D materials, a defect, an inserted relaxed cell, deformation | `mpirun -n N python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D` (`-h` lists everything; `do_BathGenerator` has four cases: three that make a bath and one that converts a configuration to POSCAR) | numpy, psutil, mpi4py with an MPI |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files CCEX reads, from a VASP OUTCAR or a QE output | `python3 main.py VASP\|QE POSCAR DEFECT --A OUT --Q OUT [--spin S] [--rangeformat vertex] ...` (`-h`; `--spin` required for QE `--A`); `examples/01..05` each have a `do` | numpy, pandas, psutil |
| `IsotopeTable/` | `src/input/isotope_table.cpp`, the 154 nuclear isotopes (spin, gyromagnetic ratio) CCEX uses for a species no gyrofile names, from the EasySpin database in `SpinTensorReader/` | `python3 example/tools/IsotopeTable/make_isotope_table.py` from the repository root | python3 |
| `Runner/` | the job folders of the tutorials: one per value of a variable (order, ppm, rbath, ...) from the templates in `MkCCEJobs/mkjobs.sngapp/`, which were driven by a tutorial's `domkjobs` and `opt/` files (removed on 2026-09-23) | `example/tools/Runner/domkjobs MKE` (`-h` lists the options), with `SourceFileName` and `SourceFileDir` set in the environment to the input template; until 2026-09-23 a copy sat in each tutorial calculation and was run there as `./domkjobs MKE` | idm Slurm |
| `AutoConvergence/` | the automatic convergence: `main.out -f input.json` with an `auto_convergence` block (or `mpirun -n N main.out -f input.json` since 2026-10-10: rank 0 runs the runner and every CCEX run gets `mpirun -n N`; the block's `launcher`, when given, comes first) finds the time window, rdip, rbath and the configurations and runs the production on the baths you prepare (it says which baths to add when they are too small or too few); `main.out --auto_status -f input.json` shows where it is. `CCEX_PYTHON` picks the interpreter, `CCEX_AUTOCONV_DIR` the tool folder; it may be a shared installation (owned by you, root or the owner of `main.out`), but neither it nor any folder above it may be writable by others (other users, or a group other than the owner's own private group; `chmod go-w`). Also `python -m autoconv analyze` for existing runs (`AutoConvergence/README.md`) | `main.out -f input.json`; tests `python -m unittest discover -s tests` | numpy, scipy, matplotlib |
| `InputKeys/` | an input in the keys of 2026-10-13 from one written before (`method` cCCE/gCCE with `cluster_mode`, the `hamiltonian`, `pcce`, `bundle`, `pulse` blocks, ...; the table is [Old and new names](../../docs/user/InputFile.md#OldNames)); comments stay, a file already in the new keys is left alone | `python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites in place and lists what changed; `--check FILE...` lists only (exit 1 while old keys remain); `--stdout FILE` prints the new text | python3 |

## BathGenerator

A newer version of this generator lives at
<https://github.com/CQML-SKKU/BathGeneratorSpinCluster>. This copy is the one the
tutorials' baths were made with and is kept for that reason.

The isotopes are drawn without a fixed seed, so two runs give the same sites with
different isotopes. On this WSL the anaconda `mpi4py` aborts at interpreter exit
("double free or corruption") even for a bare `from mpi4py import MPI`; the bath
file is complete before that and the tool is not at fault.

## SpinTensorReader

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`),
`main.py` the command line. `examples/` are five worked cases with their inputs and
the tensor files they produced; `do` in each reproduces them (`04.Diamond_N2p/do`
takes two arguments, the JT index and the case, e.g. `bash do 4 1-2`). On
2026-09-23 the package's internal imports were made package-relative, and every
example was checked to reproduce its tracked `Afile`/`Qfile` byte for byte.

**Tensor format v3 (2026-10-13).** Every `Afile` starts with the tag (`VASP` / `QUANTUMESPRESSO`) and the provenance
line `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`,
and holds the hyperfine of the spin-S defect, which CCEX reads x 1 (`docs/user/InputFile_General.md#tag_hyperfine_file`).
VASP: the OUTCAR numbers (VASP divides by the total spin), S from the OUTCAR magnetization (`--spin`, if given, must
agree). QE: GIPAW's `hyperfine.f90` uses the bare spin density with the S = 1/2 prefactor, so `--spin S` is required
and the file holds `hyperfine.out` x 1/(2S). `--version v1` is refused: it wrote the OUTCAR x 2, which is the defect's
value only for S = 1 and only when CCEX then divides by 2S; with a `VASP` line CCEX read it 2x (the h-BN tutorial file,
and example 03's MoS2 `Afile`). `--fidentify` is gone (the lines are always written). `Qfile`s (EFG, independent of the
spin) are unchanged. The examples were re-run that day: 01 (QE, S = 1) nearest ¹⁴N A_zz 50.4 MHz (PBEsol); 02 (VASP)
48.8 MHz without the core term (`Afile_nocore`, = the removed `Afile_v1` / 2) and 35.4 MHz with it (the default
`--addhfcore True`); 03 and 05 S = 1/2; 04 v2 only (the v1 outputs removed). Every `Qfile`, and every v2 `Afile` apart
from its two new lines, is byte-identical.

**Wide values (2026-10-13).** The `Afile`/`Qfile` fields are `%10.5f` between tabs, and CCEX's tensor reader split a
line at spaces only: a value of 10 characters or more (-123.45678, 1234.56789) ran into the field before it, and the
row was dropped (CCEX then stopped with "doesn't exist in A-file") or read shifted. The writers now put a space in front
of every value, and CCEX splits at any blank, so old files with such values are read too (`runningScripts_Exspin/`'s
NVHn `Afile`s have an iso of -141 MHz). `tests/check_wide_fields.py` (run by the `tensor_format` stage of
`test/code_verification/check_input_forms`) writes values up to -12345.6 with both writers, in the new and in the old
format, and checks that CCEX reads them as written. Examples 01-05 were re-run: every `Afile`/`Qfile` has the same
numbers as before (only the added spaces). The unused copies `v2_vasp/`, `v2_qe/` (not imported anywhere) were removed.
Python 3.12 prints SyntaxWarnings for a few regular expressions; they are harmless.

`runningScripts_Exspin/` holds two example runs (N2p, NVHn) with their OUTCARs and
a `do_summary` written for the author's own machine (absolute paths).

## Runner

Templates in `MkCCEJobs/mkjobs.sngapp/`: `srcjob_ConfigArray`,
`srcjob{Nstate}_ConfigArray` and `srcjob{Nstate}_VariableArray` are current
(intel-23.2 modules, the 64core partition); `srcjob` and `srcjob_VariableArray`
are legacy (icc18, 24core) and say so at their top. The run scripts `srcdo` writes
call `QosChange -j job -q <qos> -n <nodes> -b` and then `SbatchPath --array=... job`,
and `domkjobs RunL` calls `SbatchPath largejob`: the wrappers idm uses now (SbatchPath
records the job path). Until 2026-10-02 they called `shortcut/QosChanger` and `sbatch`
directly. `shortcut/` still holds the old `QosChanger` and `CdMkDir`.
`alphaFix`/`betaFix` are the qubit's `alphams`/`betams` (in `Qubit`, or in the qubit of `Qubit.qubit`); until
2026-10-14 `srcfin_Qubit` wrote them at the top of the input, where nothing reads them, and its qubit was always -1/0.
