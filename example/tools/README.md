# example/tools/

Programs that prepare what CCEX reads, plus the automatic convergence runner and the input-key tools.
[Tools](../../docs/user/Tools.md) describes each in full.

| Tool | Makes | Run | Needs |
|---|---|---|---|
| `BathGenerator/` | a bath file (positions and isotopes) from a POSCAR unit cell, with options for 2D materials, a defect, an inserted relaxed cell and deformation | `mpirun -n N python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D` (`-h` lists the options; `do_BathGenerator` has four worked cases) | numpy, psutil, mpi4py |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files, from a VASP OUTCAR or a Quantum ESPRESSO output | `python3 main.py VASP\|QE POSCAR DEFECT --A OUT --Q OUT [--spin S] [--rangeformat vertex]` (`--spin` is required for QE `--A`); each of `examples/01..05` has a `do` script | numpy, pandas, psutil |
| `IsotopeTable/` | `src/input/isotope_table.cpp`, the built-in table of 154 nuclear isotopes (spin, gyromagnetic ratio) used for species that no gyrofile names | `python3 example/tools/IsotopeTable/make_isotope_table.py` from the repository root | python3 |
| `StateFiles/` | the `bath_state_file`, `defect_axis_file` and `defect_state_file` files an input names, drawn reproducibly from `configuration.seed` with the same rule the automatic convergence uses | `python3 example/tools/StateFiles/make_state_files.py INPUT [--force]` | numpy |
| `AutoConvergence/` | the automatic convergence: an input with an `auto_convergence` block finds the time window, rdip, rbath and the number of configurations, then runs the production | `mpirun -n N main.out -f input.json`; see [its README](AutoConvergence/README.md) | numpy, scipy, matplotlib |
| `InputKeys/` | `convert_input_keys.py` rewrites an input written with older keys (comments are kept); `set_input_keys.py` writes a copy of an input with some keys changed | `python3 example/tools/InputKeys/convert_input_keys.py FILE...` (`--check` only lists); `python3 example/tools/InputKeys/set_input_keys.py IN OUT KEY=VALUE ...` | python3 |

The tests of the Python tools: `python3 -m unittest discover -s example/tools/InputKeys/tests`, the same for
`StateFiles/tests`, and `python -m pytest -q tests` in `AutoConvergence/`.

## BathGenerator

A newer version of this generator is at <https://github.com/CQML-SKKU/BathGeneratorSpinCluster>. This copy is the one
the tutorial baths were made with. The isotopes are drawn without a fixed seed, so two runs give the same sites with
different isotopes. With some `mpi4py` builds (seen with anaconda on WSL) the interpreter aborts at exit with
"double free or corruption"; the bath file is already complete by then.

## SpinTensorReader

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`) and `main.py` the command line. `examples/` holds five
worked cases with their inputs and the tensor files they produce; the `do` script in each reproduces them
(`04.Diamond_N2p/do` takes the Jahn-Teller index and the case, for example `bash do 4 1-2`).

Every `Afile` is in tensor format v3. Its first line is `VASP` or `QUANTUMESPRESSO`, the second
`# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ..., S = ...`, and the values are the
hyperfine of the spin-S defect, which CCEX reads as written
([hyperfine file](../../docs/user/InputFile_General.md#tag_hyperfine_file)). For VASP the values are the OUTCAR's and S
comes from its magnetization (a `--spin` given must agree). Quantum ESPRESSO's GIPAW writes the S = 1/2 coupling, so
`--spin S` is required and the file holds `hyperfine.out` × 1/(2S). The older format (`--version v1`) is refused.
`Qfile`s (electric-field gradients) do not depend on the spin.

Every value is written with a space in front of it, so values of 10 characters or more (such as −123.45678) stay in
their own field; `tests/check_wide_fields.py` checks that CCEX reads them as written.

`runningScripts_Exspin/` holds two example runs (N2p, NVHn) with their OUTCARs; its `do_summary` uses absolute paths of
the author's machine.
