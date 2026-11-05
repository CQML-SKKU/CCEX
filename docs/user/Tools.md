# Tools

The tools in `example/tools/` prepare the inputs of CCEX and run the automatic convergence. Each section below covers
one tool.

1. [BathGenerator](#tag_tool_bathgenerator): makes a bath file, the sites and their isotopes, from a POSCAR unit cell.
   Needs numpy, psutil and mpi4py.
2. [SpinTensorReader](#tag_tool_spintensorreader): makes the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files
   from a VASP OUTCAR or a Quantum ESPRESSO output. Needs numpy, pandas and psutil.
3. [IsotopeTable](#tag_tool_isotopetable): regenerates the built-in isotope table of CCEX.
4. [InputKeys](#tag_tool_inputkeys): converts old inputs to the current keys, and writes a copy of an input with some
   keys changed.
5. [StateFiles](#tag_tool_statefiles): writes the state and axis files an input names. Needs numpy.
6. [AutoConvergence](#tag_tool_autoconvergence): the runner of the automatic convergence, which `main.out` starts.
   Needs Python 3.10 or newer with numpy, scipy and matplotlib.

<a id="tag_tool_bathgenerator"></a>
## 1. BathGenerator

BathGenerator repeats the unit cell to the requested size and draws the isotope of every site with the probabilities
in a `prob` file. It has options for a 2D material, a defect, an inserted relaxed cell and a deformation. It sets no
random seed itself, so two runs give the same sites with different isotopes; `make_natab_baths.sh` (below) fixes the
seed.

```
mpirun -n 4 python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D
```

<a id="tag_prob_file"></a>
**The `prob` file** has three lines for every isotope that can sit on a site:

1. its name, the element symbol with the mass number;
2. its nuclear spin;
3. its abundance in %.

The isotopes of one element are drawn together. The element is the name without its digits. Every atom name of the
POSCAR except the `--Defect` one must have an isotope in the file. When the abundances of an element add up to less than 100 %, the rest of its sites hold no spin and
are left out of the bath. More than 100 % is refused. The natural-abundance diamond file of the tutorials
(`example/tutorial/basic/01_first_run/bath/genbath/prob_natab`) is

```
13C
0.5
1.1
```

so 1.1 % of the carbon sites hold a ¹³C. For h-BN (`example/tutorial/advanced/06_hBN/bath/prob_hBN`) the file lists ¹⁰B (3,
19.9), ¹¹B (3/2, 80.1), ¹⁴N (1, 99.632) and ¹⁵N (1/2, 0.368).

| Option | Meaning |
|---|---|
| `-o OUT` | the output bath file; the defect position goes to `OUT_defect` |
| `-i POSCAR` | the unit cell of the pure crystal (or, with `--convert`, a bath file) |
| `-p PROB` | the `prob` file |
| `--size L` | the edge of the bath cell in Å (a cube; with `--opt2D` the z extent is the unit cell's) |
| `--Defect NAME` | the atom name of the defect site: of the sites with that name, the one closest to `--Defect_center` is the defect, the centre of the bath, written to `OUT_defect` |
| `--Defect_center x y z` | the fractional point (each between 0 and 1) the defect site is chosen near (default 0.5 0.5 0.5) |
| `-I POSCAR` | a relaxed supercell with the defect, inserted around it |
| `-D POSCAR ... --axis AX ...` | deformed cells and their axes (`x`, `y`, `z`: plane shape; `xy`, `xz`, `yz`: line shape) |
| `--opt2D` | a two-dimensional bath (one layer, as h-BN) |
| `--convert` | convert a bath file (`-i`) to a POSCAR (`-o`) to look at it |
| `--test` | also write a POSCAR of the bath, to check it |

`-h` lists every option. `do_BathGenerator` has four worked cases: three make a bath, one converts a bath to POSCAR.

The baths of the tutorials were made with this copy. `example/tutorial/make_natab_baths.sh N OUTDIR [SIZE [RMIN]]`
calls it for natural-abundance ¹³C around an NV. It makes N configurations with the seeds 1 ... N, so the same call
gives the same files, and with RMIN it leaves out the ¹³C closer than RMIN Å to the NV. The basic tutorials use it.
A newer generator is at <https://github.com/CQML-SKKU/BathGeneratorSpinCluster>.

<a id="tag_tool_spintensorreader"></a>
## 2. SpinTensorReader

SpinTensorReader reads the hyperfine and electric-field-gradient tensors a DFT code wrote for every nucleus around the
defect. It writes them in the layout described under [hyperfine.file](InputFile_General.md#tag_hyperfine_file) and
[quadrupole.file](InputFile_General.md#tag_quadrupole_file): the species lines, the range the file covers, and one row
per atom with its position relative to the defect and the nine components.

```
python3 main.py VASP POSCAR DEFECT --A OUTCAR --Q OUTCAR [--spin S] [--rangeformat vertex]
python3 main.py QE   POSCAR DEFECT --A hyperfine.out --spin S --Q efg.out [--rangeformat vertex]
```

`--rangeformat` chooses how a file states the region it covers, which is the DFT supercell around the defect:

| Value | Writes |
|---|---|
| `vertex` (default) | the eight corners of the supercell relative to the defect (lines `v1` to `v8`, Å) |
| `coord` | the smallest and largest relative coordinates along x, y and z (lines `MinDif[A]`, `MaxDif[A]`) |

A bath spin inside that region must have a row in the file
([hyperfine.missing_row](InputFile_General.md#tag_hyperfine_missing_row)).

Every hyperfine file (`Afile`) is in tensor format v3:

- line 1 is `VASP` or `QUANTUMESPRESSO`;
- line 2 is `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ..., S = ...`;
- the numbers are the hyperfine of the spin-S defect, which CCEX reads as written.

The two DFT codes differ:

- VASP's OUTCAR is already per S, since VASP divides by the total spin. The numbers are the OUTCAR's, and S comes from
  the OUTCAR's total magnetization. A `--spin` given must agree with it.
- Quantum ESPRESSO's GIPAW `hyperfine.out` is the S = 1/2 coupling. `--spin S` is required, and the file holds the
  coupling × 1/(2S).

The quadrupole file (`Qfile`) does not depend on the spin. The spin is read, and `--spin` checked, only when an `Afile`
is made, so `--Q` alone also takes an OUTCAR without magnetization.

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`), and `main.py` the command line (`-h`). `examples/01..05`
are five worked cases with their inputs and the tensor files they produced, each with a `do` that reproduces it. The
tensors of the h-BN tutorials came from this tool, and
[Advanced 06: V_B⁻ in h-BN with DFT tensors](../../example/tutorial/advanced/06_hBN/README.md) runs it.

<a id="tag_tool_isotopetable"></a>
## 3. IsotopeTable

`make_isotope_table.py` writes `src/input/isotope_table.cpp`: the nuclear isotopes, with spin and gyromagnetic ratio in
rad/ms/G. CCEX uses them for a bath species the gyrofile does not name, or for every species when the input has no
gyrofile ([gyrofile](InputFile_General.md#tag_gyrofile)).

The source is the EasySpin isotope database in
`example/tools/SpinTensorReader/spHamilReader/spin_database/isotopedata.txt`, with
$\gamma = g_n \mu_N / \hbar$. To regenerate the table after the database changes, run this from the repository root
and rebuild CCEX:

```
python3 example/tools/IsotopeTable/make_isotope_table.py
```

<a id="tag_tool_inputkeys"></a>
## 4. InputKeys

Two scripts for input files:

```
python3 example/tools/InputKeys/convert_input_keys.py FILE...         # rewrite inputs written with old keys, in place
python3 example/tools/InputKeys/convert_input_keys.py --check FILE... # list only; exit 1 while old keys remain
python3 example/tools/InputKeys/set_input_keys.py IN OUT KEY=VALUE ...
```

`convert_input_keys.py` reads `.json` and `.jsonc` inputs, keeps the comments (`//`, `#` or `!`) and lists what it
changed. What it maps is the table of [old names](InputFile_OldNames.md).

`set_input_keys.py` writes a copy of `IN` with keys set:

- `KEY` is a key or a dotted path into a block (`configuration.s_range=[1,5]`; the block has to exist in `IN`). It is
  matched without regard to case.
- `VALUE` is read as JSON when it is JSON (`50`, `true`, `["bath_1"]`), and as a string otherwise.
- `null` removes the key (`configuration.s_range=null`: the ensemble).
- `OUT` is plain JSON; the comments are not kept.

The tutorials use it to run one input with other values, for example in the scans of
[Basic 03](../../example/tutorial/basic/03_convergence/README.md).

<a id="tag_tool_statefiles"></a>
## 5. StateFiles

`make_state_files.py` writes the state and axis files a CCEX input names. Sampled states and Defect axes then come from
files made by a fixed rule, in place of CCEX's own random numbers. The files can also serve other codes and other runs.

```
python3 example/tools/StateFiles/make_state_files.py INPUT [--force]
```

It reads `bathfile`, `k_range`, `s_range`, `seed`, `path` and `gyrofile` from the `configuration` block, and the
Defect species from `Defect.species`. For every configuration k and state s it writes the files of the keys below that
the input names ([file formats](InputFile_General.md#tag_configuration_files)). A bath file whose name in the key is
`""` is skipped.

| Key | Each line after the first (one per line of the bath file) | Drawn with numpy |
|---|---|---|
| `bath_state_file` | m of the spin, uniform over its 2S+1 levels (S from the gyrofile, else the built-in isotope table) | `default_rng([seed, k, s])` |
| `defect_axis_file` | the axis state, uniform over 1..`n_axis_states` of its Defect species; 0 for another species (CCEX reads an axis file only when every spin is a Defect species with axis states) | `default_rng([seed, k, 0])` |
| `defect_state_file` | the m of each on-site spin of its Defect species (uniform over 2I+1, `spins`), space-separated; 0 for another species | `default_rng([seed, k, s, 1])` |

The first line of each file is the first field of its bath file's first line. `seed` is the number given, its entry for k when it
is a list, or 1 when absent; without `k_range`, k is 0. A file that exists is overwritten only with `--force`.

The automatic convergence makes the files it needs with the same functions, so the two give the same states.
[Basic 02](../../example/tutorial/basic/02_bath_and_states/README.md) uses this tool in section 2.

Tests: `python3 -m unittest discover -s example/tools/StateFiles/tests`. One test runs `bin/main.out` on the files and
checks that CCEX reads every one.

<a id="tag_tool_autoconvergence"></a>
## 6. AutoConvergence

AutoConvergence is the runner of the automatic convergence. `main.out` starts it when the input has an
`auto_convergence` block ([Automatic convergence tags](InputFile_AutoConvergence.md), [how it works](../methods/auto-convergence.md)). `CCEX_PYTHON` picks the Python
interpreter.

Its own commands compute nothing. Run them from `example/tools/AutoConvergence`, or with it on `PYTHONPATH`:

<a id="tag_ac_commands"></a>

```sh
python3 -m autoconv status INPUT.json                       # where a run is
python3 -m autoconv analyze RUN_DIR [RUN_DIR ...] [--config INPUT] [--report-dir DIR]
python3 -m autoconv judge --setting RUN_DIR ... --reference RUN_DIR ... [--report-dir DIR]
```

- `status` prints how the run stands, then `progress.txt`.
- `analyze` judges existing runs (one `RUN_DIR` per configuration) and gives the decay time and the sampling verdict.
- `judge` compares runs of one setting with reference runs of the same configurations, as the rdip and rbath steps
  do, with the default threshold 0.02.
- The `--config` of `analyze` is the input with the `auto_convergence` block, read as the runner reads it.

The state and axis files it makes for a run are written by the functions of [StateFiles](#tag_tool_statefiles).
