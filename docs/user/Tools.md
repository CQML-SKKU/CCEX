# Tools

What makes the inputs CCEX reads: three tools in `example/tools/`. A fourth,
`AutoConvergence/`, is the automatic convergence: an input with an `auto_convergence`
block run as `main.out -f input.json` (without mpirun) chooses the time window, rdip,
rbath and the configurations itself; see `example/tools/AutoConvergence/README.md`.

| Tool | Makes | Needs |
|---|---|---|
| `BathGenerator/` | a bath file, site positions and isotopes, from a POSCAR unit cell | numpy, psutil, mpi4py with an MPI |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files, from a VASP OUTCAR or a Quantum ESPRESSO output | numpy, pandas, psutil |
| `Runner/` | the job folders of the tutorials on the cluster, one per value of a variable | the idm Slurm cluster |

## BathGenerator

Repeats the unit cell to the requested size and draws the isotope of every site
with the probabilities in a `prob` file; options for a 2D material, a defect, an
inserted relaxed cell and a deformation. The isotopes are drawn without a fixed
seed, so two runs give the same sites with different isotopes.

```
mpirun -n 4 python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D
```

`-h` lists everything; `do_BathGenerator` has three worked cases. The baths of
the tutorials were made with this copy. A newer generator lives at
<https://github.com/CQML-SKKU/BathGeneratorSpinCluster>.

## SpinTensorReader

Reads the hyperfine and electric-field-gradient tensors a DFT code wrote for
every nucleus around the defect and writes them in the layout the reference
describes under [hyperfine.file](InputFile_General.md#tag_hyperfine_file) and
[quadrupole.file](InputFile_General.md#tag_quadrupole_file): the species lines,
the range the file covers and one row per atom with its position relative to
the defect and the nine components.

```
python3 main.py VASP POSCAR DEFECT --A OUTCAR --Q OUTCAR [--spin S] [--rangeformat vertex]
python3 main.py QE   POSCAR DEFECT --A hyperfine.out --spin S --Q efg.out [--rangeformat vertex]
```

Since 2026-10-13 every hyperfine file (`Afile`) is in tensor format v3: line 1 `VASP` or
`QUANTUMESPRESSO`, line 2 `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ..., S = ...`,
and the numbers are the hyperfine of the spin-S defect, which CCEX reads as written. VASP's OUTCAR is already per
S (VASP divides by the total spin), so the numbers are the OUTCAR's and S comes from the OUTCAR's total magnetization
(a `--spin` given must agree). Quantum ESPRESSO's GIPAW `hyperfine.out` is the S = 1/2 coupling, so `--spin S` is
required and the file holds it x 1/(2S). `--version v1` (the OUTCAR x 2, which CCEX read right only untagged and
for S = 1) is refused and `--fidentify` is gone. The quadrupole file (`Qfile`) does not depend on the spin and is
unchanged; the spin is read (and `--spin` checked) only when an `Afile` is made, so `--Q` alone takes an OUTCAR with
magnetization 0 or none (an EFG run need not be spin polarized; refused until review 20261013-180049 T1).

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`), `main.py` the
command line (`-h`). `examples/01..05` are five worked cases with their inputs
and the tensor files they produced, each with a `do` that reproduces it; the
tensors of the h-BN tutorial (07) came from it.

## Runner

Makes the Slurm job folders the tutorials were recorded with: one per value of
a variable (`order`, ppm, `rbath`, ...) from the templates in
`MkCCEJobs/mkjobs.sngapp/`, driven by a tutorial's `domkjobs` and `opt/` files
(`./domkjobs -h`). It is written for the group's idm cluster and its module
names; on another machine the tutorials' `run.sh` do the same job by hand.
