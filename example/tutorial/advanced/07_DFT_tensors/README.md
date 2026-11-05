# Advanced 07: DFT hyperfine and quadrupole tensors

Folder: `example/tutorial/advanced/07_DFT_tensors/`

This topic takes the hyperfine and quadrupole tensors from DFT instead of the point dipole of the positions.
`example/tools/SpinTensorReader` reads a VASP OUTCAR (or a Quantum ESPRESSO output) and writes two files. `Afile` holds
the hyperfine tensor of each atom of the cell (MHz, tensor format v3), and `Qfile` the electric-field gradient of each
atom (for spins I > 1/2). An input names them in `hamiltonian.common.hyperfine` (`dipole` dft, `contact` true,
`file`) and `hamiltonian.common.quadrupole` (`file`).

The rows are matched to the bath spins by position from the defect, so the bath has to be built on the lattice of the
DFT cell. A bath spin inside the cell without a row stops the run, and spins outside the cell take the point dipole
(`missing_row` point_dipole extends that to every missing row). The cell of `makeTensors_hBN` is not the one of the
tutorial bath, so `cCCE2_Hahn_tensors` uses the files kept with
[Advanced 06: V_B⁻ in hexagonal BN](../06_hBN/README.md).

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

1. [makeTensors_hBN](#dft_make): how the tensor files are made from a VASP run.
2. [cCCE2_Hahn_tensors](#dft_effect): what the DFT tensors change.

<a id="dft_make"></a>
## 1. makeTensors_hBN: making the tensor files

- **What it shows:** how the tensor files are made. `run.sh` runs `example/tools/SpinTensorReader/main.py` on the VASP
  run of V_B⁻ kept in `example/tools/SpinTensorReader/examples/02.hBN_VB_VASP`, with its relaxed cell
  `VASP_CONTCAR_flat`, the defect named `VB`, and `VASP_OUTCAR_flat` for both `--A` and `--Q`. It checks that
  `results/Afile` and `results/Qfile` are the example's `Afile_v2` and `Qfile_v2`, byte for byte.
- **What the files hold:**
  - The first line names the DFT code (`VASP`).
  - The second line, which CCEX reads, says what the values are
    (`# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = VASP OUTCAR, S = 1`).
  - Then come the gyromagnetic g-factors and the corners `v1` .. `v8` of the cell (Å from the defect, the range of the
    rows).
  - Last comes one row per atom, with the position from the defect, the isotropic (contact) term and the nine
    components, in MHz. Here there are 53 B and 54 N, and the three ¹⁴N next to the vacancy have a contact term of
    121 MHz.
  - `Qfile` has the same layout with the field gradient (Hartree/Bohr²).
- **How to run:** `sh run.sh`. It takes under a second and needs python3 with numpy, pandas and psutil, but not CCEX.
- **What to change next:** your own OUTCAR. `--addhfcore False` leaves out the core term of the contact interaction,
  and for Quantum ESPRESSO give `--spin S` ([Tools](../../../../docs/user/Tools.md)).

<a id="dft_effect"></a>
## 2. cCCE2_Hahn_tensors: what the tensors change

- **What it shows:** a cCCE2 Hahn echo at 3 T over five configurations (`k_range`) with the hyperfine and quadrupole
  tensors (`ccein_hf3_qd2.jsonc`), with the quadrupole only (the hyperfine a point dipole, `hf0_qd2`) and with the
  hyperfine only (`hf3_qd0`).
- **Keys that matter:** `hamiltonian.common.hyperfine` {`dipole` dft, `contact` true, `file`},
  `hamiltonian.common.quadrupole` {`structure` plain, `file`}, `bfield` 30000 G, `rbath` 20 Å, `rdip` 10 Å, `deltat`
  3 µs × `nstep` 50 and `k_range` [1, 5].
- **How to run:** `NP=4 bash run.sh`. It takes 8 s and 0.07 GB per rank.
- **Expected result:** `expected.png`. The T2 is 56.5 µs with both tensors, 56.6 µs with the quadrupole only and
  34.0 µs with the hyperfine only.
- **What to change next:** `k_range` [1, 50], or `rbath` 30 to check convergence.

![cCCE2_Hahn_tensors: expected result](cCCE2_Hahn_tensors/expected.png)

Inputs: [`ccein_hf0_qd2.jsonc`](cCCE2_Hahn_tensors/ccein_hf0_qd2.jsonc), [`ccein_hf3_qd0.jsonc`](cCCE2_Hahn_tensors/ccein_hf3_qd0.jsonc), [`ccein_hf3_qd2.jsonc`](cCCE2_Hahn_tensors/ccein_hf3_qd2.jsonc)
