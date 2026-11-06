# check_master_equation/

Checks the master-equation methods (`method` meCCE and megCCE, see `docs/methods/master-equation.md`) against their
limits at rate 0, an outside code (PyCCE) and closed forms. The main system is the `adding_MECCE` branch's test: an
S = 1/2 electron qubit in an electron-spin bath of 1e16 cm⁻³, B = 300 G along z, order 2, Hahn echo, ensemble.
`run.sh` writes every case from one template.

```
bash test/code_verification/check_master_equation/run.sh     # about 85 s; bin/main.out or CCEX_BIN, CCEX_MPI; 1 rank
```

## What it checks

`compare.py A B --tol T` takes max \|A − B\| on the shared time grid of `_noDiv`. It exits 2 when a file cannot be
read or the grids differ. The jump rate is 1.5915e-4 MHz (1 rad/ms).

| Case | Compared with | Tolerance |
|---|---|---|
| `test_superoperator_unit.cpp`: superoperators and jump operators | unit test | no failure |
| meCCE at rate 0: Hahn (also `_wiDiv`), CPMG 2 and 3, a pure state of `p1_cce2_single`, `natab_knight_shift` | cCCE (`knight` true for the last) | 1e-6 |
| meCCE at rate 0, Ramsey | cCCE | 1e-8 |
| meCCE with bath jumps | the branch output in `reference/` | 1e-6 |
| meCCE with bath jumps | PyCCE | 1e-7 |
| meCCE, one bath spin on the z axis, FID, two rates | PRL 132, 250401, Eq. 6 | 1e-8 |
| megCCE at rate 0, also two NV qubits (S5) | gCCE | 1e-6 |
| gCCE, megCCE with bath jumps, megCCE with bath and qubit `z` jumps | conj(PyCCE) | 1e-7 |
| megCCE order 0, qubit `z` or `-` jump | exp(−g t/2) | 1e-9 |

The PyCCE tolerance is about five times the floor the input sets (bath coordinates to 5 decimals, `H_BAR` in
the ninth digit). PyCCE and CCEX are complex conjugates after an odd number of π pulses, hence `conj`. Only the
branch's meCCE output is a reference, because its gCCE keeps the qubit gyromagnetic ratio and deltat in float.

## Files

| File | Contents |
|---|---|
| `run.sh` | Writes the cases from the template and runs the comparisons. |
| `compare.py` | Max \|A − B\| on a shared grid (`--conj-b`, `--subset-b`, `--scale-exp RATE`). |
| `test_superoperator_unit.cpp` | The unit test of `src/calculation/superoperator.cpp`. |
| `inputs/` | The branch's qubit, bath and gyro files, and `bath_1e_z` (one electron at (0, 0, 500 Å)). |
| `reference/branch-2c63eac52/` | The branch's meCCE output, with `MD5SUMS`. |
| `reference/pycce/` | PyCCE coherences of the same system, how they were made (`README.md`), `MD5SUMS`. |
