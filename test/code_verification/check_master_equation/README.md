# check_master_equation/ — ME-CCE and ME-gCCE (`method` meCCE, megCCE)

    bash test/code_verification/check_master_equation/run.sh     # bin/main.out, or CCEX_BIN; CCEX_MPI as for the others

What the master-equation methods (docs/methods/master-equation.md) must give: their limits at rate 0, an external
reference (PyCCE) and closed forms. One run of the suite is a unit test and 32 CCEX runs on 1 rank (about 85 s on WSL, 50 s
of it the pure-state meCCE run);
it decides with the exit status (0 = every check passed).

## What runs

1. `test_superoperator_unit.cpp`: the Liouville-space superoperators and jump operators (`src/calculation/superoperator.cpp`).
2. The adding_MECCE branch's test system, written by `run.sh` (`mk`) from one template in the current keys: an S=1/2
   electron qubit (`inputs/qubit_espin`), an electron-spin bath of 1e16 cm^-3 (`inputs/bath_espin`, `inputs/gyro_espin`),
   B = 300 G along z, order 2, rbath 500 Å (three spins, one pair), rdip 300 Å, Hahn echo, ensemble (`nstate` 0),
   deltat 0.04 ms × 50 unless the case says otherwise. Jump rate `RATE` = 1.5915e-4 MHz = 1 rad/ms (the branch's value).
3. `cce_knight`, `me_r0_knight`: `check_numerics/cases/natab_knight_shift` as it is and with method meCCE.
4. `bath_1e_z`: one electron at (0, 0, 500 Å), on the qubit's z axis, for PRL 132, 250401 Eq. 6.
5. `cce_single`, `me_r0_single`: `check_numerics/cases/p1_cce2_single` (NV, P1 1 ppm, CCE2, Hahn) as cCCE and as
   meCCE on its first pure bath state (`bath_state_file` key and `-s <prefix> -N 1`; the log must read `<prefix>1`).
6. `gcce_2nv`, `meg_r0_2nv`: `check_refusals`' S5 system (two NV with their ZFS as `intmap` tensors, two 13C,
   `example/tutorial/03_Diamond_multiNV_natab/NV2_13C2`, order 1, one state, the state file as in 5) as gCCE and as
   megCCE, quantity coherence, deltat 1e-4 ms × 20 (S5's 1e-7 keeps the curve within 7e-7 of 1; here |L| reaches 0.43).

| case | method | differences from the template |
|---|---|---|
| `cce`, `cce_ramsey` | cCCE | `cce_ramsey`: `pulse.count` 0 |
| `me_r0`, `me_r0_ramsey` | meCCE | no jump operators; `me_r0_ramsey`: `pulse.count` 0 |
| `cce_cpmg2`, `cce_cpmg3` | cCCE | `pulse.count` 2, 3 (CPMG) |
| `me_r0_cpmg2`, `me_r0_cpmg3` | meCCE | no jump operators, `pulse.count` 2, 3 |
| `me_bath` | meCCE | S+ and S- on the bath spins `e` at `RATE` |
| `me_bath_d01` | meCCE | as `me_bath`, deltat 0.1 × 20 (PyCCE's grid) |
| `me1_eq6_g1`, `me1_eq6_g10` | meCCE | `bath_1e_z`, order 1, rbath 1000, FID, deltat 0.05 × 41, S+ and S- at 1 and 10 rad/ms |
| `gcce` | gCCE | deltat 0.1 × 20 |
| `meg_r0` | megCCE | deltat 0.1 × 20, no jump operators |
| `meg_bath` | megCCE | deltat 0.1 × 20, S+ and S- on `e` at `RATE` |
| `meg_bath_qz` | megCCE | as `meg_bath`, and a qubit jump (`bath_name` `qubit`) `z` at `RATE` |
| `meg0_qz`, `meg0_qminus` | megCCE | order 0 (the 0th cluster alone), one qubit jump (`bath_name` `qubit`) `z` or `-` at `RATE` |

## How it decides (`compare.py`)

`compare.py A B --tol T` is max |A − B| over the shared time grid; `--conj-b` compares with conj(B), `--subset-b` lets B
carry PyCCE's extra end point (t = 2 ms), `--scale-exp RATE` multiplies B by exp(−π·10³·RATE·t). It exits 2 when a file
cannot be read or the grids differ, so a run that did not finish is not a pass. Each check compares `_noDiv` unless it
says otherwise.

| check | tol | measured max \|diff\| | to compare with |
|---|---|---|---|
| superoperator unit | — | FAIL=0 | |
| meCCE rate 0 == CCE (Hahn, noDiv) | 1e-6 | 7.000e-10 | branch summary 5.2e-9 (`test/verification/summary.txt`) |
| meCCE rate 0 == CCE (Hahn, wiDiv) | 1e-6 | 7.000e-10 | |
| meCCE rate 0 == CCE (Ramsey, complex) | 1e-8 | 8.062e-10 | |
| meCCE rate 0 == CCE (CPMG 2, complex) | 1e-6 | 7.000e-10 | |
| meCCE rate 0 == CCE (CPMG 3, complex) | 1e-6 | 1.100e-09 | |
| meCCE rate 0 == cCCE (pure bath state, p1_cce2_single, complex) | 1e-6 | 4.509e-08 | the cCCE with `evolution_update` exact (instead of its default recursive): 1.0e-10 from the default and 4.5e-8 from meCCE |
| meCCE rate 0 == cCCE with `hamiltonian.cCCE.knight` true (natab, Hahn) | 1e-6 | 0.000e+00 | knight moves this case by 0.97 |
| meCCE bath jumps == branch 2c63eac52 | 1e-6 | 1.910e-08 | |
| meCCE bath jumps == PyCCE (S=1/2) | 1e-6 | 3.578e-09 | branch measured 1.25e-8 |
| meCCE 1 spin FID == PRL2024 Eq. 6 (2g < a) | 1e-8 | 2.421e-10 | |
| meCCE 1 spin FID == PRL2024 Eq. 6 (2g > a) | 1e-8 | 2.691e-10 | |
| megCCE rate 0 == gCCE | 1e-6 | 4.428e-09 | branch summary 1.2e-8 (`test/verification/2_verify_megcce/summary.txt`) |
| megCCE rate 0 == gCCE (two NV qubits, S5 system) | 1e-6 | 1.000e-10 | |
| gCCE == conj(PyCCE) (S=1/2) | 1e-6 | 1.886e-08 | model 1.8e-8 (the CCEX algorithm in double precision, numpy) |
| megCCE bath jumps == conj(PyCCE) (S=1/2) | 1e-6 | 5.408e-09 | model 5.1e-9 |
| megCCE bath + qubit z jumps == conj(PyCCE) (S=1/2) | 1e-6 | 4.183e-09 | the bath-jump check above |
| megCCE order 0, qubit `z` == exp(−g t/2) | 1e-9 | 2.654e-10 | |
| megCCE order 0, qubit `-` == exp(−g t/2) | 1e-9 | 2.654e-10 | |

Measured on WSL (g++ `wsl-gnu`, `CCEX_BLAS=mkl`), 2026-10-13 (the CPMG, pure-state and two-qubit rows 2026-10-14). Where the checks can tell a wrong answer apart:

- `gCCE == conj(PyCCE)` is the first comparison of the existing gCCE with an outside code. Without the conjugation the
  difference is 1.66 (max |Im| of the gCCE curve 0.83); for megCCE 0.38. PyCCE measures ⟨Rα|ρ|Rβ⟩ and CCEX ⟨α|ρ|β⟩,
  complex conjugates after an odd number of π pulses (`reference/pycce/README.md`).
  What is left with it, 1.886e-08, is the input and not the compiler or data type (numpy model of this row, 2026-10-14):
  `inputs/bath_espin` keeps PyCCE's `random_bath(seed=2)` coordinates to 5 decimals (2.2e-8 from the bath PyCCE used),
  and `H_BAR` is 1.054571729 (`include/common/utilities.h:72`) against PyCCE's `HBAR_MU0_O4PI` 1.05457172 (1.5e-8);
  both: 1.85e-8.
  One ulp of the electron gyromagnetic ratio moves the curve by 3.2e-9 (Zeeman phase about 1e7 rad at 300 G, t = 2 ms):
  the floor between two double-precision codes on this system.
- `megCCE bath jumps`: the same run without jumps (= gCCE) is 0.69 from conj(PyCCE).
- `megCCE bath + qubit z jumps`: the qubit jump in the clusters with bath spins, C_q ⊗ 1_B, which order 0 (bath
  dimension 1) cannot tell from 1_B ⊗ C_q; put the other way round for the test, this check failed with 5.213e-01 while
  the order-0 checks passed. Without the qubit jump (`meg_bath`) the run is 0.22 from this reference; `meg_bath` times
  exp(−g t/2) is 0.046 from it (in the gCCE the qubit dephasing does not factor out on this resonant system).
- order 0: a 0th cluster without the qubit jump (the branch before 6022d7998, a `nspin > 0` condition on the
  dissipator) gives 1, 0.62 from exp(−g t/2); put back for the test, both checks failed with 6.247e-01. S=1/2 with B ∥ z:
  `z` and `-` (amplitude damping) both damp ρ_αβ at g/2 before and after the π pulse, so L = exp(−g t/2), g = 1 rad/ms.
- tol 1e-9 for order 0: two files of 10 decimals multiplied (the single-file tolerance of `check_analytic` is 2e-10).

Why the branch's own gCCE and ME-gCCE outputs are not references (only its meCCE output is):
the branch keeps the qubit gyromagnetic ratio and deltat in float (`qubit.cpp:249`, `general.h:42`, fixed in this code),
and this system is resonant (qubit and bath electrons at the same Larmor frequency, so the result depends on ω_L t mod 2π),
so its gCCE and ME-gCCE curves are 0.28 and 0.052 from this code's; meCCE uses no qubit Zeeman term and agrees.

## Reference data

- `reference/branch-2c63eac52/ccex_mecce_noDiv`: the branch's meCCE output (`test/verification/1_verify_mecce`).
- `reference/pycce/pycce_{mecce,gcce,megcce}.dat`: PyCCE coherences of the same system on linspace(0, 2, 21)
  (`reference/pycce/README.md`); `pycce_megcce_qz.dat` (the qubit `z` jump added) was generated with PyCCE 1.1.1 on
  2026-10-13, the command in that README.
- `inputs/{bath_espin,qubit_espin,gyro_espin}`: the branch's `test/mecce_test/comparison_megcce_S05/` files.
  Every copied file's md5 against the branch blob is in the `MD5SUMS` beside it (checked when copied), and that of
  `pycce_megcce_qz.dat` as generated.
