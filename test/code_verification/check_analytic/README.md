# check_analytic/

Compares CCEX with closed forms and exact answers rather than with another CCEX run, so it catches errors that every
CCEX-to-CCEX comparison would share (a wrong constant, a wrong time axis). It is small and machine-independent, and CI
runs it.

```
bash test/code_verification/check_analytic/run.sh      # about 15 s; bin/main.out or CCEX_BIN, CCEX_MPI
```

## What it checks

| Case | Compared with | Tolerance |
|---|---|---|
| One ¹³C beside an NV, free induction and Hahn echo, order 1 | the closed-form coherence (`analytic.py`) | 2e-10 (the 10 written decimals) |
| The same with a gyromagnetic ratio 1 % off | must fail (a self-test of the check) | — |
| `hfmedi` at zero field: degenerate qubit levels coupled or not | refused when coupled, finite otherwise; an energy shift must not change the decision | — |
| `hamiltonian.cCCE.delta_A`, one ¹³C 2.55 Å from the NV | gCCE order 1 (exact for one spin) | 1e-2 at 5 G, 5e-2 with a tilted field |
| `bath_bath_mean_field` "off", two spins | the product of each spin alone | 1.5e-10 |
| gCCE with a non-axial ZFS (E = 100 MHz), as written and turned by 30° | the exact 6 × 6 answer (`gcce_zfs_exact.py`) | 2e-8 (the GHz qubit phase amplifies rounding) |

The exit status is 0 only when every check passes.

## Files

| File | Contents |
|---|---|
| `run.sh` | runs the cases and prints one PASS or FAIL per check |
| `analytic.py` | the closed forms of the free-induction and Hahn-echo coherence (in its docstring) and the comparison |
| `gcce_zfs_exact.py` | writes the non-axial ZFS input and its exact answer |
| `inputs/` | the qubit, bath and gyro files of the cases |
