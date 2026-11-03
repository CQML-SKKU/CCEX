# check_analytic/: CCEX against a closed form

    bash test/code_verification/check_analytic/run.sh     # bin/main.out, or CCEX_BIN; CCEX_MPI as for the others

Every other suite under `code_verification/` compares CCEX with CCEX (another input form, another frame, a recorded
run), so an error common to both sides passes there: a wrong time scale in the free evolution, a wrong constant in
the dipole formula. This one does not. It is small and machine-independent, so the CI workflow runs it.

Exit status 0 only if all fourteen checks pass. Everything written goes to `_work/`.

## One ¹³C, closed form

One ¹³C beside the NV (`inputs/`: qubit at (10, 10, 10) Å, the ¹³C at (13, 11.5, 14) Å, B = 500 G along z), order 1,
which is exact for one bath spin, ensemble average over its state, 400 points to 0.2 ms. `analytic.py` writes the
point-dipole hyperfine out itself and compares the free-induction and Hahn-echo coherence with their closed forms
(its docstring has them). The tolerance, 2e-10, is the ten decimals CCEX writes. A self-test runs the echo with a
gyro 1 % off and requires the check to fail.

`analytic.py` exits 1 on a numerical mismatch and 2 when the output is missing, cannot be read or does not hold the
requested grid (400 points, k * 0.0005 ms). The self-test counts only a clean CCEX run that `analytic.py` rejects
with 1. The inputs use the `Qubit.file` form with the branch keys `Qubit.alpha_state`/`Qubit.beta_state`
`{"ms": [m]}` and the ZFS as a `Qubit.intmap` self entry (CCEX refuses the old `Qubit.zfs`).

Measured on WSL (Intel MPI, Make and CMake builds): deviation 5.0e-11 for both pulse counts. The 1 % wrong gyro
gives 0.37. Reading the time column as tau instead of the total time gives 0.37 (echo) and 1.96 (free induction).

## `hfmedi` at zero field

More runs, not closed forms:

- With an axial ZFS the levels +1 and -1 are degenerate but not coupled by S. The run has to finish with finite values.
- Without ZFS the levels 0 and +1 are degenerate and coupled. The run has to be refused.
- The same degenerate levels with a zero hyperfine (gyro 0) have no coupling. The run has to finish with finite
  values (the check must look at the coupling, not at S alone).
- ZFS `diag(0,0,1)` and the same plus 1e13 I (an energy shift only) have to get the same decision, both accepted and
  finite (the degeneracy tolerance follows the spread of the levels, not the largest level).

## `hamiltonian.cCCE.delta_A`

The second-order δA(m_s)·I of each bath spin, against gCCE order 1, which keeps the whole S·A·I and is exact for one
spin: one ¹³C 2.55 Å from the NV (`inputs/bath_dA`), ensemble, 400 points to 1 ms. The Hahn echo is judged by
max |L − L_gCCE| and the free induction by max ||L| − |L_gCCE|| (the constant second-order shift of the qubit
frequency is not in the cCCE).

| field | without delta_A (echo / FID) | with delta_A | tolerance |
|---|---|---|---|
| B = 5 G along z | 0.81 / 0.65 | 4.8e-3 / 5.3e-3 | 1e-2 |
| B = (2, 0, 5) G | knight alone 0.56 / 0.81 | knight with delta_A 2.0e-2 / 2.9e-2 | 5e-2 |

The deviation without delta_A has to be at least ten times the one with it (27 to 170 measured; four checks). What
is left is the qubit Zeeman energy left out of the denominators (as for knight) and, with a field across the axis,
knight's own higher orders.

## `hamiltonian.common.bath_bath_mean_field` "off"

Two electron-like spins off the axis (`inputs/bath_off`, one fixed state each, order 1) give the product of each spin
alone in its bath (1.5e-10, the written decimals), while the default "z" differs by 0.44 (one check).

## gCCE with a non-axial ZFS

An S = 1 qubit (D = 2.87 GHz, E = 100 MHz) and one ¹³C (`gcce_zfs_exact.py`, which writes the input and the exact
6 × 6 answer), B = 5 G, ensemble, free induction, as written and turned by 30° about z. The turned run keeps the
m<sub>s</sub> states, i.e. another preparation phase relative to the ZFS axes, and its exact |L| differs by up to
3.6e-5. CCEX has to follow the exact answer in both (5.0e-10 and 7.6e-10). The tolerance is 2e-8 because the GHz
qubit phase amplifies rounding. Two checks.
