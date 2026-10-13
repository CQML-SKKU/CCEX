# check_deer/ — DEER bath pulses (`pulse.bath`)

    bash test/code_verification/check_deer/run.sh     # bin/main.out, or CCEX_BIN; CCEX_MPI as for the others

What `pulse.bath` (docs/methods/deer.md) must do, checked against answers that are known rather than a recorded run,
plus the paths a bath pulse has to reach. One run of the suite is a unit test and 30 CCEX runs on 1 rank (two on
`check_numerics/cases/p1_manual_hahn`, seventeen on 12 spins at order 1 with 41 points, seven at order 2 or 3, four that
must stop before any cluster; about 90 s on WSL); it decides with the exit status (0 = every check passed).

## What runs

1. `test_merge_unit.cpp`: the qubit/bath sequence merge (`src/physical/pulse_bath.cpp`), each pulse on exactly one edge.
2. `pulse.bath: {"count": 0}` on `check_numerics/cases/p1_manual_hahn` gives the output bytes of the input without the
   block.
3. The six cases of the adding_MECCE branch (`example/z.Tests/DEER`, tip 2c63eac52) and five pure-state cases,
   written by `run.sh` from one template (the branch input in the current keys) and the differences below, then judged
   by `check_deer.py`:

| case | `pulse.bath` | state | evolution |
|---|---|---|---|
| `nopulse` | none | ensemble (`nstate` 0) | default (matrix, recursive) |
| `onres` | X 180 at 0.5 of the sequence (on the qubit pi), tone 79.975539 MHz | ensemble | default |
| `offres` | as `onres`, tone 84.98 MHz (5 MHz from every P1 line) | ensemble | default |
| `multi_hit` | as `onres`, tones [79.975539, 200] MHz | ensemble | default |
| `multi_miss` | as `onres`, tones [200, 300] MHz | ensemble | default |
| `nearpulse` | X 180 at 0.505 (0.005 after the qubit pi), tone 79.975539 MHz | ensemble | default |
| `path_matrix` | as `onres` | one pure state (`inputs/state_1`) | matrix, exact |
| `path_matrix_recursive` | as `onres` | pure | matrix, recursive |
| `path_vector` | as `onres` | pure | vector, exact |
| `path_recursive` | as `onres` | pure | vector, recursive |
| `path_nopump` | none | pure | vector, exact |

   Task C2 adds megCCE, `cluster_mode` pcce and multi-spin clusters (the same template; `seed` 3528 in every run, which
   also fixes the pcce partition search):

| case | from | change |
|---|---|---|
| `meg_onres` | `onres` | `method` megCCE, no jump operators (rate 0; `evolution_update` exact, the ME default) |
| `meg_nearpulse` | `nearpulse` | `method` megCCE, rate 0 (the bath pi on an edge of its own) |
| `multi_hit_rev` | `onres` | tones [200, 79.975539] MHz (the resonant one second) |
| `pg_onres` | `onres` | `cluster_mode` pcce, `pcce.partition_size` 1 (inside pgCCE) |
| `pmeg_onres` | `meg_onres` | the same (inside pmegCCE) |
| `o2`, `meg_o2` | `onres`, `meg_onres` | `order` 2 on the pure state `inputs/state_1` (66 pairs) |
| `o2_rev`, `meg_o2_rev` | `o2`, `meg_o2` | the bath file, the axis file and the state file with the 12 spins in reverse order (written by `run.sh`) |
| `o3` | `o2` | `order` 3 (216 triples) |
| `o3_rot` | `o3` | the three files with the spins rotated by 5 (the last five first) |
| `onres_expm` | `onres` | `propagator` expm (the log: evolution matrix, `evolution_update` exact) — final review |
| `pg_o2` | `o2` | `cluster_mode` pcce, `pcce.partition_size` 1 (pcce clusters of two spins with the pump) — final review |

   The pure-state runs pass the state file as the `bath_state_file` key and as `-s <prefix> -N 1`, and `run.sh` checks
   that the log names the file read (`Read StateFile : <prefix>1`).
4. Refusals (exit status and message) and warnings:

| case | input | expected |
|---|---|---|
| `bad_bspin` | `onres` with `pulse.bath.spin` 1 (alphams 1, betams 0) on P1 (S = 1/2) | `pulse.bath.spin 1 differs from the spin 0.5` (before any cluster, not a kron size error) |
| `bad_defect` | `onres` with `pulse.bath.defect` "P2" | `pulse.bath.defect "P2" names no bath spin` |
| `bad_ms_equal` | `onres` with betams 0.5 = alphams | `... must be two different levels` |
| `bad_ms_range` | `onres` with alphams 1.5 for S = 1/2 | `... must be two different levels` |
| `offres`, `multi_miss` | tones that reach no P1 line | run, with `Warning: pulse.bath drives no spin` |

## How it decides (`check_deer.py`)

At order 1 the bath spins do not interact and the echo is a product over the spins the pump drives:

    V(t) = prod_j cos(a_j t min(f, 1-f)),   a_j = 2 pi 52.04 MHz nm^3 (1 - 3 cos^2 theta_j) / r_j^3

with f the fraction of the total time t at which the bath pi lands (f = 0.5: a_j t / 2). Five of the twelve spins are
on the pumped line (configuration 1 of `inputs/avaax_DEERtest`).

| check | tolerance | source |
|---|---|---|
| `nopulse`: \|V − 1\| | 1e-5 | branch `check_deer.py` (`TOL_REFOCUS`, the propagator's precision with a 2.88 GHz ZFS) |
| `offres` = `nopulse`, `multi_miss` = `nopulse` | 1e-12 | controller ruling A3 ulp-drift: a merged sequence (Manual, fractions as edge differences) against the generated one; the branch had 0 |
| `multi_hit` = `onres` | 0 | branch (`TOL_EXACT`): the same merged sequence and the same driven spins |
| `onres` vs closed form at f = 0.5, `nearpulse` at f = 0.505 (real part) | 2e-4 | branch (`TOL_CLOSED`) |
| `nearpulse` ≠ `nopulse` | > 1e-3 | branch: a pump applied twice cancels and would fall onto `nopulse` |
| `_wiDiv` = `_noDiv` (the six ensemble cases, bytes) | — | branch: at order 1 there is no division |
| vector = matrix (exact), recursive = exact (vector and matrix), complex | 1e-9 | `evolution_update` bound of CHANGELOG 2.0.0 |
| `path_vector` ≠ `path_nopump` | > 1e-3 | the pump acts on the vector path (its block product assumes kron(U, 1)) |
| `meg_onres` = `onres`, `meg_nearpulse` = `nearpulse`; `meg_o2` = `o2` (`_noDiv`, `_wiDiv`) | 1e-6 | megCCE at rate 0 is the gCCE: the bound of `check_master_equation` (megCCE r0 vs gCCE) |
| `multi_hit_rev` = `onres` | 0 | as `multi_hit`: every tone is tried, not the first one only |
| `pg_onres` = `onres`, `pmeg_onres` = `meg_onres` | 1e-9 | controller ruling R9: pcce with bath pulses is open only while partition_size 1 gives the standard result |
| `onres_expm` = `onres` | 1e-9 | final review ruling: the expm and eigen propagators of one Hamiltonian, rounding only (`onres_expm` − `nopulse` is 1.34: the pump acts on the expm path) |
| `pg_o2` = `o2` (`_noDiv`, `_wiDiv`) | 1e-8 | as the reorderings below: partition_size 1 gives the standard pairs, the order of the two spins inside a pair may differ (rounding only) |
| `o2` ≠ `path_vector` (`_wiDiv`) | > 1e-3 | the pairs act (otherwise the reordering checks below would compare order-1 values) |
| `o2_rev` = `o2`, `meg_o2_rev` = `meg_o2`, `o3_rot` = `o3` (`_noDiv`, `_wiDiv`) | 1e-8 | rounding only (the 2.88 GHz ZFS phase through 79 or 295 cluster factors); the factor over the measured values is arbitrary |

Everything is judged on `_noDiv` (`<case>_noDiv`, `<case>_state1_noDiv` for the pure states), the whole complex value
except where the closed form (real) is the reference; the order-2 and order-3 comparisons also on `_wiDiv`.

Why both reorderings: the spins of a cluster (`cluster_mode` standard) are kept in the order of their index in the bath
file. Reversing the file
reverses every pair, so a $U_b$ that puts a spin's rotation on the other spin of the pair (only the first spin driven,
say) changes the curve; but a $U_b$ built in reverse kron order gives the same wrong operator in both runs (in each
run spin $a$'s rotation lands on $b$ and $b$'s on $a$). A rotation of the file turns some triples (those that wrap
around) into rotated, not reversed, order, and that $U_b$ then moves the curve. Measured with mutants of
`bathPulseOperator` (Task C2): reverse kron order → only `o3_rot` fails (2.0); first spin of a cluster only → `o2_rev`,
`meg_o2_rev` (1.0) and `o3_rot` (2.9) fail; megCCE without $U_b$ → `meg_onres` (1.3) and `meg_o2` (2.0) fail.

## Changes from the branch test

- One merge tolerance, 1e-9 (the branch matched qubit edges to 1e-12 and bath pulses to 0.01; `nearpulse` exists for
  that bug), and a pulse is never placed twice (`test_merge_unit.cpp`, Task A1).
- `_noDiv` is judged (the branch read `_wiDiv`); the two are compared byte for byte separately.
- The pure-state path checks, the refusals and the off-resonant warning are new.
- Merged against unmerged runs are compared within 1e-12 instead of bit for bit (ruling above); the measured
  difference is 0.
- `multi_hit` uses the branch's tone 79.975539 MHz. The plan's 79.98 MHz is 28.0 rad·kHz from the P1 line, outside the
  10 rad·kHz window: the pump drove nothing and `multi_hit` came out byte-identical to `nopulse` (measured, 2026-10-13).
  `offres` keeps the plan's 84.98 MHz (the branch had 84.975539; both are 5 MHz off every line).
- The input in the current keys (`-I`, `-a`, `-N 0` of the branch `run_test.sh` written into the input):

| branch `ccein.js` | here |
|---|---|
| `avaaxfile` / `-a` | `Defect.axis_file` |
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` |
| Defect `dfname`, `naddspin`, `navaax` | `name`, `n_nuclear_spins`, `n_configurations` |
| Defect `apprx` | removed (it changed nothing) |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| `Qubit.nqubit` | removed (the length of `Qubit.qubit`) |
| `npulse` | `pulse.count` |
| `addsubclus` 0 | `add_subclusters` false |
| `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift`, `bspin`, `balphams`, `bbetams` | `pulse.bath.count`, `.sequence`, `.defect`, `.energy_shift`, `.spin`, `.alphams`, `.betams` |
| top-level `alphams`, `betams` | removed (the same values are in the `Qubit.qubit` entry) |
| `nk` [0, 0], `hf_readmode` 0, `savemode` normal | left out (the defaults) |
| `-I bath_DEERtest_12spin.txt`, `-N 0`, `-S Random` | `bathfile`, `nstate` 0, `seed` 3528 (the input's) |

No key of the template is reported as not read.

## Inputs (md5 against the branch, 2026-10-13)

| file | md5 | from |
|---|---|---|
| `inputs/bath_DEERtest_12spin.txt` | 58cc77aaba9e91da23b9db04832a6d30 | `origin/adding_MECCE:example/z.Tests/DEER/` (same) |
| `inputs/avaax_DEERtest` | 4dca4e0bc9eeafc1681042ff4cc5fb44 | same (same) |
| `../bath_data/single_data/DiaP1_gyro` | 1132615a81ffc6c7e56c4093be91de01 | = `example/Diamond_NV_P1Bath/bath/DiaP1_gyro` of 503fe203a and of the branch tip |
| `inputs/state_1` | e636e210df8e06b8f439c3c4d512fc1d | written here (plan Task A3 Step 1): 12 spins, m = ±1/2, one fixed pure state for the path cases |

## Measured (WSL, MKL build, 1 rank, 2026-10-13)

| check | value |
|---|---|
| `nopulse` max \|V − 1\| | 8.76e-07 |
| `offres` − `nopulse`, `multi_miss` − `nopulse`, `multi_hit` − `onres` | 0 |
| `onres` − closed form (f = 0.5) | 5.33e-05 |
| `nearpulse` − closed form (f = 0.505) | 5.34e-05 |
| `nearpulse` − `nopulse` | 1.35 |
| vector − matrix, recursive − exact (vector), recursive − exact (matrix) | 0 (10 written decimals) |
| `path_vector` − `path_nopump` | 2.00 |
| `meg_onres` − `onres`, `meg_nearpulse` − `nearpulse`; `meg_o2` − `o2` | 1.0e-10, 1.0e-10 (the last written decimal); 8.6e-10 |
| `multi_hit_rev` − `onres`; `pg_onres` − `onres`; `pmeg_onres` − `meg_onres` | 0; 0; 0 |
| `o2` − `path_vector` (`_wiDiv`) | 1.17 |
| `o2_rev` − `o2`; `meg_o2_rev` − `meg_o2`; `o3_rot` − `o3` | 6.4e-10; 1.4e-10; 5.1e-10 |
| `onres_expm` − `onres`; `pg_o2` − `o2` (2026-10-14) | 1.0e-10 (the last written decimal); 3.2e-10 |

Not a check, for the record: the pure-state `path_vector` curve is `exp(-i 2 pi t sum_j a_j m_j)` over the driven spins
within 8.6e-4 (5.8e-4 with 52.0412 instead of 52.04 MHz nm^3; the phase reaches 12 rad), and `path_nopump` is 1 within
8.6e-7.
