# check_deer/

Checks the DEER bath pulses (`pulse.bath`, see `docs/methods/deer.md`) against known answers rather than a recorded
run, and checks that a bath pulse reaches every evolution path. The cases are 12 P1 spins at order 1 to 3, written by
`run.sh` from one template.

```
bash test/code_verification/check_deer/run.sh     # about 90 s; bin/main.out or CCEX_BIN, CCEX_MPI; 1 rank
```

## What it checks

At order 1 the echo is V(t) = ∏ⱼ cos(aⱼ t min(f, 1−f)) over the driven spins, aⱼ = 2π · 52.04 MHz nm³ (1 − 3cos²θⱼ) / rⱼ³,
with f the fraction of the time at which the bath pi lands.

| Case | Compared with | Tolerance |
|---|---|---|
| `test_merge_unit.cpp`: the qubit and bath sequence merge | each pulse on one edge, never twice | 1e-9 |
| `pulse.bath` count 0 on `../check_numerics/cases/p1_manual_hahn` | the input without the block | same bytes |
| `nopulse` | 1 | 1e-5 |
| `offres`, `multi_miss` (no line hit, must warn) | `nopulse` | 1e-12 |
| `multi_hit`, `multi_hit_rev` | `onres` | 0 |
| `onres` (f = 0.5), `nearpulse` (f = 0.505), real part | the closed form | 2e-4 |
| `nearpulse`, `path_vector`, `o2` (`_wiDiv`) | must differ from `nopulse`, `path_nopump`, `path_vector` | > 1e-3 |
| `_wiDiv` of the ensemble cases | their `_noDiv` | same bytes |
| `path_*` (one pure state): vector, matrix, exact, recursive | each other | 1e-9 |
| `meg_onres`, `meg_nearpulse`, `meg_o2` (megCCE, rate 0) | `onres`, `nearpulse`, `o2` | 1e-6 |
| `pg_onres`, `pmeg_onres` (pcce, `partition_size` 1), `onres_expm` | `onres`, `meg_onres`, `onres` | 1e-9 |
| `pg_o2`; `o2_rev`, `meg_o2_rev` (spins reversed), `o3_rot` (spins rotated by 5) | `o2`; `o2`, `meg_o2`, `o3` | 1e-8 |
| `units_onres`, `bw_wide` (bandwidth 5 kHz) | `onres`, and `bw_default` must drive nothing and warn | same bytes |
| `bad_bspin`, `bad_defect`, `bad_ms_equal`, `bad_ms_range`, `bad_bw` | refused before any cluster, with their message | — |

Everything is judged on `_noDiv`, and the order-2 and order-3 comparisons also on `_wiDiv`. Every run uses `seed` 3528.
The pure-state runs read `inputs/state_1` through `configuration.bath_state_file` with `nstate` 1, and the log must
name that file. Reversing and rotating the spin order catch a bath rotation put on the wrong spin of a cluster.
The exit status is 0 only when every check passes.

## Files

| File | Contents |
|---|---|
| `run.sh` | Writes the cases from the template and runs them. |
| `check_deer.py` | The closed form and the comparisons. |
| `test_merge_unit.cpp` | The unit test of the merge in `src/physical/pulse_bath.cpp`. |
| `inputs/` | The 12-spin bath and axes of the `adding_MECCE` branch's DEER test, and `state_1`. |
