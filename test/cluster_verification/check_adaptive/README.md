# check_adaptive/: what a change to src/clustering/algo_adaptive/ has to pass

    sh test/cluster_verification/check_adaptive/verify.sh unit        # build the driver (cpp/bin/ccex-adaptive) and its unit tests, run them
    sh test/cluster_verification/check_adaptive/verify.sh regress     # the merge policy against 33 recorded runs and their golden answers
    sh test/cluster_verification/check_adaptive/verify.sh inprocess   # the driver from an "adaptive" block inside bin/main.out against the separate one (11 cases), then probe_rotation.py
    sh test/cluster_verification/check_adaptive/verify.sh sanitize    # unit and regress again, built with ASan and UBSan
    sh test/cluster_verification/check_adaptive/verify.sh fuzz [SEED [COUNT]]   # random JSON through the reader and the writer

`regress` and `inprocess` need `bin/main.out` and `mpirun`. Run them from a shell where Intel MPI's `mpirun` is first
on `PATH`. `regress.py` treats any nonzero exit after a timeout alike, so the timeout fault case fault-probe-timeout also passes
under anaconda's or Ubuntu's `mpirun`.

## How `regress` decides

Every decision is held exactly to `regression/golden/`: exit status, accepted patches, rounds, policy choices, and
every string, integer and boolean of the summaries. A float is held to 1e-9 relative.

When `bin/main.out` was built with MKL (the build turns it on when it finds MKL; the run logs of the case say
`MKL for large clusters : on`), a float may instead move by 1e-5 of the golden value, and never less than 1e-12
absolute. The absolute floor covers a golden value of 0 and round-off-sized values such as regroup's
max_abs_difference_from_prediction (1.0e-15 in golden against 8.9e-16 with MKL). Decisions stay exact. The reason:
MKL changes the last digits of the clusters it takes (zgemm in the cCCE ensemble Hahn echo's recursive path from bath
dimension 21, zheevd from dimension 63: the bath for cCCE, qubit × bath for gCCE), and the policies' error numbers are
differences of nearly equal coherences. Measured: up to 4.4e-6 relative (WSL, MKL 2026.0: 8 of the 33 cases,
scaled_error and local_before/after, with the same decisions). The golden answers are the ones recorded without MKL;
`test_adaptive_tools.py` (`RegressVerdict`, run by `unit`) holds the rule.

## Files

| here | is |
|---|---|
| `cpp/` | the unit tests and the Makefile that builds them and the driver from `src/clustering/algo_adaptive/`; `regroup_replay` (built, not run by `unit`) is the stored-probe replay of method D below |
| `regression/` | `cases.json` (33 cases: 14 merge, cache and fault requests, among them fault-probe-timeout, a merge run whose first probe hangs, and the 19 regroup cases below), `inputs/`, `golden/` (the answers, compared field by field except the measured times) |
| `reference/` | the JSON round-trip cases and the number-text boundaries the reader and writer are held to |
| `regress.py`, `inprocess_compare.py`, `json_expect.py`, `json_roundtrip.py`, `compare_summaries.py`, `config_verdicts.py`, `test_adaptive_tools.py`, `probe_rotation.py` | the comparisons (`probe_rotation.py`: the probe step across ranks and states, with its memory budget, run by `inprocess`) |
| `fuzz_json.py`, `fault_mpiexec.py`, `hang.py` | the fuzzer and the fault injector (a fake `mpiexec` that fails, hangs, ignores SIGTERM or drops a response file on a chosen step) |

`cpp/bin*/` and `cpp/out/` are what a run writes; both are ignored.

## Method D (`--regroup`)

The unit test `cpp/test_regroup.cpp` checks every function of `policy_regroup` on its own (NaN, inf, zero and ties
included). It also runs whole rounds against a bath whose expansion is known in closed form, so that each of these
is seen to act: the three acceptance conditions, the budget refusal, the expansion-limit refusal, the revisit block,
a rollback, the rollback bound and the change limit.

For the factors and except paths it also checks that:

- the ranking from point values read back from a factors response is the ranking from the curves (NaN, inf, 0, ties,
  no positive culprit);
- a factors or except response with a wrong file, row or header is refused;
- the term difference X and the new terms N_c rebuild every candidate's exponent table on a sparse graph (singletons
  outside the region that change power, cancelled and new terms);
- P_K x N_c is the whole expansion within 1e-12. The except path reassembles in long double as CCEX does:
  approximate_from from L0 agrees with a __float128 product and with approximate_weights within 1e-13, is non-finite
  at the same points, and with singleton powers near -1000 stays below 1e-15 of the __float128 product where double
  errs by ~4e-13;
- rounds through the four paths decide alike on the closed-form bath (complete and sparse graphs);
- an except step that does not reproduce the run ends the round with an error;
- the except path drops its regions' curves.

The regression cases, all `--regroup`, run with `--atol 1e-3 --rtol 0` (D does not read them):

| case | is |
|---|---|
| `regroup-p1` | six P1 sites, three bundles of two; one split, then physical (exit 0) |
| `regroup-p1-sparse` | eight P1 sites with rdip 45 Å, so not every pair is connected: a site moves between bundles, neighbours are checked, a revisit is blocked, then no candidate (exit 2) |
| `regroup-p1-limit` | `regroup-p1` with `bundle.max_clusters` 6: the three splits are refused before CCEX is asked and the rest are judged, no candidate (exit 2) |
| `regroup-p1-budget` | `regroup-p1` with `--max-total-work 5e6`: the region's direct solution does not fit, `budget_limited` (exit 2) |
| `regroup-c13-physical` | the 13C input of `merge-cce`: physical from the start, nothing changes (exit 0) |
| `regroup-fault-probe-drop` | `regroup-p1` whose first probe loses its response: `failed` (exit 1) |
| `regroup-refuse-gcce`, `-order1`, `-auto`, `-prefer` | the scope refusals, before anything runs (exit 1, 1, 1, 2) |
| `regroup-p1-factors` | `regroup-p1` with `--regroup-factors factors`: the same decisions, the culprit from one factors step |
| `regroup-p1-probecap`, `regroup-p1-factors-probecap` | `regroup-p1` with `--max-probes 5`: the probe path cannot rank the culprits (budget_limited before any), the factors path ranks them and is refused at the region's direct solutions (both exit 2) |
| `regroup-fault-factors-drop` | `regroup-p1-factors` whose factors step loses rank 0's response: `failed` (exit 1) |
| `regroup-p1-except`, `regroup-p1-sparse-except` | `regroup-p1` and `regroup-p1-sparse` with `--regroup-whole except`: the same decisions, (c) from one except step per region |
| `regroup-fault-except-drop` | `regroup-p1-except` whose except step loses rank 0's excepted curves: `failed` (exit 1) |
| `regroup-refuse-factors-alone`, `regroup-refuse-except-contributions` | `--regroup-factors` without `--regroup` (exit 2); `--regroup-whole except` with `cluster_contributions`, which CCEX's except step refuses (exit 1) |

The P1 inputs come from the check_input_forms P1 10 ppm input
(`test/code_verification/check_input_forms/inputs/Diamond_NV_P1_10ppm_cce`): the six (eight) sites nearest the NV,
their positions pulled halfway toward the NV (`bath_p1_6`, `bath_p1_8`), the 15N dropped and the P1 hyperfine given
as the axis detunings the stored P1 subsystems used, five states, 30 points.

The factors and except paths (`--regroup-factors factors`, `--regroup-whole except`) are compared with the probe
path on the regroup cases, on 1 and 3 ranks.

### Replay of stored probes (a study, not a standing check)

The replay of the stored P1 probes needs the stored data (`~/tmp/ccex-p1-cluster`). `rebundle2.py trace` writes the
prototype's decisions with `REBUNDLE_NOREVISIT=1`, `cpp/bin/regroup_replay` makes the C++ decisions on the same
probes, and `regroup_parity.py PY.jsonl CPP.jsonl` compares them. `regroup_replay --factors --except` takes the paths
of the two options on the same probes.
