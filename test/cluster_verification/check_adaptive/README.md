# check_adaptive/

What a change to `src/clustering/algo_adaptive/` has to pass. The C++ driver (`cpp/bin/ccex-adaptive`) runs the
refinement policies, launching `bin/main.out` for every step, and its decisions are held to recorded golden answers.

```
sh test/cluster_verification/check_adaptive/verify.sh unit        # build the driver and its unit tests, run them
sh test/cluster_verification/check_adaptive/verify.sh regress     # the policies against 33 recorded runs and their golden answers
sh test/cluster_verification/check_adaptive/verify.sh inprocess   # an "adaptive" block inside bin/main.out against the separate driver, then probe_rotation.py
sh test/cluster_verification/check_adaptive/verify.sh sanitize    # unit and regress again, built with ASan and UBSan
sh test/cluster_verification/check_adaptive/verify.sh fuzz [SEED [COUNT]]   # random JSON through the reader and the writer
```

`regress` and `inprocess` need `bin/main.out` and `mpirun`. Run them from a shell where Intel MPI's `mpirun` is first
on `PATH`.

## What it checks

| Case | Compared with | Tolerance |
|---|---|---|
| `cpp/test_*.cpp`: policy functions, JSON reader and writer, NaN, inf, zero, ties | unit tests | — |
| `cpp/test_regroup.cpp`: method D (`--regroup`) on a bath with a closed-form expansion | the closed form | P_K × N_c 1e-12, except path 1e-13 of `__float128`, 1e-15 with powers near −1000 |
| `regress`: 14 merge, cache and fault cases, 19 regroup cases (`--atol 1e-3 --rtol 0`) | `regression/golden/` | decisions exact (exit status, patches, rounds, every string, integer and boolean), floats 1e-9 relative, times not compared |
| The same with an MKL `bin/main.out` (log: `MKL for large clusters : on`) | `regression/golden/` | decisions exact, floats 1e-5 relative, at least 1e-12 absolute |
| `inprocess`: 11 cases, and `probe_rotation.py` across ranks and states | the separate driver | — |
| regroup factors and except paths, 1 and 3 ranks | the probe path | the same decisions |

Any nonzero exit after a timeout counts alike, so `fault-probe-timeout` also passes under another `mpirun`.
`test_adaptive_tools.py` (run by `unit`) holds the MKL rule.

## Updating the golden answers

The golden answers change only on purpose. Run `python3 regress.py WORKDIR --update` after the difference has been
reviewed, and commit the new `golden/` in the same commit as the change, with the reason. They are recorded without
MKL.

## Files

| File | Contents |
|---|---|
| `verify.sh` | Runs `unit`, `regress`, `inprocess`, `sanitize` or `fuzz`. |
| `cpp/` | The unit tests and the Makefile that builds them and the driver. `regroup_replay` is built but not run by `unit`. |
| `regression/` | `cases.json`, `inputs/`, `golden/`. |
| `reference/` | The JSON round-trip cases and the number-text boundaries of the reader and writer. |
| `regress.py`, `inprocess_compare.py`, `json_expect.py`, `json_roundtrip.py`, `compare_summaries.py`, `config_verdicts.py`, `test_adaptive_tools.py`, `probe_rotation.py` | The comparisons. |
| `fuzz_json.py`, `fault_mpiexec.py`, `hang.py` | The fuzzer and a fake `mpiexec` that fails, hangs or drops a response. |
| `regroup_parity.py` | A study, not a standing check. Method D against its prototype on stored probes. |

`cpp/bin*/` and `cpp/out/` are what a run writes. Both are ignored.
