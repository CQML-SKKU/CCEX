# cluster_verification/

Checks of the clustering methods. Bundle clustering builds the bundles once (`cluster_mode` bundle,
`src/clustering/bundle/`). Adaptive refinement rebuilds them from measured error (the `adaptive` block,
`src/clustering/algo_adaptive/`). pCCE partitions the bath into centres (`src/clustering/cluster_pcce.cpp`).

```
sh test/cluster_verification/check_adaptive/verify.sh unit|regress|inprocess|sanitize|fuzz
make cpu-test                                   # CMake: --target cpu-test; the bundle and partition cache unit tests
python3 test/cluster_verification/check_pcce/check_<name>.py ...   # each pCCE check on its own
```

## What it checks

| Folder | A change to | Run |
|---|---|---|
| `check_adaptive/` | `src/clustering/algo_adaptive/`: unit tests, 33 recorded runs with golden answers, the in-process driver against the separate one, sanitizers, fuzzers. | `verify.sh unit`, `regress`, `inprocess` (and `sanitize`) |
| `check_bundle/` | `src/clustering/bundle/`: the three bundle unit tests in `check_pcce/` and the recorded runs of `check_adaptive/`, which use the bundle partition. | `make cpu-test`, then `check_adaptive/verify.sh regress` and `inprocess` |
| `check_pcce/` | `src/clustering/cluster_pcce.cpp` or the partition cache: eight `check_*.py` and seven `test_*.cpp`. | each script with the options its row in `check_pcce/README.md` lists |

## Files

| Folder | Contents |
|---|---|
| `check_adaptive/` | The adaptive driver's checks and golden answers. |
| `check_bundle/` | Bundle inputs and the records of the 2026-09 experiments. No runnable check. |
| `check_pcce/` | The pCCE partition, cache and reproducibility checks, their benchmarks and the bundle unit tests. |
