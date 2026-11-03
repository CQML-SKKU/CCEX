# test/cluster_verification/

The clustering methods of this branch. (internal notes)
names them: **bundle clustering** builds the bundles once (`cluster_mode=bundle`,
`src/clustering/bundle/`); **adaptive refinement** rebuilds them from measured
error (the `adaptive` block, `src/clustering/algo_adaptive/`); **pCCE** partitions the bath
into centres.

| folder | what a change has to pass | run |
|---|---|---|
| `check_adaptive/` | a change to `src/clustering/algo_adaptive/`: the unit tests, the 33 recorded runs with their golden answers, the in-process driver against the separate one, the sanitizers, the fuzzers | `sh check_adaptive/verify.sh unit\|regress\|inprocess\|sanitize\|fuzz` |
| `check_bundle/` | a change to `src/clustering/bundle/`: the three bundle unit tests of `check_pcce/` (`test_adaptive_partition.cpp`, `test_adaptive_connections.cpp`, `test_adaptive_high_order.cpp`) and the recorded runs of `check_adaptive/`, which use the bundle partition. `check_bundle/` keeps the records of the 2026-09-12..16 experiments; their scripts, which no longer ran, were removed on 2026-11-03 | `make cpu-test` (CMake: `--target cpu-test`) runs the three since 2026-10-02 (D31), then `sh check_adaptive/verify.sh regress\|inprocess` |
| `check_pcce/` | a change to `src/clustering/cluster_pcce.cpp` or the partition cache: eight `check_*.py`, seven `test_*.cpp` | each `python3 check_pcce/check_<name>.py ...` on its own, with the options its row in `check_pcce/README.md` lists (they differ: `--cpu`, `--fixtures`, `--baseline`); the `test_*.cpp` compile lines are in its README |

`check_bundle/` and `check_pcce/` still hold the records and tools of the
2026-09-12..16 experiments beside the checks; their READMEs say which is which.
