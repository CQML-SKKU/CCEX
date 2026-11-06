# check_pcce/

Checks of the pCCE partition, its cache and its reproducibility (`src/clustering/cluster_pcce.cpp`). Each script runs
on its own, with its own options. The folder also holds the three bundle unit tests and the partition benchmarks.

```
python3 test/cluster_verification/check_pcce/check_<name>.py --output DIR [...]   # options per script below
make cpu-test                                                                    # CMake: --target cpu-test; four unit tests of this folder
```

## What it checks

| Script | Checks | Needs |
|---|---|---|
| `check_reproducibility.py --cpu BIN --output DIR` | The global trial count, rank-independent partitions and bath states. | CPU |
| `check_partition.py`, `check_partition_search.py`, `check_partition_cycle.py [--fixtures DIR]`, each `--output DIR` | The partition unit checks, differential assignment checks, exact final partitions with and without cycle skipping (they build `test_partition_checks.cpp`, `test_partition_search.cpp`, `test_partition_cycle.cpp`, and `--sanitize` rebuilds with AddressSanitizer). | a build |
| `check_partition_cache.py --cpu BIN --output DIR` | The production cache: MPI, sampled states, corrupt files. | CPU (CUDA variant optional) |
| `check_partition_cycle_integration.py --fixtures DIR --baseline DIR --output DIR` | Cycling on real baths, MPI statistics, the legacy RNG, cache migration. | a baseline run |
| `check_edge_cases.py --cpu BIN --output DIR` | The empty-order regressions against full-CCE and first-order references. | CPU |
| `check_baths.py --output DIR` | On non-PL6 fixtures, pCCE CPU against CUDA and the old CPU path. | `../../code_verification/bath_data/`, `../../code_verification/check_input_forms/inputs` |
| `recheck_baths.py --original DIR --fixtures DIR --output DIR [--cpu BIN] [--reference-results DIR]` | Rechecks an earlier `check_baths.py` output and runs only the missing references. With `--reference-results`, it does so without running CCEX. | an earlier `check_baths.py` output |
| `test_partition_log.py`, `test_profile_parser.py` | The log parser here and the profile parser of `../../gpu_verification/check_gpu/benchmark.py`. | `python3 -m unittest` |
| `make cpu-test`: `test_adaptive_partition.cpp`, `test_adaptive_connections.cpp`, `test_adaptive_high_order.cpp` | The bundle partition, connections and high-order expansion. | a build |
| `make cpu-test`: `test_partition_cache_format.cpp` | The cache file codec: binary64 roundtrip, every byte flipped, every truncation, trailing data, a changed key, a missing file (also for an ASan/UBSan build). | a build |

The unit-check drivers find the build through `unit_build.py`, which asks `../../unit_build_env.sh` (the `unit-env`
target of the `do_compile.sh` Makefile). For a CMake build set `CCEX_BUILD_DIR=<folder>`, and `CCEX_BUILD_CONFIG=<config>`
for a multi-config generator. The tests get the build's own flags.

## Files

| File | Contents |
|---|---|
| `check_*.py`, `recheck_baths.py` | The checks above. |
| `test_partition_checks.cpp`, `test_partition_search.cpp`, `test_partition_cycle.cpp` | Built by the partition checks. |
| `test_partition_cache_format.cpp`, `test_adaptive_*.cpp` | The unit tests of `make cpu-test`. |
| `test_partition_log.py`, `test_profile_parser.py` | Python unit tests. |
| `unit_build.py` | Finds the build's compiler, flags and objects. |
| `benchmark_baths.py`, `benchmark_partition_search.py`, `benchmark_partition_cache.py` | Time the partition search and the cache. |
| `sweep_support.py`, `inventory_sweep.py`, `partition_log.py` | Helpers of the benchmarks. |
| `product-convergence-20260912.json` | A record of the 2026-09 PL6 work. |
