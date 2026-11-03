# check_pcce/ — pCCE: the partition, its cache, and reproducibility

| script | checks | needs |
|---|---|---|
| `check_reproducibility.py --cpu BIN --output DIR` | the global trial count, rank-independent partitions and bath states | CPU |
| `check_partition.py --output DIR`, `check_partition_search.py --output DIR`, `check_partition_cycle.py --output DIR [--fixtures DIR]` | the partition unit checks (they build `test_partition_checks.cpp`, `test_partition_search.cpp`, `test_partition_cycle.cpp` against the objects of the current build; `--sanitize` rebuilds the partition sources with AddressSanitizer), the differential assignment checks, exact final partitions with and without cycle skipping | a build: `unit_build.py` asks `../../unit_build_env.sh`, which brings the build up to date and prints its compiler, flags, libraries and objects as the build writes them (the `unit-env` target of the Makefile of `do_compile.sh`; or, with `CCEX_BUILD_DIR=<a CMake build folder>` of any generator, its `unit-build.env`, `CCEX_BUILD_CONFIG=<config>` for a multi-config one; D26). The tests get the build's own flags, not a `-march=native` of their own |
| `check_partition_cache.py --cpu BIN --output DIR` | the production cache: MPI, sampled states, corrupt files | CPU (CUDA variant optional) |
| `check_partition_cycle_integration.py --fixtures DIR --baseline DIR --output DIR` | cycling on real baths, MPI statistics, the legacy RNG, cache migration | a baseline run |
| `check_edge_cases.py --cpu BIN --output DIR` | the empty-order regressions against full-CCE and first-order references | CPU |
| `check_baths.py --output DIR` | non-PL6 fixtures: pCCE CPU against CUDA and the old CPU path | `../../code_verification/bath_data/`, `check_input_forms/inputs` |
| `recheck_baths.py --original DIR --fixtures DIR --output DIR [--cpu BIN] [--reference-results DIR]` | rechecks the results of an earlier `check_baths.py` run (left unchanged) and runs only the physical references it lacks; with `--reference-results` it rechecks a downloaded run without running CCEX | an earlier `check_baths.py` output |
| `test_partition_log.py`, `test_profile_parser.py` | unit tests of the log parser here and of the profile parser in `../../gpu_verification/check_gpu/benchmark.py` | `python3 -m unittest` |

On 2026-09-24 `check_reproducibility.py`, `check_edge_cases.py` and the three
unit-check drivers passed on the WSL build. `test_adaptive_partition.cpp`,
`test_adaptive_connections.cpp` and `test_adaptive_high_order.cpp` (the bundle partition, connections and
high-order expansion) are built and run by `make cpu-test` / the CMake `cpu-test` target since 2026-10-02
(D31). `test_partition_cache_format.cpp`
(the cache file codec: binary64 roundtrip, every byte flipped and every truncation, trailing data, a changed key,
a missing file) runs in `cpu-test` as well since 2026-10-02 (D32), with a file in `obj/` or the CMake build folder
that it removes; it is also meant for an ASan/UBSan build. `benchmark_*.py` time the partition search and the cache;
`sweep_support.py`, `inventory_sweep.py` and `partition_log.py` are their helpers.

Removed on 2026-09-24 (last in commit a212af9e2: `git show a212af9e2:test/cluster_verification/check_pcce/<name>`):
the sweep and analysis tools of the 2026-09-11..12 PL6 work (`run_sweep.py`,
`prepare_sweep.py`, `instrument_product.py`, `run_product_trace.py`,
`analyze_*.py`, `finalize_product_convergence.py`, `check_product_convergence.py`,
`product_trace.h`, `replay_product.cpp`, `probe_product_mpi.cpp`,
`cluster_inventory.cpp`, `Makefile.inventory`) and 22 dated records
(`grid-*`, `benchmark-*`, `pilot-*`, `validation-*`, `performance-*`, 20 MB);
`product-convergence-20260912.json` stays as a record.
