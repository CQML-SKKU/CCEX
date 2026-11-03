# check_pcce/: pCCE partition, cache and reproducibility

| script | checks | needs |
|---|---|---|
| `check_reproducibility.py --cpu BIN --output DIR` | the global trial count, rank-independent partitions and bath states | CPU |
| `check_partition.py --output DIR`, `check_partition_search.py --output DIR`, `check_partition_cycle.py --output DIR [--fixtures DIR]` | the partition unit checks, the differential assignment checks, exact final partitions with and without cycle skipping. They build `test_partition_checks.cpp`, `test_partition_search.cpp`, `test_partition_cycle.cpp` against the objects of the current build; `--sanitize` rebuilds the partition sources with AddressSanitizer | a build (see below) |
| `check_partition_cache.py --cpu BIN --output DIR` | the production cache: MPI, sampled states, corrupt files | CPU (CUDA variant optional) |
| `check_partition_cycle_integration.py --fixtures DIR --baseline DIR --output DIR` | cycling on real baths, MPI statistics, the legacy RNG, cache migration | a baseline run |
| `check_edge_cases.py --cpu BIN --output DIR` | the empty-order regressions against full-CCE and first-order references | CPU |
| `check_baths.py --output DIR` | non-PL6 fixtures: pCCE CPU against CUDA and the old CPU path | `../../code_verification/bath_data/`, `check_input_forms/inputs` |
| `recheck_baths.py --original DIR --fixtures DIR --output DIR [--cpu BIN] [--reference-results DIR]` | rechecks the results of an earlier `check_baths.py` run (left unchanged) and runs only the physical references it lacks; with `--reference-results` it rechecks a downloaded run without running CCEX | an earlier `check_baths.py` output |
| `test_partition_log.py`, `test_profile_parser.py` | unit tests of the log parser here and of the profile parser in `../../gpu_verification/check_gpu/benchmark.py` | `python3 -m unittest` |

The unit-check drivers find the build through `unit_build.py`, which asks `../../unit_build_env.sh`. That script
brings the build up to date and prints its compiler, flags, libraries and objects (the `unit-env` target of the
Makefile of `do_compile.sh`). For a CMake build folder of any generator set `CCEX_BUILD_DIR=<folder>` (it reads
`unit-build.env`), and `CCEX_BUILD_CONFIG=<config>` for a multi-config one. The tests get the build's own flags, not
a `-march=native` of their own.

`make cpu-test` / the CMake `cpu-test` target builds and runs four tests of this folder:

- `test_adaptive_partition.cpp`, `test_adaptive_connections.cpp`, `test_adaptive_high_order.cpp`: the bundle
  partition, connections and high-order expansion.
- `test_partition_cache_format.cpp`: the cache file codec (binary64 roundtrip, every byte flipped and every
  truncation, trailing data, a changed key, a missing file). It writes a file in `obj/` or the CMake build folder and
  removes it. It is also meant for an ASan/UBSan build.

`benchmark_*.py` time the partition search and the cache; `sweep_support.py`, `inventory_sweep.py` and
`partition_log.py` are their helpers. `product-convergence-20260912.json` is a record of the 2026-09 PL6 work.

The sweep and analysis tools of that work and 22 dated records were removed. They are in
`git show a212af9e2:test/cluster_verification/check_pcce/<name>`: `run_sweep.py`, `prepare_sweep.py`,
`instrument_product.py`, `run_product_trace.py`, `analyze_*.py`, `finalize_product_convergence.py`,
`check_product_convergence.py`, `product_trace.h`, `replay_product.cpp`, `probe_product_mpi.cpp`,
`cluster_inventory.cpp`, `Makefile.inventory`, and the records `grid-*`, `benchmark-*`, `pilot-*`, `validation-*`,
`performance-*`.
