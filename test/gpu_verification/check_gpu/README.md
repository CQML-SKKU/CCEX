# check_gpu/

Checks that the CUDA backend gives the CPU's numbers. The unit tests check the kernels, batching and host-side
bookkeeping; the `check_*.py` scripts run CPU and CUDA builds on the same inputs. Most need `bin/main_gpu.out` and a
GPU, so they run in a job on idm4-gpu.

```
make numeric-test      # CPU only
make cpu-test          # CPU and mpirun (MPIRUN=... picks another); also in CI
make gpu-test          # CUDA build and the GPU, inside a job
# CMake: cmake --build build --target numeric-test | cpu-test | gpu-test
```

`bash do_compile.sh gpu` builds only `bin/main.out` and `bin/main_gpu.out`, not the tests. Build the tests through
these targets, which use the build's own flags (and `-UNDEBUG` for `cpu-test`).

## What it checks

| Case | Compared with | Tolerance |
|---|---|---|
| `numeric-test`: `test_scaled_product`, `test_hamiltonian_check`, `test_small_cluster_host` (Hamiltonians of scale 1e±200 included) | the CPU | right, or back to the CPU |
| `cpu-test`: `test_cluster_log`, `test_review_controls`, `test_mpi_seed`, `test_host_lifetimes` (1 and 4 ranks) | unit tests | — |
| `gpu-test`: the numeric tests, `test_gcce_*`, `test_cce_cuda`, `test_review_failure`, `test_review_memory` | the CPU | unit tests |
| `test_small_cluster_cuda`: the packed kernel (`gpu_options.packed_max_dim`) on `data/small_cluster/` | CPU coherences | 1e-8 + 2e-6 \|CPU\| per point; both lanes the same bits |
| `test_warp_cluster_cuda`: the warp kernel (`gpu_options.warp_max_dim`) on `data/warp_cluster/` | CPU coherences | 1e-8 + 2e-6 \|CPU\|; both lanes the same bits |
| `check_cce_pulses.py --cpu BIN --cuda BIN --output DIR` | cCCE with every pulse count | CPU |
| `check_acceleration.py`, `check_pcce.py`, `check_options.py` | gCCE4 accuracy and batching, pCCE, `execution.backend` and `execution.gpu_options` | CPU; `check_options.py` 1e-8 + 2e-6 \|CPU\| |
| `check_auto_k_range.py --cuda BIN --output DIR [--mpirun M]`: `backend` auto with `k_range` on two ranks, one opening the card | both configurations written | within 180 s |
| `check_warp_spread.py REFINED WARP` | `main.out` with `gcce_refinement` on and `evolution_update` exact, in the verdict region | 1e-8 + 2e-6 \|CPU\| per value |
| `check_warp_spread.py --spread ICPX GCC WARP` | the spread of two CPU builds | \|warp − icpx\| ≤ max(3 \|gcc − icpx\|, 2e-6 \|CPU\|, 1e-12) |
| `check_packed_emulation.py --output DIR` (about 15 s, CPU build): the packed batch on 3 ranks, emulated card | the CPU | 1e-8 + 2e-6 \|CPU\|; 6 MiB window budget; auto with no card same bytes |
| `check_seed_states.py --cpu BIN --output DIR` | seeds and sampled states across ranks | CPU build |
| `check_review_controls.py --cpu BIN`, `check_review_supplement.py`, `check_review_baseline.py --baseline BIN --results DIR` | the controls, supplement and baseline of the 2026-09-12 review | `--before-cuda` / `--baseline` binaries |
| `check_host_leaks.py --binary BIN --baseline BIN --output DIR [--cuda BIN]` | host memory over repeated runs | CPU build |

`check_warp_spread.py` fails on nan, inf, unreadable files, different grids or nothing judged. In `check_options.py`, `cce3-ready-packed-warp` needs packed and
warp counts above 0 and more than two windows through the feeder; the `*-oom-late-*` cases set
`CCEX_TEST_GPU_OOM_BYTES=8192` and the `[gpu-memory]` line must show allocations above 0 and `cached_bytes=0`. Hidden
cards are `CUDA_VISIBLE_DEVICES=-1`, and `=99` only off WSL (its driver crashes in `cuInit`). A `check_seed_states.py`
failure in `benchmark.parse_profiles` can be two ranks' interleaved `[profile]` lines; rerun it.

## Files

| File | Contents |
|---|---|
| `test_*.cpp` | the unit tests (`test_small_cluster_host` is `test_small_cluster_cuda.cpp` built with `-DSMALL_CLUSTER_HOST_ONLY`) |
| `small_cluster_host_emulation.cpp` | stands in for the card in `check_packed_emulation.py` |
| `check_*.py` | the CPU against CUDA checks above |
| `test_check_warp_spread.py` | the unit test of `check_warp_spread.py` |
| `benchmark.py` | the helper the checks and `../gpu_benchmarks/` import |
| `examples/` | the P1 and pCCE inputs |
| `data/` | the gCCE accuracy files, `small_cluster/` and `warp_cluster/` dumps, each with `make_dumps.py` (see `data/README.md`) |
| `validation-*.json` | recorded results of the checks |
