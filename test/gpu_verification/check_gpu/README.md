# check_gpu/ — the CUDA backend gives the CPU's numbers

| script | checks | run |
|---|---|---|
| `test_*.cpp` (14) | unit tests of the kernels, the batching, the scaled product, the Hamiltonian check before the CUDA cCCE path, the seed sharing and the host-side lifetimes | `make numeric-test`, `cpu-test`, `gpu-test` (not `bash do_compile.sh gpu`, which builds `bin/main.out` and `bin/main_gpu.out` only): see "Unit tests" below |
| `check_cce_pulses.py --cpu BIN --cuda BIN --output DIR` | cCCE with every pulse count, CPU against CUDA | CUDA build |
| `check_acceleration.py`, `check_pcce.py`, `check_options.py` | gCCE4 accuracy and batching, pCCE, `backend` and the `gpu_options` block, CPU against CUDA. `check_options.py` also (review 20261013-221715, ④): X13 `cce3-ready-packed-warp` — three states and nstep 3000, so states 2–3 run after the card is up and an order fills several 131-cluster windows; packed and warp counts must both be > 0 and more than two windows must pass the feeder (the earlier `cce3-auto-packed-warp` may legitimately finish on the CPU); X12 `*-oom-late-*` — `CCEX_TEST_GPU_OOM_BYTES=8192` fails a cluster after its first allocations, and the `[gpu-memory]` line must show allocations > 0 and `cached_bytes=0`. The hidden-card cases run with `CUDA_VISIBLE_DEVICES=-1` everywhere (CCEX answers "no GPU" without the driver) and `=99` only off WSL: WSL's CUDA driver 577.00 crashes inside `cuInit` for an index it does not have, with or without CCEX (A55) | CUDA build |
| `check_warp_spread.py REFINED WARP` (`--spread ICPX GCC WARP`) | the warp kernel's whole production curves (A50b (iii), user decision 2026-10-13): against `main.out` with `gcce_refinement` on (and `evolution_update` exact), in the A40 verdict region (outside \|CPU\| > 1 and between), 1e-8 + 2e-6 \|CPU\| per value. Why the refined curve: every code's double eigenvalues are rounded at an ulp of \|H\| (~2.4e7 on P1 gCCE2) and the product amplifies it to 1e-6 - 1e-5 relative; two refined CPU builds agree to 2.5e-7 (CR60), and the warp kernel is as close to the refined curve as the production CPU ((internal notes)). `--spread`: A50 (a), \|warp − icpx\| ≤ max(3 \|gcc − icpx\|, 2e-6 \|CPU\|, 1e-12). Strict reading (review 20261013-180049 G1): time column and every complex value; nan/inf, an empty file, a line that does not parse, a missing file, different times or row/value counts, or nothing judged FAIL; `test_check_warp_spread.py` (CI python job) checks those. `test_warp_cluster_cuda` and `check_options.py` keep 1e-8 + 2e-6: per cluster the warp kernel is 5-9 orders inside it | output folders of server runs |
| `check_packed_emulation.py --output DIR` | the packed batch (`gpu_options.packed_max_dim`) end to end on 3 ranks with `small_cluster_host_emulation.cpp` for the card: a card that cannot be opened (auto: all on the CPU, byte-identical; cpu with an explicit value: refused; fast, late and without clusters), the 6 MiB window budget, the batch off without clusters or beyond the budget, the batch on within 1e-8 + 2e-6 \|CPU\| (review 20261012-233058-452278 P3-01, P3-04); about 15 s | CPU build is enough (it links the build's objects, `CCEX_BUILD_DIR` for CMake) |
| `check_seed_states.py --cpu BIN --output DIR` | seeds and sampled bath states are reproduced across ranks | CPU build is enough |
| `check_review_controls.py --cpu BIN --output DIR` | the controls the 2026-09-12 review asked for | CPU build is enough |
| `check_review_supplement.py --cpu BIN --cuda BIN --before-cuda BIN --output DIR` | the review's supplement | CUDA build, and the CUDA binary from before the review |
| `check_review_baseline.py --baseline BIN --results DIR` | the review's baseline: reruns the CPU inputs of an earlier `check_acceleration`/`check_cce_pulses`/`check_pcce` results folder with a frozen pre-review binary | a CPU binary (`--baseline`) and that results folder |
| `check_host_leaks.py --binary BIN --baseline BIN --output DIR [--cuda BIN]` | host memory over repeated runs | CPU build is enough (the current and a baseline binary); `--cuda` adds the CUDA variants |

## Unit tests: which command builds and runs which

After `bash do_compile.sh <preset>` has written the Makefile (or a CMake configure):

| tests | built and run by | needs |
|---|---|---|
| `test_scaled_product`, `test_hamiltonian_check`, `test_small_cluster_host` (`test_small_cluster_cuda.cpp` with `-DSMALL_CLUSTER_HOST_ONLY`: the packed arithmetic on the 11 dumps, the dump reader on cut copies, Hamiltonians of scale 1e±200 that must come out right or go back to the CPU) | `make numeric-test` (CMake: `cmake --build build --target numeric-test`) | a CPU |
| `test_cluster_log`, `test_review_controls`, `test_mpi_seed` (1 and 4 ranks), `test_host_lifetimes` (1 and 4 ranks) | `make cpu-test` (CMake: `--target cpu-test`; also in the CI workflow) | a CPU and `mpirun` (`make cpu-test MPIRUN=...` picks another; CMake uses the launcher FindMPI found) |
| the numeric two and `test_gcce_cuda`, `test_cce_cuda`, `test_gcce_accuracy` (twice: built-in case and `data/gcce_scale_separation_three_spins.txt`), `test_gcce_refinement` (the parallel refinement at dim 192, 15, 2 and 1), `test_gcce_batch`, `test_review_failure`, `test_review_memory` (also X12: after an out-of-memory part way, `gpuOutOfMemoryRecovery` empties the cache, a work released afterwards is freed, and caching resumes), `test_small_cluster_cuda` (the packed kernel of `gpu_options.packed_max_dim`, A26: the same arithmetic on the host and on the card against the CPU coherences of the 11 files in `data/small_cluster/`, 1e-8 + 2e-6 |CPU| per point, both lanes the same bits; `data/small_cluster/make_dumps.py` writes the files again with a CPU build, byte for byte; `small_cluster_host_emulation.cpp` stands in for the card so the node feeder and shared-memory windows can be run end to end without CUDA, as its header says), `test_warp_cluster_cuda` (the warp kernel of `gpu_options.warp_max_dim`, A40: through the batch API on the CPU-dumped clusters of the 15 files in `data/warp_cluster/` — every compiled dimension, gCCE 6–27 and cCCE 8, 16, single state and ensemble, coherence and dm, Hahn, CPMG, Ramsey, 64 production P1 clusters; gCCE dimension 12 Hahn goes through the X11 pulse-in-eigenbasis variant — 1e-8 + 2e-6 \|CPU\| per value, both lanes the same bits; a file cut short anywhere (header, plan, body, one byte) is refused, and every file is also read cut in memory at those places (review 20261013-180049 G2, as `test_small_cluster_cuda`); `data/warp_cluster/make_dumps.py` writes them with a CPU build; the host emulation has no warp kernel) | `make gpu-test` (CMake: `--target gpu-test`) | a CUDA build and the GPU, so inside a job; `server/build.sbatch` builds `test_cce_cuda` and `test_gcce_cuda`, and `server/measure.sbatch` runs those two only |

`cpu-test` compiles the four with the build's own flags and `-UNDEBUG` for their asserts;
`test_host_lifetimes` links the objects of `bin/main.out` but `main.o` (CMake compiles those
sources again for it). Until 2026-10-02 these four and the two CUDA review tests had no target
and were built by hand (D17); built without `-march=native` against objects built with it,
`test_host_lifetimes` crashed (exit 139) on 2026-10-01.

On 2026-09-24 `check_review_controls.py --cpu` passed on the CPU build;
`check_seed_states.py` wrote `adaptive_cluster_max_dim` in its adaptive cases until
2026-09-24; it wrote `bundle_max_dim` from then on (`bundle.max_dim` since 2026-10-13) and its nine cases pass on the CPU build.
One run in two failed in `benchmark.parse_profiles` because MPI interleaved two
ranks' `[profile]` lines in the p1 case; that is the output forwarding, not the
code, and a rerun passed.

`benchmark.py` is the helper the checks and `../gpu_benchmarks/` import: it runs a
build on the four standard calculations of `code_verification/check_input_forms/inputs`
and reads what it wrote. `examples/` holds the P1 and pCCE inputs the checks use, `data/` the gCCE accuracy
test's files and, in `data/small_cluster/`, the 11 `.bin` dumps of `test_small_cluster_cuda` with `make_dumps.py`. The seven `validation-*.json` files are the
results recorded on idm4-gpu and the RTX 3050 on 2026-09-08..11 that (internal notes)
cites; the other records of 2026-09-12 and the probe that led to `scaled_product.h`
were removed on 2026-09-24 (git history up to 8f2bcbd75).
