# check_gpu/: the CUDA backend gives the CPU's numbers

| script | checks | run |
|---|---|---|
| `test_*.cpp` (14) | unit tests of the kernels, the batching, the scaled product, the Hamiltonian check before the CUDA cCCE path, the seed sharing and the host-side lifetimes | `make numeric-test`, `cpu-test`, `gpu-test` (not `bash do_compile.sh gpu`, which builds only `bin/main.out` and `bin/main_gpu.out`). See "Unit tests" below. |
| `check_cce_pulses.py --cpu BIN --cuda BIN --output DIR` | cCCE with every pulse count, CPU against CUDA | CUDA build |
| `check_auto_k_range.py --cuda BIN --output DIR [--mpirun M]` | `backend` `auto` with `configuration.k_range` on two ranks, only one of which opens the card. Both configurations must be written and the run must end within 180 s. | CUDA build and a GPU (an RTX 3050 is enough, a few seconds) |
| `check_acceleration.py`, `check_pcce.py`, `check_options.py` | gCCE4 accuracy and batching, pCCE, `execution.backend` and the `execution.gpu_options` block, CPU against CUDA. `check_options.py` also runs `cce3-ready-packed-warp`: three states and nstep 3000, so states 2 and 3 run after the card is up and an order fills several 131-cluster windows. Packed and warp counts must both be > 0 and more than two windows must pass the feeder (`cce3-auto-packed-warp` may legitimately finish on the CPU). The `*-oom-late-*` cases set `CCEX_TEST_GPU_OOM_BYTES=8192` so a cluster fails after its first allocations; the `[gpu-memory]` line must show allocations > 0 and `cached_bytes=0`. The hidden-card cases run with `CUDA_VISIBLE_DEVICES=-1` everywhere (CCEX answers "no GPU" without the driver) and with `=99` only off WSL, because WSL's CUDA driver 577.00 crashes inside `cuInit` for an index it does not have, with or without CCEX. | CUDA build |
| `check_warp_spread.py REFINED WARP` (`--spread ICPX GCC WARP`) | the warp kernel's whole production curves against `main.out` with `gcce_refinement` on (and `evolution_update` exact), in the verdict region (outside \|CPU\| > 1 and between), 1e-8 + 2e-6 \|CPU\| per value. The refined curve is the reference because every code rounds its double eigenvalues at an ulp of \|H\| (~2.4e7 on P1 gCCE2) and the product amplifies that to 1e-6 to 1e-5 relative; two refined CPU builds agree to 2.5e-7 (CR60), and the warp kernel is as close to the refined curve as the production CPU. `--spread`: \|warp − icpx\| ≤ max(3 \|gcc − icpx\|, 2e-6 \|CPU\|, 1e-12). It reads the time column and every complex value. FAIL on nan/inf, an empty file, a line that does not parse, a missing file, different times or row/value counts, or nothing judged; `test_check_warp_spread.py` (CI python job) checks those. `test_warp_cluster_cuda` and `check_options.py` keep 1e-8 + 2e-6: per cluster the warp kernel is 5 to 9 orders inside it. | output folders of server runs |
| `check_packed_emulation.py --output DIR` | the packed batch (`execution.gpu_options.packed_max_dim`) end to end on 3 ranks, with `small_cluster_host_emulation.cpp` standing in for the card: a card that cannot be opened (auto: all on the CPU, byte-identical; cpu with an explicit value: refused; fast, late and without clusters), the 6 MiB window budget, the batch off without clusters or beyond the budget, the batch on within 1e-8 + 2e-6 \|CPU\|. About 15 s. | a CPU build is enough (it links the build's objects; `CCEX_BUILD_DIR` for CMake) |
| `check_seed_states.py --cpu BIN --output DIR` | seeds and sampled bath states are reproduced across ranks | a CPU build is enough |
| `check_review_controls.py --cpu BIN --output DIR` | the controls the 2026-09-12 review asked for | a CPU build is enough |
| `check_review_supplement.py --cpu BIN --cuda BIN --before-cuda BIN --output DIR` | the review's supplement | CUDA build, and the CUDA binary from before the review |
| `check_review_baseline.py --baseline BIN --results DIR` | the review's baseline: reruns the CPU inputs of an earlier `check_acceleration`/`check_cce_pulses`/`check_pcce` results folder with a frozen pre-review binary | a CPU binary (`--baseline`) and that results folder |
| `check_host_leaks.py --binary BIN --baseline BIN --output DIR [--cuda BIN]` | host memory over repeated runs | a CPU build is enough (the current and a baseline binary); `--cuda` adds the CUDA variants |

## Unit tests: which command builds and runs which

After `bash do_compile.sh <preset>` has written the Makefile (or after a CMake configure):

| tests | built and run by | needs |
|---|---|---|
| `test_scaled_product`, `test_hamiltonian_check`, `test_small_cluster_host` (`test_small_cluster_cuda.cpp` with `-DSMALL_CLUSTER_HOST_ONLY`: the packed arithmetic on the 11 dumps, the dump reader on cut copies, Hamiltonians of scale 1e±200 that must come out right or go back to the CPU) | `make numeric-test` (CMake: `cmake --build build --target numeric-test`) | a CPU |
| `test_cluster_log`, `test_review_controls`, `test_mpi_seed` (1 and 4 ranks), `test_host_lifetimes` (1 and 4 ranks) | `make cpu-test` (CMake: `--target cpu-test`; also in the CI workflow) | a CPU and `mpirun` (`make cpu-test MPIRUN=...` picks another; CMake uses the launcher FindMPI found) |
| the numeric two and `test_gcce_cuda`, `test_cce_cuda`, `test_gcce_accuracy` (twice: the built-in case and `data/gcce_scale_separation_three_spins.txt`), `test_gcce_refinement` (the parallel refinement at dim 192, 15, 2 and 1), `test_gcce_batch`, `test_review_failure`, `test_review_memory` (after an out-of-memory part way, `gpuOutOfMemoryRecovery` empties the cache, a work released afterwards is freed, and caching resumes), `test_small_cluster_cuda`, `test_warp_cluster_cuda` (both below) | `make gpu-test` (CMake: `--target gpu-test`) | a CUDA build and the GPU, so inside a job. `server/build.sbatch` builds `test_cce_cuda` and `test_gcce_cuda`, and `server/measure.sbatch` runs only those two. |

- `test_small_cluster_cuda`: the packed kernel of `execution.gpu_options.packed_max_dim`. It runs the same arithmetic on the host and on the card against the CPU coherences of the 11 files in `data/small_cluster/`, 1e-8 + 2e-6 \|CPU\| per point, both lanes the same bits. `data/small_cluster/make_dumps.py` writes the files again with a CPU build, byte for byte. `small_cluster_host_emulation.cpp` stands in for the card so the node feeder and shared-memory windows can be run end to end without CUDA (see its header).
- `test_warp_cluster_cuda`: the warp kernel of `execution.gpu_options.warp_max_dim`, through the batch API on the CPU-dumped clusters of the 15 files in `data/warp_cluster/`. It covers every compiled dimension (gCCE 6 to 27, cCCE 8 and 16), single state and ensemble, coherence and dm, Hahn, CPMG, Ramsey, and 64 production P1 clusters; gCCE dimension 12 Hahn goes through the pulse-in-eigenbasis variant. Tolerance 1e-8 + 2e-6 \|CPU\| per value, both lanes the same bits. A file cut short anywhere (header, plan, body, one byte) is refused, and every file is also read cut in memory at those places, as in `test_small_cluster_cuda`. `data/warp_cluster/make_dumps.py` writes the files with a CPU build. The host emulation has no warp kernel.

`cpu-test` compiles its four tests with the build's own flags and `-UNDEBUG` for their asserts.
`test_host_lifetimes` links the objects of `bin/main.out` except `main.o` (CMake compiles those
sources again for it). Build these tests through the targets: built by hand without `-march=native`
against objects built with it, `test_host_lifetimes` crashed (exit 139).

## Notes

- `check_review_controls.py --cpu` passed on the CPU build (2026-09-24).
- `check_seed_states.py` writes `bundle.max_dim` in its adaptive cases; its nine cases pass on the CPU build.
  One run in two failed in `benchmark.parse_profiles` because MPI interleaved two ranks' `[profile]` lines in
  the p1 case. That is the output forwarding, not the code; a rerun passed.
- `benchmark.py` is the helper the checks and `../gpu_benchmarks/` import: it runs a build on the four standard
  calculations of `code_verification/check_input_forms/inputs` and reads what it wrote.
- `examples/` holds the P1 and pCCE inputs the checks use. `data/` holds the gCCE accuracy test's files,
  the 11 `.bin` dumps of `test_small_cluster_cuda` in `data/small_cluster/` with `make_dumps.py`, and the
  `data/warp_cluster/` files.
- The seven `validation-*.json` files are the results recorded on idm4-gpu and the RTX 3050 on 2026-09-08..11.
