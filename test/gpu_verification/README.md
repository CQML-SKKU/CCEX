# gpu_verification/

Checks of the CUDA backend (`src/calculation/gpu/`, `backend` in the input). The unit tests and the CPU against CUDA
comparisons are in `check_gpu/`, the timings in `gpu_benchmarks/`, and the Slurm jobs that run them on idm4-gpu in
`server/`.

```
make gpu-test                                   # inside a GPU job; CMake: --target gpu-test
SbatchPath test/gpu_verification/server/build.sbatch   # on idm; see server/README.md
```

## What it checks

| Folder | Checks | Needs |
|---|---|---|
| `check_gpu/` | Unit tests of the kernels, batching and host bookkeeping. `check_*.py` that run CPU and CUDA builds on the same inputs (pulses, gCCE4 accuracy and batching, pCCE, seeds and sampled states, options, host leaks). | a CUDA build (`bin/main_gpu.out`) and a GPU for most, CPU builds for `check_review_controls.py --cpu`, `check_seed_states.py --cpu`, `check_host_leaks.py` without `--cuda`, `check_review_baseline.py` and `check_packed_emulation.py` |
| `gpu_benchmarks/` | Timing, not correctness: CPU against CUDA, rank scaling, the PL6 pCCE case. | a CUDA build and an idle GPU |
| `server/` | The build, CCE2 measurement and high-order jobs, and the `min_dim` and `gpu_sharing` measurement jobs. | idm4-gpu |

## Files

| Folder | Contents |
|---|---|
| `check_gpu/` | The unit tests, the comparison scripts, their inputs (`examples/`, `data/`) and recorded results. |
| `gpu_benchmarks/` | The timing scripts. |
| `server/` | The Slurm jobs and their environment. |
