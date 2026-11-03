# gpu_benchmarks/

Measures how much faster the CUDA backend is, not whether it is right. The scripts time CPU against CUDA runs and
import `../check_gpu/benchmark.py`, the helper the checks share. All need a CUDA build and an idle GPU.

```
SbatchPath test/gpu_verification/server/min_dim.sbatch       # on idm; runs benchmark_min_dim.py on idm4-gpu
SbatchPath test/gpu_verification/server/gpu_sharing.sbatch   # runs benchmark_gpu_sharing.py
python3 test/gpu_verification/gpu_benchmarks/test_backend_check.py   # the backend check, no GPU; also in CI
```

The other scripts have no job.

## What it checks

| Script | Times |
|---|---|
| `benchmark_acceleration.py` | gCCE4 clusters by dimension, CPU against CUDA |
| `benchmark_pcce_pl6.py`, `benchmark_pl6_acceleration.py` | the PL6 pCCE case; `benchmark_pcce_pl6.py` reads coherence files through `example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py` |
| `benchmark_rank_scaling.py` | wall time against the number of MPI ranks |
| `benchmark_min_dim.py` | CPU against forced CUDA per cluster dimension on cut P1 inputs, where `backend=auto` should switch (`gpu_options.min_dim`); one rank per GPU in eight of its ten conditions, eight ranks in `cce_ranks8` and `gcce_ranks8` |
| `benchmark_gpu_sharing.py` | the same with 1 to 32 ranks sharing one GPU |

A `benchmark_min_dim.py` run whose profile shows another backend than asked stops, except `cpu-batch` in a
`backend=cpu` run; `test_backend_check.py` checks this.

## Files

| File | Contents |
|---|---|
| `benchmark_*.py` | the timing scripts above |
| `test_backend_check.py` | the unit test of the backend check in `benchmark_min_dim.py` |
