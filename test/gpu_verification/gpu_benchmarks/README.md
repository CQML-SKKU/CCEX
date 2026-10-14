# gpu_benchmarks/ — how much faster, not whether right

| script | times |
|---|---|
| `benchmark_acceleration.py` | gCCE4 clusters by dimension, CPU against CUDA |
| `benchmark_pcce_pl6.py`, `benchmark_pl6_acceleration.py` | the PL6 pCCE case, where the GPU helps most |
| `benchmark_rank_scaling.py` | wall time against the number of MPI ranks |
| `benchmark_min_dim.py` | CPU against forced CUDA per cluster dimension on cut P1 inputs, one rank per GPU in eight of its ten conditions and eight ranks in `cce_ranks8`, `gcce_ranks8`: where `backend=auto` should switch (`gpu_options.min_dim`); a run whose profile shows another backend than asked stops, except `cpu-batch` (A30) in a `backend=cpu` run (review 20261013-180049 G3; `test_backend_check.py`, CI python job) |
| `benchmark_gpu_sharing.py` | the same with 1 to 32 ranks sharing one GPU: how that switch point depends on the rank count |

All want a CUDA build and an idle GPU (of the jobs that ran them on idm4-gpu, only
`../server/min_dim.sbatch` and `gpu_sharing.sbatch` remain; the rank-scaling and PL6 jobs were removed on 2026-10-01). They import `../check_gpu/benchmark.py`, the helper every check shares
(it runs a build on the four standard calculations and reads what it wrote);
`benchmark_pcce_pl6.py` reads coherence files through
`example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py`.
