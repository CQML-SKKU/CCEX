# gpu_benchmarks/: how much faster, not whether right

| script | times |
|---|---|
| `benchmark_acceleration.py` | gCCE4 clusters by dimension, CPU against CUDA |
| `benchmark_pcce_pl6.py`, `benchmark_pl6_acceleration.py` | the PL6 pCCE case, where the GPU helps most |
| `benchmark_rank_scaling.py` | wall time against the number of MPI ranks |
| `benchmark_min_dim.py` | CPU against forced CUDA per cluster dimension on cut P1 inputs: where `backend=auto` should switch (`gpu_options.min_dim`). One rank per GPU in eight of its ten conditions, eight ranks in `cce_ranks8` and `gcce_ranks8`. A run whose profile shows another backend than asked stops, except `cpu-batch` in a `backend=cpu` run (`test_backend_check.py`, CI python job, checks this). |
| `benchmark_gpu_sharing.py` | the same with 1 to 32 ranks sharing one GPU: how that switch point depends on the rank count |

All need a CUDA build and an idle GPU. On idm4-gpu, `../server/min_dim.sbatch` and `gpu_sharing.sbatch` run
the last two; the others have no job. They import `../check_gpu/benchmark.py`, the helper every check shares.
`benchmark_pcce_pl6.py` reads coherence files through
`example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py`.
