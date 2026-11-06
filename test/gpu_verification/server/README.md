# server/

Slurm jobs that build CCEX on idm4-gpu and compare the CPU, forced CUDA and automatic backends on the eight inputs
`check_input_forms/inputs/Diamond_*/inputfiles/ccein_nstate{0,1}.json` (NV + ¹³C and NV + P1, CCE and gCCE,
ensemble and single state). The jobs change only backend, profiling, output, order and `savemode`.
`validate_inputs.py` checks the rest before running.

```bash
mkdir server-results                                   # on idm, in an isolated copy of the source with Eigen, uthash and the inputs
SbatchPath test/gpu_verification/server/build.sbatch   # prints "Submitted batch job N"
SbatchPath --dependency=afterok:N test/gpu_verification/server/measure.sbatch
python3 test/gpu_verification/server/test_run_orders.py   # the runner itself, no GPU or Slurm, under a second; also in CI
bash test/gpu_verification/server/test_mps.sh             # the MPS start and stop, no GPU
```

## How to submit

- idm4-gpu is the `gpu` partition of idm's Slurm, so submit on idm with `SbatchPath`, never `sbatch` directly.
- The jobs carry their own `--qos` line. Do not run `QosChange` on them, because it would also rewrite their
  partition, `--out`, task and CPU lines.
- Every job has `--export=NONE`. `environment.sh` starts from `module purge` and loads GCC 11.5.0, Intel oneAPI
  2023.2 (icpx, Intel MPI, MKL), CUDA 13.2 and Python 3.13, so `mpirun` is the Intel MPI of the build.
  The card has compute capability 12.0 (`CUDA_ARCH=120`).
- `mps.sh` starts CUDA MPS for the measurement, order and `min_dim` jobs, and stops it keeping the job's exit
  status. `CCEX_MPS=0` leaves it off.
- The jobs write in the source directory on the shared home. A large run may work in `/scratch/$USER/<job>` (local,
  no backup) and copy what it keeps back. The node has no outside network.

## What it checks

| Job | Resources | Does | Passes when |
|---|---|---|---|
| `build.sbatch` | 2 CPUs, 8 GiB, 20 min, no GPU | `bash do_compile.sh gpu` builds `bin/main.out`, `bin/main_gpu.out` and `test_cce_cuda`, `test_gcce_cuda`. Both binaries run with `backend=cpu` on the eight cases. | The two binaries agree. |
| `measure.sbatch` | 1 CPU, 1 GPU, 8 GiB, 10 min | The two CUDA unit tests. CPU against forced CUDA and `auto`, three repetitions, median wall time, one rank. | Comparisons pass, and forced CUDA shows GPU computation. |
| `order.sbatch ORDER` | whole node exclusive: 32 ranks, 1 GPU, `--mem=0`, `gpu-long` (72 h) | `run_orders.py`: CPU, forced CUDA and `auto` at order 3, 4 or 5 on the eight cases, one run per backend. | Every comparison passes, and a non-finite CPU result is rejected first. |
| `min_dim.sbatch` | 8 ranks, exclusive, `gpu-long` | `../gpu_benchmarks/benchmark_min_dim.py`: CPU against forced CUDA per cluster dimension, ten conditions. | Timing, no verdict. |
| `gpu_sharing.sbatch` | 32 ranks, `gpu-short` | `../gpu_benchmarks/benchmark_gpu_sharing.py`: the same with 1, 4, 8, 16 and 32 ranks on one GPU, without MPS. | Timing, no verdict. |

`auto` keeps `min_dim` 256, so gCCE on these fixtures stays on the CPU. For cCCE, `auto` sends bath dimension 2 and
4 to the packed kernel only from 40000 such clusters per rank per state (20000 where the CPU batch does not take the
run), never for the ensemble, so these fixtures leave it off. The `[gpu-packed]` line says why.

Submit order 3 with `afterok` on a passed `measure`, then 4 after 3 and 5 after 4, with `--kill-on-invalid-dep=yes`.
Time batches start at four to bound GPU memory per rank. The driver samples node RAM and GPU memory every two seconds and stops a case
when free RAM drops below 8 GiB or GPU use exceeds 90 percent. Any error, guard, timeout or failed comparison stops
the progression.

Outputs go to `server-results/<job>-JOB.log` and `server-results/<job>-JOB/` (`summary.json`, `environment.json`,
every input, log and raw output). `source-manifest.json` holds the hashes of the transferred files. The records of
past runs are in git history.

## Files

| File | Contents |
|---|---|
| `build.sbatch`, `measure.sbatch`, `order.sbatch`, `min_dim.sbatch`, `gpu_sharing.sbatch` | The jobs. |
| `environment.sh` | `module purge` and the modules the jobs use. |
| `mps.sh`, `test_mps.sh` | CUDA MPS for a job, and its test without a GPU. |
| `mpi_launch.sh` | `mpirun` bound to cores. |
| `run_orders.py`, `test_run_orders.py` | The high-order runner and its test. |
| `validate_inputs.py` | Checks the eight inputs before a run. |
| `record_environment.py` | Saves the allocated environment beside the results. |
