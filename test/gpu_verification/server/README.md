# CCE2 validation on `ssh gpu`

These Slurm jobs target `idm4-gpu`, using its installed GCC 11.5.0, Intel oneAPI
2023.2 (`intel/oneapi-2023.2`: icpx, Intel MPI 2021.10, MKL), CUDA 13.2 and Python 3.13.10 modules.
The NVIDIA RTX PRO 6000 Blackwell Server Edition reports compute capability 12.0, so the CUDA build uses
`CUDA_ARCH=120`. CPU and CUDA executables are both built on that server.

The inputs are the eight existing files under
`test/code_verification/check_input_forms/inputs/Diamond_*/inputfiles/ccein_nstate{0,1}.json`:

| Input family | Methods | States | Order | rbath | rdip | nstep | deltat |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| NV–13C | CCE, gCCE | ensemble, single | 2 | 30 Å | 8 Å | 30 | 0.08 ms |
| NV–P1 | CCE, gCCE | ensemble, single | 2 | 200 Å | 100 Å | 40 | 0.0003 ms |

Both families keep `pulse.count` 1 and `seed=1`. The original qubit, gyro, bath and
defect data are copied with the source snapshot. The benchmark changes backend,
profiling, output paths and `savemode=coherence`; it does not select fewer bath sites
or raise the order. `validate_inputs.py` checks the expected order, radii, time
point count, state count, pulse count, seed and method before execution.

`test_run_orders.py` checks the runner itself (the physics is passed through unchanged, failures and limits are
reported, memory and time limits stop a run) with mocked samplers and small child processes. It needs no GPU and
no Slurm and takes under a second: `python3 test/gpu_verification/server/test_run_orders.py`. CI runs it; run it
before submitting a changed runner.

## Submitting

Create an isolated source directory on the server, including the bundled Eigen
and uthash headers and all eight input files with their referenced data. Create
`server-results/` before submission so Slurm can open its log files. Submit from
the source directory with `SbatchPath`, the wrapper these servers use (it calls
`sbatch` with the same arguments and records the job in `~/JOBpath`; do not call
`sbatch` directly). idm4-gpu is the `gpu` partition of idm's Slurm, so submit on idm.
Job ids from before 2026-10-13 (640 to 698 below) are idm4-gpu's own, recorded in `~/JOBpath.idm4-gpu`.

```bash
mkdir server-results
SbatchPath test/gpu_verification/server/build.sbatch          # prints "Submitted batch job N"
SbatchPath --dependency=afterok:N test/gpu_verification/server/measure.sbatch
```

- The jobs carry their own `--qos` line, so do not run `QosChange` on them: it would also rewrite their
  partition, `--out`, task and CPU lines (for a GPU job it writes one task with 32 CPUs).
- Every job has `--export=NONE`, so idm's login modules (intel, anaconda) do not reach the GPU node;
  `environment.sh` starts from `module purge` and loads the modules above, so the jobs' `mpirun` is the
  Intel MPI the binaries are built with.
- The GPU node has its own `/scratch` (a local NVMe, `/scratch/$USER`; not seen from idm, no RAID or backup,
  no automatic deletion). These jobs build and write in the submitted source directory on the shared home,
  which is enough for their small outputs. A large run can work in `/scratch/$USER/<job>` inside the job and
  copy what it keeps back to home before it ends (never leave the only copy there).
- The node has no outside network (packages come through idm; see the group's server notes).
- `mps.sh`, sourced after `environment.sh` by the measurement, order and `min_dim` jobs, starts CUDA MPS for
  the job, so the ranks share the GPU at once instead of taking turns (GPU runs of 32 ranks took 0.49 to 0.62
  of the time in job 698, curves byte-identical). It stops MPS when the job ends and keeps the job's exit
  status even when the stop fails. `bash test_mps.sh` checks this without a GPU, with a stand-in control
  command. `CCEX_MPS=0` leaves MPS off. `gpu_sharing.sbatch` measures the sharing itself and does not use it.

## The build job

`build.sbatch` requests two CPUs, 8 GiB RAM and 20 minutes. It requests no GPU;
the scheduler warning about GPU invisibility is expected for this stage. It
builds `bin/main.out` and `bin/main_gpu.out` with `bash do_compile.sh gpu` (the
preset for this server: `mpiicpc` with icpx, MKL and CUDA 13.2 from modules,
`CUDA_ARCH=120`) and the two CUDA unit tests the measurement job runs,
then compares both executables with `backend=cpu` on the eight cases. These are
CPU compatibility checks, not GPU measurements. The other unit tests and the commands that build them
are in `../check_gpu/README.md`.

## The measurement job

`measure.sbatch` requests one CPU, one GPU, 8 GiB RAM and ten minutes. It runs
the independent cCCE and gCCE CUDA tests and compares CPU with forced CUDA and
with automatic selection. Each comparison uses three repetitions per backend,
alternating execution order and reporting the median whole-process wall time.
The CPU references are rerun inside the same allocation. Thread counts are
fixed at one; each executable uses one MPI rank.

Automatic selection keeps the default dimension threshold of 256. The CCE2 fixtures are
smaller than that, so gCCE `auto` profiles are expected to show CPU computation.
With method cCCE (coherence), `auto` sends cCCE clusters of bath dimension 2 and 4 to the packed GPU kernel
by default (`execution.gpu_options.packed_max_dim` 4, `Config_resolvePackedMaxDim` in `src/input/general.cpp`),
but only where it beats the CPU batch, which the GPU binary also runs: one state from 40000 clusters of bath
dimension 2 and 4 per rank per state, never for the ensemble (a run the CPU batch does not take keeps the
20000-per-rank threshold). The small fixtures therefore leave the packed kernel off unless `packed_max_dim` is
given, and their `[gpu-packed]` log line says why; such a run writes the bytes of `main.out`.
Forced-CUDA profiles must show GPU computation.

## Last runs

- 2026-10-02, clean copy of commit 37f841284:
  - `build` (job 640): `main.out` and `main_gpu.out` agree with `backend=cpu` on the eight cases (difference 0).
  - `measure` (job 641): the CUDA unit tests pass (39 PASS lines, no FAIL). CPU against forced CUDA differs by
    at most 1e-10, against automatic selection by 0. These cases are below the 256 dimension threshold, so
    forced CUDA takes 2 to 5 times the CPU time.
  - `order.sbatch 3` (job 644, 4 minutes): passed. CPU against forced CUDA differs by 0 on five of the eight
    cases and at most 2.2e-6 (P1 CCE, ensemble). With 32 ranks on one GPU and clusters of dimension 192 at
    most, forced CUDA takes 11 to 45 times the CPU time.
- Orders 4 and 5 have not been run.
- `make gpu-test` (not one of these jobs) passed 84/0 in job 692 (2026-10-12).
- Job 17391 (2026-10-13, commit a72bba1b9, 32 ranks, MPS), first on idm's `gpu` partition: the card visible,
  `make gpu-test` 185/0, `check_gpu/check_options.py` 95/0, P1 CR0 cCCE2 `main_gpu.out` against `main.out`
  within 2.2e-10. Job 17401 ran `build.sbatch` with `--export=NONE`: the build and the eight `backend=cpu`
  comparisons passed (max error 0).
- Jobs 17444 and 17447 (2026-10-14) used the `/scratch` way: `make gpu-test` 201/0, `check_options.py` 110 PASS,
  and P1 gCCE2 warp runs old against new.
- Job 17708 (2026-10-21, CUDA 13.2 and Intel MPI, in a two-CPU step): `build.sbatch` passed in 788 s, at most
  4.9 GB of memory with the page cache, so its 20 minutes and 8 GiB still hold; `make gpu-test` passed 203/0.

Use `squeue -j JOB` for a queued or running job and
`sacct -j JOB --format=JobID,State,ExitCode,Elapsed,MaxRSS` for completion and
accounting. A submitted or pending job is not a completed GPU validation.

## Results

Job output stays under the isolated source directory:

- `server-results/build-JOB.log`: build and CPU compatibility check output.
- `server-results/build-JOB-cpu-checks/summary.json`: CPU comparisons.
- `server-results/measure-JOB.log`: numerical checks and timing summaries.
- `server-results/measure-JOB/environment.json`: allocated environment and binary hashes.
- `server-results/measure-JOB/cuda/summary.json`: CPU versus forced CUDA.
- `server-results/measure-JOB/auto/summary.json`: CPU versus automatic selection.

Each benchmark directory also keeps every generated input, full run log,
profile and raw `_noDiv`/`_wiDiv` output. A `source-manifest.json` in the source
directory records hashes of the transferred files, including the original inputs.

## Measurement jobs for `execution.gpu_options.min_dim`

Both build with `do_compile.sh gpu` themselves.

- `min_dim.sbatch` (8 ranks, gpu-long): `gpu_benchmarks/benchmark_min_dim.py`, CPU against forced CUDA per
  cluster dimension for ten conditions; `server-results/min-dim-JOB.log` and
  `min-dim-JOB/conditions/summary.json` (jobs 648, 659).
- `gpu_sharing.sbatch` (32 ranks, gpu-short): `gpu_benchmarks/benchmark_gpu_sharing.py`, the same comparison
  with 1, 4, 8, 16 and 32 ranks on the one GPU, enough clusters per rank to keep them busy;
  `server-results/sharing-JOB.log` and `sharing-JOB/conditions/summary.json`.

## Whole-server orders 3, 4 and 5

`order.sbatch` takes one order as its positional argument. It requests the entire
node exclusively: 32 MPI ranks, one CPU per rank, one GPU and all schedulable RAM
(`--mem=0`). Existing jobs keep their allocations; each new job waits for Slurm.
All 32 CUDA ranks share the one GPU. This measures that configuration's actual
performance; reserving the entire server does not guarantee better GPU scaling.

Submit order 3 with an `afterok` dependency on the successful CCE2 measurement,
then order 4 after order 3, then order 5 after order 4. Use
`--kill-on-invalid-dep=yes` so a failed prerequisite cancels its dependent job.
Each order uses the account's allowed `gpu-long` QoS with a 72-hour Slurm limit
(the default `gpu-short` QoS permits only 24 hours). The first high-order comparison uses one
run per backend and case, not the CCE2 suite's three repetitions.

`run_orders.py` keeps all eight original input conditions and changes only order
and execution/output settings. Both CPU and GPU use 32 ranks. It runs CPU, forced
CUDA and automatic selection on each case, reusing the same CPU reference for
the two comparisons. The first order-3 job first checks a 32-rank CPU/CUDA launch
using the small gCCE2 P1 ensemble fixture. Every comparison must pass before the
job can complete successfully and release the next dependency.

Time batches start at four to bound per-rank GPU storage. For example, a P1 gCCE5
ensemble cluster has dimension 3072; its two main time-batch buffers occupy
1.125 GiB per rank at batch four, or 36 GiB across 32 ranks. Eigensystems, initial
states, solver workspaces and CUDA contexts need additional memory. No change
is made to the bath, `rbath`, `rdip`, pulse sequence, time grid or seed.

The driver records whole-node RAM and GPU memory every two seconds and aborts
its case if available host RAM drops below 8 GiB or GPU usage exceeds 90 percent.
Samples are approximate observations, not exact per-process peaks or a guarantee
against an allocation failure between samples. A calculation error, resource
guard failure, timeout or failed numerical comparison stops the progression.

Outputs are under `server-results/orderORDER-JOB/comparisons/`. `summary.json`
records progress, per-case times, errors, sampled memory peaks and actual summed
cluster counts by Hilbert-space dimension. Each run keeps its generated input,
full log and memory samples. A CPU result containing non-finite values is rejected
before running the two GPU-executable variants.
