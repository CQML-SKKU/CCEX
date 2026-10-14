# CCE2 validation on `ssh gpu`

These Slurm jobs target `idm4-gpu`, using its installed GCC 11.5.0, Open MPI
4.1.5, CUDA 12.8 and Python 3.13.10 modules. The NVIDIA RTX PRO 6000 Blackwell
Server Edition reports compute capability 12.0, so the CUDA build uses
`CUDA_ARCH=120`. CPU and CUDA executables are both built on that server.

The inputs are the eight existing files under
`test/code_verification/check_input_forms/inputs/Diamond_*/inputfiles/ccein_nstate{0,1}.json`:

| Input family | Methods | States | Order | rbath | rdip | nstep | deltat |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| NV–13C | CCE, gCCE | ensemble, single | 2 | 30 Å | 8 Å | 30 | 0.08 ms |
| NV–P1 | CCE, gCCE | ensemble, single | 2 | 200 Å | 100 Å | 40 | 0.0003 ms |

Both families retain `pulse.count` 1 and `seed=1`. The original qubit, gyro, bath and
defect data are copied with the source snapshot. The benchmark changes backend,
profiling, output paths and `savemode=coherence`; it does not select fewer bath sites
or raise the order. `validate_inputs.py` checks the expected order, radii, time
point count, state count, pulse count, seed and method before execution.
`test_run_orders.py` checks the runner itself (the physics is passed through unchanged, failures and limits are
reported, memory and time limits stop a run) with mocked samplers and small child processes: no GPU and no Slurm,
under a second, `python3 test/gpu_verification/server/test_run_orders.py`. CI runs it; run it before submitting a
changed runner.

## Jobs

Create an isolated source directory on the server, including the bundled Eigen
and uthash headers and all eight input files with their referenced data. Create
`server-results/` before submission so Slurm can open its log files. Submit from
the source directory with `SbatchPath`, the wrapper these servers use (it calls
`sbatch` with the same arguments and records the job in `~/JOBpath`;
`sbatch` is not called directly). Since 2026-10-13 idm4-gpu is the `gpu` partition of idm's Slurm: submit on idm
(job ids are idm's; jobs before that are in `~/JOBpath.idm4-gpu`). The jobs carry their own `--qos` line, so
`QosChange` is not run on them: it would also rewrite their partition, `--out`,
task and CPU lines (for a GPU job it writes one task with 32 CPUs). Every job has `--export=NONE`, so idm's login
modules (intel, anaconda) do not reach the GPU node; `environment.sh` starts from `module purge`. The GPU node has its own
`/scratch` (a local NVMe, `/scratch/$USER`; not seen from idm, no RAID or backup, no automatic deletion): these jobs
still build and write in the submitted source directory on the shared home, which their outputs are small enough for; a
large run can work in `/scratch/$USER/<job>` inside the job and copy what it keeps back to home before it ends (never
leave the only copy there). The node has no outside network (packages come through idm; see the group's server notes). Job 17391
(2026-10-13, commit a72bba1b9, 32 ranks, MPS) checked the new setup: the card visible, `make gpu-test` 185/0,
`check_gpu/check_options.py` 95/0, P1 CR0 cCCE2 `main_gpu.out` against `main.out` within 2.2e-10. Job 17401 ran `build.sbatch` with `--export=NONE`: the
build and the eight `backend=cpu` comparisons passed (max error 0).
Jobs 17444 and 17447 (2026-10-14, review 20261013-221715 group 4) used the `/scratch` way: `make gpu-test` 201/0, `check_options.py` 110 PASS, and
P1 gCCE2 warp runs old against new ((internal notes)).

```bash
mkdir server-results
SbatchPath test/gpu_verification/server/build.sbatch          # prints "Submitted batch job N"
SbatchPath --dependency=afterok:N test/gpu_verification/server/measure.sbatch
```

The build job requests two CPUs, 8 GiB RAM and 20 minutes. It requests no GPU;
the scheduler warning about GPU invisibility is expected for this stage. It
builds `bin/main.out` and `bin/main_gpu.out` with `bash do_compile.sh gpu` (the
preset for this server: Open MPI and CUDA 12.8 from modules, icpx when oneAPI is
installed, `CUDA_ARCH=120`) and the two CUDA unit tests the measurement job runs,
then compares both executables with `backend=cpu` on the eight cases. These are
CPU compatibility checks, not GPU measurements. `environment.sh` loads the
modules and, when oneAPI is installed, its runtime libraries, with Open MPI kept
first on `PATH`. `mps.sh`, sourced after it by the measurement, order and `min_dim` jobs, starts CUDA MPS for the job
(the ranks then share the GPU at once instead of taking turns: GPU runs of 32 ranks took 0.49-0.62 of the time in
job 698, curves byte-identical) and stops it when the job ends, keeping the job's exit status even when the stop
fails (`bash test_mps.sh` checks it without a GPU, with a stand-in control command); `CCEX_MPS=0` leaves it off. `gpu_sharing.sbatch`
measures the sharing itself and does not use it.

The measurement job requests one CPU, one GPU, 8 GiB RAM and ten minutes. It runs
the independent cCCE and gCCE CUDA tests and compares CPU with forced CUDA and
with automatic selection. Each comparison uses three repetitions per backend,
alternating execution order and reporting the median whole-process wall time.
The CPU references are rerun inside this same allocation. Thread counts are
fixed at one; each executable uses one MPI rank.

Automatic selection retains the default dimension threshold of 256. The current
CCE2 fixtures are smaller than that threshold, so gCCE `auto` profiles are expected
to show CPU computation. Since 2026-10-12 (A26) `auto` with method cCCE (coherence)
sends cCCE clusters of bath dimension 2 and 4 to the packed GPU kernel by default
(`gpu_options.packed_max_dim` 4, `Config_resolvePackedMaxDim` in `src/input/general.cpp`).
Since A48 (2026-10-13) that default is used only where it beats the CPU batch (A30), which the
GPU binary now also runs: one state from 40000 clusters of bath dimension 2 and 4 per rank per
state, never for the ensemble (a run the CPU batch does not take keeps the 20000-per-rank
threshold). The small fixtures therefore leave the packed kernel off unless `packed_max_dim` is
given, and their `[gpu-packed]` log line says why; such a run writes the bytes of `main.out`.
Forced-CUDA profiles must show GPU computation.

## Status of these jobs (2026-10-02)

`build.sbatch`, `measure.sbatch` and `order.sbatch` were moved to the current
build (`do_compile.sh gpu`, `bin/main.out`, `bin/main_gpu.out`) on 2026-10-01.
Until then they called `Makefile.accel` and `bin/main.cpu.out`/`bin/main.cuda.out`,
which no longer exist, and `run_orders.py` looked for `benchmark.py` in the wrong
folder.

They ran on idm4-gpu on 2026-10-02 from a clean copy of commit 37f841284 and passed:
- `build` (job 640): `main.out` and `main_gpu.out` agree with `backend=cpu` on the eight
  cases (difference 0).
- `measure` (job 641): the CUDA unit tests pass (39 PASS lines, no FAIL); CPU against forced
  CUDA differs by at most 1e-10, against automatic selection by 0. These cases are below the
  256 dimension threshold, so forced CUDA takes 2-5 times the CPU time.
- `order.sbatch 3` (job 644, 4 minutes): status passed; CPU against forced CUDA differs by 0 on
  five of the eight cases and at most 2.2e-6 (P1 CCE, ensemble). With 32 ranks on one GPU and
  clusters of dimension 192 at most, forced CUDA takes 11-45 times the CPU time.

Orders 4 and 5 have not been run.

`make gpu-test` (not one of these jobs) passed 84/0 in idm4-gpu job 692 on 2026-10-12.

## Removed jobs

The other 19 jobs here were one-off experiments of 2026-09-08..16, and all of them
called the removed build interface: the pCCE build, measurement, review and bath
checks (`pcce_*`), the partition cache, search, cycle and reference benchmarks,
the review fixes and supplement (`review_fixes`, `review_supplement`), host leaks,
rank scaling, the PL6 acceleration benchmark and the three bundle-clustering
validations (`adaptive_*validation`). They were removed on 2026-10-01 with
`paths.example.sh`, which only they read; git history up to `ad358ed0c` has them.
The unit tests some of them compiled by hand are described, with the commands
that build them now, in `../check_gpu/README.md`. Earlier, the jobs of the
2026-09-12..16 experiments (28 files) were removed on 2026-09-24; git history up
to 8f2bcbd75 has them.

## Results and status

Job output remains under the isolated source directory:

- `server-results/build-JOB.log`: build and CPU compatibility check output.
- `server-results/build-JOB-cpu-checks/summary.json`: CPU comparisons.
- `server-results/measure-JOB.log`: numerical checks and timing summaries.
- `server-results/measure-JOB/environment.json`: allocated environment and binary hashes.
- `server-results/measure-JOB/cuda/summary.json`: CPU versus forced CUDA.
- `server-results/measure-JOB/auto/summary.json`: CPU versus automatic selection.

Two measurement jobs for `gpu_options.min_dim` (A3 of (internal notes)) build with `do_compile.sh gpu` themselves:

- `min_dim.sbatch` (8 ranks, gpu-long): `gpu_benchmarks/benchmark_min_dim.py`, CPU against forced CUDA per cluster dimension for ten
  conditions; `server-results/min-dim-JOB.log` and `min-dim-JOB/conditions/summary.json` (jobs 648, 659).
- `gpu_sharing.sbatch` (32 ranks, gpu-short): `gpu_benchmarks/benchmark_gpu_sharing.py`, the same comparison with 1, 4, 8, 16 and 32
  ranks on the one GPU, enough clusters per rank to keep them busy; `server-results/sharing-JOB.log` and `sharing-JOB/conditions/summary.json`.

Each benchmark directory also retains every generated input, full run log,
profile and raw `_noDiv`/`_wiDiv` output. A `source-manifest.json` in the source
directory records hashes of the transferred files, including the original inputs.

Use `squeue -j JOB` for a queued/running job and
`sacct -j JOB --format=JobID,State,ExitCode,Elapsed,MaxRSS` for completion and
accounting. A submitted or pending job is not a completed GPU validation.

## Deferred whole-server orders 3, 4 and 5

`order.sbatch` takes one order as its positional argument. It requests the entire
node exclusively: 32 MPI ranks, one CPU per rank, one GPU and all schedulable RAM
(`--mem=0`). Existing jobs retain their allocations; each new job waits for Slurm.
All 32 CUDA ranks share the one GPU. This measures that configuration's actual
performance; reserving the entire server does not guarantee better GPU scaling.

Submit order 3 with an `afterok` dependency on the successful CCE2 measurement,
then order 4 after order 3, then order 5 after order 4. Use
`--kill-on-invalid-dep=yes` so a failed prerequisite cancels its dependent job.
Each order uses the account's allowed `gpu-long` QoS with a 72-hour Slurm limit
(the default `gpu-short` QoS permits only 24 hours). The first high-order comparison uses one
run per backend and case, rather than the CCE2 suite's three repetitions.

`run_orders.py` keeps all eight original input conditions and changes only order
and execution/output settings. Both CPU and GPU use 32 ranks. It runs CPU, forced
CUDA and automatic selection on each case, reusing the same CPU reference for
the two comparisons. The first order-3 job first checks a 32-rank CPU/CUDA launch
using the small gCCE2 P1 ensemble fixture. Every comparison must pass before the
job can complete successfully and release the next dependency.

Time batches start at four to bound per-rank GPU storage. For example, a P1 gCCE5
ensemble cluster has dimension 3072; its two main time-batch buffers occupy
1.125 GiB per rank at batch four, or 36 GiB across 32 ranks. Eigensystems, initial
states, solver workspaces and CUDA contexts require additional memory. No change
is made to the bath, `rbath`, `rdip`, pulse sequence, time grid or seed.

The driver records whole-node RAM and GPU memory every two seconds and aborts
its case if available host RAM drops below 8 GiB or GPU usage exceeds 90 percent.
Samples are approximate observations, not exact per-process peaks or a guarantee
against an allocation failure between samples. A calculation error, resource
guard failure, timeout or failed numerical comparison stops the progression.

Outputs are under `server-results/orderORDER-JOB/comparisons/`. `summary.json`
records progress, per-case times, errors, sampled memory peaks and actual summed
cluster counts by Hilbert-space dimension. Each run retains its generated input,
full log and memory samples. A CPU result containing non-finite values is rejected
before running the two GPU-executable variants.
