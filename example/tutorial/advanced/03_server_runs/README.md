# Advanced 03: Runs on a cluster (Slurm)

Folder: `example/tutorial/advanced/03_server_runs/`

A run too long for a workstation goes to a cluster as a batch job. This folder shows it with Slurm (`sbatch`) on the
production of [Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md)'s `cCCE2_Hahn_states`,
5 bath configurations × 50 sampled states, 250 curves. That run is still small (about 20 s on 4 ranks), so the
scripts can be tried anywhere, and the same steps carry over to runs of hours. Add what your cluster needs to the
`#SBATCH` lines (`--partition`, `--account`) and to the `module load` line (the compiler and MPI CCEX was built with).

1. [cCCE2_Hahn_slurm](#sr_slurm): sizing a run, one job, and a job array.

<a id="sr_slurm"></a>
## 1. cCCE2_Hahn_slurm: sizing a run, one job, and a job array

- **What it shows:** three steps.
  1. `estimate.sh` runs one state of one configuration on the spot and prints three lines of the log. They are the
     clusters per order, the memory estimate per rank (0.05 GB, against the memory available per rank) and the time of
     the state (0.08 s on 4 ranks). The whole run is then about 250 × 0.08 s.
  2. `job.sbatch` runs everything in one job, `mpirun -n $SLURM_NTASKS main.out -f ccein.jsonc -o results`, and then
     `analyse.sh`, which does the error correction, the average of the 250 curves and the T2 fit as in
     [Advanced 02: Analysis of the output](../02_analysis/README.md).
  3. `job_array.sbatch` splits the same run by configuration. Array task k writes an input with `k_range` [k, k] and
     `seed` [k] (`example/tools/InputKeys/set_input_keys.py`) and runs it, and the tasks may run at the same time.
     Their files have the names and the values of the one-job run, since the same states are drawn, so `bash
     analyse.sh` afterwards gives the same average.
- **Keys and settings that matter:**
  - `configuration.seed`, one per configuration, so that a configuration run alone draws the same states.
  - `output.file`, one output name for every task.
  - `--ntasks`, the MPI ranks, over which the clusters are shared.
  - `$SLURM_SUBMIT_DIR`. Slurm runs a copy of the script, so the script changes to the folder it was submitted from.
- **How to run:** `NP=4 bash estimate.sh`, then `sbatch job.sbatch`. Or `sbatch job_array.sbatch`, and `bash
  analyse.sh` when all five tasks are done. `CCEX_BIN` names another binary and `CCEX_ROOT` the repository.
- **Expected result:** `expected.png`, from the scripts run with the Slurm variables set by hand on 4 ranks. The
  analyzer's fit gives T2 263 µs and n 0.96, and `plot_coherence.py` 241 µs (the last time |L| ≥ 1/e). The one-job
  and the array run gave the same 250 files.
- **What to change next:** split the states as well, since `s_range` [26, 50] computes the states 26 .. 50 of the run
  [1, 50]. Use `srun` in place of `mpirun` where the cluster asks for it.

![cCCE2_Hahn_slurm: expected result](cCCE2_Hahn_slurm/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_slurm/ccein.jsonc)
