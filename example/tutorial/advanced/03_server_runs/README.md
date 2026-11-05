# Advanced 03: Runs on a cluster (Slurm)

Folder: `example/tutorial/advanced/03_server_runs/`

A run too long for a workstation goes to a cluster as a batch job. This folder shows it with Slurm (`sbatch`) on the
system of [Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md)'s `cCCE2_Hahn_states`, with
5 bath configurations × 50 sampled states, 250 curves. That run is still small (about 20 s on 4 ranks), so the
scripts can be tried anywhere, and the same steps carry over to runs of hours. Add what your cluster needs to the
`#SBATCH` lines (`--partition`, `--account`) and to the `module load` line (the compiler and MPI CCEX was built with).

This page has one case:

1. [Sizing a run, one job, and a job array](#sr_slurm) (`cCCE2_Hahn_slurm/`).

<a id="sr_slurm"></a>
## 1. Sizing a run, one job, and a job array

Folder: `cCCE2_Hahn_slurm/`

<a id="a03_input"></a>
### 1.1 The input

The input is the P1 bath of Basic 06 with 5 configurations of 50 states, one seed per configuration and a fixed
output name:

```jsonc
"configuration" : {"gyrofile": "../../../basic/06_defect_bath/bath/DiaP1_gyro",
                   "bathfile": ["../../../basic/06_defect_bath/bath/14N_1ppm/bath_DiaP1_1ppm_{k}"],
                   "k_range": [1, 5], "s_range": [1, 50], "seed": [1, 2, 3, 4, 5]},
"output"        : {"file": "ccein"}
```

- `configuration.seed` has one entry per configuration, so a configuration run alone draws the same states.
- `output.file` gives every task the same output name.

Three scripts in the folder use it:

- `estimate.sh` runs one state of one configuration on the spot and prints three lines of the log.
- `job.sbatch` runs everything in one job, `mpirun -n $SLURM_NTASKS main.out -f ccein.jsonc -o results`, and then
  `analyse.sh`, which does the error correction, the average of the 250 curves and the T2 fit as in
  [Advanced 02: Analysis of the output](../02_analysis/README.md).
- `job_array.sbatch` splits the same run by configuration. Array task k writes an input with `k_range` [k, k] and
  `seed` [k] (`example/tools/InputKeys/set_input_keys.py`) and runs it, and the tasks may run at the same time.

In the job scripts `--ntasks` sets the MPI ranks, over which the clusters are shared. Slurm runs a copy of the
script, so each script first changes to `$SLURM_SUBMIT_DIR`, the folder it was submitted from.

<a id="a03_run"></a>
### 1.2 Running it

```sh
NP=4 bash estimate.sh
sbatch job.sbatch                     # or: sbatch job_array.sbatch, then bash analyse.sh when all five tasks are done
```

`estimate.sh` prints the clusters per order, the memory estimate per rank (0.05 GB, against the memory available per
rank) and the time of the state (0.08 s on 4 ranks). The whole run is then about 250 × 0.08 s. `CCEX_BIN` names
another binary and `CCEX_ROOT` the repository.

<a id="a03_result"></a>
### 1.3 The result

![cCCE2_Hahn_slurm: expected result](cCCE2_Hahn_slurm/expected.png)

The figure comes from the scripts run with the Slurm variables set by hand on 4 ranks. The analyzer's fit gives T2
263 µs and n 0.96, and `plot_coherence.py` 241 µs (the last time |L| ≥ 1/e). The one-job and the array run gave the
same 250 files, since the same states are drawn, so `bash analyse.sh` gives the same average after either.

Next, split the states as well. `s_range` [26, 50] computes the states 26 .. 50 of the run [1, 50]. Use `srun` in
place of `mpirun` where the cluster asks for it.

Input: [`ccein.jsonc`](cCCE2_Hahn_slurm/ccein.jsonc)
