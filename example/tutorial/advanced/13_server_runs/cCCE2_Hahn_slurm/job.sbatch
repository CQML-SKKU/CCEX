#!/bin/bash
#SBATCH --job-name=ccex_p1
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
# advanced/13_server_runs/cCCE2_Hahn_slurm: the whole run (5 configurations x 50 states) in one Slurm job.
# Submit from this folder: sbatch job.sbatch. Add what your cluster needs: --partition, --account, the modules.
set -e
cd "$SLURM_SUBMIT_DIR"                       # Slurm runs a copy of this script; the folder is where sbatch was called
# module load gcc openmpi                    # the compiler and MPI that CCEX was built with
ROOT=${CCEX_ROOT:-$(cd ../../../../.. && pwd)}
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
mkdir -p results
mpirun -n "$SLURM_NTASKS" "$CCEX" -f ccein.jsonc -o results > results/ccein.log   # or: srun "$CCEX" ...
bash analyse.sh "$ROOT"
