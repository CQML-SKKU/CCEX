#!/bin/bash
#SBATCH --job-name=ccex_p1_k
#SBATCH --array=1-5
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
# advanced/13_server_runs/cCCE2_Hahn_slurm: the same run split by configuration, one array task per k = 1 .. 5.
# Each task writes results/ccein_k<k>_state<s>_..., the files job.sbatch writes; the tasks may run at the same time.
# Submit from this folder: sbatch job_array.sbatch; when all five are done: bash analyse.sh
set -e
cd "$SLURM_SUBMIT_DIR"
# module load gcc openmpi
ROOT=${CCEX_ROOT:-$(cd ../../../../.. && pwd)}
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
K=$SLURM_ARRAY_TASK_ID
mkdir -p results/inputs
# configuration K alone, with its own seed (the K-th of ccein.jsonc's list); output.file stays "ccein"
python3 "$ROOT/example/tools/InputKeys/set_input_keys.py" ccein.jsonc results/inputs/ccein_k$K.json \
    "configuration.k_range=[$K, $K]" "configuration.seed=[$K]"
mpirun -n "$SLURM_NTASKS" "$CCEX" -f results/inputs/ccein_k$K.json -o results > results/ccein_k$K.log
