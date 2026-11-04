#!/bin/sh
# advanced/13_server_runs/cCCE2_Hahn_slurm: one state of one configuration on this machine, to size the job first.
# Prints the memory estimate per rank and the time of the state from the log. Usage: NP=4 bash estimate.sh
set -e
ROOT=$(cd ../../../../.. && pwd)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results/estimate
python3 "$ROOT/example/tools/InputKeys/set_input_keys.py" ccein.jsonc results/estimate/one.json \
    "configuration.k_range=[1, 1]" "configuration.s_range=[1, 1]" "configuration.seed=[1]"
mpirun -n $NP "$CCEX" -f results/estimate/one.json -o results/estimate > results/estimate/one.log
grep -E "Clusters per order|Memory estimate|State done" results/estimate/one.log
