#!/bin/sh
# advanced/11_cluster_shapes/pCCE2_Hahn_states: pCCE2 of the NV in 1 ppm P1 with partitions of 1, 2 and 3 spins
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
SIZES=${SIZES:-1 2 3}      # partition_size values to run; 4 takes much longer
for k in $SIZES; do
  python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/sK$k.json cluster.pcce.partition_size=$k
  mpirun -n $NP $CCEX -f results/sK$k.json -o results > results/sK$k.log
done
set -- ; for k in $SIZES; do set -- "$@" "partition_size $k=results/sK${k}_state*_wiDiv"; done
python3 $PLOT results/coherence.png --unit us "$@"
