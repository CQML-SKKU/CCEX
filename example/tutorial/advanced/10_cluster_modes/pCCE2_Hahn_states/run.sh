#!/bin/sh
# advanced/10_cluster_modes/pCCE2_Hahn_states: pCCE2 of the NV in 1 ppm P1 with partitions of 1, 2 and 3 spins, and 3 with the connections of spins
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
# partition_size 3 with the connections of spins instead of centroids
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.jsonc results/sK3_spins.json cluster.pcce.partition_size=3 \
    cluster.pcce.connections=spins
mpirun -n $NP $CCEX -f results/sK3_spins.json -o results > results/sK3_spins.log
set -- ; for k in $SIZES; do
  if [ $k = 1 ]; then l="standard CCE2 (partition_size 1)"; else l="pCCE, $k spins per partition"; fi
  set -- "$@" "$l=results/sK${k}_state*_wiDiv"
done
python3 $PLOT results/coherence.png --unit us "$@" "pCCE, 3 spins per partition, connections spins=results/sK3_spins_state*_wiDiv"
