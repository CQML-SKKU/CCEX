#!/bin/sh
# basic/02_bath_and_states/cCCE2_Hahn_builtinIsotopes: the 13C gamma from a gyro file against the built-in isotope table
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in gyrofile builtin; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
# where the 13C spin and gamma came from (the run log says it per species)
grep -H -E 'built-in isotope table|Read GyroFile' results/ccein_builtin.log results/ccein_gyrofile.log || true
# the gammas differ by 2e-7 (relative): over 2 ms at 500 G about 1e-3 rad of 13C precession, so the curves differ by about 2e-5
python3 $PLOT results/coherence.png "gyro file=results/ccein_gyrofile_wiDiv" "built-in table=results/ccein_builtin_wiDiv"
