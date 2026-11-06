#!/bin/sh
# advanced/01_auto_convergence/cCCE2_Hahn_states: ../cCCE2_Hahn_ensemble with sampled bath states (configuration.s_range):
# the runner also tests whether test_states states per configuration are enough (step 4)
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
NBATH=${NBATH:-40}
mkdir -p results
# the bath configurations results/bath/config_1..NBATH (110 A cell, rbath up to 53 A; about 1 min for 40 the first time, NP ranks)
sh $ROOT/example/tutorial/make_natab_baths.sh $NBATH results/bath
# rank 0 runs the runner (python3 with numpy, scipy, matplotlib), which starts every CCEX run with mpirun -n NP.
# Exit status 0: converged; 3: unresolved (the summary says what is missing); calling again continues.
mpirun -n $NP $CCEX -f ccein.jsonc > results/process_auto || echo "automatic convergence ended with status $?"
cat results/auto/result/summary.txt
cp results/auto/result/6_converged_coherence.png results/coherence.png     # the runner's figure of the converged curve
