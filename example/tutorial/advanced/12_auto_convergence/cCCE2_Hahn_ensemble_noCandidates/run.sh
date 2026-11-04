#!/bin/sh
# advanced/12_auto_convergence/cCCE2_Hahn_ensemble_noCandidates: ../cCCE2_Hahn_ensemble without the rdip and rbath candidates: the time
# window, rdip, rbath (from their automatic starts) and the number of bath configurations found by main.out itself, then the production
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
NBATH=${NBATH:-40}
mkdir -p results
# the bath configurations results/bath/config_1..NBATH in 140 A cells (rbath up to 68 A: the automatic start, 33.6 A,
# raised twice by 1.25 needs 65.5 A, more than the 53 A of ../cCCE2_Hahn_ensemble's 110 A cells; about 10 s each the first time)
sh $ROOT/example/tutorial/make_natab_baths.sh $NBATH results/bath 140
# rank 0 runs the runner (python3 with numpy, scipy, matplotlib), which starts every CCEX run with mpirun -n NP.
# Exit status 0: converged; 3: unresolved (the summary says what is missing); calling again continues.
mpirun -n $NP $CCEX -f ccein.jsonc > results/process_auto || echo "automatic convergence ended with status $?"
cat results/auto/result/summary.txt
cp results/auto/result/6_converged_coherence.png results/coherence.png     # the runner's figure of the converged curve
