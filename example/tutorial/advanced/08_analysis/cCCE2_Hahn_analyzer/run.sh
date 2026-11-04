#!/bin/sh
# advanced/08_analysis/cCCE2_Hahn_analyzer: the analyzer of example/analysis on 20 curves of a P1 bath: error correction,
# the average over configurations and states, the T2 fit
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
A=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.jsonc -o results > results/ccein.log
K="1 2"; S=$(seq 1 10)
# 1. error correction: each _wiDiv repaired against its _noDiv where they differ by the tolerance (-err_tol, default 2)
python3 $A -d results/ -fi ccein_kCONF_stateST_wiDiv -err_d results/ -err_fi ccein_kCONF_stateST_noDiv \
    -vn CONF ST -v $K -v $S -uc ms -err -ntn -fo results/ccein_kCONF_stateST_ErrCorr > results/analyzer.log
# 2. the average of the 20 corrected curves (as complex numbers) and 3. the fit exp(-(t/T2)^n) of the average.
# Keep -uc ms in both calls: the files the analyzer writes have their times in the -uc unit, and a second call reads
# them as ms.
python3 $A -d results/ -fi ccein_kCONF_stateST_ErrCorr -vn CONF ST -v $K -v $S -uc ms -ea -ntn \
    -fo results/ccein_Avg -fit_w -fit_space 5 >> results/analyzer.log
echo "N, T2 (ms), n: $(cat results/ccein_Avg_params)"
# how much the correction changed (nothing here: no curve's division blew up)
python3 - <<PY
import glob, sys
sys.path.insert(0, "$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1")
from coherence import read_coherence
n = 0
for f in sorted(glob.glob("results/ccein_k*_state*_wiDiv")):
    _, w = read_coherence(f); _, c = read_coherence(f.replace("_wiDiv", "_ErrCorr"))
    n += int((abs(w - c) > 1e-9).sum())
print("points changed by the error correction:", n)
PY
python3 $ROOT/example/tutorial/plot_coherence.py results/coherence.png --unit us \
    "average (analyzer)=results/ccein_Avg_conv20" "k1 state 1=results/ccein_k1_state1_ErrCorr"
