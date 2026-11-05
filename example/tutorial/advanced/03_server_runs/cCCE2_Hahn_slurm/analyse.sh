#!/bin/sh
# advanced/03_server_runs/cCCE2_Hahn_slurm: the 250 curves of results/ corrected, averaged and fitted (example/analysis),
# then plotted. Usage: bash analyse.sh [repository root]
set -e
ROOT=${1:-$(cd ../../../../.. && pwd)}
A=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
K="1 2 3 4 5"; S=$(seq 1 50)
# -uc ms in both calls: the files the analyzer writes keep the time in the -uc unit, so mixed units would compound
python3 "$A" -d results/ -fi ccein_kCONF_stateST_wiDiv -err_d results/ -err_fi ccein_kCONF_stateST_noDiv \
    -vn CONF ST -v $K -v $S -uc ms -err -ntn -fo results/ccein_kCONF_stateST_ErrCorr > results/analyse.log
python3 "$A" -d results/ -fi ccein_kCONF_stateST_ErrCorr -vn CONF ST -v $K -v $S -uc ms -ea -ntn \
    -fo results/ccein_Avg -fit_w -fit_space 5 >> results/analyse.log
echo "N, T2 (ms), n: $(cat results/ccein_Avg_params)"
MPLBACKEND=Agg python3 "$ROOT/example/tutorial/plot_coherence.py" results/coherence.png --unit us \
    "mean of 250=results/ccein_Avg_conv250" "k1 state 1=results/ccein_k1_state1_wiDiv"
