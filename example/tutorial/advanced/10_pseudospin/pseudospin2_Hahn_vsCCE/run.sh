#!/bin/sh
# advanced/10_pseudospin/pseudospin2_Hahn_vsCCE: V_B- in h-BN at 3 T, the pairs' part of the Hahn echo by method pseudospin
# (channels grouped and pairwise) times cCCE1 against cCCE2
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in cce2 cce1 grouped pairwise; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
# pseudospin gives the pairs alone; times cCCE1 (the single spins' ESEEM) it is comparable with cCCE2
python3 - <<'PY'
rows = lambda f: [l.split('\t') for l in open(f) if l.strip()]
for name in ('grouped', 'pairwise'):
    with open(f'results/{name}_x_cce1_wiDiv', 'w') as out:
        for a, b in zip(rows(f'results/ccein_{name}_wiDiv'), rows('results/ccein_cce1_wiDiv')):
            q = complex(a[1].strip()) * complex(b[1].strip())
            out.write(f"{a[0]}\t{q.real:+.10f}{q.imag:+.10f}j\t\n")
PY
python3 $PLOT results/coherence.png --unit us "cCCE2=results/ccein_cce2_wiDiv" \
    "pseudospin (grouped) x cCCE1=results/grouped_x_cce1_wiDiv" \
    "pseudospin (pairwise) x cCCE1=results/pairwise_x_cce1_wiDiv" \
    "pseudospin (grouped), the pairs alone=results/ccein_grouped_wiDiv"
