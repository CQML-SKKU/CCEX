#!/bin/sh
# advanced/13_pseudospin/pseudospin2_Hahn_vsCCE: V_B- in h-BN at 3 T, the pairs' part of the Hahn echo by method pseudospin
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
# the numbers (T2 and the largest difference from cCCE2) as the other tutorials print them
python3 $PLOT results/coherence_numbers.png --unit us "cCCE2=results/ccein_cce2_wiDiv" \
    "pseudospin (grouped) x cCCE1=results/grouped_x_cce1_wiDiv" \
    "pseudospin (pairwise) x cCCE1=results/pairwise_x_cce1_wiDiv" \
    "pseudospin (grouped), the pairs alone=results/ccein_grouped_wiDiv"
# the figure: the curves (left) and their difference from cCCE2 (right), since they lie within 4e-2 of each other
python3 - <<'PY2'
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def L(f):
    t, v = [], []
    for l in open(f):
        p = l.split()
        if len(p) >= 2 and not l.startswith('#'): t.append(float(p[0])); v.append(abs(complex(p[1].replace('+-', '-'))))
    return np.array(t) * 1e3, np.array(v)
runs = [('cCCE2', 'results/ccein_cce2_wiDiv', dict(lw=3, color='C0')),
        ('pseudospin (grouped) x cCCE1', 'results/grouped_x_cce1_wiDiv', dict(lw=1.2, ls='--', color='C1')),
        ('pseudospin (pairwise) x cCCE1', 'results/pairwise_x_cce1_wiDiv', dict(lw=1.2, ls=':', color='C2')),
        ('pseudospin (grouped), the pairs alone', 'results/ccein_grouped_wiDiv', dict(lw=1.2, ls='-.', color='C3'))]
t0, ref = L(runs[0][1])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
for name, f, st in runs:
    t, v = L(f); ax[0].plot(t, v, label=name, **st)
    if name != 'cCCE2': ax[1].plot(t, v - ref, label=name, **st)
ax[0].set_xlabel('time (us)'); ax[0].set_ylabel('|L(t)|'); ax[0].legend()
ax[1].axhline(0, color='C0', lw=3); ax[1].set_xlabel('time (us)'); ax[1].set_ylabel('|L| - |L(cCCE2)|'); ax[1].legend()
fig.tight_layout(); fig.savefig('results/coherence.png', dpi=100)
print('wrote results/coherence.png')
PY2
