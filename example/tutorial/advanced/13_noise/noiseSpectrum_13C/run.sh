#!/bin/sh
# advanced/13_noise/noiseSpectrum_13C: method noise with compute spectrum on basic/01_first_run's 13C bath: C(t), S(f) and the
# coherence, against the coherence of compute coherence
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in spectrum coherence; do
  mpirun -n $NP $CCEX -f ccein_$name.jsonc -o results > results/ccein_$name.log
done
grep -h "hamiltonian.noise.fft:" results/ccein_spectrum.log || true    # the grid the run chose
python3 $PLOT results/coherence_only.png "compute coherence=results/ccein_coherence_wiDiv" \
    "compute spectrum=results/ccein_spectrum_wiDiv"
# one figure: the two coherences, C(t) over its first 10 us, S(f) up to 1.5 MHz
python3 - <<'PY'
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
L = lambda f: np.array([[float(l.split('\t')[0]), abs(complex(l.split('\t')[1].strip()))] for l in open(f) if l.strip()])
c = np.loadtxt('results/ccein_spectrum_noise_correlation')
s = np.loadtxt('results/ccein_spectrum_noise_spectrum')
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for name, style in (('coherence', '-'), ('spectrum', '--')):
    d = L(f'results/ccein_{name}_wiDiv'); ax[0].plot(d[:, 0], d[:, 1], style, label=f'compute {name}')
ax[0].set_xlabel('time (ms)'); ax[0].set_ylabel('|L(t)|'); ax[0].legend()
m = c[:, 0] <= 0.01
ax[1].plot(c[m, 0] * 1e3, c[m, 1]); ax[1].set_xlabel('t (us)'); ax[1].set_ylabel('C(t) (kHz$^2$)')
m = (s[:, 0] > 0) & (s[:, 0] < 1500)
ax[2].semilogy(s[m, 0], np.abs(s[m, 1])); ax[2].set_xlabel('f (kHz)'); ax[2].set_ylabel('|S(f)| (kHz$^2$/kHz)')
fig.tight_layout(); fig.savefig('results/coherence.png', dpi=100)
print('wrote results/coherence.png')
PY
