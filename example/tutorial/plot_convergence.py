"""plot_convergence.py: a convergence scan in one figure, the curves (left) and how far each is from the last (right).

    python3 plot_convergence.py OUT.png --name NAME [--unit ms|us] [--tol 0.01] VALUE=GLOB [VALUE=GLOB ...]
    python3 plot_convergence.py OUT.png --name NAME [--unit ms|us] [--tol 0.01] --each PATTERN VALUE [VALUE ...]
    python3 plot_convergence.py OUT.png --name configurations [--unit ms|us] [--tol 0.01] --first GLOB N [N ...]

Each VALUE=GLOB is one value of the scanned input key (a number, in the order given): the files GLOB matches are
averaged as complex numbers, as in plot_coherence.py. The last value is the reference. --each PATTERN is the same with
the glob PATTERN, {v} replaced by each VALUE. --first takes the files of GLOB
in the order of their k (..._k<k>_...) and averages the first N of them for each N, the last N being the reference.

Printed per value: the number of files, the largest |L| difference from the reference and from the previous value.
The right panel draws the difference from the reference against the value, on a log scale, with the line --tol
(default 0.01). Needs numpy and matplotlib.
"""
import glob
import re
import sys

import numpy as np

from plot_coherence import read


def mean_of(files):
    data = [read(f) for f in files]
    return data[0][0], np.abs(np.nanmean(np.array([d[1] for d in data]), axis=0))


def main(argv):
    out, name, unit, tol, first, each, pairs = None, 'value', 'ms', 0.01, None, None, []
    it = iter(argv)
    for a in it:
        if a == '--name':
            name = next(it)
        elif a == '--unit':
            unit = next(it)
        elif a == '--tol':
            tol = float(next(it))
        elif a == '--first':
            first = next(it)
        elif a == '--each':
            each = next(it)
        elif out is None:
            out = a
        elif each is not None:
            pairs.append(f'{a}={each.replace("{v}", a)}')
        else:
            pairs.append(a)
    if out is None or not pairs:
        sys.exit(__doc__)
    runs = []
    if first is not None:
        files = sorted(glob.glob(first), key=lambda f: int(re.search(r'_k(\d+)_', f).group(1)))
        for n in map(int, pairs):
            if n > len(files):
                sys.exit(f'plot_convergence.py: {n} files asked, {len(files)} match {first}')
            runs.append((n, files[:n]))
    else:
        for p in pairs:
            v, pattern = p.split('=', 1)
            files = sorted(glob.glob(pattern))
            if not files:
                sys.exit(f'plot_convergence.py: no file matches {pattern}')
            runs.append((float(v), files))
    curves = [(v, *mean_of(f), len(f)) for v, f in runs]
    ref = curves[-1][2]
    scale = {'ms': 1.0, 'us': 1e3}[unit]
    diffs = []
    for i, (v, t, a, n) in enumerate(curves):
        d_ref = float(np.max(np.abs(a - ref)))
        d_prev = float(np.max(np.abs(a - curves[i - 1][2]))) if i else float('nan')
        diffs.append(d_ref)
        print(f'{name} {v:g}: {n} file(s), max |L - L(last)| = {d_ref:.3g}, from the previous = {d_prev:.3g}')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(curves)))
    for (v, t, a, n), c in zip(curves, colors):
        left.plot(t * scale, a, color=c, lw=1.0, label=f'{name} {v:g}')
    left.set_xlabel(f'time ({unit})')
    left.set_ylabel('|L(t)|')
    left.set_ylim(bottom=0)
    left.legend(fontsize=7, ncol=2)
    x = [c[0] for c in curves[:-1]]
    right.semilogy(x, [max(d, 1e-6) for d in diffs[:-1]], 'o-', color='k')
    right.axhline(tol, color='r', ls='--', lw=1, label=f'{tol:g}')
    right.set_xlabel(name)
    if all(float(v).is_integer() for v in x):
        right.set_xticks(x)
    right.set_ylabel(f'max |L - L({name} {curves[-1][0]:g})|')
    right.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print('wrote', out)


if __name__ == '__main__':
    main(sys.argv[1:])
