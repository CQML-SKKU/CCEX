"""plot_coherence.py: the |L(t)| of CCEX's _wiDiv files in one figure, and one number per curve (the tutorials' expected.png).

    python3 plot_coherence.py OUT.png [--unit ms|us|ns] [--title TEXT] [--free-y] [--log-time] LABEL=GLOB [LABEL=GLOB ...]

Each LABEL=GLOB is one curve: the files GLOB matches (several bath configurations or sampled states) are averaged as
complex numbers, then |.| is plotted. A density-matrix file (quantity dm) gives the off-diagonal element that is largest
at t = 0, divided by its value there (or, when the qubit starts in one level, that level's population).

Printed per curve: the number of files; T2 as the last time |L| is at or above 1/e (linearly interpolated to the next
point; the fast ESEEM dips of a 13C bath do not cut it short), or |L| at the last time when the curve never falls below
1/e; from the second curve on, the largest |L| difference from the first. Points that are not finite are left out of
the average and counted; a curve above 1.2 (a diverged division) is reported and the plot cut at 1.2.
--log-time draws the time axis on a log scale (t = 0 left out), for curves whose decay times differ tenfold.
Needs numpy and matplotlib.
"""
import glob
import sys

import numpy as np


def read(path):
    """(t, values): t in ms; values complex, one per time (coherence) or one matrix per time (dm)."""
    times, blocks = [], []
    for line in open(path):
        if not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        row = [complex(x.strip().replace(' ', '')) for x in parts[1:] if x.strip()]
        if parts[0].strip():
            times.append(float(parts[0]))
            blocks.append([row])
        else:
            blocks[-1].append(row)
    t = np.array(times)
    if all(len(b) == 1 and len(b[0]) == 1 for b in blocks):
        return t, np.array([b[0][0] for b in blocks])
    flat = [[x for row in b for x in row] for b in blocks]   # a D x D matrix per time, written a few values per line
    D = int(round(len(flat[0]) ** 0.5))
    rho = np.array(flat).reshape(len(flat), D, D)          # (time, D, D), row-major
    off = np.abs(rho[0]) * (1 - np.eye(rho.shape[1]))
    i, j = np.unravel_index(np.argmax(off), off.shape)
    if off[i, j] < 1e-12:                                    # no coherence at t = 0: the population of the occupied level
        i = j = int(np.argmax(np.abs(np.diag(rho[0]))))
    return t, rho[:, i, j] / rho[0, i, j]


def main(argv):
    out, unit, title, series, free_y, log_time = None, 'ms', None, [], False, False
    it = iter(argv)
    for a in it:
        if a == '--unit':
            unit = next(it)
        elif a == '--free-y':          # the y range from the data (a population near 1), not from 0
            free_y = True
        elif a == '--log-time':        # the time axis on a log scale; t = 0 cannot be drawn there
            log_time = True
        elif a == '--title':
            title = next(it)
        elif out is None:
            out = a
        else:
            label, pattern = a.rsplit('=', 1)
            series.append((label, pattern))
    if out is None or not series:
        sys.exit(__doc__)
    scale = {'ms': 1.0, 'us': 1e3, 'ns': 1e6}[unit]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    first, clip = None, False
    for label, pattern in series:
        files = sorted(glob.glob(pattern))
        if not files:
            sys.exit(f'plot_coherence.py: no file matches {pattern}')
        data = [read(f) for f in files]
        t = data[0][0]
        L = np.array([d[1] for d in data])
        bad = ~np.isfinite(L)
        L = np.where(bad, np.nan, L)
        mean = np.abs(np.nanmean(L, axis=0))
        above = np.nonzero(mean >= 1 / np.e)[0]
        if len(above) and above[-1] < len(t) - 1:
            k = above[-1] + 1
            t1e = t[k - 1] + (t[k] - t[k - 1]) * (mean[k - 1] - 1 / np.e) / (mean[k - 1] - mean[k])
            number = f'T2 {t1e * scale:.4g} {unit}'
        else:
            number = f'|L| = {mean[-1]:.4f} at {t[-1] * scale:.4g} {unit} (no 1/e crossing)'
        note = f', {bad.sum()} non-finite points left out' if bad.any() else ''
        if np.nanmax(mean) > 1.2:
            note += f', |L| up to {np.nanmax(mean):.3g} (> 1: the divided curve diverged; the plot is cut at 1.2)'
            clip = True
        diff = ''
        if first is None:
            first = (label, t, mean)
        elif len(t) == len(first[1]):
            diff = f', max |L - L({first[0]})| = {np.nanmax(np.abs(mean - first[2])):.3e}'
        print(f'{label}: {len(files)} file(s), {number}{diff}{note}')
        keep = t > 0 if log_time else slice(None)
        ax.plot((t * scale)[keep], mean[keep], lw=0.7, label=f'{label} ({len(files)})' if len(files) > 1 else label)
    if log_time:
        ax.set_xscale('log')
    ax.set_xlabel(f'time ({unit})')
    ax.set_ylabel('|L(t)|')
    if clip:
        ax.set_ylim(0, 1.2)
    elif not free_y:
        ax.set_ylim(bottom=0)
    if title:
        ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print('wrote', out)


if __name__ == '__main__':
    main(sys.argv[1:])
