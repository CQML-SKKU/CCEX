"""plot_coherence.py: the |L(t)| of CCEX's _wiDiv files in one figure, and one number per curve (the tutorials' expected.png).

    python3 plot_coherence.py OUT.png [--unit ms|us|ns] [--title TEXT] [--free-y] [--log-time] [--fit] [--bold-first] [--members] [--dashed] [--spread] LABEL=GLOB ...

Each LABEL=GLOB is one curve: the files GLOB matches (several bath configurations or sampled states) are averaged as
complex numbers, then |.| is plotted. A density-matrix file (quantity dm) gives the off-diagonal element that is largest
at t = 0, divided by its value there (or, when the qubit starts in one level, that level's population).

Printed per curve: the number of files; T2 as the last time |L| is at or above 1/e (linearly interpolated to the next
point; the fast ESEEM dips of a 13C bath do not cut it short), or |L| at the last time when the curve never falls below
1/e; from the second curve on, the largest |L| difference from the first. Points that are not finite are left out of
the average and counted; a curve above 1.2 (a diverged division) is reported and the plot cut at 1.2.
--log-time draws the time axis on a log scale (t = 0 left out), for curves whose decay times differ tenfold.
--spread prints, for each averaged curve, the T2 of every file it averages: their mean, standard deviation and
standard error (the deviation / sqrt of the number of files).
--members draws, under each averaged curve, every file it averages as a thin grey line (one legend entry).
--dashed draws the first curve thick and the others thin, dashed or dotted on top of it, so that curves that lie on
top of each other stay visible.
--bold-first draws the first curve thick and on top, the others thin and faint (a mean against its members).
--fit fits exp(-(t/T2)^n) to the upper envelope of each curve (the largest |L| in every window of 1/100 of the
time axis, so the ESEEM dips are passed over; every such point with 0 < |L| < 1), as a straight line of
log(-log|L|) against log t; it draws the fit (dashed, over the whole time axis) and the fitted points, and prints T2 and n.
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


def t1e_of(t, a):
    """The last time |L| >= 1/e (linear between points), or None when it never falls below."""
    above = np.nonzero(a >= 1 / np.e)[0]
    if not len(above) or above[-1] >= len(t) - 1:
        return None
    k = above[-1] + 1
    return t[k - 1] + (t[k] - t[k - 1]) * (a[k - 1] - 1 / np.e) / (a[k - 1] - a[k])


def fit_t2(t, a):
    """(T2, n, t of the fitted points, |L| there): exp(-(t/T2)^n) through the whole upper envelope."""
    w = max(1, len(t) // 100)
    idx = [i + int(np.argmax(a[i:i + w])) for i in range(0, len(t), w)]   # the largest |L| of each window
    idx = [i for i in idx if t[i] > 0 and 0 < a[i] < 1]                  # log(-log|L|) needs 0 < |L| < 1
    if len(idx) < 3:
        return None
    x, y = np.log(t[idx]), np.log(-np.log(a[idx]))
    n, c = np.polyfit(x, y, 1)                                           # y = n log t - n log T2
    return np.exp(-c / n), n, t[idx], a[idx]


def main(argv):
    out, unit, title, series, free_y, log_time, fit, bold_first, members, dashed, spread = (None, 'ms', None, [], False, False,
                                                                                         False, False, False, False, False)
    it = iter(argv)
    for a in it:
        if a == '--unit':
            unit = next(it)
        elif a == '--free-y':          # the y range from the data (a population near 1), not from 0
            free_y = True
        elif a == '--log-time':        # the time axis on a log scale; t = 0 cannot be drawn there
            log_time = True
        elif a == '--spread':          # the T2 of every averaged file: mean, deviation, standard error
            spread = True
        elif a == '--dashed':          # the second curve dashed, the third dotted, ...
            dashed = True
        elif a == '--members':         # every averaged file as a thin grey line
            members = True
        elif a == '--bold-first':      # the first curve thick, the others thin and faint
            bold_first = True
        elif a == '--fit':             # fit exp(-(t/T2)^n) to the upper envelope and draw it
            fit = True
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
        if spread and len(files) > 1:
            each = [t1e_of(t, np.abs(x)) for x in L]
            each = np.array([x for x in each if x is not None])
            if len(each) > 1:
                sd = each.std(ddof=1)
                number += (f'; T2 of the {len(each)} files: mean {each.mean() * scale:.4g} {unit}, deviation {sd * scale:.3g},'
                           f' standard error {sd / np.sqrt(len(each)) * scale:.3g}')
        f = fit_t2(t, mean) if fit else None
        if f:
            number += f'; fit T2 {f[0] * scale:.4g} {unit}, n {f[1]:.3g} ({len(f[2])} envelope points)'
        print(f'{label}: {len(files)} file(s), {number}{diff}{note}')
        keep = t > 0 if log_time else slice(None)
        style = dict(lw=0.7)
        if bold_first:
            style = dict(lw=1.6, color='k', zorder=3) if first[0] == label and first[1] is t else dict(lw=0.6, alpha=0.45)
        if members and len(files) > 1:
            for i, one in enumerate(np.abs(L)):
                ax.plot((t * scale)[keep], one[keep], color='0.7', lw=0.5, zorder=1,
                        label=f'each of the {len(files)}' if i == 0 else None)
            style = dict(lw=1.6, zorder=3)
        if dashed:
            n = [s[0] for s in series].index(label)
            style.update(ls=['-', (0, (5, 3)), (0, (2, 2)), '-.'][n % 4], lw=1.3 if n else 3.0, zorder=2 + n)
        line, = ax.plot((t * scale)[keep], mean[keep], label=f'{label} ({len(files)})' if len(files) > 1 else label, **style)
        if f:
            T2, n, tf, af = f
            tt = np.linspace(0, t[-1], 400)                           # the fit over the whole axis; the points show where it was fitted
            ax.plot(tt * scale, np.exp(-(tt / T2) ** n), '--', color='k', lw=1.2,
                    label=f'fit exp(-(t/T2)^n): T2 {T2 * scale:.3g} {unit}, n {n:.2f}')
            ax.plot(tf * scale, af, 'o', ms=3, mfc='none', color='k', lw=0, label='fitted points (envelope)')
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
