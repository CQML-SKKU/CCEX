"""compare.py A B --tol T [--scale-exp RATE_MHZ] [--conj-b] [--subset-b]:
max |A - B' exp(-pi 1e3 RATE t)| over the shared grid, B' = conj(B) with --conj-b.
--subset-b lets B carry exactly one extra trailing row (PyCCE linspace keeps the end point).
Exit 0 within T, 1 above, 2 when a file cannot be read or the grids differ (so a broken run is not a pass)."""
import sys
import numpy as np


def opt(name, cast, default):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


def read(p):
    try:
        rows = [l.split() for l in open(p) if l.strip() and not l.lstrip().startswith('#')]
        t = np.array([float(r[0]) for r in rows])
        if len(rows[0]) >= 3 and 'j' not in rows[0][1]:          # PyCCE: t Re Im
            return t, np.array([float(r[1]) + 1j * float(r[2]) for r in rows])
        return t, np.array([complex(r[1].replace('+-', '-')) for r in rows])   # CCEX: t complex
    except (OSError, ValueError, IndexError) as e:
        print(f'  cannot read {p}: {e}'); sys.exit(2)


a, b, tol, rate = sys.argv[1], sys.argv[2], opt('--tol', float, None), opt('--scale-exp', float, 0.0)
ta, La = read(a); tb, Lb = read(b)
if '--subset-b' in sys.argv and len(tb) == len(ta) + 1:
    tb, Lb = tb[:-1], Lb[:-1]
if '--conj-b' in sys.argv:
    Lb = Lb.conj()
if tol is None or len(ta) != len(tb) or np.max(np.abs(ta - tb)) > 1e-9 or not (np.all(np.isfinite(La)) and np.all(np.isfinite(Lb))):
    print(f'  grids differ or values not finite ({len(ta)} vs {len(tb)} rows)'); sys.exit(2)
dev = float(np.max(np.abs(La - Lb * np.exp(-np.pi * 1e3 * rate * ta))))
print(f'    max |diff| = {dev:.3e} (tol {tol:g})')
sys.exit(0 if dev <= tol else 1)
