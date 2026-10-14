#!/usr/bin/env python3
"""The warp kernel's whole production curves judged against the CPU (A50b (iii), user decision 2026-10-13).

  python3 check_warp_spread.py REFINED WARP [--prefix out]                      (default: against the refined CPU)
  python3 check_warp_spread.py --spread ICPX GCC WARP [--prefix out] [--k 3]    (A50 (a): against two CPU builds' spread)

Default. REFINED is the output folder of main.out with gcce_refinement on and evolution_update exact (the eigenvalues
beyond double: two CPU builds then agree to 2.5e-7 on P1 gCCE2 CR60, against 2.5e-6 / 1.0e-5 without it, A50b job 17410),
WARP that of main_gpu.out with the warp kernel on the same input. Every value of every PREFIX_*Div file is judged in the
A40 verdict region - per curve (file and complex position), the values outside the diverged stretch (|REFINED| > 1 and
everything between the first and last such value) - with the GPU tolerance |WARP - REFINED| <= 1e-8 + 2e-6 |REFINED|.
Why the refined curve and not the production CPU: on these baths the double-precision eigenvalues of every code (CPU and
warp alike) are rounded at an ulp of |H| (~2.4e7), which the product amplifies to 1e-6 - 1e-5 relative; against the
refined curve the warp kernel is as close as the production CPU (A50b). The diverged stretch is not judged (A40).
--spread (A50 (a), kept as an option): ICPX, GCC = main.out of the usual build and of a gcc build (do_compile.sh generic,
CXX=mpicxx), both exact; every value |WARP - ICPX| <= max(k |GCC - ICPX|, rfloor |ICPX|, afloor).

Strict reading (review 20261013-180049 G1): the time column and every complex value of every line (a line without a time
continues the previous one's matrix, as benchmark.py reads them); FAIL on a value that is not finite, an empty file, a
line that does not parse, a missing file, different times, different row or value counts, no files, or no value judged.
One line per kind (wiDiv, noDiv) and a verdict; exit 1 on FAIL. Tests: test_check_warp_spread.py (CI python job).
"""
import argparse, math, os, re, sys

NONFINITE = re.compile(r'nan|inf', re.IGNORECASE)


class Bad(Exception):
    pass


def read(path):
    """[(time, [complex values])] of an output file; Bad on anything that is not a clean, finite, non-empty curve."""
    rows = []
    try:
        text = open(path).read()
    except OSError as e:
        raise Bad(f'cannot read {path}: {e}')
    for number, line in enumerate(text.splitlines(), 1):
        parts = line.split()
        if not parts: continue
        if NONFINITE.search(line): raise Bad(f'not finite: {path}:{number}')
        try:
            values = [complex(x) for x in parts]
        except ValueError:
            raise Bad(f'does not parse: {path}:{number}')
        if 'j' in parts[0]:
            if not rows: raise Bad(f'values before the first time: {path}:{number}')
            rows[-1][1].extend(values)
            continue
        if values[0].imag != 0 or len(values) < 2: raise Bad(f'no time and value: {path}:{number}')
        rows.append((values[0].real, values[1:]))
    if not rows: raise Bad(f'empty: {path}')
    width = len(rows[0][1])
    if any(len(v) != width for _, v in rows): raise Bad(f'rows of different length: {path}')
    return rows


def same_grid(a, b, name):
    if len(a) != len(b): raise Bad(f'{name}: {len(a)} against {len(b)} rows')
    for (ta, va), (tb, vb) in zip(a, b):
        if not math.isclose(ta, tb, rel_tol=1e-12, abs_tol=1e-15): raise Bad(f'{name}: time {ta} against {tb}')
        if len(va) != len(vb): raise Bad(f'{name}: {len(va)} against {len(vb)} values at time {ta}')


def judge(kind, ref_dir, warp_dir, prefix, allowed, gcc_dir=None, region=True):
    """(ok, line) for one kind. allowed(u, v_gcc) -> the bound for |warp - u|."""
    names = sorted(f for f in os.listdir(ref_dir) if f.startswith(prefix + '_') and f.endswith('_' + kind))
    judged = left = out = 0; worst = 0.0; at = ''
    try:
        if not names: raise Bad(f'no {prefix}_*_{kind} file in {ref_dir}')
        for name in names:
            others = [warp_dir] + ([gcc_dir] if gcc_dir else [])
            for d in others:
                if not os.path.exists(os.path.join(d, name)): raise Bad(f'missing {os.path.join(d, name)}')
            A = read(os.path.join(ref_dir, name)); W = read(os.path.join(warp_dir, name))
            same_grid(A, W, name)
            G = None
            if gcc_dir:
                G = read(os.path.join(gcc_dir, name)); same_grid(A, G, name)
            for k in range(len(A[0][1])):
                col = [v[k] for _, v in A]
                big = [i for i, u in enumerate(col) if abs(u) > 1] if region else []
                for i, u in enumerate(col):
                    if big and big[0] <= i <= big[-1]:
                        left += 1; continue
                    d = abs(W[i][1][k] - u)
                    bound = allowed(u, G[i][1][k] if G else None)
                    judged += 1
                    if d > bound: out += 1
                    if d / bound > worst: worst, at = d / bound, f'{name} t={A[i][0]:g} |CPU| {abs(u):.2e} diff {d:.2e}'
        if judged == 0: raise Bad('no value judged')
    except Bad as e:
        return False, f'FAIL {kind}: {e}'
    ok = out == 0
    return ok, (f"{'PASS' if ok else 'FAIL'} {kind}: files {len(names)}, values judged {judged}, left out (diverged stretch) {left}, "
                f'outside {out}, largest diff/allowed {worst:.3f}' + (f' ({at})' if at else ''))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('dirs', nargs='+', help='REFINED WARP, or with --spread ICPX GCC WARP')
    p.add_argument('--spread', action='store_true', help='A50 (a): against the spread of two CPU builds')
    p.add_argument('--prefix', default='out')
    p.add_argument('--k', type=float, default=3.0)
    p.add_argument('--rfloor', type=float, default=2e-6)
    p.add_argument('--afloor', type=float, default=1e-12)
    a = p.parse_args(argv)
    if len(a.dirs) != (3 if a.spread else 2): p.error('two folders (REFINED WARP), or three with --spread (ICPX GCC WARP)')
    ok = True
    for kind in ('wiDiv', 'noDiv'):
        if a.spread:
            icpx, gcc, warp = a.dirs
            good, line = judge(kind, icpx, warp, a.prefix, lambda u, v: max(a.k * abs(v - u), a.rfloor * abs(u), a.afloor),
                               gcc_dir=gcc, region=False)
        else:
            ref, warp = a.dirs
            good, line = judge(kind, ref, warp, a.prefix, lambda u, v: 1e-8 + 2e-6 * abs(u))
        ok &= good
        print(line)
    print(('warp spread verdict (k %g, floor max(%g |CPU|, %g)): ' % (a.k, a.rfloor, a.afloor) if a.spread else
           'warp verdict against the refined CPU (A40 region, 1e-8 + 2e-6 |CPU|): ') + ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
