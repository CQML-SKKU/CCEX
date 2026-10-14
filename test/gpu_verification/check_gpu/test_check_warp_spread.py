#!/usr/bin/env python3
"""check_warp_spread.py on small made-up output folders (review 20261013-180049 G1): the cases that must pass and the
broken ones that must fail - an empty file, nan only, different times, different row or value counts, a missing file, a
line that does not parse, nothing judged. No binary, no GPU: python3 test_check_warp_spread.py (CI python job)."""
import contextlib, io, os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_warp_spread as cws

TIMES = [0.0, 0.0006, 0.0012, 0.0018]


def curve(values, times=TIMES):
    return ''.join(f'{t:15g}\t  {v.real:+.10f}{v.imag:+.10f}j\t\n' for t, v in zip(times, values))


class Folders:
    """REF, WARP (and GCC) folders with out_state1_{wiDiv,noDiv}; overrides: {(folder, kind): text}"""
    def __init__(self, base, ref, warp, gcc=None, overrides=None):
        self.tmp = tempfile.TemporaryDirectory()
        self.dirs = {}
        for name, vals in (('ref', ref), ('warp', warp), ('gcc', gcc)):
            if vals is None: continue
            d = os.path.join(self.tmp.name, name); os.makedirs(d); self.dirs[name] = d
            for kind in ('wiDiv', 'noDiv'):
                text = (overrides or {}).get((name, kind), curve(vals))
                if text is not None:
                    open(os.path.join(d, f'out_state1_{kind}'), 'w').write(text)

    def run(self, *extra):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            order = ('ref', 'gcc', 'warp') if '--spread' in extra else ('ref', 'warp')
            code = cws.main([*extra, *(self.dirs[n] for n in order)])
        return code, out.getvalue()


REF = [1, 0.9 - 0.01j, 0.8 - 0.04j, 0.7 - 0.05j]
NEAR = [v * (1 + 1e-9) for v in REF]


class Refined(unittest.TestCase):
    def check(self, code, out, ok):
        self.assertEqual(code, 0 if ok else 1, out)

    def test_pass(self):
        self.check(*Folders(None, REF, NEAR).run(), True)

    def test_outside(self):
        self.check(*Folders(None, REF, [v * (1 + 1e-4) for v in REF]).run(), False)

    def test_diverged_stretch_left_out(self):
        ref = [1, 5, 50, 0.5]           # |ref| > 1 at 1 and 2: those left out; 0 and 3 judged
        warp = [1, 6, 40, 0.5]
        code, out = Folders(None, ref, warp).run()
        self.check(code, out, True); self.assertIn('left out (diverged stretch) 2', out)

    def test_nothing_judged(self):
        ref = [2, 5, 50, 3]
        self.check(*Folders(None, ref, ref).run(), False)

    def test_empty_file(self):
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): ''}).run(), False)

    def test_nan_only(self):
        nan = ''.join(f'{t:15g}\t  +nan+nanj\t\n' for t in TIMES)
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): nan}).run(), False)
        self.check(*Folders(None, REF, REF, overrides={('ref', 'noDiv'): nan}).run(), False)

    def test_other_times(self):
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): curve(REF, [0, 0.001, 0.002, 0.003])}).run(), False)

    def test_row_count(self):
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): curve(REF[:3])}).run(), False)

    def test_value_count(self):
        two = ''.join(f'{t:15g}\t  {v.real:+.10f}{v.imag:+.10f}j\t  +0.1+0.0j\n' for t, v in zip(TIMES, REF))
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): two}).run(), False)

    def test_missing_file(self):
        self.check(*Folders(None, REF, REF, overrides={('warp', 'noDiv'): None}).run(), False)

    def test_does_not_parse(self):
        self.check(*Folders(None, REF, REF, overrides={('warp', 'wiDiv'): curve(REF) + '0.0024 garbage\n'}).run(), False)

    def test_no_files(self):
        f = Folders(None, REF, REF, overrides={(n, k): None for n in ('ref', 'warp') for k in ('wiDiv', 'noDiv')})
        self.check(*f.run(), False)


class Spread(unittest.TestCase):
    def test_pass_and_fail(self):
        gcc = [v * (1 + 1e-6) for v in REF]
        code, out = Folders(None, REF, [v * (1 + 2e-6) for v in REF], gcc).run('--spread')
        self.assertEqual(code, 0, out)
        code, out = Folders(None, REF, [v * (1 + 1e-4) for v in REF], gcc).run('--spread')
        self.assertEqual(code, 1, out)

    def test_broken_gcc(self):
        code, out = Folders(None, REF, REF, REF, overrides={('gcc', 'wiDiv'): ''}).run('--spread')
        self.assertEqual(code, 1, out)


if __name__ == '__main__':
    unittest.main()
