"""Tests of the T2 sampling criterion, the truncation judgement and the report files (2026-10-07 decisions).

    cd example/tools/AutoConvergence && python -m unittest discover -s tests -v

Synthetic cases have known answers; the last class re-derives numbers recorded in
(internal notes) (14-t2-criterion-6pct.txt, 18-rdip-rejudge.txt) from the study data
when it is present and is skipped otherwise.
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from autoconv import config, t2, truncation, window  # noqa: E402
from autoconv.__main__ import analyze, check_pairs, judge_runs, main  # noqa: E402
from test_autoconv import write_run  # noqa: E402

T = np.linspace(0, 15e-3, 151)          # ms, as CCEX writes


def sample(rng, C, S, T2=3e-3, n=1.2, spread=0.1, noise=0.02, mod=None):
    """C configurations x S states: per configuration T2 * (1 + spread * N), per state additive complex noise."""
    out = []
    for _ in range(C):
        Tc = T2 * (1 + spread * rng.normal())
        base = t2.model(T, Tc, n) * (1 if mod is None else mod)
        out.append(np.array([base + noise * (rng.normal(size=len(T)) + 1j * rng.normal(size=len(T))) * base for _ in range(S)]))
    return out


class TestT2(unittest.TestCase):
    def test_fit_recovers_the_model(self):
        T2, n, k = t2.fit(T, t2.model(T, 2.9e-3, 1.2))
        self.assertAlmostEqual(T2, 2.9e-3, delta=1e-9)
        self.assertAlmostEqual(n, 1.2, delta=1e-6)
        self.assertGreaterEqual(len(k), 3)                               # the fitted points (upper envelope)

    def test_linearized_se_against_bootstrap_and_truth(self):
        rng = np.random.default_rng(1)
        data = sample(rng, 20, 10)
        comp = t2.components(T, data)
        lin, boot = t2.se_rel(comp, 20, 10), t2.bootstrap_se_rel(T, data, n_boot=300)
        self.assertLess(abs(lin / boot - 1), 0.3)
        # truth: T2 of the mean over repeated samples
        reps = [t2.fit(T, np.array(sample(rng, 20, 10)).mean((0, 1)))[0] for _ in range(60)]
        self.assertLess(abs(lin / (np.std(reps) / np.mean(reps)) - 1), 0.4)

    def test_envelope_fit_of_a_strongly_modulated_curve(self):     # NV-13C ESEEM (2026-10-07)
        t = np.linspace(0, 3, 301)
        env = t2.model(t, 1.2, 2.2)
        y = env * (0.55 + 0.45 * np.cos(2 * np.pi * t / 0.031))   # dips to 0.1 x every few points
        T, n, idx = t2.fit(t, y + 0j)
        # the sampled points rarely sit on the modulation's crests, so the upper points lie a little below the true
        # envelope and T2 comes out low: about 6 % here (NV-13C: 1.144 ms with the near 13C against 1.204 ms without)
        self.assertLess(abs(T / 1.2 - 1), 0.10)
        self.assertLess(abs(n / 2.2 - 1), 0.25)                       # n is pulled down too: 18 % here, 12 % on NV-13C
        with self.assertRaises(ValueError):
            t2.fit(t, y + 0j, points="all")

    def test_fit_outside_its_range_is_a_failure(self):
        with self.assertRaisesRegex(ValueError, "range"):
            t2.fit(T, t2.model(T, 3e-3, 6.0) + 0j)                        # n = 6, above N_MAX
        with self.assertRaises(ValueError):
            t2.fit(T, np.ones(len(T)) * 0.01 + 0j)                        # no point above the floor

    def test_failed_bootstrap_resamples_are_left_out_up_to_five_percent(self):   # user 2026-10-10 (P1 CR0 cCCE2 single)
        t = np.linspace(0, 5, 51)
        good = np.exp(-(t / 1.5) ** 2)[None, None, :]
        bad = np.where(t == 0, 1.0, 0.0)[None, None, :]                     # below the floor at once: its fit fails
        def draws(every):
            calls = iter(range(10 ** 6))
            return lambda r: bad if next(calls) % every == every - 1 else good * (1 + 0.01 * r.normal())
        info = {}
        areas = t2._boot_areas(t, None, draws(50), 400, np.random.default_rng(1), 0.05, "envelope", (0.3, 4.0), info)
        self.assertEqual((info["bootstrap_failed"], len(areas)), (8, 392))      # 2 %: left out
        with self.assertRaisesRegex(ValueError, "40 of 400 bootstrap resamples could not be fitted"):
            t2._boot_areas(t, None, draws(10), 400, np.random.default_rng(1), 0.05, "envelope", (0.3, 4.0), {})   # 10 %

    def test_allocation_is_minimal_and_bounded(self):
        comp = dict(T2=1.0, n=1.0, A=1.0, B=0.01, W=0.04, C=10, S=10)
        a = t2.allocate(comp, 0.05, 1.0, 1.0)
        for r in a["table"]:
            self.assertLessEqual(r["se_rel"], 0.05 + 1e-12)
            self.assertGreater(t2.se_rel(comp, r["C"] - 1, r["S"]), 0.05 - 1e-12) if r["C"] > 2 else None
        tight = t2.allocate(comp, 0.05, 1.0, 1.0, max_configurations=3)
        self.assertTrue(tight["status"].startswith("infeasible"))
        self.assertEqual(tight["best"]["C"], 3)
        more = t2.allocate(comp, 0.05, 1.0, 1.0, correction=1.2)
        self.assertGreater(more["best"]["C"] * more["best"]["S"], a["best"]["C"] * a["best"]["S"] - 1)


class TestTruncation(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(2)

    def pair(self, factor=1.0, mod_ref=None, mod_set=None, C=12):
        ref, st = [], []
        for _ in range(C):
            Tc = 3e-3 * (1 + 0.1 * self.rng.normal())
            r = t2.model(T, Tc, 1.2) * (1 if mod_ref is None else mod_ref)
            s = t2.model(T, Tc * factor, 1.2) * (1 if mod_set is None else mod_set)
            ref.append(r); st.append(s)
        return np.array(st, complex), np.array(ref, complex)

    def test_identical_passes(self):
        s, r = self.pair()
        j = truncation.judge(T, s, r, 12.6e-3)
        self.assertEqual((j["route"], j["verdict"]), ("fit", "pass"))
        self.assertEqual(j["fit_max"], 0.0)

    def test_small_and_large_t2_shifts(self):
        s, r = self.pair(1.005)
        self.assertEqual(truncation.judge(T, s, r, 12.6e-3)["verdict"], "pass")
        s, r = self.pair(1.10)
        j = truncation.judge(T, s, r, 12.6e-3)
        self.assertEqual(j["verdict"], "fail")
        self.assertGreater(j["fit_max"], 0.01)

    def test_common_modulation_keeps_the_fit_route(self):
        mod = 1 + 0.08 * np.cos(2 * np.pi * T / 1.2e-3)
        s, r = self.pair(1.005, mod, mod)
        j = truncation.judge(T, s, r, 12.6e-3)
        self.assertEqual((j["route"], j["verdict"]), ("fit", "pass"))

    def test_a_failing_fit_goes_to_the_smoothed_route(self):
        s, r = self.pair()
        flat = np.full_like(s, 0.01)                     # nothing above the fit floor: no fit is possible
        j = truncation.judge(T, flat, flat.copy(), 12.6e-3, fit_only=False)
        self.assertEqual((j["route"], j["verdict"], j["fit_max"]), ("smoothed", "pass", None))
        # fit lines only (default): a failed fit is judged on the raw area since 2026-10-12 (equal flat curves pass)
        self.assertEqual(truncation.judge(T, flat, flat.copy(), 12.6e-3)["route"], "raw-area")
        self.assertIn("fit_error", j)

    def test_raw_area_runs_to_t_d_exactly(self):          # review 20261012-233058-452278 P4-02 (and P4-01: no np.trapezoid)
        t = np.array([0.0, 1.0, 2.0])
        y = np.array([1.0, 1.0, 0.0])
        self.assertAlmostEqual(truncation.raw_area(t, y, 1.5), 1.0 + 0.5 * (1.0 + 0.5) / 2, places=15)  # partial interval in
        self.assertAlmostEqual(truncation.raw_area(t, y, 2.0), 1.5, places=15)
        self.assertTrue(np.isnan(truncation.raw_area(t, y, 2.5)))                              # beyond the data: not judged
        # the review's counterexample: the last partial interval carries the difference (10 %, a fail at 2 %)
        tt = np.linspace(0, 1, 6)
        T2 = .9 / np.log(100) ** (1 / 1.2)
        bb = np.exp(-(tt / T2) ** 1.2)
        aa = bb.copy(); aa[1] = .01; aa[2] += bb[1] - .01; aa[-1] = 1
        self.assertGreater(abs(truncation.raw_area(tt, aa, .9) / truncation.raw_area(tt, bb, .9) - 1), 0.02)

    def test_a_failed_fit_on_one_side_is_judged_on_the_raw_area(self):     # user 2026-10-12, plan 2.1
        _, r = self.pair()
        A_ref = 3e-3 * truncation.gamma(1 + 1 / 1.2)
        same = np.array([t2.model(T, A_ref / truncation.gamma(1 + 1 / 6), 6.0)] * len(r), complex)   # n 6: outside 0.3-4
        j = truncation.judge(T, same, r, 12.6e-3)
        m = T <= 12.6e-3
        raw = truncation.raw_area(T, np.abs(same.mean(0)), 12.6e-3) / truncation.raw_area(T, np.abs(r.mean(0)), 12.6e-3) - 1
        self.assertEqual((j["route"], j["verdict"], j["area_rel"]), ("raw-area", "pass", None))
        self.assertAlmostEqual(j["raw_area_rel"], raw, places=12)
        self.assertLessEqual(abs(raw), 0.02)
        j = truncation.judge(T, same * 1.05, r, 12.6e-3)                  # 5 % more raw area: fails with the same 2 %
        self.assertEqual((j["route"], j["verdict"]), ("raw-area", "fail"))
        j = truncation.judge(T, r, same * 1.05, 12.6e-3)                  # the reference side failing: the same
        self.assertEqual((j["route"], j["verdict"]), ("raw-area", "fail"))
        self.assertEqual(truncation.judge(T, same * 1.05, r, 12.6e-3, rel_tol=0.1)["verdict"], "pass")

    def test_fitted_comparisons_keep_the_fit_verdict(self):
        for f, verdict in ((1.005, "pass"), (1.10, "fail")):
            s, r = self.pair(f)
            j = truncation.judge(T, s, r, 12.6e-3)
            self.assertEqual((j["route"], j["verdict"]), ("fit", verdict))
            self.assertAlmostEqual(abs(j["area_rel"]), abs(truncation.fit_area(j["T2"], j["n"]) / truncation.fit_area(j["T2_ref"], j["n_ref"]) - 1))
            self.assertNotIn("raw_area_rel", j)

    def test_a_non_finite_curve_stays_unjudged(self):    # (the runner leaves such configurations out before judging)
        s, r = self.pair()
        s[:, 20:] = np.nan                                # no fit, and no finite raw area
        with np.errstate(invalid="ignore"):
            j = truncation.judge(T, s, r, 12.6e-3)
        self.assertEqual((j["route"], j["verdict"]), ("fit", "fit-failure"))
        self.assertNotIn("raw_area_rel", j)

    def test_modulation_that_differs_goes_to_the_smoothed_route(self):
        s, r = self.pair(1.0, 1 + 0.08 * np.cos(2 * np.pi * T / 4e-3), 1 + 0.08 * np.cos(2 * np.pi * T / 4e-3 + 1.5))
        j = truncation.judge(T, s, r, 12.6e-3, fit_only=False)
        self.assertEqual(j["route"], "smoothed")
        self.assertIn(j["verdict"], ("pass", "fail"))


class TestCommandLine(unittest.TestCase):
    def test_analyze_and_judge_write_reports(self):
        rng = np.random.default_rng(3)
        with tempfile.TemporaryDirectory() as tmp:
            runs = []
            for c, x in enumerate(sample(rng, 8, 4)):
                d = Path(tmp) / f"cfg{c}"
                write_run(d, T, x)
                write_run(d, T, x * 1.0, prefix="out_rdip200")       # a ladder rung equal to the main output
                (d / "input.json").write_text(json.dumps(dict(bathfile=[f"bath_{c}"], seed=1, rdip=400, rdip_ladder=[200])))
                runs.append(str(d))
            rep = analyze(runs, config.parse_auto_convergence({}), cost_config=7.0, cost_state=2.6,
                          max_configurations=50, report_dir=f"{tmp}/rep")
            self.assertIn(rep["sampling"]["status"], ("t2-met", "sampling-unresolved"))
            self.assertIn("T2", rep["sampling"]["t2"])
            self.assertTrue((Path(tmp) / "rep/sampling.png").stat().st_size > 0)
            self.assertIn("SE(A)/A", (Path(tmp) / "rep/sampling.txt").read_text())
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = main(["judge", "--setting", *runs, "--reference", *runs, "--setting-prefix", "out_rdip200",
                             "--name", "rdip200", "--report-dir", f"{tmp}/rep"])
            self.assertEqual(code, 0)
            j = json.loads(out.getvalue())
            self.assertEqual((j["verdict"], j["fit_max"]), ("pass", 0.0))
            self.assertTrue((Path(tmp) / "rep/truncation_rdip200.png").stat().st_size > 0)
            with self.assertRaises(ValueError):
                judge_runs(runs, runs[:-1], config.parse_auto_convergence({}))

    def test_pairs_must_be_the_same_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            def run(name, **inp):
                d = Path(tmp) / name; d.mkdir()
                if inp: (d / "input.json").write_text(json.dumps(inp))
                return str(d)
            a = run("a", configuration=dict(bathfile=["b1"], bath_state_file=["s1"], seed=1), rdip=200, output=dict(file="x"))
            b = run("b", configuration=dict(bathfile=["b1"], bath_state_file=["s1"], seed=1), rdip=400, output=dict(file="y"))      # only setting keys differ
            self.assertEqual(check_pairs([a], [b]), [])
            c = run("c", configuration=dict(bathfile=["b2"], bath_state_file=["s1"], seed=1), rdip=400)                  # another configuration
            with self.assertRaisesRegex(ValueError, "configuration.bathfile"):
                check_pairs([a], [c])
            e = run("e", configuration=dict(bathfile=["b1"], bath_state_file=["s1"], seed=2), rdip=400)                  # another seed
            with self.assertRaisesRegex(ValueError, "configuration.seed"):
                check_pairs([a], [e])
            n = run("n")                                                                    # no input.json
            with self.assertRaisesRegex(ValueError, "unverified"):
                check_pairs([a], [n])
            self.assertEqual(len(check_pairs([a], [n], unverified=True)), 1)


PILOT = Path.home() / "tmp/ccex-autoconv-pilot/attempt-16985/runs"
PHASEB = Path.home() / "tmp/ccex-auto-trunc/0c-res/B"


@unittest.skipUnless(PILOT.is_dir() and PHASEB.is_dir(), "study data not present")
class TestRecordedNumbers(unittest.TestCase):
    def test_t2_sampling_matches_14(self):           # 14-t2-criterion-6pct.txt
        from autoconv.curves import read_runs
        runs = sorted(str(p) for p in PILOT.iterdir() if p.name.endswith("a"))
        times, data, _ = read_runs(runs)
        comp = t2.components(times, data, points="all", quantity="T2")   # 14: every point fitted, the SE of T2
        self.assertAlmostEqual(comp["T2"] * 1e3, 2.3231, delta=1e-3)
        self.assertAlmostEqual(comp["B"] * 1e6, 0.6757, delta=2e-3)
        self.assertAlmostEqual(comp["W"] * 1e6, 1.199, delta=3e-3)
        self.assertAlmostEqual(t2.se_rel(comp, 20, 10), 0.0859, delta=5e-4)
        a = t2.allocate(comp, 0.06, 7.03, 2.62, s_candidates=(2, 3, 5, 10, 20, 30), max_configurations=50)
        self.assertEqual((a["best"]["C"], a["best"]["S"]), (48, 5))

    def test_rdip_judgement_matches_18(self):        # 18-rdip-rejudge.txt: rdip 251.359 vs 490.935, 57 configurations
        from autoconv.curves import read_runs
        cfgs = sorted(p for p in PHASEB.iterdir() if p.is_dir() and p.name != "cfg64")
        self.assertEqual(check_pairs([str(c / "r251.359") for c in cfgs], [str(c / "r490.935") for c in cfgs]), [])
        with self.assertRaisesRegex(ValueError, "bathfile"):         # shifted by one: another configuration
            check_pairs([str(c / "r251.359") for c in cfgs[1:]], [str(c / "r490.935") for c in cfgs[:-1]])
        times, ds, _ = read_runs([str(c / "r251.359") for c in cfgs])
        _, dr, _ = read_runs([str(c / "r490.935") for c in cfgs])
        j = truncation.judge(times, np.array([x.mean(0) for x in ds]), np.array([x.mean(0) for x in dr]), 12.6e-3, points="all",
                             fit_only=False)
        self.assertEqual((j["route"], j["verdict"]), ("fit", "pass"))
        self.assertAlmostEqual(j["fit_max"], 1.35e-3, delta=1e-5)
        self.assertAlmostEqual(j["fit_area"], 6.82e-4, delta=1e-6)
        self.assertAlmostEqual(j["dT2_rel"], 0.00274, delta=1e-5)


class TestWindowSpikes(unittest.TestCase):          # rule 2 (2026-10-07): a running median ignores one or two points
    def setUp(self):
        self.t = np.linspace(0, 10, 501)
        self.base = [np.array([np.exp(-self.t / 1.0) + 0j for _ in range(3)]) for _ in range(4)]

    def test_single_point_spike_is_ignored_and_reported(self):
        data = [x.copy() for x in self.base]
        data[0][1, 400] = 50.0                      # one state, one point, long after the decay
        w = window.choose_window(self.t, data, 0.01, 1.2)
        self.assertEqual(w["status"], "resolved")
        self.assertIn(self.t[400], w["ignored_spikes"])
        self.assertLess(w["decay_time"], self.t[400])
        off = window.choose_window(self.t, data, 0.01, 1.2, median_points=1)
        self.assertGreater(off["decay_time"], self.t[400])   # without the median the spike holds the envelope up

    def test_a_revival_over_several_points_is_kept(self):
        data = [x.copy() for x in self.base]
        for x in data:
            x[:, 380:395] += 0.05                   # 15 points back above the threshold
        w = window.choose_window(self.t, data, 0.01, 1.2)
        self.assertGreater(w["decay_time"], self.t[390])


class TestFitWindow(unittest.TestCase):              # window from the T2 fit line (2026-10-07)
    def test_td_is_where_the_fit_reaches_the_threshold(self):
        t = np.linspace(0, 10, 501)
        data = [np.array([t2.model(t, 2.0, 1.5) + 0j for _ in range(3)]) for _ in range(3)]
        self.assertEqual(window.fit_window(t, data, 0.01, 1.2)["status"], "needs-extension")   # < 2.5 T_d: revivals unseen
        t = np.linspace(0, 20, 1001)
        data = [np.array([t2.model(t, 2.0, 1.5) + 0j for _ in range(3)]) for _ in range(3)]
        w = window.fit_window(t, data, 0.01, 1.2)
        self.assertEqual(w["status"], "resolved")
        self.assertAlmostEqual(w["decay_time"], 2.0 * np.log(100) ** (1 / 1.5), delta=1e-3)
        self.assertEqual(w["late_rise"], [])

    def test_short_window_and_late_rise(self):
        t = np.linspace(0, 1, 51)
        data = [np.array([t2.model(t, 2.0, 1.5) + 0j for _ in range(3)]) for _ in range(3)]
        self.assertEqual(window.fit_window(t, data, 0.01, 1.2)["status"], "needs-extension")   # not decayed: nothing to fit
        t = np.linspace(0, 10, 501)
        data = [np.array([t2.model(t, 1.0, 2.0) + 0.03 * (t > 6) + 0j for _ in range(3)]) for _ in range(3)]
        w = window.fit_window(t, data, 0.01, 1.2)
        self.assertGreater(len(w["late_rise"]), 0)                                              # a tail the fit has not

    def test_a_revival_lengthens_the_window(self):
        t = np.linspace(0, 10, 1001)
        bump = 0.05 * np.exp(-((t - 5.0) / 0.3) ** 2)                 # a revival around 5 (T_d about 2.1)
        data = [np.array([t2.model(t, 1.0, 2.0) + bump + 0j for _ in range(3)]) for _ in range(3)]
        w = window.fit_window(t, data, 0.01, 1.2)
        self.assertTrue(w["revival"])
        self.assertAlmostEqual(w["recommended_window"], 1.2 * w["revival_end"])
        self.assertGreater(w["revival_end"], 5.0)



class TestEnvelopeStep(unittest.TestCase):
    """steps() (the analyze driver's envelope) and the peak fit points of the runner (user 2026-10-10)."""

    def test_steps_is_the_span_over_the_grid(self):
        t = np.linspace(0, 1, 101)                                      # dt 0.01
        self.assertEqual(t2.steps(t, 0.04), 4)
        self.assertEqual(t2.steps(np.linspace(0, 1, 201), 0.04), 8)     # half the step: twice the points ahead
        self.assertEqual(t2.steps(t, 0.035), 3)                         # floor(h / dt)
        with self.assertRaises(ValueError):
            t2.steps(t, 0.005)                                          # the grid coarser than the span

    def test_steps_reads_rounded_output_times(self):           # review 20261010-160920 R1: CCEX writes %15g
        dt = 0.0123456789
        grid = lambda n: np.array([float(f"{i * dt:g}") for i in range(n)])
        span = (grid(401)[-1] - grid(401)[0]) / 400
        self.assertEqual(t2.steps(grid(201), span), 1)                 # the same step, fewer points
        self.assertEqual(t2.steps(grid(201), 8 * span), 8)
        half = np.array([float(f"{i * dt / 2:g}") for i in range(801)])
        self.assertEqual(t2.steps(half, span), 2)                      # half the step: twice the points ahead
        with self.assertRaises(ValueError):
            t2.steps(grid(201), 0.999 * span)                          # a grid really coarser than the span

    def test_smooth_curve_takes_all_points(self):                  # user 2026-10-10: fewer than 3 peaks, all points
        t = np.linspace(0, 6, 401)
        y = t2.model(t, 1.5, 2.0)
        self.assertEqual(list(t2.peaks(y)), [0])
        self.assertEqual(list(t2.peak_idx(y)), list(range(t2.fit_points(y))))

    def test_modulated_curve_fits_the_peaks(self):                 # user 2026-10-10: t = 0 and the peaks, to the floor
        t = np.linspace(0, 6, 1501)
        y = t2.model(t, 1.5, 2.0) * (1 - 0.3 * (1 - np.cos(2 * np.pi * t / 0.1)))
        p = t2.peak_idx(y)
        self.assertEqual(p[0], 0)
        self.assertTrue(all(y[i] > y[i - 1] and y[i] >= y[i + 1] for i in p[1:]))
        self.assertTrue((y[p] >= 0.05).all() and y[p[-1] + 1:][y[p[-1] + 1:] >= 0.05].size < 0.1 * len(t))   # cut at the floor
        T2, n, idx = t2.fit(t, y + 0j, points="peaks")
        self.assertEqual(list(idx), list(p))
        self.assertAlmostEqual(T2, 1.5, delta=0.03); self.assertAlmostEqual(n, 2.0, delta=0.1)

    def test_window_reports_the_fit_points(self):
        t = np.linspace(0, 6, 401)
        data = [np.array([t2.model(t, 1.5, 2.0) + 0j for _ in range(2)]) for _ in range(3)]
        w = window.fit_window(t, data, 0.01, 1.2, look=0, points="peaks")
        self.assertEqual((w["status"], w["fit_points"]), ("resolved", "all"))           # no peaks: all points
        y = t2.model(t, 1.5, 2.0) * (1 - 0.3 * (1 - np.cos(2 * np.pi * t / 0.1)))
        w = window.fit_window(t, [np.array([y + 0j])] * 3, 0.01, 1.2, look=0, points="peaks")
        self.assertEqual((w["status"], w["fit_points"]), ("resolved", "peaks"))
        self.assertEqual(w["fit_npoints"], len(t2.peak_idx(y, w["fit_floor"])))

    def test_smoothing_without_the_square_matrix(self):
        t = np.linspace(0, 1, 2500)
        x = np.sin(9 * t) + 0.1 * np.cos(70 * t)
        k = np.exp(-0.5 * ((t[:, None] - t[None, :]) / 0.03) ** 2)
        # the same sum in another order: equal to rounding, not bit for bit (numpy 2 and other BLAS builds order it
        # differently)
        np.testing.assert_allclose(truncation._smooth(t, 0.03, x), (k / k.sum(1, keepdims=True)) @ x, rtol=1e-12, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
