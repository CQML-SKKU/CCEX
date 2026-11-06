"""Automatic starts of the automatic convergence (autoconv/autostart.py and the runner without candi_rdip, candi_rbath, pilot, pool and
deltat; plan "자동 시작값", 2026-10-08), against fake_ccex.py.

    cd example/tools/AutoConvergence && python -m unittest tests.test_autostart -v
"""
import json
import math
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

from autoconv import autostart                       # noqa: E402
from test_runner import PY, call, nest, summary      # noqa: E402

# spins in a 100 A box around a qubit at its centre: rho 4.5e-3 /A^3 puts 300 spins inside about 25 A; the gyro 3.44 rad/ms/G
# puts the coupling clock near the fake curve's decay (T_d about 3.2 ms)
NSPIN, GYRO = 4500, 3.44


def make_auto(tmp, ids=tuple(range(1, 31)), nspin=NSPIN, flat=False, species=("13C",), gyro=GYRO, base=None, **block):
    tmp = Path(tmp)
    data = tmp / "data"; data.mkdir()
    rng = np.random.default_rng(7)
    for k in ids:
        xyz = rng.uniform(0, 100, (nspin, 3))
        if flat:
            xyz[:, 2] = 50.0
        names = [species[i % len(species)] for i in range(nspin)]
        (data / f"bath_{k}").write_text(f"{nspin + 1:.5f} 0 0 0\n" + "".join(f"{x:.5f} {y:.5f} {z:.5f} {n}\n" for (x, y, z), n in zip(xyz, names)))
    (data / "notes_bath.txt").write_text("not a bath\n")          # not matched by the pattern
    (data / "qubit").write_text("50 50 50\n")
    (data / "gyro").write_text(f"13C 0.5 {gyro}\n14N 1.0 1.93\n")
    fake = tmp / "fake_main.out"
    fake.write_text(f"#!/bin/sh\nexec {PY} {HERE / 'fake_ccex.py'} \"$@\"\n"); fake.chmod(0o755)
    blk = dict(workdir="auto", configurations=dict(bathfile=str(data / "bath_{k}")), threshold_config=0.06, min_configurations=8)
    blk.update(block)
    conf = dict(gyrofile=str(data / "gyro"))       # the configuration block (A64; the bath pattern in it A65); the ensemble
    inp = dict(method="cCCE", quantity="coherence", order=2, Qubit=dict(qubit=[dict(file=str(data / "qubit"))]),
               configuration=conf, bfield=500, auto_convergence=nest(blk, conf), **(base or {}))
    (tmp / "input.json").write_text(json.dumps(inp, indent=1))
    return tmp / "input.json", fake


class TestAutostartParts(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_a_tool_on_its_own_finds_the_table_by_the_binary(self):     # review 6548a9 R3
        from autoconv import autostart
        src = Path(autostart.__file__).read_text()
        g = dict(__file__="/opt/AutoConvergence/autoconv/autostart.py", __name__="autostart_shallow")
        exec(compile(src, "autostart.py", "exec"), g)                    # a tool fewer than four folders deep imports
        self.assertFalse(g["ISOTOPES"].is_file())
        repo = Path(autostart.__file__).resolve().parents[4]
        self.assertEqual(g["find_isotope_table"](repo / "bin" / "main.out"), repo / "src" / "input" / "isotope_table.cpp")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_discover_reads_the_pattern_literally_and_sorts_numerically(self):
        for name in ("b_2", "b_10", "b_1", "b_01", "b_x", "b_3.bak", "c_4"):
            (self.tmp / name).write_text("")
        (self.tmp / "b_5").mkdir()
        self.assertEqual(autostart.discover("b_{k}", self.tmp), [1, 2, 10])
        with self.assertRaises(ValueError):
            autostart.discover("dir{k}/b", self.tmp)

    @staticmethod
    def _cfg(xyz, q, names=None):
        d = np.linalg.norm(xyz - q, axis=1)
        return (names if names is not None else ["13C"] * len(d), d)

    def test_counted_starts_in_bulk(self):        # plan "자동 시작값 — 2D·섞인 bath까지, 단순판" S1-S3
        rng = np.random.default_rng(3)
        q = np.zeros(3)
        rho = 2e-3
        configs = [self._cfg(rng.uniform(-60, 60, (int(rho * 120 ** 3), 3)), q) for _ in range(6)]
        props = {"13C": (0.5, 6.728)}
        sd = autostart.count_seeds(configs, 0.0, 60.0, props, autostart.GAMMA_ELECTRON, 1.0, 0.0, n_rdip=5, n_rbath=3)
        b = (3 * autostart.N_SPINS / (4 * math.pi * rho)) ** (1 / 3)
        self.assertAlmostEqual(sd["b_n"] / b, 1.0, delta=0.05)                 # the old 3D radius, from counts
        self.assertEqual(len(sd["rbath_candidates"]), 3)
        self.assertEqual(len(sd["rdip_rungs"]) + 1, 5)
        dense = [self._cfg(rng.uniform(-60, 60, (int(8 * rho * 120 ** 3), 3)), q) for _ in range(6)]
        sd8 = autostart.count_seeds(dense, 0.0, 60.0, props, autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertAlmostEqual(sd8["b_n"] / sd["b_n"], 0.5, delta=0.05)          # density x 8: half the radius
        self.assertAlmostEqual(sd8["t_seed"] / sd["t_seed"], 1 / 8, delta=0.02)  # and an eighth of the clock

    def test_counted_starts_for_a_layer_and_a_mixture(self):
        rng = np.random.default_rng(4)
        sigma = 0.05                                   # spins per A^2 in the plane z = 0, the qubit in it
        pts = lambda n: np.column_stack([rng.uniform(-60, 60, (n, 2)), np.zeros(n)])
        configs = [self._cfg(pts(int(sigma * 120 ** 2)), np.zeros(3)) for _ in range(6)]
        sd = autostart.count_seeds(configs, 0.0, 60.0, {"13C": (0.5, 6.728)}, autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertAlmostEqual(sd["b_n"] / math.sqrt(autostart.N_SPINS / (math.pi * sigma)), 1.0, delta=0.05)
        # two species: their rates add; a species absent from the start sphere is named and left out of the clock
        mixed = []
        for _ in range(6):
            xyz = rng.uniform(-60, 60, (3000, 3))
            names = ["13C"] * 2998 + ["E", "E"]
            xyz[-2:] = [[55, 0, 0], [0, 56, 0]]          # two electrons far out
            mixed.append(self._cfg(xyz, np.zeros(3), names))
        props = {"13C": (0.5, 6.728), "E": (0.5, autostart.GAMMA_ELECTRON)}
        sm = autostart.count_seeds(mixed, 0.0, 60.0, props, autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertAlmostEqual(sm["omega"], sum(v["omega"] for v in sm["species"].values()))
        self.assertEqual(sm["species"]["E"]["count"], 0.0)
        self.assertTrue(any("species E has no spin within the start" in x for x in sm["notes"]))
        self.assertGreater(sm["r_seed"], (sm["species"]["13C"]["d_ss"] * sm["t_seed"] / autostart.PAIR_PHASE) ** (1 / 3))  # E pairs set rdip

    def test_plane_coverage(self):                   # plan S4: sqrt(a^2 + h^2)
        xyz = np.column_stack([np.random.default_rng(5).uniform(-20, 20, (4000, 2)), np.zeros(4000)])
        xyz[:2, :2] = [[-20, -20], [20, 20]]
        self.assertAlmostEqual(autostart.plane_coverage(xyz, np.array([0.0, 0.0, 30.0])), math.sqrt(20 ** 2 + 30 ** 2), places=6)
        self.assertAlmostEqual(autostart.plane_coverage(xyz, np.array([5.0, 0.0, 0.0])), 15.0, places=6)

    def test_counted_start_edge_cases(self):       # review 20261008-autoconv-counted-starts #3, #9, #10
        # a slab: the height is to its nearest occupied layer, not to its mean
        rng = np.random.default_rng(8)
        slab = np.column_stack([rng.uniform(-20, 20, (4000, 2)), rng.uniform(-3, 3, 4000)])
        slab[:2] = [[-20, -20, -3], [20, 20, 3]]
        self.assertAlmostEqual(autostart.plane_coverage(slab, np.array([0.0, 0.0, 30.0])), math.hypot(20, 27), places=6)
        # a tied shell that reaches 300 stays inside the rounded start
        d = np.full(300, 10.0004)
        sd = autostart.count_seeds([(["13C"] * 300, d)] * 2, 0.0, 20.0, {"13C": (0.5, 6.728)}, autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertGreaterEqual(sd["b_start"], 10.0004)
        self.assertEqual(sd["species"]["13C"]["count"], 300.0)
        # a species beyond the coverage is still named, with the reason
        names = ["13C"] * 300 + ["E"]
        dist = np.concatenate([np.full(300, 10.0), [60.0]])
        sm = autostart.count_seeds([(names, dist)] * 2, 0.0, 50.0, {"13C": (0.5, 6.728), "E": (0.5, autostart.GAMMA_ELECTRON)},
                                   autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertIn("E", sm["species"])
        self.assertTrue(any("beyond what the bath files cover" in x for x in sm["notes"]))

    def test_seed_boundaries(self):          # review 20261008-autoconv-autostart #7, #9, #11
        rng = np.random.default_rng(6)
        few = [self._cfg(rng.uniform(-10, 10, (50, 3)), np.zeros(3)) for _ in range(3)]
        props = {"13C": (0.5, 6.728)}
        clipped = autostart.count_seeds(few, 0.0, 10.0006, props, autostart.GAMMA_ELECTRON, 1.0, 0.0)
        self.assertIsNone(clipped["b_n"])
        self.assertLessEqual(clipped["rbath_candidates"][-1], 10.0006)          # rounded down, inside the coverage
        self.assertTrue(any("limited by the source bath" in x for x in clipped["notes"]))
        many = [self._cfg(rng.uniform(-60, 60, (3000, 3)), np.zeros(3)) for _ in range(3)]
        cut = 10.187097945369365
        sd = autostart.count_seeds(many, 0.0, 60.0, props, autostart.GAMMA_ELECTRON, 1.0, cut)
        self.assertTrue(all(np.float32(r) > np.float32(cut) for r in sd["rdip_rungs"]))
        for p, dm in (({"13C": (0.5, 0.0)}, 1.0), (props, 0.0)):
            with self.assertRaisesRegex(ValueError, "no coupling clock"):
                autostart.count_seeds(many, 0.0, 60.0, p, autostart.GAMMA_ELECTRON, dm, 0.0)
        empty = autostart.count_seeds(many, 0.0, 60.0, props, autostart.GAMMA_ELECTRON, 1.0, 1e3)
        self.assertEqual(empty["rdip_rungs"], [])                               # not an error until it is used (#8)

    def test_gyros_from_the_file_then_the_table(self):
        g = self.tmp / "gyro"
        g.write_text("13C 0.5 9.0\n")
        out = autostart.gyros(g, {"13C", "14N"})
        self.assertEqual(out["13C"], (0.5, 9.0, "gyrofile"))
        self.assertEqual(out["14N"][2], "built-in isotope table")
        self.assertAlmostEqual(out["14N"][1], 1.93377963, places=6)


class TestAutostartRuns(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_k_range_chooses_the_configurations(self):      # A65: the configuration block's k_range, else the files found
        from autoconv.runner import Runner
        inp, fake = make_auto(self.tmp)
        j = json.loads(inp.read_text())
        j["configuration"].update(path="data", bathfile=["bath_{k}"], k_range=[5, 26])   # path in front of the block's names
        inp.write_text(json.dumps(j))
        r = Runner(inp, fake)
        r.resolve()
        self.assertEqual(sorted(r.pilot + r.pool), list(range(5, 27)))
        self.assertEqual(r.patterns["bathfile"], "data/bath_{k}")
        child, _ = r.child(7, 0, 50, 30.0, 20.0, ())
        self.assertEqual(child["configuration"]["bathfile"], [str(inp.parent / "data/bath_7")])
        self.assertFalse({"k_range", "path"} & set(child["configuration"]))   # one configuration, absolute names
        self.assertEqual(child["output"], {"file": "out"})
        j["configuration"]["k_range"] = [5, 40]                  # a k without its bath file
        inp.write_text(json.dumps(j))
        r = Runner(inp, fake)
        with self.assertRaisesRegex(Exception, "missing files"):
            r.resolve()

    def test_bath_files_alone_converge(self):
        inp, fake = make_auto(self.tmp, num_candi_rdip=5, num_candi_rbath=3)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(sorted(s["found_automatically"]), ["deltat", "rbath.candidates", "rdip.candidates", "sampling.production", "sampling.test"])
        self.assertEqual((len(s["test"]), len(s["test"]) + len(s["production_available"])), (10, 30))
        start = (inp.parent / "auto/result/summary.txt").read_text()          # section 0 and the window (summary.txt, 2026-10-10)
        for word in ("30 found", "13C spin 0.5", "300 spins at", "rdip candidates", "estimate", "150 points"):
            self.assertIn(word, start)
        self.assertGreaterEqual(len(s["rdip"]["values"]), 5)             # the five first candidates, raised above when needed
        self.assertGreaterEqual(len(s["rbath"]["values"]), 4)            # three first, and the confirming step
        self.assertEqual(s["rdip"]["need"], 2)
        run = json.loads(next((inp.parent / "auto/runs/production").glob("*/input.json")).read_text())
        self.assertEqual(run["nstep"], 150)          # the default (101 until 2026-10-10: too few for the envelope step)
        self.assertAlmostEqual(run["deltat"] * 149, s["window"]["recommended_window"], places=9)

    def test_final_rbath_needs_two_rdip_passes(self):   # review 20261008-autoconv-autostart #2 and #4
        inp, fake = make_auto(self.tmp, num_candi_rdip=5, num_candi_rbath=3)
        code, out = call(inp, fake, dict(FAKE_KINK="30"))    # at rbath >= 30 the third largest rdip is off by 10 %
        self.assertIn(code, (0, 3), out)                       # converged or an honest truncation-unresolved, never a crash
        s = summary(inp)
        first = s["rdip_at_final_rbath"][0]
        self.assertEqual(first["top_neighbours"], ["pass", "fail"])         # one pass at the top is not enough
        self.assertIn("reconverged", first)
        text = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("rdip rechecked at rbath", text)
        # step 2 of summary.txt names the final rdip when the recheck at the larger rbath raised it (2026-11-05)
        final, chosen = (s.get("truncation") or {}).get("chosen", {}).get("rdip"), s["rdip"].get("chosen")
        self.assertEqual("final rdip" in text, final is not None and chosen is not None and final != chosen, text)

    def test_summary_names_the_final_rdip_when_the_recheck_raised_it(self):   # 2026-11-05: step 2 read only "chosen"
        inp, fake = make_auto(self.tmp, num_candi_rdip=5, num_candi_rbath=3)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        from autoconv.runner import Runner
        r = Runner(inp, fake)
        r.resolve()
        lines = lambda: [l for l in r._summary_lines(s, "converged", None, None) if "final rdip" in l]
        self.assertEqual(lines(), [])                                      # the recheck kept the rdip of step 2
        s["truncation"]["chosen"]["rdip"] = s["rdip"]["chosen"] * 1.25     # as if the recheck at the larger rbath raised it
        self.assertEqual(lines(), [f"  final rdip {s['rdip']['chosen'] * 1.25:g}: {s['rdip']['chosen']:g} did not hold "
                                   "at the larger rbath (step 3)"])

    def test_the_window_is_found_again_when_the_decay_time_moves(self):   # #3 (a full step after rounding) and #12
        inp, fake = make_auto(self.tmp, num_candi_rdip=5, num_candi_rbath=3)
        code, out = call(inp, fake, dict(FAKE_SHAPE="30:3"))   # same area, another decay time from rbath 30 up
        s = summary(inp)
        # the area rule cannot see a shape change, so a small rbath may still be chosen; its production curve then does not
        # decay in the new window and the run says so (window-unresolved at the production) instead of passing
        self.assertTrue(code == 0 or (code == 3 and s["stage"] == "production"), out)
        self.assertTrue(any(r.get("moved") for r in s["window"]["rounds"]), s["window"]["rounds"])
        self.assertIn("window_2", [json.loads(x)["id"] for x in (inp.parent / "auto/decisions.jsonl").read_text().splitlines() if x.strip()])
        self.assertIn("found again at rbath", (inp.parent / "auto/result/summary.txt").read_text())

    def test_window_options_with_deltat_found(self):     # #6 a start far too long, #10 revival factor 0 (#5 envelope: the method key was removed)
        for i, (kw, gyro, check) in enumerate((({}, 0.5, "coarse_grid"),
                                               (dict(revival_check_factor=0.0), GYRO, None))):
            d = Path(self.tmp) / f"c{i}"; d.mkdir()
            inp, fake = make_auto(d, gyro=gyro, nspin=1500, **kw)
            code, out = call(inp, fake)
            self.assertEqual(code, 0, (kw, out))
            if check:
                self.assertTrue(any(r.get(check) for r in summary(inp)["window"]["rounds"]), summary(inp)["window"])

    def test_an_unconfirmed_rbath_is_reported_as_advice(self):   # review 20261008-autoconv-confirmation-rule
        inp, fake = make_auto(self.tmp, nspin=500, num_candi_rbath=2)   # 300 spins lie beyond the files: [40, 50], then 62.5 > 50
        code, out = call(inp, fake)
        s = summary(inp)
        self.assertEqual((code, s["status"], s["rbath"]["chosen"]), (3, "truncation-unresolved", None), out)
        self.assertEqual(s["rbath"]["advisory"], s["rbath"]["values"][0])
        txt = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("rbath advisory (not confirmed)", txt)
        self.assertIn("confirming needs rbath", txt)
        own = next(x for x in s["rbath"]["vs_reference"] if x["value"] == s["rbath"]["advisory"])
        self.assertIn(f"({100 * own['rel']:.2f} %)", txt)                       # its own difference (review #8)
        self.assertFalse((inp.parent / "auto/runs/production").exists())      # no production from an advisory

    def test_a_given_rdip_is_not_blocked_by_the_seed(self):   # #8: inner rdip above the seed's ladder, rdip given
        inp, fake = make_auto(self.tmp, nspin=1500, base=dict(rdip=[15.0, 30.0]), candi_rdip=[20, 30], candi_rbath=[20, 25])
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)

    def test_given_lists_win_and_a_list_with_its_count_is_refused(self):
        inp, fake = make_auto(self.tmp, candi_rbath=[20, 25], candi_rdip=[10, 20, 30])
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(sorted(s["found_automatically"]), ["deltat", "sampling.production", "sampling.test"])
        self.assertEqual(s["rdip"]["need"], 1)
        self.assertEqual(s["rdip"]["values"][:3], [10.0, 20.0, 30.0])
        j = json.loads(inp.read_text()); j["auto_convergence"]["rdip"]["num_candidates"] = 4; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("not both", out)

    def test_a_flat_bath_is_taken_as_a_layer_and_a_mixture_runs(self):   # plan S4, S1-S3
        d = Path(self.tmp) / "flat"; d.mkdir()
        inp, fake = make_auto(d, nspin=2000, flat=True)
        code, out = call(inp, fake)
        self.assertIn(code, (0, 3), out)                          # taken as a layer (user 2026-10-08), not refused
        start = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("bath: a layer", start)
        d = Path(self.tmp) / "mixed"; d.mkdir()
        inp, fake = make_auto(d, species=("13C", "14N"))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        start = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("14N spin", start)
        self.assertIn("uncalibrated start", start)

    def test_too_few_configurations_are_reported_and_missing_state_files_refused(self):
        d = Path(self.tmp) / "few"; d.mkdir()
        inp, fake = make_auto(d, ids=tuple(range(1, 10)), nspin=300)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)                   # computed and reported short (review 69f9b6 R1), not refused
        self.assertIn("the minimum is 20", out)
        d = Path(self.tmp) / "fewer"; d.mkdir()
        inp, fake = make_auto(d, ids=(1, 2, 3), nspin=300)
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)                   # the pilot needs 2 and the production 2
        self.assertIn("3 configurations found", out)
        d = Path(self.tmp) / "states"; d.mkdir()
        inp, fake = make_auto(d, nspin=300)
        j = json.loads(inp.read_text())
        j["configuration"]["s_range"] = [1, 2]
        j["configuration"]["bath_state_file"] = [str(d / "data/cfg{k}_state")]
        inp.write_text(json.dumps(j))
        for k in range(1, 31):
            if k != 5:                       # configuration 5 lacks its state files
                for s in (1, 2):
                    (d / f"data/cfg{k}_state{s}").write_text("300\n" + "0.5\n" * 300)
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("missing files", out)
        self.assertIn("5: ", out)


if __name__ == "__main__":
    unittest.main()
