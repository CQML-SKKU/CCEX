"""Tests of the --auto_conv runner (autoconv/runner.py) against fake_ccex.py, whose answer is fixed beforehand:
the last two candidates are compared (converged?) and then each against the largest (2026-10-07): with rdip rungs 20/30/40
under rdip_max 50, T2 differs from 50 by 4.0/1.1/0.24 % (the fit-line maxima are about 0.48 x those: 0.019 fail, 0.005 and
0.001 pass) -> rdip 30; with rbath 10/30/40, by 5.6/0.16 % -> rbath 30. The runner must pick rbath 30, rdip 30.

    cd example/tools/AutoConvergence && python -m unittest tests.test_runner -v
"""
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

PY = sys.executable


def safe(env):
    """python -m autoconv as main.cpp starts it: the tool on PYTHONPATH and PYTHONSAFEPATH=1, so the current folder is not put
    first on sys.path (review 20261014 tools 5; Python 3.11+, and an older Python ignores the variable, where -P would stop it)."""
    old = env.get("PYTHONPATH")
    return dict(env, PYTHONSAFEPATH="1", PYTHONPATH=str(HERE.parent) + (os.pathsep + old if old else ""))


def nest(flat):
    """The nested auto_convergence block (time_window / rdip / rbath / configurations) from the runner's flat names; names
    the runner does not know stay at the top (the refusal cases)."""
    from autoconv.runner import WHERE
    out = {}
    for k, v in flat.items():
        if k in WHERE:
            sec, key = WHERE[k].split(".")
            out.setdefault(sec, {})[key] = v
        elif k == "configurations" and isinstance(v, dict):
            out.setdefault("configurations", {}).update({dict(pilot="test", pool="production").get(x, x): y for x, y in v.items()})
        else:
            out[k] = v
    return out


def make_case(tmp, pilot=(1, 2, 3, 4), pool=tuple(range(10, 40)), nspin=200, box=100.0, nstep=99, **block):
    tmp = Path(tmp)
    data = tmp / "data"; data.mkdir()
    rng = np.random.default_rng(5)
    for k in list(pilot) + list(pool):
        xyz = rng.uniform(0, box, (nspin, 3))          # the bath files cover rbath up to about box / 2 around the qubit
        (data / f"bath_{k}").write_text(f"{nspin + 1:.5f} 0 0 0\n" + "".join(f"{x:.5f} {y:.5f} {z:.5f} 13C\n" for x, y, z in xyz))
        for s in range(1, 31):
            (data / f"cfg{k}_state{s}").write_text(f"{nspin}\n" + "".join(f"{v:.5f}\n" for v in rng.choice([0.5, -0.5], nspin)))
    (data / "qubit").write_text(f"{box / 2} {box / 2} {box / 2}\n")
    fake = tmp / "fake_main.out"
    fake.write_text(f"#!/bin/sh\nexec {PY} {HERE / 'fake_ccex.py'} \"$@\"\n"); fake.chmod(0o755)
    blk = dict(workdir="auto", configurations=dict(bathfile=str(data / "bath_{k}"), statefile=str(data / "cfg{k}_state"),
                                                                pilot=list(pilot), pool=list(pool)),
               candi_rdip=[20, 30, 40, 50], candi_rbath=[10, 30, 40],
               threshold_config=0.06, min_configurations=2)   # the targets these cases were built for
    blk.update(block)
    inp = dict(method="cCCE", quantity="coherence", order=2, Qubit=dict(file=str(data / "qubit")), gyrofile=str(data / "qubit"),
               deltat=0.02, nstep=nstep, nstate=4, bfield=500, auto_convergence=nest(blk))
    (tmp / "input.json").write_text(json.dumps(inp, indent=1))
    return tmp / "input.json", fake


def call(inp, fake, env=None):
    e = dict(os.environ, **(env or {}))
    for v in ("PMI_RANK", "PMI_SIZE", "OMPI_COMM_WORLD_RANK", "MPI_LOCALRANKID"):
        if not (env and v in env): e.pop(v, None)
    p = subprocess.run([PY, "-m", "autoconv", "run", str(inp), "--binary", str(fake)], cwd=HERE.parent, env=safe(e),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def summary(inp):
    return json.loads((inp.parent / "auto/result/json/summary.json").read_text())


def count(path):
    return len(Path(path).read_text().splitlines()) if Path(path).exists() else 0


class TestRunner(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_converges_to_the_known_answer_and_restarts_without_new_runs(self):
        inp, fake = make_case(self.tmp)
        cnt = Path(self.tmp) / "count"
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual((s["truncation"]["chosen"]["rbath"], s["truncation"]["chosen"]["rdip"]), (30.0, 30.0))
        self.assertGreater(s["window"]["pilot_nstep"], 40)                  # the trial window had to grow
        self.assertEqual((s["rdip"]["estimate"]["status"], s["rdip"]["estimate"]["candidate"]), ("accepted", 30.0))   # trend + one check
        self.assertLessEqual(s["production"]["se"], 0.06)
        self.assertTrue((inp.parent / "auto/result/coherence_mean.txt").is_file())
        prod = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("rbath 30 Å, rdip 30 Å", prod)
        self.assertIn("6. converged coherence", prod)
        self.assertIn("sampling error", prod)
        self.assertIn("threshold_rdip 2 % / threshold_rbath 2 %", s["truncation"]["rule"])
        n1 = count(cnt)
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual((code, count(cnt)), (0, n1), out)                    # nothing recomputed
        # a changed output is detected and only that run is redone
        victim = next((inp.parent / "auto/runs/production").iterdir()) / "out_state1_wiDiv"
        victim.write_text(victim.read_text().replace("+", "-", 3))
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual((code, count(cnt)), (0, n1 + 1), out)

    def test_nstate_0_is_the_ensemble(self):
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"].pop("statefile"); j["nstate"] = 0
        inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual((s["production"]["states"], s["sampling_pilot"]["allocation"]["best"]["S"]), (0, 1))
        run = next((inp.parent / "auto/runs/production").iterdir())
        child = json.loads((run / "input.json").read_text())
        self.assertEqual((child["nstate"], "bath_state_file" in child), (0, False))
        self.assertTrue((run / "out_wiDiv").is_file())
        txt = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("configurations x ensemble", txt)
        self.assertNotIn("4. states", txt)                               # no state test for the ensemble
        self.assertFalse((inp.parent / "auto/result/4_states.png").exists())
        self.assertTrue((inp.parent / "auto/result/5_configurations.png").exists())

    def test_nstate_and_statefile_must_agree(self):
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text())
        for nstate, drop, msg in ((0, False, "does not read state files"), (4, True, "exstatefile needs configurations.statefile")):
            k = dict(j, nstate=nstate)       # nstate > 0 without statefile is allowed now (states made per run, 2026-10-09)
            k["auto_convergence"] = dict(j["auto_convergence"], configurations={x: v for x, v in j["auto_convergence"]["configurations"].items()
                                                                                  if not (drop and x == "statefile")})
            if drop:
                k["auto_convergence"]["configurations"]["exstatefile"] = "data/ex{k}_state"
            inp.write_text(json.dumps(k))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out)
            self.assertIn(msg, out)

    def test_core_hours_count_the_threads(self):              # review 10-07 A9
        inp, fake = make_case(self.tmp, threshold_rdip=0.05, threshold_rbath=0.05)
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="4"))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertIn("threshold_rdip 5 % / threshold_rbath 5 %", s["truncation"]["rule"])
        done = [json.loads(p.read_text()) for p in (inp.parent / "auto/runs").rglob("done.json")]
        self.assertAlmostEqual(s["core_hours"], sum(d["seconds"] * d["ranks"] * 4 for d in done) / 3600, delta=1e-9)

    def test_a_restart_counts_reused_runs_at_their_own_cost(self):   # review 20261008-autoconv-d42 #2
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="1"))
        self.assertEqual(code, 0, out)
        first = summary(inp)["core_hours"]
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="4"))     # every run reused: still what it cost at one thread
        self.assertEqual(code, 0, out)
        self.assertAlmostEqual(summary(inp)["core_hours"], first, delta=1e-12)

    def test_core_hours_unknown_when_the_counts_cannot_be_read(self):   # not refused: core-hours are only reported
        from autoconv.runner import Runner
        for launcher, n in ((["mpirun", "-n", "4"], 4), (["mpirun", "-np", "8"], 8), (["srun", "--ntasks=4"], 4),
                            (["srun", "--ntasks", "3"], 3), (["srun", "-n2"], 2), ([], 1), (["mpirun"], None), (["srun", "-N", "2"], None)):
            self.assertEqual(Runner._ranks(launcher), n)
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="4,2"))
        self.assertEqual(code, 0, out)
        self.assertIsNone(summary(inp)["core_hours"])
        self.assertIn("core-hours: unknown", (inp.parent / "auto/result/summary.txt").read_text())

    def test_sizes_are_counted_as_ccex_selects_them(self):    # review 20261008-autoconv-d42 #3: the inner rbath and rdip
        from autoconv.runner import Runner
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j.update(rbath=[15, 30], rdip=[8, 20]); inp.write_text(json.dumps(j))
        r = Runner(inp, fake)
        xyz = np.loadtxt(inp.parent / "data/bath_1", skiprows=1, usecols=(0, 1, 2))
        d = np.linalg.norm(xyz - 50, axis=1)
        kept = xyz[(d <= 30) & (d >= 15)]
        self.assertEqual(len(r.inside(1, 30)), len(kept))
        pd = np.linalg.norm(kept[:, None] - kept[None], axis=2)[np.triu_indices(len(kept), 1)]
        self.assertEqual(r.pair_count(kept, 20), int(((pd < 20) & (pd > 8)).sum()))

    def test_nodiv_false_runs_on_wiDiv_alone(self):
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["savemode"] = {"mode": "coherence", "nodiv": False}; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        run = next((inp.parent / "auto/runs/production").iterdir())
        self.assertFalse(any(p.name.endswith("_noDiv") for p in run.iterdir()))

    def test_a_given_window_is_used_as_it_is(self):
        inp, fake = make_case(self.tmp, converge_time_window=False, set_time_window=2.0)              # shorter than this curve's decay (about 3.4 ms)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual((s["window"]["status"], s["window"]["pilot_nstep"], s["production"]["nstep"]), ("given", 99, 99))
        self.assertTrue(any("given window" in x for x in s["window"]["warnings"]))   # not decayed: warned, not stopped
        txt = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("window given", txt)
        self.assertIn("WARNINGS\n  no T2 fit inside the given window", txt)

    def test_set_time_window_takes_a_unit(self):
        inp, fake = make_case(self.tmp, converge_time_window=False, set_time_window=2000, set_time_window_unit="us")   # the same 2 ms as the test above
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertEqual((summary(inp)["window"]["status"], summary(inp)["window"]["pilot_nstep"]), ("given", 99))

    def test_only_rdip_and_configurations_converge(self):     # user 2026-10-08: converge true/false per section
        inp, fake = make_case(self.tmp, converge_time_window=False, set_time_window=6.0, converge_rbath=False, value_rbath=30.0)
        j = json.loads(inp.read_text()); j["auto_convergence"]["rbath"].pop("candidates"); inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(s["converge"], dict(time_window=False, rdip=True, rbath=False, configurations=True))
        self.assertEqual((s["truncation"]["chosen"]["rbath"], s["rbath"]["values"]), (30.0, [30.0]))
        self.assertIsNotNone(s["truncation"]["chosen"]["rdip"])
        self.assertFalse((inp.parent / "auto/runs/rbath").exists())          # no rbath runs
        self.assertTrue((inp.parent / "auto/result/2_rdip.png").is_file())
        self.assertFalse((inp.parent / "auto/result/3_rbath.png").exists())
        self.assertIn("allocation", s["sampling_pilot"])

    def test_t2_fit_settings_reach_every_fit(self):     # t2_fit section (user 2026-10-08)
        inp, fake = make_case(self.tmp, fit_points="all", fit_floor=0.1)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)                     # the fake decays with n 1.3: inside the default n range
        d = Path(self.tmp) / "narrow"; d.mkdir()
        inp, fake = make_case(d, fit_n_range=[0.3, 1.0])   # n 1.3 is outside: every fit fails, the window cannot be found
        code, out = call(inp, fake)
        self.assertEqual((code, summary(inp)["status"]), (3, "window-unresolved"), out)
        self.assertIn("the last: T2 fit", summary(inp)["reason"])      # the last window's own reason, not "not decayed"

    def test_a_failed_fit_is_judged_on_the_raw_area(self):     # user 2026-10-12, plan 2.1
        # n 1.3 below rbath 30 is outside n_range [1.4, 4]: rbath 10 has no fit and is judged on the raw area (fails),
        # rbath 30 and 40 (n 1.6, area kept) on the fit as before
        inp, fake = make_case(self.tmp, fit_n_range=[1.4, 4.0])
        code, out = call(inp, fake, dict(FAKE_SHAPE="30:1.6"))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        v = {x["value"]: (x["verdict"], x.get("basis")) for x in s["rbath"]["vs_reference"]}
        self.assertEqual(v, {10.0: ("fail", "raw-area"), 30.0: ("pass", None)})
        self.assertEqual(s["truncation"]["chosen"]["rbath"], 30.0)
        dec = [json.loads(x) for x in (inp.parent / "auto/decisions.jsonl").read_text().splitlines() if '"truncation"' in x][0]
        self.assertIn(dict(parameter="rbath", lo=10.0, hi=40.0, verdict="fail"), dec["raw_area"])
        self.assertIn("judged on the raw area of |L| on [0, T_d] (a T2 fit failed): 10 vs 30, 10 vs 40",
                      (inp.parent / "auto/result/summary.txt").read_text())
        self.assertNotIn("basis", json.dumps(s["rdip"]))                 # every rdip fit succeeded: unchanged

    def test_progress_interruption_status_and_fresh(self):   # plan "끊김 대비·진행 보기" R1-R4 (2026-10-08)
        inp, fake = make_case(self.tmp)
        cnt = Path(self.tmp) / "count"
        status = lambda: subprocess.run([PY, "-m", "autoconv", "status", str(inp)], cwd=HERE.parent, env=safe(os.environ),
                                        capture_output=True, text=True).stdout
        self.assertIn("not started", status())
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        prog = (inp.parent / "auto/result/progress.txt").read_text()
        for word in ("status: converged", "finished CCEX runs: pilot", "chosen truncation: rbath 30", "estimated configurations"):
            self.assertIn(word, prog)
        self.assertIn("5 production", json.loads((inp.parent / "auto/result/json/progress.json").read_text())["steps_done"])   # review #7
        self.assertIn("state: converged", status())
        # a job cut while running: the result says running, nobody holds the lock
        res = inp.parent / "auto/result/json/summary.json"
        res.write_text(json.dumps(dict(status="running", reason=None)))
        self.assertIn("state: interrupted", status())
        with open(inp.parent / "auto/lock", "w") as f:              # a live runner holds the lock
            fcntl.flock(f, fcntl.LOCK_EX)
            self.assertIn("state: running", status())
        n1 = count(cnt)
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual((code, count(cnt)), (0, n1), out)            # continued: nothing computed again
        self.assertIn("was interrupted", out)
        self.assertIn("the run before was interrupted", (inp.parent / "auto/result/progress.txt").read_text())
        # resume fresh: the earlier work kept aside, everything computed again
        j = json.loads(inp.read_text()); j["auto_convergence"]["resume"] = "fresh"; inp.write_text(json.dumps(j))
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        self.assertTrue((inp.parent / "auto.old1/runs").is_dir())
        self.assertEqual(count(cnt), 2 * n1)
        # status: a held lock without any summary is running; a damaged summary is reported, not a crash (review #5)
        shutil.rmtree(inp.parent / "auto/result")
        with open(inp.parent / "auto/lock", "w") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            self.assertIn("state: running", status())
            # a refused second call does not overwrite the live runner's summary (review #6)
            (inp.parent / "auto/result/json").mkdir(parents=True)
            (inp.parent / "auto/result/json/summary.json").write_text(json.dumps(dict(status="running", reason=None)))
            k = json.loads(inp.read_text()); k["auto_convergence"]["unknown_key"] = 1; inp.write_text(json.dumps(k))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out)
            self.assertEqual(summary(inp)["status"], "running")
        (inp.parent / "auto/result/json/summary.json").write_text("{")
        self.assertIn("state: unknown", status())
        j["auto_convergence"]["resume"] = "later"; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("resume must be continue or fresh", out)

    def test_state_files_are_counted(self):              # states_available removed 2026-10-08
        inp, fake = make_case(self.tmp, converge_config=False, number_config=5, states_config=20)
        for s_ in range(13, 31):
            (inp.parent / f"data/cfg25_state{s_}").unlink()       # configuration 25 of the pool holds 12 states
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("state files of every pool configuration hold (12)", out)

    def test_fixed_rdip_and_fixed_configurations(self):
        inp, fake = make_case(self.tmp, converge_rdip=False, value_rdip=40.0, converge_config=False, number_config=5, states_config=2)
        j = json.loads(inp.read_text()); j["auto_convergence"]["rdip"].pop("candidates"); inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual((s["rdip"]["values"], s["truncation"]["chosen"]["rdip"]), ([40.0], 40.0))
        self.assertEqual((len(s["production"]["configurations"]), s["production"]["states"], s["production"]["batches"]), (5, 2, []))
        self.assertEqual(s["sampling_pilot"], dict(fixed=dict(configurations=5, states=2)))
        self.assertTrue((inp.parent / "auto/result/5_configurations.png").exists())       # the production growing (2026-10-10)
        runs = [json.loads(p.read_text()) for p in (inp.parent / "auto/runs").rglob("input.json")]
        self.assertTrue(all(r["rdip"] == 40.0 and "rdip_ladder" not in r for r in runs))

    def test_a_plateau_raises_rbath(self):
        inp, fake = make_case(self.tmp, candi_rbath=[10, 20])
        code, out = call(inp, fake, dict(FAKE_PLATEAU_BELOW="24"))
        s = summary(inp)
        self.assertEqual([r["rbath"] for r in s["rbath_raised_for_the_window"]], [25.0], out)   # 20 plateaus, 25 decays
        self.assertIn("rbath raised to 25", (inp.parent / "auto/result/summary.txt").read_text())

    def test_a_plateau_is_checked_in_a_longer_window(self):     # review 2026-10-08 R1
        inp, fake = make_case(self.tmp, nstep=299)      # 5.96 ms; the curve drops at 7 ms
        code, out = call(inp, fake, dict(FAKE_LATE_DROP="7"))
        runs = [json.loads((d / "input.json").read_text()) for d in (inp.parent / "auto/runs/pilot").iterdir()]
        self.assertEqual(sorted({(r["nstep"], round(r["deltat"] * (r["nstep"] - 1), 6), r["rbath"]) for r in runs})[:2],
                         [(299, 5.96, 40.0), (299, 8.94, 40.0)], out)   # a longer window on the same nstep, same bath

    def test_a_plateau_at_the_bath_limit_stops(self):
        inp, fake = make_case(self.tmp, candi_rbath=[10, 20], box=43.0)   # the files cover about 20 A: 25 is beyond
        code, out = call(inp, fake, dict(FAKE_PLATEAU_BELOW="24"))
        self.assertEqual(code, 3, out)
        self.assertIn("make larger baths", summary(inp)["reason"])
        self.assertIn("covering at least", summary(inp)["reason"])          # what to make (user 2026-10-09)

    def test_rbath_is_raised_until_the_last_two_agree(self):
        inp, fake = make_case(self.tmp, candi_rbath=[10, 20])
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(s["rbath"]["raised"], [25.0], out)     # two values, no trend; 10 vs 20 fails, 20 vs 25 agrees
        self.assertEqual(s["truncation"]["chosen"]["rbath"], 20.0)
        self.assertTrue((inp.parent / "auto/result/3_rbath.png").is_file())

    def test_rbath_at_the_bath_limit_asks_for_larger_baths(self):
        inp, fake = make_case(self.tmp, candi_rbath=[10, 20], box=44.0)   # the files cover about 22 A: 25 would be beyond
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertEqual((s["stage"], s["rbath"]["raised"]), ("rbath", []))
        self.assertIn("make larger baths", s["reason"])
        self.assertIn("covering at least", s["reason"])

    def test_rbath_beyond_the_bath_files_is_refused(self):
        inp, fake = make_case(self.tmp, candi_rbath=[10, 60])
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("beyond what the bath files cover", out)

    def test_window_unresolved_when_the_curve_never_decays(self):
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(FAKE_LATE_DROP="1e9"))      # stays near 0.3: longer window, larger bath, then stop
        self.assertEqual(code, 3, out)
        self.assertEqual(summary(inp)["status"], "window-unresolved")
        s = summary(inp)                                             # every window tried, drawn and listed (user 2026-10-10)
        self.assertGreaterEqual(len(s["window"]["tried"]), 2)
        self.assertTrue((inp.parent / "auto/result/1_window.png").is_file())
        self.assertIn("1_window.png shows the curves of every window tried", s["reason"])
        self.assertIn("windows tried:", (inp.parent / "auto/result/summary.txt").read_text())

    def test_rdip_is_raised_until_the_last_two_agree(self):
        inp, fake = make_case(self.tmp, candi_rdip=[5, 10, 20])
        code, out = call(inp, fake)
        s = summary(inp)
        self.assertEqual(s["rdip"]["raised"], [25.0])   # no trend (10.2 % after 12.4 %); 10 vs 20 fails, 20 vs 25 (1.9 % < 2 %) agrees
        self.assertEqual(s["truncation"]["chosen"]["rdip"], 20.0, out)    # 20 vs the raised 25: 1.9 % (within 2 %)
        code, out = call(inp, fake)                                         # a restart repeats nothing and decides the same
        self.assertEqual(summary(inp)["rdip"]["raised"], s["rdip"]["raised"])

    def test_rdip_cannot_go_beyond_every_pair(self):
        inp, fake = make_case(self.tmp, candi_rdip=[5, 10, 20], candi_rbath=[5, 10])
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)        # rdip 20 = 2 x rbath 10 already holds every pair, 10 vs 20 fails: not converged
        self.assertEqual((summary(inp)["status"], summary(inp)["stage"], summary(inp)["rdip"]["raised"]), ("truncation-unresolved", "rdip", []))

    def test_a_failing_larger_rung_blocks_the_smaller_ones(self):
        # rung 40 is pushed off (non-monotone): 30 vs 40 and 40 vs 50 fail although 30 alone would be close to 50/60; the
        # choice must start above the last failing comparison -> 50
        inp, fake = make_case(self.tmp, candi_rdip=[30, 40, 50, 60])
        code, out = call(inp, fake, dict(FAKE_NONMONO="40"))
        self.assertEqual(code, 0, out)
        v = {p["value"]: p["verdict"] for p in summary(inp)["rdip"]["vs_reference"]}   # each against the largest
        self.assertEqual(v[40.0], "fail")
        self.assertEqual(summary(inp)["truncation"]["chosen"]["rdip"], 50.0)

    def test_a_nonfinite_rung_is_never_judged(self):          # review 10-07 A8
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(FAKE_NAN_RUNG="30"))
        s = summary(inp)
        v = {p["value"]: p["verdict"] for p in s["rdip"]["vs_reference"]}
        self.assertEqual(v[30.0], "nonfinite", out)
        self.assertNotEqual(s.get("truncation", {}).get("chosen", {}).get("rdip"), 30.0)
        self.assertIn("rdip_at_final_rbath", s)                 # C3: rdip judged again at the final largest rbath

    def test_never_fewer_than_min_configurations(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5, min_configurations=24)   # 0.5 alone would need 2
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertEqual(len(summary(inp)["production"]["configurations"]), 30)     # 1.5 x 24, within the pool of 30
        (Path(self.tmp) / "b").mkdir()
        inp2, _ = make_case(Path(self.tmp) / "b", min_configurations=31, threshold_config=0.5)   # the pool holds 30
        code, out = call(inp2, fake)
        self.assertEqual(code, 3, out)                       # computed and reported short, not refused (review 69f9b6 R1)
        self.assertEqual(len(summary(inp2)["production"]["configurations"]), 30)
        self.assertIn("the minimum is 31", summary(inp2)["reason"])

    def test_a_changed_input_starts_a_new_plan(self):       # review 10-07 A4
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"]["threshold"] = 0.05; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertIn("decisions.1.jsonl", out)
        self.assertTrue((inp.parent / "auto/decisions.1.jsonl").is_file())
        ids = summary(inp)["production"]["configurations"]
        self.assertEqual(len(ids), len(set(ids)))

    def test_a_cut_decision_line_is_dropped(self):           # review 10-07 A5
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake)
        dec = inp.parent / "auto/decisions.jsonl"
        dec.write_text(dec.read_text() + '{"id": "extra_block", "ids": [')
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        for line in dec.read_text().splitlines():
            json.loads(line)                                    # every line whole

    def test_unhashed_file_inputs_are_refused(self):          # review 10-07 A6
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake)
        self.assertEqual(summary(inp)["status"], "converged", out)
        j = json.loads(inp.read_text()); j["pcce"] = {"cache": {"file": "cache.txt"}}; inp.write_text(json.dumps(j))   # a file not hashed
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("pcce.cache", out)
        self.assertEqual(summary(inp)["status"], "refused")       # an earlier result is not left as the present one

    def test_old_savemode_forms_are_refused(self):           # 2026-10-13: save_nodiv is savemode.nodiv, "normal" is "coherence"
        inp, fake = make_case(self.tmp)
        base = json.loads(inp.read_text())
        for edit, want in ((dict(save_nodiv=False), "save_nodiv is now savemode.nodiv"), (dict(savemode="normal"), "savemode")):
            inp.write_text(json.dumps(dict(base, **edit)))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out)
            self.assertIn(want, out)

    def test_keys_are_read_without_regard_to_case(self):     # review 20261013-180049 R4: as CCEX reads them
        from autoconv.runner import Runner, _method

        def upper(x, top=True):
            if isinstance(x, list):
                return [upper(v, False) for v in x]
            if not isinstance(x, dict):
                return x
            return {(k.upper() if top else k.title()): upper(v, False) for k, v in x.items()}

        runners = []
        for name in ("lower", "upper"):
            d = Path(self.tmp) / name; d.mkdir()
            inp, fake = make_case(d)
            j = json.loads(inp.read_text())
            j["method"] = "gCCE"
            j["hamiltonian"] = {"cCCE": {"delta_A": True}}
            if name == "upper":
                j = upper(j)
                self.assertIn("AUTO_CONVERGENCE", j)
                self.assertIn("Workdir", j["AUTO_CONVERGENCE"])
                self.assertIn("Ccce", j["HAMILTONIAN"])
            inp.write_text(json.dumps(j))
            runners.append(Runner(inp, fake))
        lower, up = runners
        self.assertEqual(_method(up.base), "gcce")
        self.assertEqual(json.dumps(up.base, sort_keys=True), json.dumps(lower.base, sort_keys=True).replace("/lower/", "/upper/"))
        self.assertEqual(json.dumps(up.block, sort_keys=True), json.dumps(lower.block, sort_keys=True).replace("/lower/", "/upper/"))
        self.assertIn("Qubit", up.base)                           # the names written with capitals keep them
        self.assertIn("delta_A", up.base["hamiltonian"]["cCCE"])

    def test_keys_that_differ_only_in_case_are_refused(self):
        inp, fake = make_case(self.tmp)
        base = json.loads(inp.read_text())
        for edit in (dict(Method="gCCE"), dict(auto_convergence=dict(base["auto_convergence"], WorkDir="other"))):
            inp.write_text(json.dumps(dict(base, **edit)))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out)
            self.assertIn("are one key to CCEX (case is not read)", out)

    def test_coherence_above_one_is_reported(self):
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(FAKE_SPIKE="1"))
        self.assertEqual(code, 0, out)
        self.assertTrue(any("|L| > 1.01" in x for x in summary(inp)["warnings"]), summary(inp)["warnings"])

    def test_the_production_counts(self):              # the registered table (plan "넉넉하게 먼저", 2026-10-09)
        from autoconv.runner import first_target, next_target
        self.assertEqual(first_target(503, 20, 64), 64)                 # P1: the pool limit wins
        self.assertEqual(first_target(57, 20, 1000), 90)                # 85.5 rounded up
        self.assertEqual(first_target(8, 20, 1000), 20)                 # at least 20
        self.assertEqual(first_target(10, 21, 22), 22)                  # min_configurations 21, a pool of 22
        self.assertEqual(first_target(10, 20, 15), 15)                  # a pool below the minimum: all of it
        self.assertEqual(next_target(64, 0.0201, 0.02, 20, 1000), (65, 85, 85))   # 21 more, at least 20 added
        self.assertEqual(next_target(60, 0.03, 0.02, 20, 64)[2], 64)            # the rest of the pool
        self.assertEqual(next_target(40, 0.04, 0.02, 20, 1000), (160, 220, 220))  # 1.5 x the 120 short

    def batches(self, inp):
        dec = [json.loads(l) for l in (inp.parent / "auto/decisions.jsonl").read_text().splitlines()]
        return [d for d in dec if d["id"].startswith("production_batch_")]

    def test_a_short_pool_is_computed_a_batch_at_a_time(self):    # plan "조금씩 늘리기" (user 2026-10-09)
        inp, fake = make_case(self.tmp, pool=tuple(range(100, 140)), threshold_config=0.001)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertEqual((s["status"], s["stage"]), ("sampling-unresolved", "production"), out)
        self.assertEqual([b["configurations"] for b in s["production"]["batches"]], [40])        # the pilot asks for far more
        self.assertEqual(s["production"]["configurations"], list(range(100, 140)))
        self.assertRegex(s["reason"], r"about \d+ needed \(\d+ more\); add \d+ more bath files")
        self.assertTrue((inp.parent / "auto/result/coherence_mean.txt").is_file())       # the result is kept
        self.assertTrue((inp.parent / "auto/result/6_converged_coherence.png").is_file())
        self.assertIn("recheck", s)
        self.assertEqual([d["id"] for d in self.batches(inp)],
                         ["production_batch_0_plan", "production_batch_0_done"])
        dec = (inp.parent / "auto/result/summary.txt").read_text()
        self.assertIn("5. configurations", dec)                          # 2, 5, 10, ..., 40 in the table
        self.assertEqual([r["configurations"] for r in s["production"]["decomposition"]], [2, *range(5, 41, 5)])
        self.assertFalse(list((inp.parent / "auto/result").glob("*decomposition*")))
        full = [r for r in s["production"]["decomposition"] if r["configurations"] == 40]
        self.assertTrue(full and all(r["se"] == s["production"]["se"] for r in full))   # the full point is the verdict

    def test_no_pass_before_the_minimum(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)            # met at once: still 20
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertEqual(len(summary(inp)["production"]["configurations"]), 20)
        (Path(self.tmp) / "b").mkdir()
        inp, fake = make_case(Path(self.tmp) / "b", pool=tuple(range(100, 140)), threshold_config=0.5, min_configurations=22)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertEqual([b["configurations"] for b in summary(inp)["production"]["batches"]], [35])   # 1.5 x 22 rounded up to 5

    def test_a_pool_below_the_minimum_is_computed_and_reported(self):
        inp, fake = make_case(self.tmp, pool=tuple(range(10, 20)), threshold_config=0.5)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertEqual(len(s["production"]["configurations"]), 10)
        self.assertIn("the minimum is 20", s["reason"])

    def found(self, **kw):
        """A case whose test and production configurations are found from the files (no lists given)."""
        inp, fake = make_case(self.tmp, **kw)
        j = json.loads(inp.read_text())
        for k in ("test", "production"):
            j["auto_convergence"]["configurations"].pop(k)
        inp.write_text(json.dumps(j))
        return inp, fake

    def test_the_test_configurations_spread_over_the_nearest_spin_distance(self):   # user 2026-10-10 (P1 CR90)
        inp, fake = self.found(pool=tuple(range(10, 40)), threshold_config=0.5)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        sp = json.loads((inp.parent / "auto/split.json").read_text())
        self.assertEqual(sp["rule"], "nearest-spin quantiles")
        data = inp.parent / "data"
        q = np.array([float(x) for x in (data / "qubit").read_text().split()])
        near = {}
        for f in data.glob("bath_*"):
            xyz = np.array([[float(v) for v in l.split()[:3]] for l in f.read_text().splitlines()[1:] if l.strip()])
            near[int(f.name.split("_")[1])] = np.linalg.norm(xyz - q, axis=1).min()
        d = sorted(near.values())
        picked = sorted(near[k] for k in sp["pilot"])
        self.assertLessEqual(picked[0], d[len(d) // 5])          # the closest fifth is represented ...
        self.assertGreaterEqual(picked[-1], d[4 * len(d) // 5])  # ... and so is the farthest fifth

    def test_files_added_continue_the_production(self):        # plan "bath 만들기 빼기" (user 2026-10-09)
        inp, fake = self.found(pool=tuple(range(10, 30)), threshold_config=0.001)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        first = summary(inp)
        self.assertIn("add ", first["reason"])
        self.assertIn("call again to continue", first["reason"])
        data = inp.parent / "data"
        rng = np.random.default_rng(11)
        for k in range(100, 120):                              # the user makes 20 more
            xyz = rng.uniform(0, 100.0, (200, 3))
            (data / f"bath_{k}").write_text("201.00000 0 0 0\n" + "".join(f"{x:.5f} {y:.5f} {z:.5f} 13C\n" for x, y, z in xyz))
            for st in range(1, 31):
                shutil.copy(data / f"cfg10_state{st}", data / f"cfg{k}_state{st}")
        cnt = Path(self.tmp) / "count"
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        again = summary(inp)
        self.assertEqual(again.get("test"), first["test"], out)
        old, new = first["production"]["configurations"], again["production"]["configurations"]
        self.assertEqual(new[:len(old)], old)                  # the production goes on after what was computed
        self.assertTrue(set(new[len(old):]) <= set(range(100, 120)) and len(new) > len(old))
        self.assertEqual(count(cnt), len(new) - len(old))      # only the new configurations run
        self.assertFalse((inp.parent / "auto/decisions.1.jsonl").exists())

    def test_a_changed_production_file_starts_a_new_plan(self):
        inp, fake = self.found(pool=tuple(range(10, 30)), threshold_config=0.001)
        call(inp, fake)
        k = summary(inp)["production"]["configurations"][0]
        f = inp.parent / f"data/bath_{k}"
        lines = f.read_text().splitlines()
        lines[1] = "1.00000 " + lines[1].split(maxsplit=1)[1]              # one spin moved: a different file
        f.write_text("\n".join(lines) + "\n")
        code, out = call(inp, fake)
        self.assertTrue((inp.parent / "auto/decisions.1.jsonl").exists(), out)

    def test_a_found_pool_below_min_configurations_is_computed(self):      # review 69f9b6 R1
        inp, fake = self.found(pool=tuple(range(10, 16)), threshold_config=0.5, min_configurations=10)   # 10 files
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertEqual((len(s["test"]), len(s["production"]["configurations"])), (2, 8))
        self.assertIn("the minimum is 20", s["reason"])

    def test_added_files_with_too_few_states_are_refused(self):            # review 69f9b6 R2
        inp, fake = self.found(pool=tuple(range(10, 30)), threshold_config=0.001)
        call(inp, fake)
        S = summary(inp)["production"]["states"]
        self.assertEqual(S, 5)                                   # test_states (user 2026-10-10)
        data = inp.parent / "data"
        xyz = np.random.default_rng(3).uniform(0, 100.0, (200, 3))
        (data / "bath_100").write_text("201.00000 0 0 0\n" + "".join(f"{x:.5f} {y:.5f} {z:.5f} 13C\n" for x, y, z in xyz))
        for st in range(1, S):                                   # fewer than the production's states
            shutil.copy(data / f"cfg10_state{st}", data / f"cfg100_state{st}")
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("missing files", out)
        self.assertIn("100", out)

    def test_state_files_made_per_run(self):              # plan "state 파일 자동" (user 2026-10-09)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"].pop("statefile"); inp.write_text(json.dumps(j))
        for f in (inp.parent / "data").glob("cfg*_state*"):
            f.unlink()                                          # none stored
        seen = Path(self.tmp) / "seen"
        code, out = call(inp, fake, dict(FAKE_STATES_SEEN=str(seen)))
        self.assertEqual(code, 0, out)
        rows = [l.split() for l in seen.read_text().splitlines()]
        self.assertTrue(rows and all(r[3] != "missing" for r in rows))          # present while CCEX runs
        by = {}
        for k, st, rb, h, step in (r for r in rows if r[4] != "states"):       # the state test's blocks hold states 6, 7, ... as 1, 2, ...
            by.setdefault((k, st), set()).add((rb, h))
        multi = [v for v in by.values() if len({rb for rb, _ in v}) > 1]
        self.assertTrue(multi)                                                 # the pilot ran at several rbath
        self.assertTrue(all(len({h for _, h in v}) == 1 for v in by.values()))  # same configuration and state: same file at every rbath
        self.assertEqual(len({h for (k, st), v in by.items() for _, h in v if st == "1"}), len({k for k, _ in by}))   # differ by k
        self.assertFalse(list((inp.parent / "auto/runs").rglob("cfg_state*")))  # removed after
        self.assertIn("states: made per run", (inp.parent / "auto/result/summary.txt").read_text())

    def test_avaax_files_made_per_run(self):               # plan "avaax 파일 자동" (user 2026-10-09)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        j = json.loads(inp.read_text()); j["Defect"] = {"species": [{"name": "13c", "n_configurations": 3}]}; inp.write_text(json.dumps(j))   # case as CCEX (6548a9 R4)
        seen = Path(self.tmp) / "seen"
        code, out = call(inp, fake, dict(FAKE_AVAAX_SEEN=str(seen)))
        self.assertEqual(code, 0, out)
        rows = [l.split() for l in seen.read_text().splitlines()]
        self.assertTrue(rows and all(r[2] != "missing" for r in rows))
        by = {}
        for k, rb, h in rows:
            by.setdefault(k, set()).add((rb, h))
        self.assertTrue(any(len({rb for rb, _ in v}) > 1 for v in by.values()))
        self.assertTrue(all(len({h for _, h in v}) == 1 for v in by.values()))     # one file per configuration at every rbath
        self.assertFalse(list((inp.parent / "auto/runs").rglob("cfg_avaax")))
        self.assertIn("avaax: made per run", (inp.parent / "auto/result/summary.txt").read_text())

    def test_made_avaax_follow_the_rule(self):
        from autoconv.runner import write_avaax
        d = Path(self.tmp)
        (d / "bath").write_text("5 0 0 0\n" + "".join(f"{i} 0 0 {'P1' if i % 2 else '13C'}\n" for i in range(40)))
        write_avaax(d / "a", d / "bath", 3, 9, {"p1": 12})        # the runner keys the defects in lower case
        write_avaax(d / "b", d / "bath", 3, 9, {"p1": 12})
        self.assertEqual((d / "a").read_text(), (d / "b").read_text())
        v = [int(x) for x in (d / "a").read_text().split()]
        self.assertEqual((v[0], len(v)), (41, 41))
        self.assertTrue(all(1 <= x <= 12 for x in v[2::2]) and all(x == 0 for x in v[1::2]))   # P1 lines 1..12, 13C lines 0

    def test_made_states_follow_the_rule(self):
        from autoconv.runner import write_states
        d = Path(self.tmp)
        (d / "bath").write_text("4 0 0 0\n0 0 0 13C\n1 0 0 14N\n2 0 0 13C\n3 0 0 14N\n")
        write_states(d / "a", d / "bath", 7, 2, 5, {"13C": 0.5, "14N": 1.0})
        write_states(d / "b", d / "bath", 7, 2, 5, {"13C": 0.5, "14N": 1.0})
        self.assertEqual((d / "a1").read_text(), (d / "b1").read_text())        # reproducible
        self.assertNotEqual((d / "a1").read_text(), (d / "a2").read_text())
        v = [float(x) for x in (d / "a1").read_text().split()[1:]]
        self.assertEqual(len(v), 4)
        self.assertTrue(v[0] in (0.5, -0.5) and v[2] in (0.5, -0.5) and v[1] in (1, 0, -1) and v[3] in (1, 0, -1))

    def test_bath_generation_is_refused(self):
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["auto_convergence"]["bath_generation"] = {"size": 50}; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("make the bath files yourself", out)

    def test_a_cut_batch_is_finished_the_same_way(self):
        inp, fake = make_case(self.tmp, pool=tuple(range(100, 140)), threshold_config=0.001)
        code, out = call(inp, fake)
        first = summary(inp)
        dec = inp.parent / "auto/decisions.jsonl"
        lines = dec.read_text().splitlines()
        cut = next(i for i, l in enumerate(lines) if json.loads(l)["id"] == "production_batch_0_plan")
        dec.write_text("\n".join(lines[:cut + 1]) + "\n")              # cut after the first block was planned
        planned = json.loads(lines[cut])["ids"]
        for d in (inp.parent / "auto/runs/production").iterdir():     # half of that batch not run yet
            if int(d.name.split("-")[0]) in planned[len(planned) // 2:]:
                shutil.rmtree(d)
        cnt = Path(self.tmp) / "count"
        code2, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code2, code, out)
        again = summary(inp)
        self.assertEqual(again["production"]["configurations"], first["production"]["configurations"])
        self.assertEqual(again["production"]["se"], first["production"]["se"])
        self.assertEqual(count(cnt), len(planned) - len(planned) // 2)    # only the runs missing

    def test_an_unstable_production_configuration_is_left_out(self):     # plan "불안정 configuration 빼기" (2026-10-10)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_BIG_IDS="10"))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual([(r["configuration"], r["where"]) for r in s["unstable"]], [(10, "production")])
        self.assertIn("|L| 1e+05", s["unstable"][0]["reason"])
        self.assertIsNotNone(s["unstable"][0]["nearest_spin"])
        self.assertEqual(s["production"]["n_used"], len(s["production"]["configurations"]) - 1)
        self.assertIn("unstable configurations left out", (inp.parent / "auto/result/summary.txt").read_text())
        mean = np.loadtxt(inp.parent / "auto/result/coherence_mean.txt")
        self.assertLess(np.abs(mean[:, 3]).max(), 1.5)                     # the blow-up is not in the mean

    def test_the_state_test_converges_and_the_production_keeps_test_states(self):   # plan "state 수렴 시험(점검)" (2026-10-10)
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        st = s["states"]
        self.assertEqual(st["values"], [2, 3, 4, 5])               # within the threshold at test_states: nothing more computed
        self.assertLessEqual(st["chosen"], 5)
        self.assertIsNone(st["warning"])
        self.assertFalse((inp.parent / "auto/runs/states").exists())
        self.assertEqual(s["production"]["states"], 5)
        self.assertEqual([r["S"] for r in s["sampling_pilot"]["allocation"]["table"]], [5])
        self.assertTrue((inp.parent / "auto/result/4_states.png").is_file())
        self.assertIn("4. states", (inp.parent / "auto/result/summary.txt").read_text())

    def test_states_that_do_not_converge_only_warn(self):
        inp, fake = make_case(self.tmp, fit_points="all")          # the states, not the fit points, are what this checks
        code, out = call(inp, fake, dict(FAKE_STATE_SPREAD="0.5"))
        s = summary(inp)
        self.assertGreater(s["states"]["values"][-1], 5)             # more states computed
        self.assertTrue(s["states"]["chosen"] is None or s["states"]["chosen"] > 5)
        self.assertTrue(any(w.startswith("state test:") for w in s["warnings"]), s["warnings"])
        self.assertEqual(s["sampling_pilot"]["allocation"]["table"][0]["S"], 5)     # the production keeps test_states

    def test_a_later_state_block_is_the_same_state(self):
        # states 6-10 computed as a block of their own read the files of states 6-10 (renumbered 1-5)
        inp, fake = make_case(self.tmp, fit_points="all")
        seen = Path(self.tmp) / "seen"
        code, out = call(inp, fake, dict(FAKE_STATES_SEEN=str(seen), FAKE_STATE_SPREAD="0.5"))    # states that need more
        data = inp.parent / "data"
        import hashlib
        rows = [l.split() for l in seen.read_text().splitlines()]
        hashes = {r[3] for r in rows if r[0] == "1"}
        for st in range(6, 16):
            self.assertIn(hashlib.sha256((data / f"cfg1_state{st}").read_bytes()).hexdigest(), hashes)

    def test_states_threshold_is_checked(self):
        inp, fake = make_case(self.tmp, states_threshold=1.5)
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("states_threshold", out)

    def test_max_configurations_limits_the_production_and_raising_it_continues(self):   # plan "최대 configuration·state 수 키"
        inp, fake = make_case(self.tmp, threshold_config=0.001, max_configurations=20)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertLessEqual(len(s["production"]["configurations"]), 20)
        self.assertIn("raise configurations.max_configurations (now 20, 30 production files there)", s["reason"])
        self.assertIn("production 20 of 30 (max_configurations)", (inp.parent / "auto/result/summary.txt").read_text())
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"]["max_configurations"] = 30; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertFalse(list((inp.parent / "auto").glob("decisions.*.jsonl")), out)     # the same plan goes on
        self.assertGreater(len(summary(inp)["production"]["configurations"]), 20)

    def test_a_lowered_max_configurations_is_refused_once_planned(self):   # review 01ecda #1
        inp, fake = make_case(self.tmp, threshold_config=0.001, max_configurations=30)
        code, out = call(inp, fake)
        self.assertEqual(code, 3, out)
        self.assertIn("raise configurations.max_configurations to about", summary(inp)["reason"])   # files and cap used up (#3)
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"]["max_configurations"] = 20; inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("below the 30 production configurations already planned", out)

    def test_max_states_does_not_bound_a_fixed_production_of_made_states(self):   # review 01ecda #2
        inp, fake = make_case(self.tmp, converge_config=False, number_config=2, states_config=20, max_states=5)
        j = json.loads(inp.read_text()); j["auto_convergence"]["configurations"].pop("statefile"); inp.write_text(json.dumps(j))
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertEqual(summary(inp)["production"]["states"], 20)

    def test_max_states_bounds_the_state_test(self):
        inp, fake = make_case(self.tmp, max_states=10, fit_points="all")
        code, out = call(inp, fake, dict(FAKE_STATE_SPREAD="0.5"))
        self.assertEqual(summary(inp)["states"]["values"][-1], 10)

    def test_max_keys_are_checked(self):
        for kw, word in ((dict(max_states=3), "max_states must be an integer >= test_states"),
                         (dict(max_configurations=1), "max_configurations must be an integer >= 2"),
                         (dict(max_configurations=5, min_configurations=10), "below min_configurations")):
            d = Path(self.tmp) / str(len(word)); d.mkdir()
            inp, fake = make_case(d, **kw)
            code, out = call(inp, fake)
            self.assertEqual(code, 2, (kw, out))
            self.assertIn(word, out)

    def _fake_mpirun(self):
        """An mpirun that writes its arguments and runs the command without MPI."""
        d = Path(self.tmp) / "bin"; d.mkdir(exist_ok=True)
        (d / "mpirun").write_text(f"#!/bin/sh\nprintf '%s\\n' \"$*\" >> {d}/calls\nshift 2\nexec \"$@\"\n"); (d / "mpirun").chmod(0o755)
        return d

    def test_the_rank_count_of_mpirun_without_a_launcher(self):   # plan "mpirun으로 부르면 그 rank 수로" (2026-10-10)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        d = self._fake_mpirun()
        code, out = call(inp, fake, dict(CCEX_AUTOCONV_NP="3", PATH=f"{d}:{os.environ['PATH']}"))
        self.assertEqual(code, 0, out)
        calls = (d / "calls").read_text().splitlines()
        self.assertTrue(calls and all(c.startswith("-n 3 ") for c in calls))
        self.assertIn("launcher: mpirun -n 3 (the mpirun that started main.out)", (inp.parent / "auto/result/summary.txt").read_text())

    def test_the_input_launcher_comes_first_and_a_bad_rank_count_is_refused(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        j = json.loads(inp.read_text()); j["auto_convergence"]["launcher"] = []; inp.write_text(json.dumps(j))
        d = self._fake_mpirun()
        code, out = call(inp, fake, dict(CCEX_AUTOCONV_NP="3", PATH=f"{d}:{os.environ['PATH']}"))
        self.assertEqual(code, 0, out)
        self.assertFalse((d / "calls").exists())                             # the input's (empty) launcher, no mpirun
        self.assertIn("launcher: none", (inp.parent / "auto/result/summary.txt").read_text())
        code, out = call(inp, fake, dict(CCEX_AUTOCONV_NP="x"))
        self.assertEqual(code, 2, out)
        self.assertIn("CCEX_AUTOCONV_NP", out)

    def test_the_production_window_is_lengthened_when_its_curve_has_not_decayed(self):   # user 2026-10-10 (P1 CR90)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_SLOW_IDS=",".join(str(k) for k in range(10, 40))))   # production 3 x slower than the test
        s = summary(inp)
        self.assertGreaterEqual(len(s["production"].get("windows", [])), 2, out)
        self.assertGreater(s["production"]["windows"][-1]["window"], s["production"]["windows"][0]["window"])
        self.assertIn("production window lengthened", (inp.parent / "auto/result/summary.txt").read_text())

    def test_nstep_is_never_changed(self):     # user 2026-10-10: nstep fixed, a given deltat sets the first window only
        inp, fake = make_case(self.tmp, threshold_config=0.5, nstep=40)        # deltat 0.02 given: the first window 0.78 ms
        code, out = call(inp, fake)
        runs = [json.loads((d / "input.json").read_text()) for step in ("pilot", "rbath", "production")
                for d in (inp.parent / "auto/runs" / step).iterdir()]
        self.assertEqual({r["nstep"] for r in runs}, {40}, out)               # the window grows by deltat, not by points
        self.assertGreater(len({r["deltat"] for r in runs}), 1)
        # a modulation on a fine nstep: the peaks (user 2026-10-10)
        d = Path(self.tmp) / "fine"; d.mkdir()
        inp, fake = make_case(d, threshold_config=0.5)
        j = json.loads(inp.read_text()); j.pop("deltat"); j["nstep"] = 3201; inp.write_text(json.dumps(j))
        code, out = call(inp, fake, dict(FAKE_MOD="0.0123:0.4"))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(s["window"]["fit_points"], "peaks")
        self.assertNotIn("half_grid", s["window"]); self.assertNotIn("half_grid", s["production"])     # no check (user 2026-10-10)
        self.assertEqual(s["production"]["nstep"], 3201)
        txt = (inp.parent / "auto/result/summary.txt").read_text()
        for word in ("fit points: the peaks", "a convergence measure", "whether they are enough is not checked"):
            self.assertIn(word, txt)
        self.assertTrue((inp.parent / "auto/result/1_window_fit.png").is_file())

    def test_the_removed_step_keys(self):          # user 2026-10-10: the fit points are the peaks
        for sec, key in ((None, "step_threshold"), ("t2_fit", "quality"), ("t2_fit", "step_limit")):
            d = Path(tempfile.mkdtemp(dir=self.tmp))
            inp, fake = make_case(d)
            j = json.loads(inp.read_text()); blk = j["auto_convergence"]
            (blk.setdefault(sec, {}) if sec else blk)[key] = 0.05
            inp.write_text(json.dumps(j))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out); self.assertIn("removed (2026-10-10)", out)
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        inp, fake = make_case(d, fit_points="envelope")
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out); self.assertIn("peaks or all", out)

    def test_review_first_pass_impl(self):                    # review 20261010-first-pass-impl R2, R8
        # R2: no rdip comparison left at all is an unresolved truncation, not an error
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        inp, fake = make_case(d, candi_rdip=[30, 80], candi_rbath=[30, 40])
        code, out = call(inp, fake, dict(FAKE_NAN_RUNG="30"))
        self.assertEqual((code, summary(inp)["status"]), (3, "truncation-unresolved"), out)
        # R8: the figure with points all given
        d = Path(tempfile.mkdtemp(dir=self.tmp))
        inp, fake = make_case(d, fit_points="all", threshold_config=0.5)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        self.assertTrue((inp.parent / "auto/result/1_window_fit.png").is_file())
        self.assertIn("fit points: all (", (inp.parent / "auto/result/summary.txt").read_text())
        self.assertIn("t2_fit.points all", (inp.parent / "auto/result/summary.txt").read_text())

    def test_a_curve_without_peaks_takes_all_points(self):    # user 2026-10-10: fewer than 3 peaks
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake)
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(s["window"]["fit_points"], "all")
        self.assertIn("fewer than 3 peaks", (inp.parent / "auto/result/summary.txt").read_text())
        self.assertTrue((inp.parent / "auto/result/1_window_fit.png").is_file())
        d = Path(self.tmp) / "noisy"; d.mkdir()
        inp, fake = make_case(d, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_FLOOR="0.3"))
        self.assertEqual((code, summary(inp)["status"]), (0, "converged"), out)

    def test_the_decay_ends_at_the_noise_floor_of_single_states(self):   # user 2026-10-10 (P1 CR90 single)
        inp, fake = make_case(self.tmp, threshold_config=0.5, fit_points="all")
        code, out = call(inp, fake, dict(FAKE_FLOOR="0.3"))
        s = summary(inp)
        self.assertNotEqual((s["status"], s.get("stage")), ("window-unresolved", "window"), out)
        self.assertGreater(s["window"]["noise_level"], 0.01)
        self.assertGreater(s["window"]["fit_floor"], 0.05 if s["window"]["noise_level"] > 0.05 else 0.0)

    def test_a_divergent_point_is_corrected_and_the_configuration_kept(self):    # plan "발산한 점 보정" (2026-10-10)
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_JUMP_IDS="10"))
        self.assertEqual(code, 0, out)
        s = summary(inp)
        self.assertEqual(s["unstable"], [])
        self.assertEqual([(r["configuration"], r["max_wiDiv"]) for r in s["corrected"]], [(10, 5.0)])
        self.assertIn("divergent points corrected", (inp.parent / "auto/result/summary.txt").read_text())
        self.assertEqual(s["production"]["n_used"], len(s["production"]["configurations"]))     # kept

    def test_a_divergent_point_without_nodiv_leaves_the_configuration_out(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        j = json.loads(inp.read_text()); j["savemode"] = {"mode": "coherence", "nodiv": False}; inp.write_text(json.dumps(j))
        code, out = call(inp, fake, dict(FAKE_BIG_IDS="10"))
        self.assertEqual(code, 0, out)
        self.assertEqual([r["configuration"] for r in summary(inp)["unstable"]], [10])

    def test_nan_in_production_is_left_out_too(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_NAN_IDS="10"))
        self.assertEqual(code, 0, out)
        self.assertEqual(summary(inp)["unstable"][0]["reason"], "nan/inf")

    def test_too_many_unstable_is_a_failure(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)                # 20 computed, 2 unstable = 10 % > 5 %
        code, out = call(inp, fake, dict(FAKE_BIG_IDS="10,11"))
        self.assertEqual(code, 3, out)
        s = summary(inp)
        self.assertEqual(s["status"], "numerical-failure")
        self.assertIn("2 of 20 production configurations unstable", s["reason"])

    def test_an_unstable_test_configuration_is_left_out(self):
        inp, fake = make_case(self.tmp, threshold_config=0.5)
        code, out = call(inp, fake, dict(FAKE_BIG_IDS="2"))
        self.assertEqual(code, 0, out)
        self.assertEqual([(r["configuration"], r["where"]) for r in summary(inp)["unstable"]], [(2, "test")])

    def test_crash_mid_output_then_resume(self):
        inp, fake = make_case(self.tmp)
        marker = Path(self.tmp) / "failed_once"
        code, out = call(inp, fake, dict(FAKE_FAIL_ONCE=str(marker)))
        self.assertEqual(code, 1, out)
        code, out = call(inp, fake, dict(FAKE_FAIL_ONCE=str(marker)))
        self.assertEqual(code, 0, out)

    def test_random_fallback_is_a_failure(self):
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(FAKE_RANDOM="1"))
        self.assertEqual(code, 1, out)
        self.assertIn("random", out)

    def test_subbath_fallback_only_with_subbath_spins(self):    # P1 inputs print it with "No sub-bath"
        msg = "The subbath states are randomly generated"
        inp, fake = make_case(self.tmp)
        code, out = call(inp, fake, dict(FAKE_RANDOM=msg))
        self.assertEqual(code, 0, out)
        (Path(self.tmp) / "b").mkdir()
        inp2, _ = make_case(Path(self.tmp) / "b")
        j = json.loads(inp2.read_text()); j["Defect"] = dict(species=[dict(name="P1", n_nuclear_spins=1)]); inp2.write_text(json.dumps(j))
        code, out = call(inp2, fake, dict(FAKE_RANDOM=msg))
        self.assertEqual(code, 1, out)

    def test_a_short_state_file_is_refused(self):
        inp, fake = make_case(self.tmp)
        f = inp.parent / "data/cfg1_state1"
        f.write_text("\n".join(f.read_text().splitlines()[:-1]) + "\n")       # one spin without a state
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("states per state file", out)

    def test_a_zero_linearized_se_gives_no_allocation(self):   # review 10-07 B3
        inp, fake = make_case(self.tmp, fit_points="all")       # the zero linearized SE, not the fit points, is what this checks
        code, out = call(inp, fake, dict(FAKE_IDENTICAL="1"))
        s = summary(inp)
        self.assertEqual((code, s["status"], s["stage"]), (3, "sampling-unresolved", "sampling"), out)
        self.assertIn("linearized SE", s["reason"])

    def test_inputs_outside_the_checks_are_refused(self):   # review 10-07 A10
        cases = [(dict(bath_shift=[[1.0, 0.0, 0.0]]), {}, "bath_shift"),
                 (dict(Qubit=dict(qubit=[dict(name="q1"), dict(name="q2")])), {}, "one qubit"),
                 ({}, dict(candi_rbath=[40]), "at least two"),
                 ({}, dict(window_method="fit"), "removed"),
                 ({}, dict(max_nstep=30), "removed"),
                 ({}, dict(states_available=30), "removed"),
                 ({}, dict(recheck=3), "removed"),
                 ({}, dict(default_nstep=40), "nstep of the input"),
                 ({}, dict(fit_floor=1.5), "t2_fit.floor"),
                 ({}, dict(fit_n_range=[2.0, 1.0]), "t2_fit.n_range"),
                 ({}, dict(t2_fit_points="all"), "now t2_fit.points"),
                 ({}, dict(threshold_config=2.0), "configurations.threshold"),
                 ({}, dict(se_target=0.02), "now configurations.threshold"),
                 ({}, dict(truncation_rel_tol=0.02), "now rdip.threshold and rbath.threshold"),
                 ({}, dict(max_spins=10), "removed"),
                 ({}, dict(rdip_max=50), "now rdip.candidates"),
                 ({}, dict(time_window=2.0), "time_window must be an object"),
                 ({}, dict(converge_time_window=False, set_time_window=2.0, set_time_window_unit="min"), "time_window.unit"),
                 ({}, dict(smoothing_points=9), "removed"),
                 ({}, dict(enabled=False), "enabled is false"),
                 ({}, dict(converge_rdip=False), "give rdip.value"),
                 ({}, dict(value_rdip=30.0), "set rdip.converge false"),
                 ({}, dict(converge_time_window=False, converge_rdip=False, converge_rbath=False, converge_config=False,
                           set_time_window=2.0, value_rdip=30.0, value_rbath=30.0, number_config=3, states_config=2), "nothing to converge"),
                 (dict(deltat=0), {}, "deltat")]
        for i, (base, block, msg) in enumerate(cases):
            d = Path(self.tmp) / f"c{i}"; d.mkdir()
            inp, fake = make_case(d, **block)
            j = json.loads(inp.read_text()); j.update(base); inp.write_text(json.dumps(j))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, (msg, out))
            self.assertIn(msg, out)

    def test_lock_and_mpi_context_refused(self):
        inp, fake = make_case(self.tmp)
        (inp.parent / "auto").mkdir()
        with open(inp.parent / "auto/lock", "w") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            code, out = call(inp, fake)
            self.assertEqual(code, 2, out)
            self.assertIn("another runner", out)
        code, out = call(inp, fake, dict(PMI_RANK="0"))
        self.assertEqual(code, 2, out)

    def test_input_refusals(self):
        cases = [dict(configurations=None), dict(candi_rdip=[10, 5]), dict(version=1), dict(unknown_key=1),
                 dict(time_window=dict(method="fit")), dict(time_window=dict(smoothing_points=9))]
        for i, bad in enumerate(cases):
            d = Path(self.tmp) / f"c{i}"; d.mkdir()
            inp, fake = make_case(d)
            j = json.loads(inp.read_text())
            for k, v in bad.items():
                if v is None: j["auto_convergence"].pop(k)
                else: j["auto_convergence"][k] = v
            inp.write_text(json.dumps(j))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, (bad, out))
        d = Path(self.tmp) / "overlap"; d.mkdir()
        inp, fake = make_case(d, pilot=(1, 2), pool=(2, 3))
        self.assertEqual(call(inp, fake)[0], 2)
        d = Path(self.tmp) / "missing"; d.mkdir()
        inp, fake = make_case(d)
        (d / "data/cfg3_state2").unlink()
        code, out = call(inp, fake)
        self.assertEqual(code, 2, out)
        self.assertIn("missing", out)
        for i, (edit, msg) in enumerate(((dict(method="DM"), "supports method"),
                                         (dict(cluster_mode="pcce", pcce=dict(partition="legacy")), "pcce.partition"),
                                         (dict(adaptive=dict(policy="merge")), "adaptive"), (dict(pcce=dict(cache=dict(file="x.cache"))), "pcce.cache"),
                                         (dict(cluster_mode="pcce", pcce=dict(cache=dict(mode="readwrite"))), "own partition cache"))):
            d = Path(self.tmp) / f"method{i}"; d.mkdir()
            inp, fake = make_case(d)
            j = json.loads(inp.read_text()); j.update(edit); inp.write_text(json.dumps(j))
            code, out = call(inp, fake)
            self.assertEqual(code, 2, (edit, out))
            self.assertIn(msg, out)

    def test_other_methods_run_every_rung_and_find_the_same_answer(self):   # plan "지원 방법 넓히기" M1, M3; A23
        for i, edit in enumerate((dict(method="gCCE"), dict(cluster_mode="pcce"), dict(order=3), dict(cluster_mode="bundle"))):
            d = Path(self.tmp) / f"m{i}"; d.mkdir()
            inp, fake = make_case(d)
            j = json.loads(inp.read_text()); j.update(edit); inp.write_text(json.dumps(j))
            code, out = call(inp, fake)
            self.assertEqual(code, 0, (edit, out))
            s = summary(inp)
            self.assertEqual((s["truncation"]["chosen"]["rbath"], s["truncation"]["chosen"]["rdip"]), (30.0, 30.0), edit)
            runs = [json.loads(p.read_text()) for p in (inp.parent / "auto/runs").rglob("input.json")]
            if "method" in edit or edit.get("cluster_mode") == "pcce":   # gCCE2 and pCCE2 have the ladder too (A23)
                self.assertTrue(any("rdip_ladder" in r for r in runs), edit)
                self.assertNotIn(20.0, {r["rdip"] for r in runs})
            else:                                   # order 3 and bundles: one run per rung
                self.assertFalse(any("rdip_ladder" in r for r in runs), edit)
                self.assertIn(20.0, {r["rdip"] for r in runs})

    def test_pcce_runs_share_a_partition_cache(self):     # A17 (1): runs differing only in rdip (or states) reuse the partition
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["cluster_mode"] = "pcce"; inp.write_text(json.dumps(j))
        cnt = Path(self.tmp) / "count"
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        runs = [(json.loads((p.parent / "input.json").read_text()), json.loads(p.read_text())) for p in (inp.parent / "auto/runs").rglob("done.json")]
        files = {}
        for cin, done in runs:
            self.assertEqual(cin["pcce"]["cache"]["mode"], "readwrite")
            self.assertEqual(cin["seed"], 1)                                # the case has no seed: 1, not the start time
            self.assertEqual(cin["pcce"]["cache"]["file"], "../../../pcce_cache/" + done["pcce_cache"]["file"])
            files.setdefault(done["pcce_cache"]["file"], set()).add((cin["bathfile"][0], cin["rbath"]))
            self.assertTrue((inp.parent / "auto/pcce_cache" / done["pcce_cache"]["file"]).is_file())
        self.assertTrue(all(len(v) == 1 for v in files.values()), files)      # one file per configuration and rbath
        self.assertEqual(len(files), len({x for v in files.values() for x in v}))
        status = [done["pcce_cache"]["status"] for _, done in runs]
        self.assertEqual(status.count("stored"), len(files))                # searched once per file, read by every other run
        self.assertGreater(status.count("hit"), 0)
        # with the ladder (A23) fewer runs differ only in rdip, but a raised top rdip still reads the partition it stored
        rdips = {}
        for cin, done in runs:
            rdips.setdefault(done["pcce_cache"]["file"], set()).add(cin["rdip"])
        self.assertTrue(any(len(v) > 1 for v in rdips.values()), rdips)
        for f in files:                                                     # the hash recorded when the file was stored
            c = inp.parent / "auto/pcce_cache" / f
            self.assertEqual(Path(f"{c}.sha256").read_text().strip(), hashlib.sha256(c.read_bytes()).hexdigest())
        n1 = count(cnt)
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))              # a restart reuses every run
        self.assertEqual((code, count(cnt)), (0, n1), out)
        c = next((inp.parent / "auto/pcce_cache").glob("*.cache"))
        c.write_text(c.read_text())                                         # the same bytes: still read
        stale = next((inp.parent / "auto/pcce_cache").glob("*.sha256"))
        stale.write_text("0" * 64 + "\n")                                  # a file whose hash is not the one stored is not read
        shutil.rmtree(inp.parent / "auto/runs")
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        stored = [json.loads(p.read_text())["pcce_cache"] for p in (inp.parent / "auto/runs").rglob("done.json")]
        self.assertEqual(sorted(x["file"] for x in stored if x["status"] == "stored"), [stale.name[:-len(".sha256")]])

    def test_partition_key_holds_what_the_partition_depends_on(self):     # A17 (1)
        from autoconv.runner import Runner
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["cluster_mode"] = "pcce"; inp.write_text(json.dumps(j))
        r = Runner(inp, fake)
        data = inp.parent / "data"
        cin = dict(method="cCCE", cluster_mode="pcce", order=2, pcce=dict(partition_size=2, trials=100), seed=1,
                   bathfile=[str(data / "bath_1")], rbath=30, rdip=20, nstate=4, bath_state_file=str(data / "cfg1_state"),
                   nstep=99, deltat=0.02, outfile="out")
        fs = dict(bathfile=data / "bath_1", state1=data / "cfg1_state1")
        extra = {"Qubit.file": str(data / "qubit")}
        key = r.partition_key(cin, fs, extra)
        for same in (dict(rdip=40), dict(nstate=8), dict(nstep=150), dict(deltat=0.05), dict(Defect=dict(state_file="x"))):
            self.assertEqual(r.partition_key(dict(cin, **same), fs, extra), key, same)
        self.assertEqual(r.partition_key(cin, dict(fs, state1=data / "cfg1_state2"), extra), key)
        for other in (dict(rbath=40), dict(pcce=dict(partition_size=3, trials=100)), dict(seed=2),
                      dict(pcce=dict(partition_size=2, trials=200)), dict(rbath=[5, 30]),
                      dict(pcce=dict(partition_size=2, trials=100, kmeans_pp=False))):
            self.assertNotEqual(r.partition_key(dict(cin, **other), fs, extra), key, other)
        self.assertNotEqual(r.partition_key(cin, dict(fs, bathfile=data / "bath_2"), extra), key)
        (data / "bath_1").write_text((data / "bath_1").read_text() + "\n")       # the same name, other contents
        self.assertNotEqual(r.partition_key(cin, fs, extra), key)
        key = r.partition_key(cin, fs, extra)
        (data / "qubit").write_text("1 2 3\n")                                 # the qubit moved: other spins inside rbath
        self.assertNotEqual(r.partition_key(cin, fs, extra), key)
        r.binary_sha = "0" * 64
        self.assertNotEqual(r.partition_key(cin, fs, extra), key)

    def test_partition_cache_rule_is_in_the_plan_fingerprint(self):     # A17 (1)
        import autoconv.runner as runner
        for shape, moves in (("pcce", True), ("standard", False)):
            d = Path(self.tmp) / shape; d.mkdir()
            inp, fake = make_case(d)
            j = json.loads(inp.read_text()); j["cluster_mode"] = shape; inp.write_text(json.dumps(j))
            r = runner.Runner(inp, fake)
            r.resolve(); r.preflight()
            fp = r.plan()
            rule = runner.PARTITION_RULE
            try:
                runner.PARTITION_RULE = dict(rule, version=rule["version"] + 1)
                self.assertEqual(r.plan() != fp, moves, shape)
            finally:
                runner.PARTITION_RULE = rule

    def test_reused_rung_runs_count_their_cost(self):   # review 4d2d3a R2
        inp, fake = make_case(self.tmp)
        j = json.loads(inp.read_text()); j["method"] = "gCCE"; inp.write_text(json.dumps(j))
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="1"))
        self.assertEqual(code, 0, out)
        first = summary(inp)["core_hours"]
        total = sum(json.loads(p.read_text())["core_hours"] for p in (inp.parent / "auto/runs").rglob("done.json"))
        self.assertAlmostEqual(first, total, delta=1e-12)               # every run once, rungs included
        code, out = call(inp, fake, dict(OMP_NUM_THREADS="1"))          # all reused: the same total
        self.assertAlmostEqual(summary(inp)["core_hours"], total, delta=1e-12)

    def test_tensor_files_are_hashed(self):          # plan M2
        inp, fake = make_case(self.tmp)
        tf = inp.parent / "data/hf.dat"
        tf.write_text("A tensor\n")
        j = json.loads(inp.read_text())
        j["hamiltonian"] = {"common": {"hyperfine": {"dipole": "dft", "contact": True, "file": "data/hf.dat"}}}
        inp.write_text(json.dumps(j))
        cnt = Path(self.tmp) / "count"
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        run = json.loads(next((inp.parent / "auto/runs/production").glob("*/input.json")).read_text())
        self.assertEqual(run["hamiltonian"]["common"]["hyperfine"]["file"], str(tf.resolve()))
        n1 = count(cnt)
        tf.write_text("another tensor\n")                  # the same name, other contents: computed again, a new plan
        code, out = call(inp, fake, dict(FAKE_COUNT=str(cnt)))
        self.assertEqual(code, 0, out)
        self.assertEqual(count(cnt), 2 * n1)
        self.assertIn("decisions.1.jsonl", out)


if __name__ == "__main__":
    unittest.main()
