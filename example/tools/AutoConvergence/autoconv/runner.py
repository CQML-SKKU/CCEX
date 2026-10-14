"""`main.out --auto_conv -f input.json`: the runner ((internal notes), "--auto_conv 첫 판 설계 v2").

    python -m autoconv run INPUT --binary /abs/path/main.out          (main.out --auto_conv calls this; not under mpirun)

The input is a CCEX input with an `auto_convergence` block (v1; example: examples/nv13c/input.json). Every calculation is an
ordinary CCEX run of a child input (the block removed) in its own folder, launched with the block's `launcher`. Steps:
window (pilot at the largest rbath and rdip with the rdip ladder on nstep points, the window grown until it resolves, the fit
points chosen), rdip
then rbath (converge(): relative area of the T2 fit, rdip judged again at the final rbath), sampling allocation of the area
from the pilot, production generously first (1.5 x the pilot's estimate), taken apart, more when short, production window and combination recheck, result (README.md). Decisions follow the current contract of the plan (window.py, t2.py,
truncation.py). A run is reused on restart only when its done.json and output hashes still match.

Exit codes: 0 converged, 3 unresolved (window, truncation or sampling; see result/summary.txt), 2 input refused,
4 stopped by a resource limit (call again to continue), 1 internal error.
"""
import copy
import fcntl
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from . import autostart
from . import report as report_mod
from . import t2 as t2mod
from . import truncation, window
from .config import KeyCaseError, read_input
from .curves import read_curve

VERSION = 1
METHODS = ("cce", "gcce", "pcce")    # Adaptive excluded (plan "지원 방법 넓히기", 2026-10-08); pcce: cCCE with cluster_mode pcce
# the input files besides the bath, as paths into the input's blocks (input keys of A47, 2026-10-13): resolved to absolute
# paths in every child input, and hashed into every run key and the plan
INPUT_FILES = (("Qubit", "file"), ("gyrofile",), ("hamiltonian", "common", "hyperfine", "file"),
               ("hamiltonian", "common", "quadrupole", "file"), ("hamiltonian", "common", "quadrupole", "file_without_defect"))
MADE_STATES = 10 ** 6     # states the runner can make (no files): no limit of their own; max_states bounds the state test (review 01ecda #2)
DEFAULTS = dict(pilot_states=10, test_number=10, test_states=5, max_configurations=None, max_states=30,
                threshold_time_window=0.01, production_window_factor=1.5,
                fit_points="peaks", fit_floor=0.05, fit_n_range=[0.3, 4.0],
                revival_check_factor=2.5,
                set_time_window=None, set_time_window_unit="ms", num_candi_rdip=7, num_candi_rbath=4,
                threshold_rdip=0.02, threshold_rbath=0.02, threshold_config=0.05, states_threshold=0.05, min_configurations=10,
                launcher=[])
KEYS = set(DEFAULTS) | {"workdir", "configurations", "candi_rdip", "candi_rbath", "value_rdip", "value_rbath",
                        "number_config", "states_config"}
# the block is nested by step (user 2026-10-08): section -> {key -> the flat name the runner uses}
SECTIONS = {
    "time_window": dict(converge="converge_time_window", value="set_time_window", unit="set_time_window_unit",
                        threshold="threshold_time_window", revival_check_factor="revival_check_factor",
                        production_factor="production_window_factor"),
    "rdip": dict(converge="converge_rdip", value="value_rdip", candidates="candi_rdip", num_candidates="num_candi_rdip",
                 threshold="threshold_rdip"),
    "rbath": dict(converge="converge_rbath", value="value_rbath", candidates="candi_rbath", num_candidates="num_candi_rbath",
                  threshold="threshold_rbath"),
    "t2_fit": dict(points="fit_points", floor="fit_floor", n_range="fit_n_range"),
    "configurations": dict(converge="converge_config", number="number_config", states="states_config", threshold="threshold_config",
                           test_number="test_number", test_states="test_states", states_threshold="states_threshold",
                           max_configurations="max_configurations", max_states="max_states",
                           min_configurations="min_configurations"),
}
FILE_KEYS = ("bathfile", "statefile", "avaaxfile", "exstatefile")
SPLIT_KEYS = dict(test="pilot", production="pool")   # configurations.test / .production (user 2026-10-08) -> the runner's names
TOP_KEYS = {"enabled", "resume", "workdir", "launcher"} | set(SECTIONS)
STEPS = ("0 start", "1 time window", "2 rdip", "3 rbath", "4 configurations", "5 production")
WHERE = {flat: f"{sec}.{key}" for sec, m in SECTIONS.items() for key, flat in m.items()}
TIME_UNITS = dict(s=1e3, ms=1.0, us=1e-3, ns=1e-6)        # to ms, the CCEX time unit
RANDOM_FALLBACK = ("bath state is randomly generated", "paxes is randomly generated")
SUBBATH_FALLBACK = "subbath states are randomly generated"   # a failure only when a defect has sub-bath spins (review 10-07 A7;
                                                             # CCEX prints it also with "No sub-bath", e.g. the P1 inputs)
MAX_GROWS = 8            # longer windows when deltat is found (1.5^8 ~ 26 x the start)
MAX_WINDOW_ROUNDS = 3    # the window found again at a grown reference (plan "자동 시작값" (8c))
COARSE_POINTS = 10       # deltat found: the mean below the fit floor (t2_fit.floor) within this many points is a grid too coarse
MAX_COARSE = 3           # cuts of the window for a grid too coarse
WALLTIME_MARGIN_S = 300.0     # a run is not started with less Slurm walltime left than 1.5 x its longest + this
DEFAULT_NSTEP = 150       # time points when the input has no nstep (101 until 2026-10-10: too few for the envelope step)
SMOOTHING_POINTS = 9      # running median of the window's plateau and revival checks
# the production (plan "넉넉하게 먼저, 나눠 보기, 모자라면 추가", user 2026-10-09; review-20261009-generous-production)
C_FLOOR = 20             # the least passed (with min_configurations) and the least added
MARGIN = 1.5             # the first block and every addition: this many times the configurations estimated
ROUND_TO = 5             # planned counts rounded up to a multiple of this (the pool limit wins)
STATE_STEP = 5           # the state test adds this many states at a time (user 2026-10-10)
STATE_RULE = "m_s uniform over the 2S+1 sublevels, numpy default_rng([seed, k, s]), one draw per bath-file line"
AVAAX_RULE = "axis uniform over 1..navaax of the defect species, numpy default_rng([seed, k, 0]), one draw per bath-file line, 0 elsewhere"
POINT_JUMP = 2.0          # |wiDiv - noDiv| at or above this: the point takes noDiv (the analyzer's error_correct, 2026-10-10)
UNSTABLE_L = 10.0         # |L| above this (or nan/inf) anywhere: the configuration is unstable and left out (user 2026-10-10)
UNSTABLE_SHARE = 0.05     # more than this share of the production unstable: a numerical failure (leaving them out could bias)
MAX_PROD_WINDOWS = 3     # a production curve not decayed inside its window: computed again on a longer one at most this often
SPLIT_RULE = "nearest-spin quantiles"   # the test configurations evenly over the nearest bath spin's distance (2026-10-10)
RECHECK = 3              # production configurations of the final combination check (fixed, user 2026-10-08)
REPLACED = dict(step_threshold="removed (2026-10-10): the fit points are the peaks", quality="removed (2026-10-10): the fit points are the peaks",
                step_limit="removed (2026-10-10): the fit points are the peaks",
                max_points="removed (2026-10-10): nstep is fixed, the grid is never made finer",
                recheck="removed (2026-10-08): the combination check always uses 3 production configurations",
                max_nstep="removed (2026-10-08): the window grows at most 8 times in every mode",
                states_available="removed (2026-10-08): the state files of each configuration are counted",
                version="removed (2026-10-08): not needed", source_radius="removed (2026-10-08): rbath is bounded by the bath files' box",
                trial_nstep="removed: nstep of the input (default 150)", default_nstep="removed (2026-10-08): nstep of the input (default 150)", decay_threshold="now time_window.threshold",
                revival_look="now time_window.revival_check_factor", margin="now time_window.production_factor",
                window_method="removed (2026-10-08): the decay time always comes from the T2 fit", t2_fit_points="now t2_fit.points",
                window_median_points="removed (2026-10-08): fixed at 9 points", smoothing_points="removed (2026-10-08): fixed at 9 points",
                method="removed (2026-10-08): the decay time always comes from the T2 fit",
                rdip_rungs="now rdip.candidates (the rungs and the largest value in one list)",
                rdip_max="now rdip.candidates (its largest value is the reference)",
                rbath_candidates="now rbath.candidates", truncation_max_tol="removed (it did not take part in the verdict)",
                truncation_area_tol="removed (it did not take part in the verdict)", cost_config="removed (2026-10-08)", cost_state="removed (2026-10-08)", margin_s="removed (2026-10-08)",
                run_timeout_s="removed (2026-10-08)", truncation_rel_tol="now rdip.threshold and rbath.threshold",
                se_target="now configurations.threshold", max_core_hours="removed (2026-10-08)",
                max_spins="removed (2026-10-08)", max_pairs="removed (2026-10-08)",
                bath_generation="removed (2026-10-09): make the bath files yourself; the run says which ones to add")
CODE = ("runner.py", "window.py", "t2.py", "truncation.py", "curves.py", "autostart.py")      # the decision code, in the plan fingerprint
# pCCE partition cache (A17 (1), user 2026-10-11): every pCCE run gets pcce.cache.mode readwrite and a file named by a key of
# its child input without the keys the partition does not read, the contents of the files it reads but the state files, and
# the binary - so runs differing only in rdip (or in states, nstep, deltat) share one partition. CCEX compares the stored
# partition's own inputs (bath spins, seed, sK, trials, iterations, kmeans_pp, build) byte for byte and searches again on
# any difference (pcce_partition_cache.cpp), so the key only keeps different partitions in different files.
PARTITION_FREE = ("rdip", "rdip_ladder", "nstate", "bath_state_file", "Defect.state_file", "nstep", "deltat")
PARTITION_RULE = dict(mode="readwrite", key="child input without " + ", ".join(PARTITION_FREE) +
                      "; contents of the bath and other files but the states; binary",
                      seed="1 when the input has none (CCEX would take the start time, and the cache holds the seed)", version=2)


class Refused(Exception):
    """The input is outside what the runner supports (exit 2)."""


class Stop(Exception):
    """A resource limit stops the runner before the next run (exit 4); calling again continues."""


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _at(d, *path, default=None):
    """d[path[0]][path[1]]...: the value at a path of keys of the input's blocks, or default."""
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


def _radius(v):
    """rbath or rdip of an input: a number (the outer radius) or [inner, outer] (A47: the inner one was rbathcut, rdipcut)."""
    if isinstance(v, list) and len(v) == 2:
        return v[0], v[1]
    return 0.0, v


def _method(base):
    """The method as the runner names it: cce, gcce or pcce (cCCE with cluster_mode pcce); the input's method and
    cluster_mode otherwise (cluster_mode bundle keeps the method's name)."""
    m = str(base.get("method", "cCCE")).lower()
    m = "cce" if m == "ccce" else m
    shape = str(base.get("cluster_mode", "standard")).lower()
    if shape == "pcce" and m == "cce":
        return "pcce"
    return m if shape in ("standard", "bundle") else f"{m} with cluster_mode {shape}"


def _nodiv(cin):
    """savemode.nodiv of a CCEX input (default true): false, CCEX writes _wiDiv alone (docs/user/InputFile_Output.md)."""
    sm = cin.get("savemode")
    return not (isinstance(sm, dict) and sm.get("nodiv", True) is False)


def _canon(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"))


def _atomic_write(path, text):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def _json(res):
    """result/json/: the files the program reads (summary.json, progress.json), apart from what people read (user 2026-10-10)."""
    d = Path(res) / "json"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _gname(r):
    """CCEX names a rung with %g of the float value (rdip_ladder.cpp)."""
    return f"{float(np.float32(r)):g}"


def unjudged_reason(name, rec, figure):
    """Why an unresolved rdip or rbath stopped, when its last comparison could not be judged (2026-10-12, job 17163): the
    curves were not finite (or a failed T2 fit left no usable raw area), so whether more would converge is unknown and larger
    values are not advised. None when the last comparison was judged (a real "fail", the raw-area fallback of a failed fit
    included: the usual reason stands)."""
    checks = (rec or {}).get("checks") or []
    if not checks or checks[-1].get("verdict") not in ("fit-failure", "nonfinite"):
        return None
    bad = [c for c in checks if c.get("verdict") in ("fit-failure", "nonfinite")]
    pairs = ", ".join(f"{c['lo']:g} vs {c['hi']:g} ({c['verdict']})" for c in bad)
    return (f"{name} could not be judged at its largest values: the T2 fit failed or the curves are not finite for {pairs}; "
            f"look at {figure} (noisy or unstable curves) before computing more")


def raw_area_pairs(rec):
    """(lo, hi, verdict) of every comparison in a converge record judged on the raw area because a T2 fit failed (user
    2026-10-12)."""
    rec = rec or {}
    est = rec.get("estimate") or {}
    ref = rec.get("reference")
    rows = [(x["lo"], x["hi"], x) for x in rec.get("checks", []) + (est.get("neighbours") or []) + ([est["check"]] if est.get("check") else [])]
    rows += [(x["value"], ref, x) for x in rec.get("vs_reference", [])]
    return {(lo, hi, x["verdict"]) for lo, hi, x in rows if x.get("basis") == "raw-area"}


def converge(values, curve, extend, limit, clip, judge, rel_tol, need=1):
    """One truncation parameter (user decisions 2026-10-07). judge(lo_curves, hi_curves) -> dict(verdict, rel). values: the candidates computed so far (sorted);
    curve(v): configurations x time of state means; extend(v): compute a new largest value; limit: the largest allowed.
    1. Trend: with three values or more whose neighbouring changes (relative area of the T2 fit) shrink by q <= 0.8
       (every neighbouring comparison judged), each value's distance from the converged one is the changes above it plus
       last x q/(1-q); when some value is within the tolerance, ONE larger value (1.25 x the largest) is computed.
    2. Then, trend or not (review 10-07 A1: the trend only proposes the next value), the largest is raised by 1.25 (up
       to limit) until the last two agree, every value is compared with the largest, and the smallest from which all
       larger ones pass is chosen. At the limit: not converged.
    need: consecutive passing neighbour comparisons at the top (2 for a start found automatically, plan "자동 시작값" (7)).
    clip: a step cut short by the limit still counts (rdip: 2 x rbath holds every pair, the exact value); without it
    (rbath) only a full 1.25 step is taken, a shorter one would agree too easily."""
    values = list(values)
    rec = dict(estimate=None, raised=[], checks=[], chosen=None)
    pair = lambda lo, hi: dict(lo=lo, hi=hi, **judge(curve(lo), curve(hi)))
    def bigger(top):
        nxt = round(top * 1.25, 3)
        if nxt > limit * (1 + 1e-9):
            nxt = round(limit, 3) if clip else top
        return nxt
    if len(values) >= 3:
        nb = [pair(lo, hi) for lo, hi in zip(values, values[1:])]
        est = dict(neighbours=nb, status="no trend")
        # review 10-07 A2: a failed comparison gives no trend; nor does one judged on the raw area (a fit failed: not the
        # quantity the trend extrapolates, user 2026-10-12 changed only the verdict)
        if any(x["rel"] is None or x.get("basis") == "raw-area" for x in nb):
            est["status"] = "no trend (a comparison failed: " + ", ".join(sorted(
                {"raw area" if x.get("basis") == "raw-area" else x["verdict"] for x in nb if x["rel"] is None or x.get("basis") == "raw-area"})) + ")"
        elif nb[-2]["rel"] > 0 and nb[-1]["rel"] / nb[-2]["rel"] <= 0.8:
            q = nb[-1]["rel"] / nb[-2]["rel"]
            tail = nb[-1]["rel"] * q / (1 - q)
            est.update(q=q, tail=tail, estimates=[dict(value=v, rel=sum(x["rel"] for x in nb[i:]) + tail)
                                                  for i, v in enumerate(values)])
            ok = [e["rel"] <= rel_tol for e in est["estimates"]]
            k = len(ok)
            while k > 0 and ok[k - 1]:
                k -= 1
            nxt = bigger(values[-1])
            if k == len(values):
                est["status"] = "no value within the tolerance by the estimate"
            elif nxt <= values[-1] * (1 + 1e-9):
                est["status"] = "no larger value allowed for the check"
            else:
                extend(nxt); values.append(nxt)
                chk = pair(values[k], nxt)
                est.update(candidate=values[k], check=chk, status="accepted" if chk["verdict"] == "pass" else "check failed")
        rec["estimate"] = est
    settled = False
    while len(values) > 1:
        last = [pair(values[-2 - i], values[-1 - i]) for i in range(min(need, len(values) - 1))]
        rec["checks"].append(last[0])
        if len(last) == need and all(x["verdict"] == "pass" for x in last):
            settled = True
            break
        nxt = bigger(values[-1])
        if nxt <= values[-1] * (1 + 1e-9):
            break
        extend(nxt); values.append(nxt); rec["raised"].append(nxt)
    rec["values"], rec["reference"] = values, values[-1]
    rec["vs_reference"] = [dict(value=v, **judge(curve(v), curve(values[-1]))) for v in values[:-1]]
    rec["need"] = need
    if settled:
        rec["chosen"] = suffix(values, rec["vs_reference"])
    elif rec["checks"] and rec["checks"][-1]["verdict"] == "pass":
        # not confirmed (fewer consecutive passes than needed): the smallest passing suffix against the unconfirmed largest
        # value is reported as advice only - chosen stays None and the run stays unresolved (review
        # 20261008-autoconv-confirmation-rule: keep need=2, report the advisory)
        rec["advisory"] = suffix(values, rec["vs_reference"])
    return rec

def write_states(prefix, bathfile, k, nstate, seed, spins, first=1):
    """State files <prefix>1..<prefix>nstate for configuration k (plan "state 파일 자동", user 2026-10-09): the spin count, then
    for every line of the bath file m_s drawn evenly from its species' sublevels with default_rng([seed, k, s]), in the file's
    order - a spin keeps its state whatever rbath or rdip a run uses, so runs of one configuration stay paired. first > 1:
    states first ... first + nstate - 1 written as <prefix>1 ... (the state test adds states a block at a time, 2026-10-10)."""
    names = [r[3] if len(r) > 3 else "" for r in (line.split() for line in Path(bathfile).read_text().splitlines()[1:]) if r]
    S = np.array([spins[n] for n in names])
    for i, st in enumerate(range(first, first + nstate), 1):
        u = np.random.default_rng([int(seed), int(k), st]).random(len(names))
        ms = S - np.floor(u * (2 * S + 1))
        Path(f"{prefix}{i}").write_text(f"{len(names)}\n" + "".join(f"{v:.5f}\n" for v in ms))


def write_avaax(path, bathfile, k, seed, navaax):
    """The avaax file of configuration k (plan "avaax 파일 자동", user 2026-10-09): the line count, then for every line of the
    bath file an axis drawn evenly from 1..navaax of its species (a defect of the input) with default_rng([seed, k, 0]), 0 for
    other species - one per configuration, the same at every rbath, rdip and state."""
    names = [r[3] if len(r) > 3 else "" for r in (line.split() for line in Path(bathfile).read_text().splitlines()[1:]) if r]
    n = np.array([navaax.get(x.lower(), 0) for x in names])      # names compared ignoring case, as CCEX does
    ax = np.where(n > 0, 1 + np.floor(np.random.default_rng([int(seed), int(k), 0]).random(len(names)) * n), 0).astype(int)
    Path(path).write_text(f"{len(names) + 1}\n" + "".join(f"{v}\n" for v in ax))


def round_up(x):
    return ROUND_TO * math.ceil(x / ROUND_TO - 1e-9)


def grow_counts(total, first):
    """The counts the computed sample is taken apart at (user 2026-10-10): first, then every 5 (5, 10, ...) up to total; every
    count when total is 10 or fewer."""
    if total <= 10:
        return list(range(first, total + 1))
    return sorted({first, *range(5, total, 5), total})


def first_target(n_pilot, cmin, pool):
    """Configurations of the first production block: MARGIN x the pilot's estimate, at least cmin, rounded up, at most the pool."""
    return min(pool, round_up(max(cmin, MARGIN * n_pilot)))


def next_target(c, r, tau, cmin, pool):
    """After c configurations measured at SE r against the target tau: (estimated need, wanted total, next total within the pool)."""
    need = max(cmin, math.ceil(c * (r / tau) ** 2 - 1e-9))
    want = round_up(c + max(C_FLOOR, MARGIN * max(0, need - c)))
    return need, want, min(pool, want)


def over_one(step, ids, data, times, bound=1.01):
    """Warnings for |L| > bound (a coherence cannot exceed 1: a CCE division by a near-zero cluster term, e.g. a strongly
    coupled spin whose own coherence passes near zero; reported, the data are left as they are; user 2026-10-07)."""
    out = []
    for k, x in zip(ids, data):
        a = np.abs(np.asarray(x))
        if (a > bound).any():
            i = np.unravel_index(np.nanargmax(a), a.shape)
            out.append(f"{step} configuration {k}: |L| > {bound:g} at {int((a > bound).sum())} points, up to {a.max():.3g} at "
                       f"{times[i[-1]]:.4g} ms (a CCE division by a near-zero cluster term)")
    return out


def suffix(values, vs_reference):
    """The smallest value from which every larger one passes against the largest."""
    ok = [x["verdict"] == "pass" for x in vs_reference]
    k = len(ok)
    while k > 0 and ok[k - 1]:
        k -= 1
    return values[k]


class Runner:
    def __init__(self, input_path, binary):
        self.input_path = Path(input_path).resolve()
        self.binary = Path(binary).resolve()
        autostart.ISOTOPES = autostart.find_isotope_table(self.binary) or autostart.ISOTOPES
        try:
            raw = read_input(self.input_path)       # keys in the names looked up below, case not read as in CCEX (R4)
        except KeyCaseError as e:
            raise Refused(str(e))
        if "auto_convergence" not in raw:
            raise Refused("the input has no auto_convergence block")
        self.block = raw.pop("auto_convergence")
        self.base = raw
        self.root = self.input_path.parent
        self._validate()
        # main.out under mpirun -n N passes N (plan "mpirun으로 부르면 그 rank 수로", 2026-10-10): the input's launcher first, then
        # mpirun -n N, else one rank
        np_env = os.environ.get("CCEX_AUTOCONV_NP")
        if np_env is not None and not re.fullmatch(r"[1-9]\d*", np_env):
            raise Refused(f"CCEX_AUTOCONV_NP = {np_env!r} is not a rank count")
        if "launcher" in self.block:
            self.launcher_from = "input"
        elif np_env is not None:
            self.o["launcher"], self.launcher_from = ["mpirun", "-n", np_env], "mpirun"
        else:
            self.launcher_from = None
        self.work = (self.root / self.block["workdir"]).resolve()
        self.work.mkdir(parents=True, exist_ok=True)
        self.decisions_path = self.work / "decisions.jsonl"
        self.binary_sha = _sha(self.binary)
        self.ranks = self._ranks(self.o["launcher"])
        self.threads = self._threads()
        self.core_hours = 0.0
        self.max_time = {}                  # step -> longest wall time (s), for the walltime check
        self.t_seed, self.npoints, self.start_warnings, self.resolved = None, None, [], {}
        self._counted = set()               # run folders whose cost is already in core_hours
        self.unstable = {}                  # configuration -> why it was left out (plan "불안정 configuration 빼기", 2026-10-10)
        self.corrected = {}                 # configuration -> divergent points that took noDiv (plan "발산한 점 보정", 2026-10-10)
        # rdip_ladder exists for order 2 of cCCE, gCCE and pCCE with standard clusters only (improvements.md A23); elsewhere
        # every rung is its own run, gathered under the ladder's names so the rest reads them alike (plan "지원 방법 넓히기" M3)
        b_ = self.base
        self.ladder = _method(b_) in ("cce", "gcce", "pcce") and b_.get("order") == 2 and \
            str(b_.get("cluster_mode", "standard")).lower() in ("standard", "pcce")
        self.pcce = _method(b_) == "pcce"
        self.step_now, self.interrupted, self.status_now = 0, None, "running"

    # ---------------------------------------------------------------- input
    @staticmethod
    def _flatten(block):
        """The nested block (time_window / rdip / rbath / configurations, user 2026-10-08) as the flat names the runner
        uses; an earlier flat name is refused with where it went."""
        if not isinstance(block, dict):
            raise Refused("auto_convergence must be an object")
        moved = {k: REPLACED.get(k) or f"now {WHERE[k]}" for k in block if k in REPLACED or k in WHERE}
        if moved:                          # user 2026-10-08: renamed, moved into a section, or removed
            raise Refused("; ".join(f"auto_convergence.{k}: {v}" for k, v in moved.items()))
        unknown = set(block) - TOP_KEYS
        if unknown:
            raise Refused(f"unknown auto_convergence keys: {sorted(unknown)} (sections: time_window, rdip, rbath, configurations)")
        on = block.get("enabled", True)
        if not isinstance(on, bool):
            raise Refused("auto_convergence.enabled must be true or false")
        if not on:
            raise Refused("auto_convergence.enabled is false: main.out -f input.json runs the ordinary calculation")
        flat = {k: v for k, v in block.items() if k not in SECTIONS and k not in ("enabled", "resume")}
        nulls = [k for k, v in flat.items() if v is None]
        if nulls:
            raise Refused(f"auto_convergence.{nulls[0]} is null: give a value or leave it out")
        for sec, names in SECTIONS.items():
            part = block.get(sec, {})
            if not isinstance(part, dict):
                raise Refused(f"auto_convergence.{sec} must be an object")
            extra = set(part) - set(names) - (set(FILE_KEYS) | set(SPLIT_KEYS) if sec == "configurations" else set())
            gone = {k: REPLACED[k] for k in extra if k in REPLACED}
            if sec == "configurations":
                gone.update({k: f"now configurations.{u}" for u, k in SPLIT_KEYS.items() if k in extra})
            if gone:
                raise Refused("; ".join(f"auto_convergence.{sec}.{k}: {v}" for k, v in gone.items()))
            if extra:
                raise Refused(f"unknown keys in auto_convergence.{sec}: {sorted(extra)} (known: {sorted(names)})")
            for k, v in part.items():
                if v is None:                  # null is not "find it": leave the key out for that (review ...-autostart)
                    raise Refused(f"auto_convergence.{sec}.{k} is null: leave it out to have it found, or give a value")
                if k in names:
                    flat[names[k]] = v
            if sec == "configurations" and sec in block:
                flat["configurations"] = {k: v for k, v in part.items() if k in FILE_KEYS}
                flat["configurations"].update({SPLIT_KEYS[k]: v for k, v in part.items() if k in SPLIT_KEYS})
        return flat

    def _validate(self):
        self.resume = self.block.get("resume", "continue") if isinstance(self.block, dict) else "continue"
        if self.resume not in ("continue", "fresh"):
            raise Refused(f"auto_convergence.resume must be continue or fresh, not {self.resume!r}")
        b = self._flatten(self.block)
        conv = {}
        for sec, flag in (("time_window", "converge_time_window"), ("rdip", "converge_rdip"), ("rbath", "converge_rbath"),
                          ("configurations", "converge_config")):
            v = b.pop(flag, True)
            if not isinstance(v, bool):
                raise Refused(f"auto_convergence.{sec}.converge must be true or false")
            conv[sec] = v
        self.converge = conv
        if not any(conv.values()):
            raise Refused("every section has converge false: nothing to converge (run CCEX itself with these values)")
        # converge false: the value is used as it is and is required; converge true: no fixed value
        for sec, val, others in (("time_window", "set_time_window", ()), ("rdip", "value_rdip", ("candi_rdip", "num_candi_rdip")),
                                 ("rbath", "value_rbath", ("candi_rbath", "num_candi_rbath")),
                                 ("configurations", "number_config", ())):
            key = WHERE.get(val, f"{sec}.number")
            if not conv[sec] and val not in b:
                raise Refused(f"auto_convergence.{sec}.converge is false: give {key} (the value to use)")
            if conv[sec] and val in b:
                raise Refused(f"auto_convergence.{key} is a fixed value: set {sec}.converge false to use it, or leave it out")
            bad = [WHERE[o] for o in others if o in b and not conv[sec]]
            if bad:
                raise Refused(f"{sec}.converge is false: {bad} are not used, leave them out")
        if "states_config" in b and conv["configurations"]:
            raise Refused("auto_convergence.configurations.states is a fixed value: set configurations.converge false to use it")
        if "workdir" not in b:
            raise Refused("auto_convergence.workdir is required")
        if "configurations" not in b:
            raise Refused("auto_convergence.configurations.bathfile is required")
        # absent keys are found automatically (plan "자동 시작값", 2026-10-08); a key given wins, null or empty is refused
        for k in ("candi_rdip", "candi_rbath"):
            if k in b and (b[k] is None or b[k] == []):
                raise Refused(f"auto_convergence.{WHERE[k]} is {b[k]!r}: leave it out to have it found, or give a value")
        self.auto = {k for k, sec in (("candi_rdip", "rdip"), ("candi_rbath", "rbath")) if k not in b and self.converge[sec]}
        self.o = dict(DEFAULTS, **{k: v for k, v in b.items() if k in DEFAULTS})
        base = self.base
        method = _method(base)
        if method not in METHODS:
            raise Refused(f"--auto_conv supports method cCCE, gCCE or cCCE with cluster_mode pcce (the adaptive block is excluded), "
                          f"not {method}")
        if str(base.get("quantity", "")).lower() != "coherence":
            raise Refused("--auto_conv supports quantity coherence only")
        order = base.get("order")
        if not (isinstance(order, int) and not isinstance(order, bool) and order >= (1 if method == "gcce" else 2)):
            raise Refused(f"order must be an integer >= {1 if method == 'gcce' else 2} for {method}")
        if "save_nodiv" in base:            # CCEX refuses it too (2026-10-13)
            raise Refused('save_nodiv is now savemode.nodiv ("savemode": {"mode": "coherence", "nodiv": false})')
        for path, bad in ((("savemode",), lambda v: str(v.get("mode", "coherence") if isinstance(v, dict) else v).lower() != "coherence"),
                          (("cluster_contributions",), lambda v: not isinstance(v, dict) or v.get("enabled", True) is not False),
                          (("max_clusters_per_order",), lambda v: any(v)), (("rdip_ladder",), lambda v: True),
                          (("adaptive",), lambda v: True), (("pcce", "cache"), lambda v: True),
                          (("pcce", "partition"), lambda v: str(v).lower() == "legacy"),   # legacy partitions depend on the MPI size
                          (("Defect", "axis_file"), lambda v: True), (("Defect", "state_file"), lambda v: True)):   # review 10-07 A6: files not hashed
            v = _at(base, *path)
            if v is not None and bad(v):
                k = ".".join(path)
                why = " (the runner gives every pCCE run its own partition cache)" if k.startswith("pcce.cache") else ""
                raise Refused(f"--auto_conv refuses {k} = {v!r} in the base input{why}")
        dt = base.get("deltat")
        if "deltat" in base and dt is None:
            raise Refused("deltat is null: leave it out to have it found from the window, or give a positive number of ms")
        if dt is None:                     # deltat found from the window: nstep points over it (user 2026-10-08)
            self.auto.add("deltat")
            self.deltat = None
        elif not (isinstance(dt, (int, float)) and not isinstance(dt, bool) and math.isfinite(dt) and dt > 0):
            raise Refused("deltat must be a positive number of ms, or left out to be found from the window")
        else:
            self.deltat = float(dt)
        n = base.get("nstep")
        if n is not None and not (isinstance(n, int) and not isinstance(n, bool) and n >= 3):
            raise Refused("nstep must be an integer >= 3")
        self.nstep0 = n if n is not None else DEFAULT_NSTEP     # time points: nstep of the input (user 2026-10-08)
        # review 10-07 A10: the preflight measures the bath around ONE qubit in the files' own coordinates
        shift = base.get("bath_shift")
        flat = [shift] if isinstance(shift, (int, float)) else (shift or [])
        while any(isinstance(x, list) for x in flat):
            flat = [y for x in flat for y in (x if isinstance(x, list) else [x])]
        if any(float(x) != 0 for x in flat):
            raise Refused("--auto_conv refuses bath_shift (its bath radius and resource checks read the files' own coordinates)")
        qubits = base.get("Qubit", {})
        if isinstance(qubits, dict) and len(qubits.get("qubit", []) or []) > 1:
            raise Refused("--auto_conv supports one qubit (its bath radius is measured from it)")
        o = self.o
        checks = (("threshold_config", lambda v: 0 < v < 1), ("states_threshold", lambda v: 0 < v < 1), ("threshold_rdip", lambda v: 0 < v < 1), ("threshold_rbath", lambda v: 0 < v < 1), ("production_window_factor", lambda v: v >= 1),
                  ("revival_check_factor", lambda v: v >= 0), ("threshold_time_window", lambda v: 0 < v < 1))
        for k, ok in checks:
            v = o[k]
            if not (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and ok(v)):
                raise Refused(f"auto_convergence.{WHERE.get(k, k)} = {v!r} is out of range")
        for k in ("num_candi_rdip", "num_candi_rbath"):      # the first candidates when the lists are left out (1.25 apart)
            if not (isinstance(o[k], int) and not isinstance(o[k], bool) and o[k] >= 2):
                raise Refused(f"auto_convergence.{WHERE[k]} must be an integer >= 2")
        for k, lst in (("num_candi_rdip", "candi_rdip"), ("num_candi_rbath", "candi_rbath")):
            if k in b and lst in b:
                raise Refused(f"give {WHERE[lst]} (the values) or {WHERE[k]} (how many to start with), not both")
        if o["fit_points"] not in ("peaks", "all"):
            raise Refused(f"t2_fit.points must be peaks or all, not {o['fit_points']!r}")
        fl, nr = o["fit_floor"], o["fit_n_range"]
        if not (isinstance(fl, (int, float)) and not isinstance(fl, bool) and 0 < fl < 1):
            raise Refused(f"t2_fit.floor must be a number between 0 and 1, not {fl!r}")
        if not (isinstance(nr, list) and len(nr) == 2 and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in nr)
                and 0 < nr[0] < nr[1]):
            raise Refused(f"t2_fit.n_range must be [low, high] with 0 < low < high, not {nr!r}")
        # one T2 fit setting for every step: the window, the rdip and rbath comparisons and the sampling (user 2026-10-08)
        # the fit points: the peaks of each curve (all points when fewer than 3), user 2026-10-10
        self.fit = dict(points=o["fit_points"], floor=float(fl), n_range=(float(nr[0]), float(nr[1])))
        u = o["set_time_window_unit"]
        if u not in TIME_UNITS:
            raise Refused(f"time_window.unit must be one of {sorted(TIME_UNITS)}, not {u!r}")
        w0 = o["set_time_window"]
        if w0 is not None:
            if not (isinstance(w0, (int, float)) and not isinstance(w0, bool) and math.isfinite(w0) and w0 > 0):
                raise Refused("time_window.value must be a positive number (in time_window.unit, default ms)")
            o["set_time_window"] = float(w0) * TIME_UNITS[u]          # kept in ms from here on
        # candi_rdip: the rdip candidates, the largest the reference (computed with the others as its ladder)
        cr = [float(x) for x in b["candi_rdip"]] if "candi_rdip" in b else None
        self.rungs, self.rdip_max = (cr[:-1], cr[-1]) if cr else (None, None)
        self.rbaths = [float(x) for x in b["candi_rbath"]] if "candi_rbath" in b else None
        for name in ("value_rdip", "value_rbath"):           # converge false: one value, used as it is
            v = b.get(name)
            if v is not None and not (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0):
                raise Refused(f"auto_convergence.{WHERE[name]} must be a positive number (Å)")
        if not self.converge["rdip"]:
            self.rungs, self.rdip_max = [], float(b["value_rdip"])
        if not self.converge["rbath"]:
            self.rbaths = [float(b["value_rbath"])]
        self.fixed_config = None
        if not self.converge["configurations"]:
            n, st = b["number_config"], b.get("states_config")
            if not (isinstance(n, int) and not isinstance(n, bool) and n >= 2):
                raise Refused("auto_convergence.configurations.number must be an integer >= 2")
            if self.base.get("nstate", 0) and not (isinstance(st, int) and not isinstance(st, bool) and st >= 1):
                raise Refused("nstate > 0: auto_convergence.configurations.states (states per configuration, >= 1) is required")
            if not self.base.get("nstate", 0) and st is not None:
                raise Refused("nstate 0 (ensemble): auto_convergence.configurations.states is not used, leave it out")
            self.fixed_config = (n, st or 0)
        c = b["configurations"]
        if not isinstance(c, dict) or "bathfile" not in c:
            raise Refused("auto_convergence.configurations.bathfile is required")
        for k in ("pilot", "pool"):
            if k in c and not (isinstance(c[k], list) and c[k]):
                raise Refused(f"configurations.{dict(pilot='test', pool='production')[k]} is {c[k]!r}: leave it out to have it found, or give the IDs")
        self.auto |= {k for k in ("pilot", "pool") if k not in c}
        # ensemble or single-sample from the CCEX nstate of the input (user decision 2026-10-07): 0 (the CCEX default) is the
        # ensemble, one curve per configuration, and only the number of configurations is chosen; > 0 samples states, the
        # value being the pilot's states per configuration, and needs configurations.statefile
        nst = self.base.get("nstate", 0)
        if not isinstance(nst, int) or isinstance(nst, bool) or nst < 0:
            raise Refused("nstate must be a non-negative integer (0: ensemble, > 0: states sampled)")
        self.ensemble = nst == 0
        if "pilot_states" in b:
            raise Refused("pilot_states is replaced by configurations.test_states (nstate of the input only chooses ensemble or states)")
        # the test's size (user 2026-10-10): test_number configurations (default 10) of test_states states (default 5); nstate
        # of the input only chooses the ensemble (0) or sampled states (> 0)
        tn, ts = self.o["test_number"], self.o["test_states"]
        if not (isinstance(tn, int) and not isinstance(tn, bool) and tn >= 2):
            raise Refused("configurations.test_number must be an integer >= 2")
        if "test_number" in b and "pilot" in c:
            raise Refused("configurations.test_number and configurations.test: give one (the list is the test configurations)")
        if self.ensemble and "test_states" in b:
            raise Refused("nstate 0 (ensemble): configurations.test_states is not used, leave it out")
        if not self.ensemble and not (isinstance(ts, int) and not isinstance(ts, bool) and ts >= 2):
            raise Refused("configurations.test_states must be an integer >= 2 (the states of a test configuration)")
        self.o["pilot_states"] = 0 if self.ensemble else ts
        # the largest numbers (plan "최대 configuration·state 수 키", user 2026-10-10): max_configurations (default every production
        # file) and max_states (default 30, the state test)
        mxc, mxs = self.o["max_configurations"], self.o["max_states"]
        if mxc is not None and not (isinstance(mxc, int) and not isinstance(mxc, bool) and mxc >= 2):
            raise Refused("configurations.max_configurations must be an integer >= 2 (leave it out to use every production file)")
        if self.ensemble and "max_states" in b:
            raise Refused("nstate 0 (ensemble): configurations.max_states is not used, leave it out")
        if not self.ensemble and not (isinstance(mxs, int) and not isinstance(mxs, bool) and mxs >= ts):
            raise Refused(f"configurations.max_states must be an integer >= test_states ({ts})")
        if self.ensemble and ("statefile" in c or "exstatefile" in c):    # the tool's pattern names (configurations.*)
            raise Refused("nstate 0 (ensemble) does not read state files: remove configurations.statefile or set nstate > 0")
        # nstate > 0 without state files: the runner makes them for each run and removes them after (plan "state 파일 자동",
        # user 2026-10-09; CCEX's own random states change with rbath and would break the paired comparisons)
        self.auto_states = not self.ensemble and "statefile" not in c
        # keyed in lower case: CCEX finds a defect by its name ignoring case (strcasecmp; review 6548a9 R4)
        self.defect_axes = {str(d_["name"]).lower(): int(d_.get("n_configurations", 0) or 0)
                            for d_ in _at(self.base, "Defect", "species", default=[])
                            if isinstance(d_, dict) and "name" in d_ and int(d_.get("n_configurations", 0) or 0) > 0}
        self.auto_avaax = {}
        if self.auto_states and "exstatefile" in c:
            raise Refused("configurations.exstatefile needs configurations.statefile (the states the runner makes cover the bath spins only)")
        mc = self.o["min_configurations"]
        if not isinstance(mc, int) or isinstance(mc, bool) or mc < 2:
            raise Refused("min_configurations must be an integer >= 2")
        if mxc is not None and mxc < mc:
            raise Refused(f"configurations.max_configurations {mxc} is below min_configurations {mc}")
        if not self.ensemble and self.o["pilot_states"] < 2:
            raise Refused("the pilot needs at least 2 states")
        self.pilot = [int(x) for x in c["pilot"]] if "pilot" in c else None
        self.pool = [int(x) for x in c["pool"]] if "pool" in c else None
        self.patterns = {k: c[k] for k in ("bathfile", "statefile", "avaaxfile", "exstatefile") if k in c}
        if not self.auto:
            self._check_lists()

    def _check_lists(self):
        """The candidate lists and the configuration split, given or found."""
        # the production configurations used: at most max_configurations of the list (which stays whole, so raising the cap
        # later keeps the plan; user 2026-10-10)
        mx = self.o["max_configurations"]
        self.cap = len(self.pool) if mx is None else min(mx, len(self.pool))
        rungs, rmax, rb = self.rungs, self.rdip_max, self.rbaths
        if (self.converge["rdip"] and not rungs) or (self.converge["rbath"] and len(rb) < 2):   # converge() compares two values
            raise Refused("rdip.candidates and rbath.candidates need at least two values each")
        for name, seq in (("rdip.candidates", rungs + [rmax]), ("rbath.candidates", rb)):
            if not seq or any(not math.isfinite(x) or x <= 0 for x in seq) or any(y <= x for x, y in zip(seq, seq[1:])):
                raise Refused(f"{name} must be positive, finite and strictly increasing")
        cut = np.float32(_radius(self.base.get("rdip"))[0] or 0.0)
        if any(np.float32(r) <= cut for r in rungs):        # rdip_ladder.cpp needs inner rdip < rung (review ...-autostart #9)
            raise Refused(f"rdip candidates at or below the inner rdip {float(cut):g}: {[r for r in rungs if np.float32(r) <= cut]}")
        names = [_gname(r) for r in rungs]
        if len(set(names)) != len(names) or _gname(rmax) in names:
            raise Refused("two rdip rungs print the same name (%g of the float value)")
        if set(self.pilot) & set(self.pool) or len(set(self.pilot)) != len(self.pilot) or len(set(self.pool)) != len(self.pool):
            raise Refused("configurations.test and configurations.production IDs must be distinct and not overlap")
        if self.fixed_config is not None:
            n, st = self.fixed_config
            if n > self.cap:
                raise Refused(f"configurations.number {n} is more than the pool holds ({len(self.pool)}"
                              + (f", max_configurations {self.o['max_configurations']})" if self.o["max_configurations"] is not None else ")"))

        if len(self.pilot) < 2:
            raise Refused("the pilot needs at least 2 configurations")

    @staticmethod
    def _shown(keys):
        """The keys as the user writes them (nested names)."""
        return [dict(pilot="configurations.test", pool="configurations.production").get(k) or WHERE.get(k, k) for k in keys]

    def _split(self, ids, seed, mc):
        """The test and production configurations not given, kept across calls in <workdir>/split.json (plan "bath 만들기 빼기",
        user 2026-10-09): while every recorded configuration is still there the split stays and the files the user added go
        after the production in a seeded order, so a short production continues on them; otherwise a new split."""
        path = self.work / "split.json"
        try:
            old = json.loads(path.read_text()) if path.is_file() else None
        except ValueError:
            old = None
        note = None
        if old is not None:
            have = set(ids)
            keep = old.get("seed") == seed and old.get("test_number", self.o["test_number"]) == self.o["test_number"] and \
                old.get("rule") == SPLIT_RULE and \
                set(old["pilot"] + old["pool"]) <= have and \
                (self.pilot is None or self.pilot == old["pilot"]) and (self.pool is None or self.pool == old["pool"])
            if keep:
                new = autostart.order([k for k in ids if k not in set(old["pilot"] + old["pool"])], seed)
                pilot, pool = old["pilot"], old["pool"] + ([] if self.pool is not None else new)
                _atomic_write(path, json.dumps(dict(old, pool=pool)))
                return pilot, pool, (f"{len(new)} configurations added since the last call: after the production" if new and self.pool is None else None)
            note = ("the configurations recorded in split.json are no longer all there (or the seed, test_number, a given list or the "
                    "split rule changed): split again")
        # the test configurations evenly over the nearest bath spin's distance to the qubit (user 2026-10-10, P1 CR90: the
        # closer the spin the faster the decay, and ten drawn at random were closer than the production); the rest in the seeded order
        near = {}
        def even(cands, n):
            ordered = autostart.order(cands, seed)
            for k in ordered:
                if k not in near:
                    near[k] = float(np.min(np.linalg.norm(self.bath_xyz(k) - self.qubit_xyz(), axis=1)))
            by = sorted(ordered, key=lambda k: near[k])            # stable: equal distances keep the seeded order
            pick = sorted({min(len(by) - 1, int((i + 0.5) * len(by) / n)) for i in range(n)})
            chosen = [by[i] for i in pick]
            return chosen, [k for k in ordered if k not in set(chosen)]
        if self.pilot is None and self.pool is None:
            p = min(self.o["test_number"], len(ids) - mc)
            if p < 2:
                if self.fixed_config is not None or len(ids) < 4:
                    raise ValueError(f"{len(ids)} configurations found: the pilot needs 2 and the pool min_configurations "
                                     f"{mc} (at least {mc + 2})")
                p = 2                       # fewer than the minimum: the pilot takes 2, the rest is computed and reported short
            pilot, pool = even(ids, p)
        else:
            keep = self.pilot or self.pool          # the list given stays as it is; the other from the rest
            rest = [k for k in ids if k not in keep]
            if self.pilot is None:
                pilot, pool = even(rest, min(self.o["test_number"], len(rest)))[0], self.pool
            else:
                pilot, pool = self.pilot, autostart.order(rest, seed)
        _atomic_write(path, json.dumps(dict(seed=seed, test_number=self.o["test_number"], rule=SPLIT_RULE, pilot=pilot, pool=pool,
                                            nearest_spin={str(k): round(v, 3) for k, v in near.items() if k in pilot})))
        return pilot, pool, note

    def resolve(self):
        """Fill the keys left out (plan "자동 시작값", 2026-10-08): the configurations from the bathfile pattern, split in a
        seeded order, and the first rbath, rdip and time window from the bath density and the gyromagnetic ratios. The
        values and their evidence go to self.resolved (summary.txt section 0); they are part of the plan fingerprint."""
        o = self.o
        rec = dict(found=sorted(self._shown(self.auto)))
        if self.pilot is None or self.pool is None:
            try:
                ids = autostart.discover(self.patterns["bathfile"], self.root)
            except ValueError as e:
                raise Refused(str(e))
            need = o["pilot_states"]
            missing = {k: [p.name for p in self.files(k, need).values() if not p.is_file()] for k in ids}
            missing = {k: v for k, v in missing.items() if v}
            if missing:                     # never left out quietly: a missing file could select the configurations
                raise Refused("configurations with missing files (give configurations.test and .production to choose them by hand): "
                              + "; ".join(f"{k}: {v[:3]}" for k, v in list(missing.items())[:5]))
            seed = self.base.get("seed", 1)
            mc = self.fixed_config[0] if self.fixed_config else o["min_configurations"]
            try:
                self.pilot, self.pool, note = self._split(ids, seed, mc)
            except ValueError as e:
                raise Refused(str(e))
            if note:
                self.start_warnings.append(note)
            rec.update(found_configurations=ids, split_seed=seed)
        need_seed = self.rbaths is None or self.rdip_max is None or \
            (self.deltat is None and o["set_time_window"] is None)
        if need_seed:
            # counted spins (plan "자동 시작값 — 2D·섞인 bath까지, 단순판", 2026-10-08): any geometry the coverage allows, any
            # number of species
            q = self.qubit_xyz()
            configs = []
            for k in self.pilot:
                sp, xyz = autostart.species_and_xyz(self.files(k, 1)["bathfile"])
                configs.append((sp, np.linalg.norm(xyz - q, axis=1)))
            names = sorted({n for sp, _ in configs for n in sp})
            r_cov = min(self.bath_radius(k) for k in self.pilot + self.pool)
            cut = float(_radius(self.base.get("rbath"))[0] or 0.0)
            gf = self.base.get("gyrofile")
            gf = None if gf is None else (Path(gf) if Path(gf).is_absolute() else self.root / gf)
            g = autostart.gyros(gf, set(names))
            unknown = [n for n in names if n not in g]
            if unknown:
                raise Refused(f"species {unknown}: neither in the gyrofile nor in the built-in isotope table"
                              + ("" if autostart.ISOTOPES.is_file() else f" (not found: {autostart.ISOTOPES}; give them in the gyrofile)"))
            qb = self.base.get("Qubit", {})
            q0 = qb.get("qubit", [{}])[0] if isinstance(qb, dict) and qb.get("qubit") else {}
            qg = float(q0.get("gyro", autostart.GAMMA_ELECTRON))
            dm = float(q0.get("alphams", _at(self.base, "Qubit", "alphams", default=1.0))) - \
                float(q0.get("betams", _at(self.base, "Qubit", "betams", default=0.0)))
            try:
                sd = autostart.count_seeds(configs, cut, r_cov, {n: g[n][:2] for n in names}, qg, dm,
                                           float(_radius(self.base.get("rdip"))[0] or 0.0), n_rdip=o["num_candi_rdip"],
                                           n_rbath=o["num_candi_rbath"])
            except ValueError as e:
                raise Refused(str(e))
            if self.rbaths is None:
                self.rbaths = sd["rbath_candidates"]
            if self.rdip_max is None:
                if not sd["rdip_rungs"]:
                    raise Refused(f"no rdip candidate above the inner rdip below the top {sd['rdip_max']:g}: give rdip.candidates")
                self.rdip_max, self.rungs = sd["rdip_max"], sd["rdip_rungs"]
            self.t_seed = sd["t_seed"]
            warn = list(sd["notes"])
            if len(names) > 1 or getattr(self, "_thin", False):
                warn.append("the clock of a layer or of several species is an uncalibrated start (K from 3D NV-13C and P1)")
            self.start_warnings += warn
            rec.update(species=sd["species"], qubit_gamma=qg, dm=dm, r_cov=r_cov,
                       **{k: v for k, v in sd.items() if k not in ("notes", "species")})
            rec["geometry"] = "a layer" if getattr(self, "_thin", False) else "bulk"
        self.npoints = self.nstep0             # every run on nstep points (user 2026-10-10: fixed, deltat only the first window)
        rec.update(pilot=self.pilot, pool=self.pool, candi_rbath=self.rbaths, candi_rdip=self.rungs + [self.rdip_max])
        self.resolved = rec
        self._check_lists()

    def _path(self, pattern, k):
        p = Path(pattern.replace("{k}", str(k)))
        return p if p.is_absolute() else (self.root / p)

    def files(self, k, nstate, first=1):
        """Every file a configuration's run reads (bath, its states first..first + nstate - 1, axes, sub-bath states)."""
        out = {"bathfile": self._path(self.patterns["bathfile"], k)}
        if not self.ensemble and not self.auto_states:
            prefix = self._path(self.patterns["statefile"], k)
            out.update({f"state{s}": Path(f"{prefix}{s}") for s in range(first, first + nstate)})
        if "avaaxfile" in self.patterns:
            out["avaaxfile"] = self._path(self.patterns["avaaxfile"], k)
        if "exstatefile" in self.patterns:
            ex = self._path(self.patterns["exstatefile"], k)
            out.update({f"exstate{s}": Path(f"{ex}{s}") for s in range(first, first + nstate)})
        return out

    def count_states(self, k):
        """State files 1, 2, ... of configuration k that exist (and its exstate files when given), counted until one is missing."""
        if self.auto_states:
            return MADE_STATES
        n = 0
        while all(p.is_file() for name, p in self.files(k, n + 1).items() if name in (f"state{n + 1}", f"exstate{n + 1}")):
            n += 1
        return n

    def preflight(self):
        """Files present, state files as long as the bath file, no two configurations with the same bath, sizes bounded."""
        seen = {}
        self.source_sha = {}
        self._state_names = set()
        for k in self.pilot + self.pool:
            fs = self.files(k, self.o["pilot_states"])
            missing = [str(p) for p in fs.values() if not p.is_file()]
            if missing:
                raise Refused(f"configuration {k}: missing files {missing[:3]}{' ...' if len(missing) > 3 else ''}")
            if not self.ensemble and k in self.pool:      # every state file of a pool configuration may be used
                fs = self.files(k, self.count_states(k))
            h = _sha(fs["bathfile"])
            self.source_sha[k] = {n: _sha(p) for n, p in sorted(fs.items())}
            if h in seen:
                raise Refused(f"configurations {seen[h]} and {k} have the same bath file content")
            seen[h] = k
            nb = int(round(float(fs["bathfile"].read_text().split()[0]))) - 1
            # CCEX reads a spin's state from the line of the state file that has the spin's line of the bath file, so a
            # state file needs at least as many entries as the bath has spins (the P1 files of the study hold one more and
            # a header with their own line count; the NV-13C example holds as many and the spin count)
            if self.auto_states or self.defect_axes:
                self._state_names |= set(autostart.species_and_xyz(fs["bathfile"])[0])
            ns = nb if self.ensemble or self.auto_states else sum(1 for line in fs["state1"].read_text().splitlines()[1:] if line.strip())
            if ns < nb:
                raise Refused(f"configuration {k}: {ns} states per state file, {nb} spins in the bath file")
        if self.auto_states:                # the spins of every species the made states need
            gf = self.base.get("gyrofile")
            gf = None if gf is None else (Path(gf) if Path(gf).is_absolute() else self.root / gf)
            g = autostart.gyros(gf, self._state_names)
            unknown = sorted(n for n in self._state_names if n not in g)
            if unknown:
                raise Refused(f"species {unknown}: no spin for the state files the runner makes (gyrofile or the built-in table"
                              + ("" if autostart.ISOTOPES.is_file() else f", not found: {autostart.ISOTOPES}") + ")")
            self.state_spins = {n: g[n][0] for n in self._state_names}
        # avaax files made per run when none are given and a defect with axes is in the baths (plan "avaax 파일 자동", 2026-10-09)
        present = {x.lower() for x in self._state_names}
        self.auto_avaax = {n: v for n, v in self.defect_axes.items() if n in present} if "avaaxfile" not in self.patterns else {}
        # states the sampling may use: the fewest over the pool (the state files are counted; states_available removed 2026-10-08)
        self.states_avail = 0 if self.ensemble else min(self.count_states(k) for k in self.pool)
        if self.fixed_config is not None and self.fixed_config[1] > self.states_avail:
            raise Refused(f"configurations.states {self.fixed_config[1]} is more than the state files of every pool configuration hold ({self.states_avail})")
        # the largest rbath the bath files allow (user decision 2026-10-07): the smallest box radius over the configurations,
        # (source_radius was removed on 2026-10-08)
        radii = [self.bath_radius(k) for k in self.pilot + self.pool]
        self.max_rbath = min(radii)
        if self.rbaths[-1] > self.max_rbath * (1 + 1e-9):
            raise Refused(f"the largest rbath candidate {self.rbaths[-1]:g} is beyond what the bath files cover ({self.max_rbath:.4g} A)")

    def bath_radius(self, k):
        """Largest rbath the bath file of configuration k covers: the distance from the qubit to the nearest face of the box
        the bath spins span; a thin direction (a layer: thin files) is no limit - sqrt(a^2 + h^2) from the in-plane
        edge distance a and the qubit's height h. Beyond it the run asks for larger baths (user 2026-10-08)."""
        xyz, q = self.bath_xyz(k), self.qubit_xyz()
        side = xyz.max(0) - xyz.min(0)
        if side.min() < autostart.FLAT_RATIO * side.max():
            self._thin = True
            return autostart.plane_coverage(xyz, q)
        return float(min((q - xyz.min(0)).min(), (xyz.max(0) - q).min()))

    def qubit_xyz(self):
        if _at(self.base, "Qubit", "file") is not None:
            p = Path(self.base["Qubit"]["file"]); p = p if p.is_absolute() else self.root / p
            return np.array([float(x) for x in p.read_text().split()[:3]])
        return np.array(self.base["Qubit"]["qubit"][0]["xyz"], float)

    def bath_xyz(self, k):
        return np.loadtxt(self.files(k, 1)["bathfile"], skiprows=1, usecols=(0, 1, 2), ndmin=2)

    def inside(self, k, rbath):
        """Bath spins CCEX keeps: inner rbath <= distance <= rbath, both read as float (bath.cpp, BathArray_readBath)."""
        xyz = self.bath_xyz(k)
        r = np.linalg.norm(xyz - self.qubit_xyz(), axis=1)
        return xyz[(r <= np.float32(rbath)) & (r >= np.float32(_radius(self.base.get("rbath"))[0]))]

    def pair_count(self, xyz, r):
        """Pairs CCEX connects: inner rdip < distance < rdip (bath.cpp, the connectivity map)."""
        lo, n = np.float32(_radius(self.base.get("rdip"))[0]), 0
        for i in range(len(xyz) - 1):
            d = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1)
            n += int(((d < np.float32(r)) & (d > lo)).sum())
        return n

    @staticmethod
    def _ranks(launcher):
        """Rank count stated by the launcher (-n N, -np N, --np N, --ntasks N, --ntasks=N, -nN); no launcher is one rank.
        A count that cannot be read is None: the core-hours (reported only) are then unknown, never guessed."""
        if not launcher:
            return 1
        found = []
        for i, a in enumerate(launcher):
            nxt = launcher[i + 1] if i + 1 < len(launcher) else ""
            if a in ("-n", "-np", "--np", "--ntasks"):
                found.append(nxt)
            elif a.startswith("--ntasks=") or a.startswith("--np="):
                found.append(a.split("=", 1)[1])
            elif re.fullmatch(r"-n\d+", a):
                found.append(a[2:])
        if len(found) != 1 or not re.fullmatch(r"[1-9]\d*", found[0]):
            return None                    # not read: the core-hours are reported as unknown
        return int(found[0])

    @staticmethod
    def _threads():
        """OpenMP threads per rank: OMP_NUM_THREADS when set (a positive integer), else 1 - the build presets do not link
        OpenMP. Any other value (a list such as "4,2", 0, text) is None: the core-hours are then unknown."""
        v = os.environ.get("OMP_NUM_THREADS")
        if v is None:
            return 1
        if not re.fullmatch(r"\s*[1-9]\d*\s*", v):
            return None                    # not read: the core-hours are reported as unknown
        return int(v)

    # ---------------------------------------------------------------- runs
    def child(self, k, nstate, nstep, rbath, rdip, ladder, first=1):
        c = copy.deepcopy(self.base)
        fs = self.files(k, nstate, first)
        # the inner radii of the base input stay ([inner, outer]; A47)
        inner_b, inner_d = (v[0] if isinstance(v, list) and len(v) == 2 else None for v in (c.get("rbath"), c.get("rdip")))
        c.update(bathfile=[str(fs["bathfile"])], nstate=nstate, nstep=nstep, deltat=self.deltat,
                 rbath=rbath if inner_b is None else [inner_b, rbath], rdip=rdip if inner_d is None else [inner_d, rdip],
                 outfile="out")
        # the configuration's files under the CCEX input keys (A47): bath_state_file, Defect.axis_file, Defect.state_file
        if self.ensemble:
            c.pop("bath_state_file", None)
        elif self.auto_states or first > 1:              # a later block of given state files: copied there as 1, 2, ...
            c["bath_state_file"] = "made/cfg_state"      # made in the run folder just before CCEX starts, removed after
        else:
            c["bath_state_file"] = str(self._path(self.patterns["statefile"], k))
        if "avaaxfile" in self.patterns:
            c.setdefault("Defect", {})["axis_file"] = str(fs["avaaxfile"])
        elif self.auto_avaax:
            c.setdefault("Defect", {})["axis_file"] = "made/cfg_avaax"
        if "exstatefile" in self.patterns:
            c.setdefault("Defect", {})["state_file"] = "made/cfg_exstate" if first > 1 else \
                str(self._path(self.patterns["exstatefile"], k))
        for path in INPUT_FILES:
            v = _at(c, *path)
            if isinstance(v, str) and not Path(v).is_absolute():
                _at(c, *path[:-1])[path[-1]] = str((self.root / v).resolve())
        if ladder and self.ladder:
            c["rdip_ladder"] = list(ladder)
        else:
            c.pop("rdip_ladder", None)
        if self.pcce:
            c.setdefault("seed", 1)                      # PARTITION_RULE: without it every run searches the partition again
        return c, fs

    def partition_key(self, cin, fs, extra):
        """The pCCE partition cache key (PARTITION_RULE): the bath spins inside rbath (bath file content, rbath with its inner radius,
        the qubit), seed, pcce.partition_size, trials, max_iterations, kmeans_pp and the binary are all in it; rdip is not."""
        free = {n: v for n, v in cin.items() if n not in PARTITION_FREE}
        if isinstance(free.get("Defect"), dict) and "state_file" in free["Defect"]:   # Defect.state_file: states too
            free["Defect"] = {n: v for n, v in free["Defect"].items() if n != "state_file"}
            if not free["Defect"]:
                free.pop("Defect")
        return hashlib.sha256(_canon(dict(rule=PARTITION_RULE, input=free,
                                          files={n: _sha(p) for n, p in fs.items() if not n.startswith(("state", "exstate"))},
                                          other={n: _sha(p) for n, p in extra.items() if Path(p).is_file()},
                                          binary=self.binary_sha)).encode()).hexdigest()

    def run(self, step, k, nstate, nstep, rbath, rdip, ladder=(), first=1):
        """Run (or reuse) one CCEX calculation; returns its folder. first > 1: states first ... first + nstate - 1 (written or
        copied to made/ as 1 ... nstate; the state test, 2026-10-10)."""
        cin, fs = self.child(k, nstate, nstep, rbath, rdip, ladder, first)
        extra = {".".join(path): str((self.root / v).resolve() if not Path(v).is_absolute() else v)
                 for path in INPUT_FILES for v in [_at(cin, *path)] if isinstance(v, str)}
        made = {}
        if self.auto_states and nstate:
            made["states"] = STATE_RULE
        if first > 1:
            made["first_state"] = first
        if self.auto_avaax:
            made.update(avaax=AVAAX_RULE, navaax=self.auto_avaax)
        if made:
            made.update(seed=int(self.base.get("seed", 1)), k=k)
        cache = None
        if self.pcce:                      # A17 (1): the partition cache, named by what the partition depends on
            ckey = self.partition_key(cin, fs, extra)
            cache = self.work / "pcce_cache" / f"{ckey}.cache"
            # relative to the run folder (runs/<step>/<name>/), so the run key does not hold where the workdir is
            cin.setdefault("pcce", {})["cache"] = dict(mode="readwrite", file=f"../../../pcce_cache/{cache.name}")
        key = hashlib.sha256(_canon(dict(input=cin, rungs=[] if self.ladder else list(ladder), files={n: _sha(p) for n, p in sorted(fs.items())},
                                         **(dict(states_made=made) if made else {}),
                                         other={n: _sha(p) for n, p in extra.items() if Path(p).is_file()},
                                         binary=self.binary_sha, launcher=self.o["launcher"])).encode()).hexdigest()
        d = self.work / "runs" / step / f"{k}-{key[:12]}"
        tags = [f"_state{s}" for s in range(1, nstate + 1)] if nstate else [""]       # nstate 0: out_wiDiv (ensemble)
        kinds = ("wiDiv", "noDiv") if _nodiv(cin) else ("wiDiv",)   # savemode.nodiv false: wiDiv only
        expected = [f"out{g}_{w}" for g in tags for w in kinds]
        expected += [f"out_rdip{_gname(r)}{g}_{w}" for r in ladder for g in tags for w in kinds]
        done = d / "done.json"
        if done.is_file():
            try:
                rec = json.loads(done.read_text())
                if rec["key"] == key and all(_sha(d / f) == h for f, h in rec["outputs"].items()):
                    # the cost it had when it ran (records before 2026-10-08 hold no threads: counted as one)
                    cost = rec.get("core_hours", None if rec.get("ranks") is None else rec["seconds"] * rec["ranks"] / 3600)
                    if d not in self._counted:     # each run counts once per call, however often it is asked for
                        self._counted.add(d)
                        self._account(step, rec["seconds"], cost)
                    if ladder and not self.ladder:   # the rungs' own runs count too (review 4d2d3a R2): reused from their records
                        for r in ladder:
                            self.run(step, k, nstate, nstep, rbath, r, ())
                    return d
            except (OSError, ValueError, KeyError):
                pass
        self._check_resources(step)
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        _atomic_write(d / "input.json", json.dumps(cin, indent=1))
        cmd = list(self.o["launcher"]) + [str(self.binary), "-f", "input.json"]
        if made:                                  # the files the runner makes for this run (removed after)
            (d / "made").mkdir()
            if "states" in made:
                write_states(d / "made" / "cfg_state", fs["bathfile"], k, nstate, made["seed"], self.state_spins, first)
            elif first > 1:                       # given state files of a later block, renumbered from 1
                for i, st in enumerate(range(first, first + nstate), 1):
                    shutil.copyfile(fs[f"state{st}"], d / "made" / f"cfg_state{i}")
                    if f"exstate{st}" in fs:
                        shutil.copyfile(fs[f"exstate{st}"], d / "made" / f"cfg_exstate{i}")
            if "avaax" in made:
                write_avaax(d / "made" / "cfg_avaax", fs["bathfile"], k, made["seed"], self.auto_avaax)
        if cache is not None:              # a cache file is read only when its hash is the one recorded when a run stored it
            cache.parent.mkdir(exist_ok=True)
            sums = Path(f"{cache}.sha256")
            if cache.is_file() and (not sums.is_file() or sums.read_text().strip() != _sha(cache)):
                cache.unlink()
                sums.unlink(missing_ok=True)
        start = time.time()
        with open(d / "process.log", "w") as log:
            p = subprocess.Popen(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
            try:
                code = p.wait()                   # no time limit of its own: the job's walltime check stops before a run
            finally:
                if p.poll() is None:              # the runner stopped (signal, error): the child goes with it
                    os.killpg(p.pid, signal.SIGKILL); p.wait()
                if made:                          # the made files are not kept (storage)
                    shutil.rmtree(d / "made", ignore_errors=True)
        seconds = time.time() - start
        cost = None if self.ranks is None or self.threads is None else seconds * self.ranks * self.threads / 3600
        self._counted.add(d)
        self._account(step, seconds, cost)
        logtext = (d / "process.log").read_text(errors="replace")
        if code != 0:
            raise RuntimeError(f"{d}: CCEX exited with {code}")
        subbath = any(int(d.get("n_nuclear_spins", 0) or 0) > 0 for d in _at(self.base, "Defect", "species", default=[])
                      if isinstance(d, dict))
        if any(m in logtext for m in RANDOM_FALLBACK) or (subbath and SUBBATH_FALLBACK in logtext):
            raise RuntimeError(f"{d}: bath states or axes fell back to random values")
        own = [f for f in expected if self.ladder or not f.startswith("out_rdip")]
        missing = [f for f in own if not (d / f).is_file()]
        if missing:
            raise RuntimeError(f"{d}: missing outputs {missing[:3]}")
        part = None
        if cache is not None:              # [partition-cache] status=hit (read) or stored (searched and written), rank 0's log
            st = re.findall(r"^\[partition-cache\] status=(hit|stored) ", logtext, re.M)
            if len(st) != 1 or not cache.is_file():
                raise RuntimeError(f"{d}: the partition cache was neither read nor stored")
            part = dict(file=cache.name, status=st[0], sha=_sha(cache))
            if st[0] == "stored":
                _atomic_write(Path(f"{cache}.sha256"), part["sha"] + "\n")
        if ladder and not self.ladder:     # every rung its own run, gathered here under the ladder's names (plan M3)
            for r in ladder:
                dr = self.run(step, k, nstate, nstep, rbath, r, ())
                for g in tags:
                    for w in kinds:
                        shutil.copyfile(dr / f"out{g}_{w}", d / f"out_rdip{_gname(r)}{g}_{w}")
        for f in expected:
            t, _ = read_curve(d / f)
            if len(t) != nstep:
                raise RuntimeError(f"{d}/{f}: {len(t)} points, expected {nstep}")
        _atomic_write(done, json.dumps(dict(key=key, seconds=seconds, ranks=self.ranks, threads=self.threads, core_hours=cost,
                                            binary=self.binary_sha, **(dict(pcce_cache=part) if part else {}),
                                            outputs={f: _sha(d / f) for f in expected}), indent=1))
        self.progress()                    # the finished-run counts
        return d

    def _account(self, step, seconds, cost):
        self.core_hours = None if self.core_hours is None or cost is None else self.core_hours + cost
        self.max_time[step] = max(self.max_time.get(step, 0.0), seconds)

    def _check_resources(self, step):
        job = os.environ.get("SLURM_JOB_ID")
        if job and shutil.which("squeue"):
            left = subprocess.run(["squeue", "-h", "-j", job, "-o", "%L"], capture_output=True, text=True).stdout.strip()
            m = re.fullmatch(r"(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+)", left)
            if m:
                d_, h, mi, s = (int(x) if x else 0 for x in m.groups())
                remaining = ((d_ * 24 + h) * 60 + mi) * 60 + s
                if remaining < 1.5 * self.max_time.get(step, 0.0) + WALLTIME_MARGIN_S:
                    raise Stop(f"{remaining} s of walltime left before a {step} run")

    def _flag(self, k, x, t, where, first=1):
        """True when configuration k is unstable in x (states x time): |L| > UNSTABLE_L or nan/inf anywhere (plan "불안정
        configuration 빼기", user 2026-10-10). The first finding is recorded with the nearest bath spin to the qubit; first: the
        number of the first state in x (a later block of the state test)."""
        x = np.asarray(x)
        finite = np.isfinite(x)
        a = np.where(finite, np.abs(np.where(finite, x, 0)), np.inf)
        if a.max() <= UNSTABLE_L:
            return False
        if k not in self.unstable:
            s_, i_ = np.unravel_index(int(np.argmax(a)), a.shape)
            try:
                near = float(np.min(np.linalg.norm(self.bath_xyz(k) - self.qubit_xyz(), axis=1)))
            except (OSError, ValueError):
                near = None
            self.unstable[k] = dict(configuration=int(k), where=where,
                                    reason="nan/inf" if not finite.all() else f"|L| {float(a[s_, i_]):.3g}",
                                    time=float(t[i_]), state=int(s_) + first if not self.ensemble else 0, nearest_spin=near)
        return True

    def _stable_means(self, dirs, prefix, nstate, t, where):
        """Per test configuration, its mean over the states; an unstable one is a row of nan, left out of every comparison
        that holds it (judge_curves keeps the configurations stable on both sides)."""
        rows = []
        for k, d in zip(self.pilot, dirs):
            x = self.curves(d, prefix, nstate)[1]
            rows.append(np.full(x.shape[1], np.nan, dtype=complex) if self._flag(k, x, t, where) else x.mean(0))
        return np.array(rows)

    def curves(self, d, prefix, nstate):
        """configurations' states x time (nstate 0, the ensemble: one row from <prefix>_wiDiv). A point where
        |wiDiv - noDiv| >= POINT_JUMP or wiDiv is nan/inf takes the noDiv value (plan "발산한 점 보정", user 2026-10-10: the
        rule of the group's analyzer, error_correct, plus nan/inf); without noDiv files (savemode.nodiv false) nothing is changed.
        Corrections are recorded per configuration (the largest count of one read)."""
        t = None; x = []
        nodiv = _nodiv(self.base)
        fixed, worst = 0, 0.0
        for g in ([f"_state{s}" for s in range(1, nstate + 1)] if nstate else [""]):
            ti, v = read_curve(d / f"{prefix}{g}_wiDiv")
            t = ti if t is None else t
            if nodiv:
                _, n = read_curve(d / f"{prefix}{g}_noDiv")
                with np.errstate(invalid="ignore"):
                    bad = ~np.isfinite(v) | (np.abs(v - n) >= POINT_JUMP)
                if bad.any():
                    fixed += int(bad.sum())
                    big = np.abs(v[bad & np.isfinite(v)])
                    worst = max(worst, float(big.max()) if big.size and big.size == int(bad.sum()) else math.inf)
                    v = np.where(bad, n, v)
            x.append(v)
        if fixed:
            try:
                k = int(Path(d).name.split("-")[0])
            except ValueError:
                k = None
            if k is not None:
                rec = self.corrected.get(k)
                if rec is None or fixed > rec["points"]:
                    self.corrected[k] = dict(configuration=k, points=fixed, max_wiDiv=worst)
        return t, np.array(x)

    def decide(self, ident, **record):
        """Append a decision; a recorded one with the same ID is returned instead (a production batch is never decided twice)."""
        old = self.decisions()
        if ident in old:
            return old[ident]
        rec = dict(id=ident, **record)
        self._drop_cut_line()
        with open(self.decisions_path, "a") as f:
            f.write(json.dumps(rec) + "\n"); f.flush(); os.fsync(f.fileno())
        return rec

    def _drop_cut_line(self):
        """A cut last line (a crash while writing) is removed, so the next decision is not glued to it (review 10-07 A5)."""
        if self.decisions_path.is_file():
            text = self.decisions_path.read_text()
            if text and not text.endswith("\n"):
                _atomic_write(self.decisions_path, text[:text.rfind("\n") + 1])

    def decisions(self):
        out = {}
        if self.decisions_path.is_file():
            for line in self.decisions_path.read_text().split("\n"):
                try:
                    r = json.loads(line)
                except ValueError:
                    continue                      # a cut last line
                out.setdefault(r["id"], r)
        return out

    # ---------------------------------------------------------------- steps
    def plan(self):
        """Decisions belong to one plan (review 10-07 A4): the input (block and base), the source files and the decision code.
        A workdir whose decisions came from another plan has them moved aside (decisions.<n>.jsonl) before anything is
        decided, so an old production batch or allocation is never mixed with a new one; the CCEX runs stay reusable (their keys
        hold everything they read)."""
        self._drop_cut_line()
        here = Path(__file__).resolve().parent
        other = {}                         # Qubit.file and gyrofile contents too (the runs hashed them, the plan did not)
        for path in INPUT_FILES:
            v = _at(self.base, *path)
            if isinstance(v, str) and (self.root / v if not Path(v).is_absolute() else Path(v)).is_file():
                other[".".join(path)] = _sha(self.root / v if not Path(v).is_absolute() else Path(v))
        files, resolved = self.source_sha, self.resolved
        if "pool" in self.auto:
            # production configurations found: files added later keep the plan - its allocation and the production batches stay
            # (review 20261008-autoconv-resume-bathgen #4, plan "bath 만들기 빼기" 2026-10-09); a changed file starts a new one
            files = {k: v for k, v in self.source_sha.items() if k in self.pilot}
            resolved = {k: v for k, v in self.resolved.items() if k not in ("pool", "found_configurations", "r_cov")}   # review 82c666 R3
        # max_configurations raised later continues the plan like files added (user 2026-10-10)
        block = copy.deepcopy(self.block)
        if isinstance(block.get("configurations"), dict):
            block["configurations"].pop("max_configurations", None)
        fp = hashlib.sha256(_canon(dict(block=block, base=self.base, files=files, other=other, binary=self.binary_sha, resolved=resolved,
                                        code={c: _sha(here / c) for c in CODE},
                                        **(dict(partition_cache=PARTITION_RULE) if self.pcce else {}))).encode()).hexdigest()
        dec = self.decisions()
        old = dec.get("plan")
        # production files seen before must be unchanged for the plan to stay (review 82c666 R2): appended ones are recorded
        seen = {}
        for k_, v_ in dec.items():
            if k_.startswith("production_seen"):
                seen.update({int(a): b for a, b in v_["sha"].items()})
        changed = [k for k, h in seen.items() if self.source_sha.get(k) != h]
        if self.decisions_path.is_file() and (old is None or old.get("fingerprint") != fp or changed):
            n = 1
            while (self.work / f"decisions.{n}.jsonl").exists():
                n += 1
            self.decisions_path.rename(self.work / f"decisions.{n}.jsonl")
            print(f"[auto] the input, the source files or the decision code changed: earlier decisions moved to decisions.{n}.jsonl",
                  file=sys.stderr)
            seen, changed = {}, []         # a new plan records every production file again (review 4d2d3a R1)
        self.decide("plan", fingerprint=fp)
        if "pool" in self.auto:
            new = {str(k): self.source_sha[k] for k in self.pool if k not in seen or changed}
            if new:
                self.decide(f"production_seen_{len(self.decisions())}", sha=new)
        return fp

    def progress(self, step=None, status=None, note=None):
        """result/progress.txt and json/progress.json: where the run is, what it has done and chosen, what comes next and how to
        continue (plan "끊김 대비·진행 보기" R3). Written at every step and after every finished CCEX run."""
        if step is not None:
            self.step_now = step
        if status is not None:
            self.status_now = status
        out = self.work / "result"
        out.mkdir(exist_ok=True)
        runs = {}
        if (self.work / "runs").is_dir():
            for d in sorted((self.work / "runs").iterdir()):
                if d.is_dir():
                    runs[d.name] = sum(1 for x in d.iterdir() if (x / "done.json").is_file())
        dec = self.decisions()
        chosen = {}
        for k in [k for k in dec if k == "window" or k.startswith("window_")]:
            chosen["window"] = dict(nstep=dec[k].get("nstep"), deltat=dec[k].get("deltat"), decay_time=dec[k].get("decay_time"))
        if "truncation" in dec:
            chosen["truncation"] = dict(rbath=dec["truncation"]["rbath"], rdip=dec["truncation"]["rdip"])
        if "allocation" in dec:
            chosen["configurations"] = dict(C=dec["allocation"]["C"], S=dec["allocation"]["S"])
        last = [v for k, v in dec.items() if re.fullmatch(r"production(_w\d+)?_batch_\d+_done", k)]   # every production window
        if last:
            chosen["production"] = dict(computed=len(last[-1]["ids"]), se=last[-1]["se"], needed=last[-1]["need_total"])
        fixed = [sec for sec, on in getattr(self, "converge", {}).items() if not on]
        again = f"main.out -f {self.input_path}"
        rec = dict(updated=time.strftime("%Y-%m-%d %H:%M:%S"), status=self.status_now, step=STEPS[self.step_now], note=note,
                   steps_done=list(STEPS[:self.step_now + (self.status_now == "converged")]), steps_left=list(STEPS[self.step_now + 1:]),
                   fixed_sections=fixed,
                   finished_runs=runs, chosen=chosen, interrupted_before=self.interrupted, continue_with=again)
        _atomic_write(_json(out) / "progress.json", json.dumps(rec, indent=1, default=float))
        rec["round_note"] = note
        lines = [f"updated: {rec['updated']}", f"status: {self.status_now}" + (f" ({note})" if note else ""),
                 f"now: step {STEPS[self.step_now]}" if self.status_now == "running" else f"last step: {STEPS[self.step_now]}",
                 "done: " + (", ".join(rec["steps_done"]) or "-"), "left: " + (", ".join(rec["steps_left"]) or "-")]
        if fixed:
            lines.append("fixed (not converged): " + ", ".join(fixed))
        lines.append("finished CCEX runs: " + (", ".join(f"{k} {v}" for k, v in runs.items()) or "none yet"))
        for k, v in chosen.items():
            if k == "configurations":         # the test's estimate, not the number used (that is the production line)
                lines.append(f"estimated configurations: {v['C']} needed at {v['S']} states")
                continue
            lines.append(f"chosen {k}: " + ", ".join(f"{a} {b:g}" if isinstance(b, (int, float)) else f"{a} {b}" for a, b in v.items()))
        if self.interrupted:
            lines.append(f"the run before was interrupted at step {self.interrupted.get('step')} ({self.interrupted.get('updated')}); "
                         "this run continued it")
        if self.status_now in ("stopped", "interrupted"):
            lines.append(f"to continue: {again}  (finished runs and decisions are reused; resume \"fresh\" starts over)")
        elif self.status_now == "running":
            lines.append(f"if the job is cut: {again} continues from here")
        _atomic_write(out / "progress.txt", "\n".join(lines) + "\n")

    def mark(self, status, reason):
        """result/json/summary.json and summary.txt say a run is going on or ended without a result, so an earlier result is
        never read as the present one (review 10-07 A9)."""
        out = self.work / "result"
        out.mkdir(exist_ok=True)
        _atomic_write(_json(out) / "summary.json", json.dumps(dict(status=status, reason=reason), indent=1))
        _atomic_write(out / "summary.txt", f"status: {status}" + (f" ({reason})" if reason else "") + "\n")
        self.progress(status=status, note=reason)

    def execute(self):
        res = self.work / "result"
        if res.is_dir():                     # an earlier run's step files (a step now fixed would leave them; review ...-autostart #12)
            for f in list(res.glob("[0-6]_*")) + [res / "coherence_mean.txt", res / "summary.json", res / "progress.json"]:   # json: before 2026-10-10
                if f.is_file():
                    f.unlink()
        self.mark("running", "starting: checking (or making) the bath files")   # a cut here is an interruption too (review #5)
        self.progress(step=0)
        self.resolve()
        self.preflight()
        self.mark("running", None)
        self.progress(step=1)
        plan_fp = self.plan()
        # CCEX writes the total evolution time: 2τ for a Hahn echo, 2Nτ for N equally spaced π pulses, t otherwise
        npulse = int(_at(self.base, "pulse", "count", default=0) or 0)
        manual = str(_at(self.base, "pulse", "name", default="")).lower() == "manual"
        report_mod.TIME_LABEL = "t" if manual else "τ" if npulse == 0 else ("2τ" if npulse == 1 else f"{2 * npulse}τ")   # Ramsey: τ (user 2026-10-10)
        o, S0, Bn, Rm = self.o, 0 if self.ensemble else self.o["pilot_states"], self.rbaths[-1], self.rdip_max
        auto_dt, dt0 = self.deltat is None, self.deltat
        # automatic starts need two passes in a row before a value counts as converged (plan "자동 시작값" (7))
        need_rdip = 2 if "candi_rdip" in self.auto else 1
        need_rbath = 2 if "candi_rbath" in self.auto else 1
        summary = dict(version=VERSION, binary=str(self.binary), binary_sha256=self.binary_sha, input=str(self.input_path),
                       test=self.pilot, production_available=self.pool, plan=plan_fp, found_automatically=sorted(self._shown(self.auto)))
        summary["warnings"] = [f"start: {x}" for x in self.start_warnings]
        # window
        given = o["set_time_window"]

        def grid(W):
            """Points and step for a window W (ms): always nstep points (user 2026-10-10: a given deltat sets the first window
            only)."""
            return self.npoints, W / (self.npoints - 1)

        if given is not None:
            nstep, self.deltat = grid(float(given))
            if not auto_dt:
                summary["warnings"].append(f"deltat {dt0:g} ms of the input is not used: the given window on nstep {nstep} points")
        elif auto_dt:
            # the start covers the revival check, at least the production window (revival_check_factor may be 0; review
            # 20261008-autoconv-autostart #10)
            nstep, self.deltat = grid(max(o["revival_check_factor"], o["production_window_factor"]) * self.t_seed)
        else:
            nstep, self.deltat = self.nstep0, dt0
        rbath_raised, window_rounds = [], []
        plateau_seen, grows, shrunk, coarse = None, 0, False, 0
        tried = []                             # every window computed: drawn when no window is found (user 2026-10-10)

        judged = {}                            # the last stable test data (the figure of a stop shows what was judged; R6)

        def window_failed(status, reason, w=None):
            if w is not None:
                summary["window"] = dict(w, tried=[dict(window=float(a[0][-1]), nstep=len(a[0]), rbath=a[2]) for a in tried])
                if judged:                     # the fit points tried on the last window (user 2026-10-10)
                    report_mod.window_fit(self.work / "result", judged["t"], judged["data"], w)
            report_mod.window_tried(self.work / "result", tried)
            return self.finish(summary, status, "window", reason + "; 1_window.png shows the curves of every window tried")
        for outer in range(MAX_WINDOW_ROUNDS):
            self.progress(step=1, note=None if outer == 0 else f"window round {outer + 1}")   # each round starts at the window (#7)
            while given is not None:           # a window given by the user: used as it is, T_d only measured (2026-10-07)
                dirs = [self.run("pilot", k, S0, nstep, Bn, Rm, self.rungs) for k in self.pilot]
                data = [self.curves(d, "out", S0)[1] for d in dirs]
                t = self.curves(dirs[0], "out", S0)[0]
                tried.append((t, data, Bn))
                keep = [not self._flag(k, x, t, "test") for k, x in zip(self.pilot, data)]
                if sum(keep) < 2:
                    return window_failed("numerical-failure", "fewer than 2 stable test configurations "
                                         f"(|L| > {UNSTABLE_L:g} or nan/inf: {sorted(self.unstable)})")
                data = [x for x, m in zip(data, keep) if m]
                stable_pilot = [k for k, m in zip(self.pilot, keep) if m]
                judged.update(t=t, data=data)
                w = self._window(t, data, look=0.0, choose=True)
                w.update(status="given", recommended_window=float(t[-1]), given=float(given), warnings=[])
                if w.get("decay_time") is None:
                    w["decay_time"] = float(t[-1])
                    w["warnings"].append(f"no T2 fit inside the given window ({w.get('reason')}); truncation judged on the whole window")
                elif w["decay_time"] > t[-1]:
                    w["warnings"].append(f"the curve has not decayed inside the given window: T_d {w['decay_time']:.4g} ms > {t[-1]:.4g} ms")
                    w["decay_time"] = float(t[-1])
                elif o["production_window_factor"] * w["decay_time"] > t[-1]:
                    w["warnings"].append(f"note: {o['production_window_factor']:g} x T_d = {o['production_window_factor'] * w['decay_time']:.4g} ms is longer than the given window")
                break
            while given is None:
                dirs = [self.run("pilot", k, S0, nstep, Bn, Rm, self.rungs) for k in self.pilot]
                data = [self.curves(d, "out", S0)[1] for d in dirs]
                t = self.curves(dirs[0], "out", S0)[0]
                tried.append((t, data, Bn))
                keep = [not self._flag(k, x, t, "test") for k, x in zip(self.pilot, data)]
                if sum(keep) < 2:
                    return window_failed("numerical-failure", "fewer than 2 stable test configurations "
                                         f"(|L| > {UNSTABLE_L:g} or nan/inf: {sorted(self.unstable)})")
                data = [x for x, m in zip(data, keep) if m]
                stable_pilot = [k for k, m in zip(self.pilot, keep) if m]
                judged.update(t=t, data=data)
                w = self._window(t, data, choose=True)
                W = float(t[-1])
                need_w = w.get("needed_window") or w.get("recommended_window")
                if w["status"] == "resolved":
                    # deltat found: a start far too long leaves few points on the decay, so the window is cut once to what
                    # the checks need (plan "자동 시작값" (5); both window methods, review 20261008-autoconv-autostart #5)
                    if not shrunk and need_w and W > 2 * need_w:
                        shrunk = True
                        nstep, self.deltat = grid(1.05 * need_w)
                        continue
                    break
                # deltat found: the mean already below the fit floor within the first points is a time step too coarse, not a
                # window too short - the window is cut to the decay instead of grown (review 20261008-autoconv-autostart #6)
                mu_abs = np.abs(np.mean([x.mean(0) for x in data], 0))
                low = np.nonzero(mu_abs[1:] < self.fit["floor"])[0]
                if len(low) and low[0] + 1 <= COARSE_POINTS and coarse < MAX_COARSE:
                    coarse += 1
                    cut = float(t[low[0] + 1])
                    nstep, self.deltat = grid(max(o["revival_check_factor"], o["production_window_factor"]) * 1.5 * cut)
                    window_rounds.append(dict(coarse_grid=True, window=W, decayed_by=cut, new_window=float(self.deltat * (nstep - 1))))
                    continue
                if w["status"] == "plateau" and plateau_seen != (Bn, nstep, self.deltat) and \
                        grows < MAX_GROWS:
                    # review 10-07 A3: a plateau is believed only when it is still there in a window 1.5 x longer
                    grows += 1
                    nstep, self.deltat = grid(1.5 * W)
                    plateau_seen = (Bn, nstep, self.deltat)
                    continue
                if w["status"] == "plateau" and not self.converge["rbath"]:
                    return window_failed("window-unresolved",
                                         f"the coherence stays near {w['plateau_level']:.3g} (a plateau) at the fixed rbath {Bn:g} A: "
                                         "a larger rbath may resolve it (rbath.converge true, or a larger rbath.value)", w)
                if w["status"] == "plateau":     # too small a bath: raise rbath by a full 1.25 step, within the bath files
                    new = round(Bn * 1.25, 3)
                    if new > self.max_rbath * (1 + 1e-9):
                        unchecked = "" if plateau_seen == (Bn, nstep, self.deltat) else \
                            " (not checked in a longer window: the pilot is at its longest window, a longer one may also resolve it)"
                        return window_failed("window-unresolved",
                                             f"the coherence stays near {w['plateau_level']:.3g} (a plateau) and rbath {Bn:g} A x 1.25 is beyond "
                                             f"what the bath files cover ({self.max_rbath:.4g} A): make larger baths, covering at least {new:.4g} A "
                                             f"around the qubit{unchecked}", w)
                    self.rbaths.append(new)
                    # the rdip top goes up with it (a small rdip also holds the curve up; plan "자동 시작값" (8a))
                    newR = min(round(Rm * 1.25, 3), round(2 * new, 3))
                    if self.converge["rdip"] and newR > Rm * (1 + 1e-9) and _gname(newR) not in [_gname(r) for r in self.rungs + [Rm]]:
                        self.rungs, Rm = self.rungs + [Rm], newR
                        self.rdip_max = Rm
                    self.decide(f"rbath_raise_{len(rbath_raised) + 1}", rbath=new, rdip_max=Rm, plateau_level=w["plateau_level"])
                    rbath_raised.append(dict(rbath=new, rdip_max=Rm, plateau_level=w["plateau_level"]))
                    Bn = new
                    continue
                # longer: 1.5 x, or what the fit asks for; at most MAX_GROWS times in either mode (max_nstep removed 2026-10-08)
                want = max(1.5 * W, w.get("needed_window") or 0.0,
                           o["production_window_factor"] * w["decay_time"] if w.get("decay_time") is not None else 0.0)
                if grows >= MAX_GROWS:
                    # the last window's own reason too: "not decayed" was said of curves that fell at once (h-BN 500 G, 2026-10-10)
                    return window_failed("window-unresolved",
                                         f"no time window found after {grows} longer windows (up to {W:.4g} ms; the last: "
                                         f"{w.get('reason') or w.get('status')})", w)
                grows += 1
                nstep, self.deltat = grid(want)
            td = w["decay_time"]
            Bw, Rw = Bn, Rm                  # the reference the window was found at
            if outer == 0:
                summary["warnings"].extend(over_one("pilot", stable_pilot, data, t))
            summary["rbath_raised_for_the_window"] = rbath_raised if given is None else []
            summary["max_rbath"] = self.max_rbath
            rec_w = self.decide("window" if outer == 0 else f"window_{outer + 1}", nstep=nstep, deltat=self.deltat, decay_time=td,
                                fit_floor=w.get("fit_floor", self.fit["floor"]), points=self.fit["points"])
            # the later fits of this run stop at the window's noise floor too (user 2026-10-10); never below t2_fit.floor
            self.fit["floor"] = max(float(self.o["fit_floor"]), float(rec_w.get("fit_floor", self.o["fit_floor"])))
            summary["window"] = dict(w, pilot_nstep=nstep, deltat=self.deltat, rounds=window_rounds)
            res_dir = self.work / "result"
            report_mod.window_fit(res_dir, t, data, w)      # always (review R8)
            report_mod.window_step(res_dir, t, data, dict(w, pilot_nstep=nstep, rbath=Bn, rbath_raised=summary["rbath_raised_for_the_window"],
                                                          max_rbath=self.max_rbath))
            sigma = td / 25
            opts = dict(sigma=sigma, smooth=False, **self.fit)
            tol_r, tol_b = o["threshold_rdip"], o["threshold_rbath"]

            def judge_curves(lo_curves, hi_curves, rel_tol):
                # configurations unstable on either side are left out of this comparison, the pairs kept (2026-10-10)
                ok = np.isfinite(lo_curves).all(1) & np.isfinite(hi_curves).all(1)
                if ok.sum() < 2:                              # review 10-07 A8: never judged
                    return dict(verdict="nonfinite", rel=None)
                lo_curves, hi_curves = lo_curves[ok], hi_curves[ok]
                try:
                    j = truncation.judge(t, lo_curves, hi_curves, td, rel_tol=rel_tol, **opts)
                    if j["route"] == "raw-area":              # a fit failed: judged on the raw area (user 2026-10-12)
                        return dict(verdict=j["verdict"], rel=abs(j["raw_area_rel"]), basis="raw-area")
                    return dict(verdict=j["verdict"], rel=None if j["area_rel"] is None else abs(j["area_rel"]))
                except (ValueError, RuntimeError):
                    return dict(verdict="fit-failure", rel=None)

            self.progress(step=2)
            # rdip, at the largest rbath, on the pilot runs (the ladder gives every rung from one run per configuration)
            state = dict(dirs=dirs, top=Rm, b=Bn)
            judge_r = lambda a, b: judge_curves(a, b, tol_r)
            judge_b = lambda a, b: judge_curves(a, b, tol_b)
            def rdip_curve(r):
                pre = "out" if r == state["top"] else f"out_rdip{_gname(r)}"
                return self._stable_means(state["dirs"], pre, S0, t, "test")
            def rdip_extend(r):
                ladder = [x for x in self.rungs + [state["top"]] if x < r]
                state["dirs"] = [self.run("pilot", k, S0, nstep, state["b"], r, ladder) for k in self.pilot]
                self.rungs, state["top"] = ladder, r
            prefix = lambda r: "out" if r == Rm else f"out_rdip{_gname(r)}"
            absmean = lambda c: np.abs(c[np.isfinite(c).all(1)].mean(0))
            if self.converge["rdip"]:
                rd = converge(self.rungs + [Rm], rdip_curve, rdip_extend, 2 * Bn, True, judge_r, tol_r, need=need_rdip)
                Rm, dirs, self.rdip_max = state["top"], state["dirs"], state["top"]
                rd["rbath"] = Bn
                report_mod.trend_step(res_dir, 2, "rdip", t, {r: absmean(rdip_curve(r)) for r in rd["values"]}, rd, rd["chosen"],
                                      w["recommended_window"], tol_r, fixed=f"rbath {Bn:g} Å")
            else:                              # rdip.converge false: the value given, used as it is
                rd = dict(estimate=None, raised=[], checks=[], values=[Rm], reference=Rm, vs_reference=[], chosen=Rm, fixed=True)
            R = rd["chosen"]
            # rbath, at the rdip chosen (runs at the largest rdip with the ladder, so the chosen rdip is in every run)
            runs = {Bn: dirs}
            rb = None
            rdip_reason = f"rdip has not converged up to {Rm:g} (2 x rbath holds every pair)"
            self.progress(step=3)
            if R is not None and not self.converge["rbath"]:     # rbath.converge false: the value given, used as it is
                rb = dict(estimate=None, raised=[], checks=[], values=[Bn], reference=Bn, vs_reference=[], chosen=Bn, fixed=True)
            elif R is not None:
                def rbath_run(b):
                    runs[b] = [self.run("rbath", k, S0, nstep, b, Rm, self.rungs) for k in self.pilot]
                for b in self.rbaths[:-1]:
                    rbath_run(b)
                def rbath_curve(b):
                    return self._stable_means(runs[b], prefix(R), S0, t, "test")
                # rbath at the rdip chosen; then rdip again at the final largest rbath (review 10-07 C3); when it does not hold
                # there, rdip is converged again at that rbath, the top raised if needed, and rbath judged once more
                # (plan "자동 시작값" (8b))
                for attempt in range(3):
                    rb = converge(list(self.rbaths), rbath_curve, rbath_run, self.max_rbath, False, judge_b, tol_b,
                                  need=need_rbath)
                    self.rbaths = rb["values"]
                    if rb["chosen"] is None:
                        break
                    if not self.converge["rdip"]:     # a fixed rdip is not judged again at the final rbath
                        rd_final = dict(chosen=R)
                        break
                    Bf = self.rbaths[-1]
                    top = runs[Bf]
                    cur = lambda r: self._stable_means(top, prefix(r), S0, t, "test")
                    rv = self.rungs + [Rm]
                    vs = [dict(value=r, **judge_r(cur(r), cur(Rm))) for r in rv[:-1]]
                    top_ok = [judge_r(cur(rv[-2 - i]), cur(rv[-1 - i])) for i in range(min(need_rdip, len(rv) - 1))]
                    settled = len(top_ok) == need_rdip and all(x["verdict"] == "pass" for x in top_ok)
                    R2 = suffix(rv, vs) if settled and vs and vs[-1]["verdict"] == "pass" else None   # review ...-autostart #2
                    rd_final = dict(rbath=Bf, reference=Rm, vs_reference=vs, top_neighbours=[x["verdict"] for x in top_ok], chosen=R2)
                    summary.setdefault("rdip_at_final_rbath", []).append(rd_final)
                    if R2 is not None and R2 <= R:
                        break
                    if R2 is None and attempt < 2:
                        state.update(dirs=top, top=Rm, b=Bf)
                        rd2 = converge(rv, rdip_curve, rdip_extend, 2 * Bf, True, judge_r, tol_r, need=need_rdip)
                        rd_final["reconverged"] = {k: v for k, v in rd2.items() if k != "values"} | dict(values=rd2["values"])
                        if rd2["chosen"] is None:
                            break
                        R2 = rd2["chosen"]
                        if state["top"] != Rm:          # a new top: every rbath run again with it in the ladder
                            Rm = self.rdip_max = state["top"]
                            runs[Bf] = state["dirs"]
                            for b in self.rbaths[:-1]:
                                rbath_run(b)
                        rd_final["chosen"] = R2
                    if R2 is None:
                        break
                    R = R2
                else:
                    rd_final["chosen"] = None
                rb["rdip"] = R
                report_mod.trend_step(res_dir, 3, "rbath", t, {b: absmean(rbath_curve(b)) for b in rb["values"]}, rb, rb["chosen"],
                                      w["recommended_window"], tol_b, fixed=f"rdip {R:g} Å")
                if rb["chosen"] is not None and rd_final["chosen"] is None:      # after the plot (review ...-autostart #4)
                    R = None
                    rdip_reason = f"rdip at the final rbath {self.rbaths[-1]:g} A: {rv[-2]:g} vs {Rm:g} does not agree"
            # the window again at the final reference when it grew by a step or more (plan "자동 시작값" (8c))
            Bf = self.rbaths[-1]
            if given is not None or not (self.converge["rdip"] and self.converge["rbath"]) or R is None or rb is None or rb["chosen"] is None or \
                    (Bf < 1.25 * Bw * (1 - 1e-4) and Rm < 1.25 * Rw * (1 - 1e-4)):     # radii carry 3 decimals (review ...-autostart #3)
                break
            fin = [self.curves(d, "out", S0)[1] for d in runs[Bf]]
            w2 = self._window(t, fin)
            td2 = w2.get("decay_time")
            moved = td2 is None or not (0.8 * td <= td2 <= 1.25 * td)
            window_rounds.append(dict(reference_rbath=Bf, reference_rdip=Rm, decay_time=td2, previous_decay_time=td, moved=moved))
            if not moved:
                break
            if outer == MAX_WINDOW_ROUNDS - 1:
                summary["window"]["rounds"] = window_rounds
                return self.finish(summary, "window-unresolved", "window",
                                   f"the decay time keeps changing with the reference ({td:.4g} -> {td2} ms at rbath {Bf:g}, rdip {Rm:g})")
            Bn = Bf
            target = max(o["revival_check_factor"], o["production_window_factor"]) * td2 if td2 is not None else 1.5 * float(t[-1])
            nstep, self.deltat = grid(target)
            plateau_seen, shrunk = None, True
        B = rb["chosen"] if rb else None
        Bn = self.rbaths[-1]
        clean = lambda r: None if r is None else {k: v for k, v in r.items() if k != "values"} | dict(values=r["values"])
        summary["rdip"], summary["rbath"] = clean(rd), clean(rb)
        summary["converge"] = dict(self.converge)
        summary["truncation"] = dict(rule=f"relative area of the T2 fit within threshold_rdip {100 * tol_r:g} % / threshold_rbath {100 * tol_b:g} %: "
                                     "trend + one check, else raise until the last "
                                     f"{'two' if max(need_rdip, need_rbath) == 1 else 'three'} agree", sigma=sigma)
        if R is None:
            return self.finish(summary, "truncation-unresolved", "rdip", unjudged_reason("rdip", rd, "2_rdip.png") or rdip_reason)
        if B is None:
            return self.finish(summary, "truncation-unresolved", "rbath", unjudged_reason("rbath", rb, "3_rbath.png") or
                               f"rbath has not converged up to {Bn:g} A (the bath files cover up to {self.max_rbath:.4g} A): make larger baths, "
                               f"covering at least {round(Bn * 1.25, 3):g} A around the qubit")
        pairs_used = sum(self.pair_count(self.inside(k, B), R) for k in self.pilot)
        # the comparisons judged on the raw area because a T2 fit failed (user 2026-10-12), recorded only when there are any
        raw = sorted({("rdip", *x) for x in raw_area_pairs(rd)} | {("rbath", *x) for x in raw_area_pairs(rb)}
                     | {("rdip", *x) for f in summary.get("rdip_at_final_rbath", [])
                        for x in raw_area_pairs(f) | raw_area_pairs(f.get("reconverged"))})
        self.decide("truncation", rbath=B, rdip=R, **({} if not raw else dict(
            raw_area=[dict(parameter=n, lo=lo, hi=hi, verdict=v) for n, lo, hi, v in raw])))
        summary["truncation"].update(chosen=dict(rbath=B, rdip=R, pairs=pairs_used))
        if self.fixed_config is not None:      # configurations.converge false: the number (and states) given
            C, S = self.fixed_config
            summary["sampling_pilot"] = dict(fixed=dict(configurations=C, states=S))
        else:
            C, S = None, None
        self.progress(step=4)
        if C is None:                          # sampling allocation from the pilot at the chosen setting
            pk = [(k, x) for k, x in ((k, self.curves(d, prefix(R), S0)[1]) for k, d in zip(self.pilot, runs[B]))
                  if not self._flag(k, x, t, "test")]
            pdata = [x for _, x in pk]
            if len(pdata) < 2:
                return self.finish(summary, "numerical-failure", "sampling", f"fewer than 2 stable test configurations ({sorted(self.unstable)})")
            if not self.ensemble:              # a check of the states, never changes the production (user 2026-10-10)
                summary["states"] = st = self._state_test(t, nstep, B, R, pk)
                report_mod.growth(res_dir, "4_states", t, st.pop("_curves"), [(p["states"], p["se"]) for p in st["points"]],
                                  o["states_threshold"], w["recommended_window"], "states",
                                  f"4. states: test {len(st['configurations'])} configurations\nrbath {B:g} Å, rdip {R:g} Å",
                                  chosen=st["chosen"], auto_y=True)
                if st["warning"]:
                    summary["warnings"].append(st["warning"])
            try:
                tc = t2mod.components(t, pdata, **self.fit)
                binfo = {}
                boot = t2mod.bootstrap_se_rel(t, pdata, **self.fit, info=binfo)
            except (ValueError, RuntimeError) as e:
                return self.finish(summary, "sampling-unresolved", "sampling", f"T2 of the pilot: {e}")
            lin = t2mod.se_rel(tc, tc["C"], tc["S"])
            if not (math.isfinite(lin) and lin > 0):
                # review 10-07 B3: a symmetric sample can make every linearized influence cancel (lin 0) while the bootstrap does
                # not; an allocation scaled from 0 would ask for nothing
                summary["sampling_pilot"] = dict(A=tc["A"], T2=tc["T2"], n=tc["n"], se_rel_linearized=lin, se_rel_bootstrap=boot)
                return self.finish(summary, "sampling-unresolved", "sampling",
                                   f"the linearized SE of the pilot is {lin:g} (bootstrap {boot:.3g}): no allocation can be scaled from it")
            corr = max(1.0, boot / lin)
            # the states stay those of the test (user 2026-10-10: the configurations converge at the initial state count)
            sc = (1,) if self.ensemble else (S0,)
            alloc = t2mod.allocate(tc, o["threshold_config"], 1.0, 1.0, s_candidates=sc,      # a configuration and a state cost alike
                                   max_configurations=self.cap, correction=corr,
                                   min_configurations=o["min_configurations"])
            summary["sampling_pilot"] = dict(A=tc["A"], T2=tc["T2"], n=tc["n"], se_rel_linearized=lin, se_rel_bootstrap=boot, correction=corr,
                                             allocation=alloc, **binfo)
            if binfo.get("bootstrap_failed"):
                summary["warnings"].append(f"test: {binfo['bootstrap_failed']} of {binfo['bootstrap_n']} bootstrap resamples could not be "
                                           "fitted and were left out of the sampling error")
            if alloc["best"] is None:
                return self.finish(summary, "sampling-unresolved", "sampling", f"allocation {alloc['status']}")
            # a pool too small for the target is computed anyway, a batch at a time (plan "조금씩 늘리기", user 2026-10-09);
            # the states stay those of the allocation, C is the pilot's estimate of the configurations needed
            S = alloc["best"]["S"]
            need = next(r["C"] for r in alloc["table"] if r["S"] == S)
            rec_alloc = self.decide("allocation", C=need, S=S)    # a recorded allocation stays the one used (82c666-b #4)
            C, S = rec_alloc["C"], rec_alloc["S"]
        if self.ensemble:
            S = 0                       # one ensemble curve per configuration (CCEX nstate 0)
        self.progress(step=5)
        # production
        if given is None:                 # 1.5 x T_d / revival end, on the same grid rule as the pilot
            nprod, self.deltat = grid(w["recommended_window"])
        else:
            nprod = nstep
        # the production window (user 2026-10-10, P1 CR90): when the production curve has not decayed inside it, the
        # production is computed again on a longer window (the same configurations first), at most MAX_PROD_WINDOWS times
        prod_windows, start = [], None
        for pr in range(MAX_PROD_WINDOWS + 1):
            tag = "" if pr == 0 else f"_w{pr + 1}"
            if pr:
                rec_p = self.decide(f"production_window{tag}", nstep=nprod, deltat=self.deltat)
                nprod, self.deltat = rec_p["nstep"], rec_p["deltat"]
            if self.fixed_config is not None:      # a number given: exactly those, its SE reported, not judged
                ids = self.pool[:C]
                prod = [self.run("production", k, S, nprod, B, R) for k in ids]
                res, batches = self._measure(prod, S, ids), []
                if res["status"] == "se":
                    res.update(status="ok", reason=None)
            else:
                ids, prod, res, batches = self._production(S, nprod, B, R, C, tag=tag, start=start)
            tp = self.curves(prod[0], "out", S)[0]
            # the unstable configurations are out of the mean, the reports and the checks below (2026-10-10)
            kept = [(k, d) for k, d in zip(ids, prod) if k not in self.unstable]
            ids_ok, prod_ok = [k for k, _ in kept], [d for _, d in kept]
            x = [self.curves(d, "out", S)[1] for d in prod_ok]
            if res["status"] not in ("ok", "se") or len(x) < 2:
                break
            # production check: the production curve must decay inside the window computed
            wp = self._window(tp, x, look=0.0)
            if wp.get("decay_time") is not None:
                wp["status"] = "resolved" if wp["decay_time"] <= tp[-1] else "needs-extension"
            prod_windows.append(dict(window=float(tp[-1]), nstep=nprod, status=wp["status"], decay_time=wp.get("decay_time")))
            if given is not None or wp["status"] == "resolved" or pr == MAX_PROD_WINDOWS:
                break
            want = max(1.5 * float(tp[-1]), wp.get("needed_window") or 0.0,
                       o["production_window_factor"] * wp["decay_time"] if wp.get("decay_time") is not None else 0.0)
            nprod, self.deltat = grid(want)
            start = len(ids)
        summary["production"] = dict(res, configurations=ids, states=S, nstep=nprod, batches=batches)
        if len(prod_windows) > 1:
            summary["production"]["windows"] = prod_windows
        if res.get("bootstrap_failed"):
            summary["warnings"].append(f"production: {res['bootstrap_failed']} of {res['bootstrap_n']} bootstrap resamples could not be "
                                       "fitted and were left out of the sampling error")
        summary["warnings"].extend(over_one("production", ids_ok, x, tp))
        mu = np.mean([v.mean(0) for v in x], 0)
        (self.work / "result").mkdir(exist_ok=True)
        _atomic_write(self.work / "result" / "coherence_mean.txt", "# t(ms)  Re  Im  |L|  (production, configurations weighted equally)\n" +
                      "".join(f"{a:.8g} {b.real:.12e} {b.imag:.12e} {abs(b):.12e}\n" for a, b in zip(tp, mu)))
        if res["status"] not in ("ok", "se") or len(x) < 2:   # a numerical failure first; a sample short of the target is reported last
            return self.finish(summary, "numerical-failure", "production", res.get("reason") or "fewer than 2 stable production configurations")
        # 5_configurations.png (user 2026-10-10): the production's configurations growing, their mean curves and measured errors
        summary["production"]["decomposition"] = dec_rows = self._decompose(tp, x, S)
        each = "" if self.ensemble else f", {S} state{'s' if S != 1 else ''} each"
        report_mod.growth(res_dir, "5_configurations", tp, {c: np.abs(np.mean([v.mean(0) for v in x[:c]], 0))
                                                      for c in [1] + [r["configurations"] for r in dec_rows]},
                          [(r["configurations"], r["se"]) for r in dec_rows], o["threshold_config"], float(tp[-1]), "configurations",
                          f"5. configurations: production{each}\nrbath {B:g} Å, rdip {R:g} Å",
                          chosen=dec_rows[-1]["configurations"] if res["status"] == "ok" else None)
        summary["production"]["window"] = wp
        report_mod.production_step(res_dir, tp, x, res.get("T2"), res.get("n"), ensemble=self.ensemble,
                                   fixed=f"rbath {B:g} Å, rdip {R:g} Å")
        if given is not None and wp["status"] != "resolved":          # a given window: a warning, not a stop
            summary.setdefault("warnings", []).append("the production curve does not decay inside the given window")
        elif wp["status"] != "resolved":
            return self.finish(summary, "window-unresolved", "production",
                               f"the production curve does not decay inside the window after {len(prod_windows)} production windows "
                               f"(up to {tp[-1]:.4g} ms)")
        # combination recheck on the first production configurations
        if (B, R) != (Bn, Rm):
            ck = ids_ok[:RECHECK]
            tdr = min(wp["decay_time"], float(tp[-1])) if wp.get("decay_time") is not None else float(tp[-1])
            refd = [self.run("recheck", k, S, nprod, Bn, Rm, [R] if R != Rm else []) for k in ck]
            lo_c = np.array([self.curves(d, "out", S)[1].mean(0) for d in prod_ok[:len(ck)]])
            hi_c = np.array([self.curves(d, "out", S)[1].mean(0) for d in refd])
            if np.isfinite(lo_c).all() and np.isfinite(hi_c).all():
                jr = truncation.judge(tp, lo_c, hi_c, tdr, rel_tol=max(tol_r, tol_b), **dict(opts, sigma=tdr / 25))
            else:                                   # review 10-07 A8: non-finite curves are never judged
                jr = dict(verdict="nonfinite", area_rel=None, _curves=None, settings=None)
            jr.pop("_curves"); jr.pop("settings")
            summary["recheck"] = dict(jr, recheck_ids=ck)
            if jr["verdict"] != "pass":     # a warning, not a stop (user 2026-10-07): the rdip and rbath errors add, and the
                # difference is reported next to the sampling SE it is compared with
                raw = jr.get("route") == "raw-area"
                rel = jr.get("raw_area_rel") if raw else jr.get("area_rel")
                summary.setdefault("warnings", []).append(
                    f"combination recheck {jr['verdict']}: rbath {B:g}, rdip {R:g} vs rbath {Bn:g}, rdip {Rm:g} on configurations {ck}"
                    + (f": {'raw area of |L| on [0, T_d] (a T2 fit failed)' if raw else 'area'} differs by {100 * abs(rel):.2f} %"
                       if rel is not None else ""))
        if res["status"] == "se":
            return self.finish(summary, "sampling-unresolved", "production", res["reason"])
        return self.finish(summary, "converged", None, None)

    def _production(self, S, nstep, B, R, n_pilot, tag="", start=None):
        """The production (plan "넉넉하게 먼저, 나눠 보기, 모자라면 추가", user 2026-10-09): first MARGIN x the pilot's estimate
        n_pilot (at least max(C_FLOOR, min_configurations), rounded up to ROUND_TO, at most the pool), judged on everything
        computed; when short, MARGIN x the shortage the measured SE implies (at least C_FLOOR) more, until the estimated SE meets
        the target or the pool is used up. Each block's configurations are recorded before they run and its verdict after, so
        a cut run finishes the same block; a pool used up goes on from the last verdict when the user has added files.
        tag: the production window round ("" the first, "_w2", ... after the window was lengthened, user 2026-10-10), its
        blocks recorded apart; start: the configurations of the first block of a later round (those computed before)."""
        tau, cmin = self.o["threshold_config"], max(C_FLOOR, self.o["min_configurations"])
        ids, prod, batches, res = [], [], [], None
        planned = sum(len(v["ids"]) for k_, v in self.decisions().items() if k_.startswith(f"production{tag}_batch_") and k_.endswith("_plan"))
        if planned > self.cap:                # lowered since (review 01ecda #1): the recorded blocks are not cut
            raise Refused(f"configurations.max_configurations {self.cap} is below the {planned} production configurations already "
                          f"planned in this workdir: raise it to {planned} or more, or start over with \"resume\": \"fresh\"")
        target = min(start, self.cap) if start else first_target(n_pilot, cmin, self.cap)
        for i in itertools.count():
            dec = self.decisions()
            plan = dec.get(f"production{tag}_batch_{i}_plan")
            if plan is None:
                plan = self.decide(f"production{tag}_batch_{i}_plan", ids=self.pool[len(ids):target], states=S, nstep=nstep, rbath=B, rdip=R,
                                   target=target, rule=f"margin {MARGIN:g}, floor {C_FLOOR}, rounded to {ROUND_TO}",
                                   **(dict(n_pilot=n_pilot) if i == 0 else {}))
            if (plan["states"], plan["nstep"], plan["rbath"], plan["rdip"]) != (S, nstep, B, R) or set(plan["ids"]) & set(ids):
                raise RuntimeError(f"production block {i + 1}: the recorded block {plan} does not continue this production")
            prod += [self.run("production", k, S, nstep, B, R) for k in plan["ids"]]
            ids += plan["ids"]
            done, res = dec.get(f"production{tag}_batch_{i}_done"), None
            if done is None:
                res = self._measure(prod, S, ids)
                se, n = res["se"], res.get("n_used", len(ids))        # the stable configurations count (2026-10-10)
                if se is None:                        # non-finite values or a failed fit: more configurations do not mend it
                    need, want, stop = None, None, "failure"
                elif se <= tau and n >= cmin:
                    need, want, stop = n, n, "met"
                else:
                    need, want, _ = next_target(n, se, tau, cmin, self.cap)
                    stop = None
                done = self.decide(f"production{tag}_batch_{i}_done", ids=list(ids), se=se, A=res.get("A"), need_total=need, want=want, stop=stop,
                                   n_used=n)
            batches.append(dict(configurations=len(done["ids"]), stable=done.get("n_used", len(done["ids"])), se=done["se"],
                                need_total=done["need_total"], want=done["want"]))
            # the wanted stable count, in computed configurations: as many more as stable ones are missing
            target = min(self.cap, len(ids) + done["want"] - done.get("n_used", len(ids))) if done["stop"] is None else 0
            if target <= len(ids):                    # met, failed, or the pool used up (files added later continue here)
                break
        if res is None:                               # the last verdict was recorded by an earlier call
            res = self._measure(prod, S, ids)
        n = res.get("n_used", len(ids))
        if res["status"] in ("ok", "se") and (res["status"] == "se" or n < cmin):
            need, want = batches[-1]["need_total"], batches[-1]["want"]
            res["status"] = "se"
            res["reason"] = (f"only {n} stable production configurations, the minimum is {cmin}" if res["se"] <= tau else
                             f"SE of the area {res['se']:.2%} > target {tau:.2%} with {n} stable configurations (the production holds "
                             f"{self.cap}): about {need} needed ({need - n} more)")
            # what to add (plan "bath 만들기 빼기": suggested, not made): the next block with its margin; files already there
            # beyond max_configurations: raise it first
            if len(self.pool) > self.cap:
                res["reason"] += (f"; raise configurations.max_configurations (now {self.cap}, {len(self.pool)} production "
                                  f"files there) to about {self.cap + want - n} and call again")
                return ids, prod, res, batches
            if self.o["max_configurations"] is not None:     # the files and the cap both used up (review 01ecda #3)
                res["reason"] += f"; raise configurations.max_configurations to about {self.cap + want - n} (or leave it out)"
            res["reason"] += (f"; add {want - n} more bath files ({self.patterns['bathfile']}"
                              + (f" and {S} state files each" if not self.ensemble and not self.auto_states else "")
                              + (") and call again to continue" if "pool" in self.auto else ") to configurations.production and call again"))
        return ids, prod, res, batches

    def _state_test(self, t, nstep, B, R, pk):
        """The states on the test configurations (plan "state 시험을 통계 오차로", user 2026-10-10): the states' share of the
        sampling error of A (the configurations kept, the states drawn again; max of the linearized sqrt(W/(C S))/A and a
        bootstrap over each configuration's states) at test_states; above states_threshold, STATE_STEP more states (each block
        computed alone at the chosen rbath and rdip) and measured again, up to max_states or the state files of the test
        configurations. The S chosen is the smallest from which every point computed is within the threshold. A check only:
        an S above test_states or none is a warning; the production keeps test_states."""
        S0, tol = self.o["pilot_states"], self.o["states_threshold"]
        top = min([self.o["max_states"]] + [self.count_states(k) for k, _ in pk])
        data = {k: x for k, x in pk}
        err = {}
        def measure(S):
            sub = [data[k][:S] for k in sorted(data)]
            try:
                tc = t2mod.components(t, sub, **self.fit)
                binfo = {}
                lin, boot = t2mod.se_rel_states(tc), t2mod.bootstrap_se_rel_states(t, sub, **self.fit, info=binfo)
                err[S] = dict(states=S, se=max(lin, boot), se_rel_linearized=lin, se_rel_bootstrap=boot, A=tc["A"], **binfo)
            except (ValueError, RuntimeError) as e:
                err[S] = dict(states=S, se=None, reason=str(e))
        for S in range(2, S0 + 1):                # the states already computed (one state: no spread to measure)
            measure(S)
        while len(data) >= 2 and (err[max(err)]["se"] is None or err[max(err)]["se"] > tol) and max(err) < top:
            s0, n = max(err), min(STATE_STEP, top - max(err))
            for k in list(data):
                x = self.curves(self.run("states", k, n, nstep, B, R, first=s0 + 1), "out", n)[1]
                if self._flag(k, x, t, "states test", first=s0 + 1):
                    data.pop(k)                   # unstable in a new block: out of every S
                else:
                    data[k] = np.concatenate([data[k], x])
            measure(s0 + n)
        values = sorted(err)
        ok = [err[v]["se"] is not None and err[v]["se"] <= tol for v in values]
        i = len(ok)
        while i > 0 and ok[i - 1]:
            i -= 1
        rec = dict(values=values, points=[err[v] for v in values], threshold=tol, max=top, configurations=sorted(data),
                   chosen=values[i] if i < len(values) else None, warning=None)
        if len(data) < 2:
            rec["warning"] = "state test: not done (fewer than 2 stable test configurations)"
        elif rec["chosen"] is None:
            rec["warning"] = (f"state test: the states' share of the sampling error is above {tol:.0%} up to {values[-1]} states; "
                              f"the production uses test_states {S0}")
        elif rec["chosen"] > S0:
            rec["warning"] = (f"state test: the states converge at {rec['chosen']}, more than test_states {S0} the production uses "
                              f"(raise configurations.test_states to {rec['chosen']} to use them)")
        last = max(values) if values else S0
        rec["_curves"] = {v: np.abs(np.mean([data[k][:v].mean(0) for k in data], 0)) for v in [*range(1, S0 + 1), *[v for v in values if v > S0]]
                          if v <= last} if data else {}
        return rec

    def _decompose(self, t, x, S):
        """The computed production growing (5_configurations, user 2026-10-10): the first c configurations, c in 2, 5, 10, ..., all
        (2, 3, ... up to 10), each fitted again from its own mean (same grid, fit and bootstrap); the last is the verdict."""
        rows = []
        for c in grow_counts(len(x), 2):
            sub = x[:c]
            try:
                tc = t2mod.components(t, sub, **self.fit)
                binfo = {}
                boot = t2mod.bootstrap_se_rel(t, sub, **self.fit, info=binfo)
                lin = t2mod.se_rel(tc, tc["C"], tc["S"])
                rows.append(dict(A=tc["A"], se_rel_linearized=lin, se_rel_bootstrap=boot, se=max(lin, boot), status="ok", configurations=c, **binfo))
            except (ValueError, RuntimeError) as e:
                rows.append(dict(A=None, se=None, status="unavailable", reason=str(e), configurations=c))
        return rows

    def _window(self, t, data, look=None, choose=False):
        """choose: the window search (the floor raised to the noise of this mean), else the run's floor."""
        o = self.o
        # the decay time from the T2 fit; the plateau and revival checks on a 9-point running median (user 2026-10-08: method
        # and smoothing_points removed; 9 ignores spikes of up to 4 points and keeps longer revivals, 6f6de5c97)
        return window.fit_window(t, data, o["threshold_time_window"], o["production_window_factor"], median_points=SMOOTHING_POINTS,
                                 look=o["revival_check_factor"] if look is None else look,
                                 raise_floor=choose, **self.fit)    # after the choice the run's floor stays (review R4)

    def _measure(self, prod, S, ids):
        t = self.curves(prod[0], "out", S)[0]
        x_all = [self.curves(d, "out", S)[1] for d in prod]
        keep = [not self._flag(k, v, t, "production") for k, v in zip(ids, x_all)]
        bad = [k for k, m in zip(ids, keep) if not m]
        x = [v for v, m in zip(x_all, keep) if m]
        base = dict(n_used=len(x), unstable=bad)
        if len(bad) > UNSTABLE_SHARE * len(ids):     # many left out could bias the result (2026-10-10)
            return dict(base, status="unstable", se=None,
                        reason=f"{len(bad)} of {len(ids)} production configurations unstable (|L| > {UNSTABLE_L:g} or nan/inf), more than "
                               f"{UNSTABLE_SHARE:.0%}: {bad[:10]}{' ...' if len(bad) > 10 else ''}")
        if len(x) < 2:
            return dict(base, status="fit", se=None, reason="fewer than 2 stable production configurations")
        try:
            tc = t2mod.components(t, x, **self.fit)
            binfo = {}
            boot = t2mod.bootstrap_se_rel(t, x, **self.fit, info=binfo)
        except (ValueError, RuntimeError) as e:
            return dict(base, status="fit", se=None, reason=f"T2 of the production: {e}")
        lin = t2mod.se_rel(tc, tc["C"], tc["S"])
        se = max(lin, boot)
        if not math.isfinite(se):
            return dict(base, status="fit", se=None, reason=f"SE of the production area is not finite (linearized {lin}, bootstrap {boot})")
        return dict(base, status="ok" if se <= self.o["threshold_config"] else "se", se=se, se_rel_linearized=lin, se_rel_bootstrap=boot,
                    A=tc["A"], T2=tc["T2"], n=tc["n"], reason=None if se <= self.o["threshold_config"] else f"SE of the area {se:.3%} > target",
                    **binfo)

    def finish(self, summary, status, stage, reason):
        summary.update(status=status, stage=stage, reason=reason, core_hours=self.core_hours)
        summary["unstable"] = sorted(self.unstable.values(), key=lambda r: r["configuration"])
        summary["corrected"] = sorted(self.corrected.values(), key=lambda r: r["configuration"])
        out = self.work / "result"
        out.mkdir(exist_ok=True)
        _atomic_write(_json(out) / "summary.json", json.dumps(summary, indent=1, default=float))
        lines = self._summary_lines(summary, status, stage, reason)
        _atomic_write(out / "summary.txt", "\n".join(lines) + "\n")
        self.progress(step=len(STEPS) - 1 if status == "converged" else None, status=status,
                      note=None if status == "converged" else f"{stage}: {reason}")
        print("\n".join(lines))
        return {"converged": 0}.get(status, 3)


    def _summary_lines(self, summary, status, stage, reason):
        """summary.txt (plan "결과 글 파일 셋으로", user 2026-10-10): the result first, the warnings, then a short section per step
        that has run (values compared on one line, the value chosen and its threshold on the next, the step's figure named)."""
        o, L = self.o, []
        pct = lambda v: "-" if v is None else f"{100 * v:.2f} %"
        head_fig = lambda text, fig: text.ljust(48) + f"({fig})"
        def table(head, xs, ys):                 # two aligned rows
            cells = [(f"{x:g}", pct(y)) for x, y in zip(xs, ys)]
            w = [max(len(c), len(d)) for c, d in cells]
            hw = max(len(head[0]), len(head[1]))
            return [f"  {head[0]:<{hw}}  " + "  ".join(f"{c:>{n}}" for (c, _), n in zip(cells, w)),
                    f"  {head[1]:<{hw}}  " + "  ".join(f"{d:>{n}}" for (_, d), n in zip(cells, w))]
        p, tr = summary.get("production") or {}, (summary.get("truncation") or {}).get("chosen")
        w = summary.get("window") or {}
        L.append(f"status: {status}" + (f" ({stage}: {reason})" if stage else ""))
        if p.get("se") is not None and p.get("T2") is not None:
            S = p.get("states")
            L.append(f"result: T2 {p['T2']:.4g} ms, n {p['n']:.3g} (A {p['A']:.4g} ms, sampling error {pct(p['se'])})")
            L.append(f"        rbath {tr['rbath']:g} Å, rdip {tr['rdip']:g} Å, {p.get('n_used', len(p['configurations']))} configurations x "
                     + ("ensemble" if S == 0 else f"{S} state" + ("s" if S != 1 else ""))
                     + (f", window {p['window']['window_end']:.3g} ms" if (p.get("window") or {}).get("window_end") else ""))
        if (self.work / "result" / "coherence_mean.txt").is_file() and p.get("configurations"):
            L.append("coherence: coherence_mean.txt")
        L.append("core-hours: unknown (the launcher's rank count or OMP_NUM_THREADS could not be read)" if self.core_hours is None
                 else f"core-hours: {self.core_hours:.3f}")
        warn = list(w.get("warnings", [])) + list(summary.get("warnings", []))
        for name in ("rdip", "rbath"):        # an unconfirmed value: advice, not a result
            r = summary.get(name) or {}
            if r.get("chosen") is None and r.get("advisory") is not None:
                own = next((x for x in r.get("vs_reference", []) if x["value"] == r["advisory"]), None)
                warn.append(f"{name} advisory (not confirmed): {r['advisory']:g} agrees with the largest computed {r['reference']:g}"
                            + (f" ({100 * own['rel']:.2f} %)" if own and own.get("rel") is not None else "")
                            + f", but {r['need']} passes in a row are needed; confirming needs {name} {round(r['reference'] * 1.25, 3):g}"
                            + (f" and the bath files cover {self.max_rbath:.4g} Å" if name == "rbath" else ""))
        if self.corrected:                    # points that took noDiv, the configurations kept (plan "발산한 점 보정")
            warn.append(f"divergent points corrected (|wiDiv - noDiv| >= {POINT_JUMP:g} or nan/inf -> noDiv): " + "; ".join(
                f"{r['configuration']} ({r['points']} points" + (f", wiDiv up to {r['max_wiDiv']:.3g}" if math.isfinite(r["max_wiDiv"]) else ", nan/inf")
                + ")" for r in sorted(self.corrected.values(), key=lambda r: r["configuration"])))
        if self.unstable:                     # left out of every average and judgement (plan "불안정 configuration 빼기")
            warn.append(f"unstable configurations left out (|L| > {UNSTABLE_L:g} or nan/inf): " + "; ".join(
                f"{r['configuration']} ({r['where']}, {r['reason']} at {r['time']:.4g} ms"
                + (f", state {r['state']}" if r["state"] else "")
                + (f", nearest bath spin {r['nearest_spin']:.3g} Å" if r["nearest_spin"] is not None else "") + ")"
                for r in sorted(self.unstable.values(), key=lambda r: r["configuration"])))
        if warn:
            L += ["", "WARNINGS"] + [f"  {x}" for x in warn]
        # 0. start
        rs = self.resolved or {}
        L += ["", "0. start"]
        if self.pilot is not None and self.pool is not None:
            found = rs.get("found_configurations")
            L.append("  configurations: " + (f"{len(found)} found ({Path(self.patterns['bathfile']).name}); " if found else "")
                     + f"test {len(self.pilot)} {self.pilot}, production {getattr(self, 'cap', len(self.pool))}"
                     + (f" of {len(self.pool)} (max_configurations)" if getattr(self, "cap", len(self.pool)) < len(self.pool) else "")
                     + ("" if found else " (given)"))
        if rs.get("r_cov") is not None:
            L.append(f"  bath: {rs.get('geometry', 'bulk')}, files cover {rs['r_cov']:.3g} Å; "
                     + ", ".join(f"{n} spin {v['spin']:g}, density {v['rho']:.3g} /Å^3" for n, v in (rs.get("species") or {}).items()))
        cand = lambda v: "[" + ", ".join(f"{x:g}" for x in v) + "]"
        if rs.get("candi_rbath"):                # the candidates at the start (larger ones computed later are in section 3)
            L.append(f"  rbath candidates {cand(rs['candi_rbath'])}" + (" (given" if "candi_rbath" not in self.auto else " (")
                     + (f"; {autostart.N_SPINS} spins at {rs['b_n']:.3g} Å)" if rs.get("b_n") is not None else ")").replace("(;", "(").replace("()", ""))
        if rs.get("candi_rdip"):
            L.append(f"  rdip candidates {cand(rs['candi_rdip'])}" + (" (given" if "candi_rdip" not in self.auto else " (")
                     + (f"; estimate {rs['r_seed']:.2g} Å)" if rs.get("r_seed") is not None else ")").replace("(;", "(").replace("()", ""))
        L.append("  states: " + ("ensemble (nstate 0)" if self.ensemble else
                                 f"made per run (seed {int(self.base.get('seed', 1))})" if self.auto_states else f"files {self.patterns.get('statefile')}"))
        if self.auto_avaax:
            L.append(f"  avaax: made per run (navaax {self.auto_avaax})")
        L.append("  launcher: " + (" ".join(self.o["launcher"]) + (" (input)" if self.launcher_from == "input" else
                                                                   " (the mpirun that started main.out)")
                                   if self.o["launcher"] else "none (one rank; give launcher or start main.out with mpirun -n N)"))
        # 1. time window
        if w:
            L += ["", head_fig("1. time window", "1_window.png")]
            if w.get("status") == "given":
                L.append(f"  window given: {w.get('recommended_window', 0):.3g} ms")
            if w.get("T2") is not None:
                L.append(f"  T2 fit: T2 {w['T2']:.4g} ms, n {w['n']:.3g}" + (f";  T_d (|L| = {w['threshold']:g}) {w['decay_time']:.3g} ms"
                                                                         if w.get("decay_time") is not None else ""))
            if w.get("fit_points") == "peaks":
                L.append(f"  fit points: the peaks ({w['fit_npoints']} points, t = 0 and the points higher than both neighbours; "
                         "1_window_fit.png)")
            elif w.get("fit_points") == "all":
                why = "t2_fit.points all" if o["fit_points"] == "all" else "fewer than 3 peaks"
                L.append(f"  fit points: all ({w['fit_npoints']} points, {why}; 1_window_fit.png)")
            if w.get("T2") is not None:
                L.append("  T2 here is a convergence measure from this fit, not the exact T2 of the system")
            L.append(f"  time points: {w.get('pilot_nstep') or (w.get('tried') or [{}])[-1].get('nstep')} per run, fixed (whether they are "
                     "enough is not checked)")
            if w.get("status") != "given" and w.get("recommended_window") is not None:
                base = (f"the revival end {w['revival_end']:.3g} ms" if w.get("revival") and w.get("revival_end") is not None else "T_d")
                L.append(f"  window = {w.get('margin', o['production_window_factor']):g} x {base} = {w['recommended_window']:.3g} ms, "
                         f"{w.get('pilot_nstep')} points")
            if w.get("status") not in ("resolved", "given"):
                L.append(f"  status: {w.get('status')}" + (f" ({w['reason']})" if w.get("reason") else ""))
            if w.get("status") != "given":
                L.append("  revival: " + (f"until {w['revival_end']:.3g} ms" if w.get("revival") else "none")
                         + (f";  slow configurations: {len(w['late_rise'])} points above {w['threshold']:g} after T_d (max {w['late_rise_max']:.2g})"
                            if w.get("late_rise") else ""))
            if w.get("tried"):
                ws = [x["window"] for x in w["tried"]]
                L.append(f"  windows tried: {len(ws)}, from {min(ws):.3g} to {max(ws):.3g} ms (each drawn in 1_window.png)")
            for r in summary.get("rbath_raised_for_the_window") or []:
                L.append(f"  plateau near {r['plateau_level']:.3g}: rbath raised to {r['rbath']:g} Å")
            for r in w.get("rounds", []):
                if r.get("coarse_grid"):
                    L.append(f"  grid too coarse: window cut to {r['new_window']:.3g} ms")
                else:
                    L.append(f"  found again at rbath {r['reference_rbath']:g}, rdip {r['reference_rdip']:g}: T_d " + ("-" if r["decay_time"] is None else f"{r['decay_time']:.3g}") + " ms -> "
                             + ("window found again" if r["moved"] else "kept"))
        # 2. rdip, 3. rbath
        for num, name, other, fig in ((2, "rdip", "rbath", "2_rdip.png"), (3, "rbath", "rdip", "3_rbath.png")):
            r = summary.get(name)
            if not r:
                continue
            tol = o[f"threshold_{name}"]
            if r.get("fixed"):
                L += ["", f"{num}. {name}", f"  {name} {r['chosen']:g} given (not converged)"]
                continue
            L += ["", head_fig(f"{num}. {name}  ({other} {r[other]:g} Å, reference {name} {r['reference']:g} Å)" if other in r
                               else f"{num}. {name}  (reference {name} {r['reference']:g} Å)", fig)]
            vs = r.get("vs_reference", [])
            L += table((name, "|ΔA/A|"), [x["value"] for x in vs], [x["rel"] for x in vs])
            L.append(f"  chosen {r['chosen']:g} (threshold {100 * tol:g} %)" if r.get("chosen") is not None else f"  not converged (threshold {100 * tol:g} %)")
            raw = sorted({(lo, hi) for lo, hi, _ in raw_area_pairs(r)})
            if raw:                              # user 2026-10-12: a failed T2 fit is judged on the raw area
                L.append("  judged on the raw area of |L| on [0, T_d] (a T2 fit failed): " + ", ".join(f"{a:g} vs {b:g}" for a, b in raw))
            if name == "rbath":
                for x in summary.get("rdip_at_final_rbath", []):
                    R = (tr or {}).get("rdip", x.get("chosen"))
                    own = next((v for v in x.get("vs_reference", []) if v["value"] == R), None)
                    L.append(f"  rdip rechecked at rbath {x['rbath']:g}: " + (f"{R:g} still passes ({pct(own['rel'])}" + (", raw area" if own.get("basis") == "raw-area" else "") + ")" if own and own["verdict"] == "pass"
                                                                            else f"chosen {x.get('chosen')}"))
        # 4. states
        st = summary.get("states")
        if st:
            L += ["", head_fig(f"4. states  (test {len(st['configurations'])} configurations)", "4_states.png")]
            L += table(("states", "σA/A"), [q["states"] for q in st["points"]], [q["se"] for q in st["points"]])
            L.append((f"  chosen {st['chosen']}" if st.get("chosen") is not None else "  not converged")
                     + f" (threshold {100 * st['threshold']:g} %); production uses test_states {o['pilot_states']}")
        # 5. configurations
        sp = summary.get("sampling_pilot") or {}
        if sp.get("fixed"):
            L += ["", "5. configurations", f"  given: {sp['fixed']['configurations']} x {sp['fixed']['states']} states (not converged)"]
        elif sp or p.get("decomposition"):
            S = p.get("states", o["pilot_states"])
            L += ["", head_fig("5. configurations  (production" + ("" if S == 0 else f", {S} states each") + ")", "5_configurations.png")]
            best = (sp.get("allocation") or {}).get("best")
            if best:
                done = [b["configurations"] for b in p.get("batches", [])]
                L.append(f"  estimate from the test: {best['C']} needed" + (f" -> {' -> '.join(map(str, done))} computed" if done else ""))
            rows = [r for r in p.get("decomposition", [])]
            if rows:
                if len(rows) > 11:              # a line that fits: about 11 counts spread over the range
                    cs = np.array([r["configurations"] for r in rows])
                    pick = sorted({int(np.argmin(np.abs(cs - v))) for v in np.linspace(cs[0], cs[-1], 11)})
                    rows = [rows[i] for i in pick]
                L += table(("configurations", "σA/A"), [r["configurations"] for r in rows], [r["se"] for r in rows])
            if p.get("se") is not None:
                L.append((f"  chosen {p.get('n_used', len(p['configurations']))}" if p.get("status") == "ok" else "  not enough")
                         + f" (threshold {100 * o['threshold_config']:g} %)")
        rc = summary.get("recheck") or {}
        if rc.get("route") == "raw-area":        # user 2026-10-12
            L += ["", f"combination recheck {rc['verdict']} on the raw area of |L| on [0, T_d] (a T2 fit failed): "
                      f"differs by {100 * abs(rc['raw_area_rel']):.2f} %"]
        # 6. converged coherence
        if p.get("T2") is not None:
            L += ["", head_fig("6. converged coherence", "6_converged_coherence.png"),
                  f"  T2 {p['T2']:.4g} ms, n {p['n']:.3g}, A {p['A']:.4g} ms, σA/A {pct(p['se'])}"]
            if p.get("windows"):
                L.append("  production window lengthened: " + " -> ".join(f"{x['window']:.3g} ms" for x in p["windows"]))
            wp = p.get("window") or {}
            if wp.get("decay_time") is not None:
                L.append(f"  window check: decay {wp['decay_time']:.3g} ms " + ("inside" if wp.get("status") == "resolved" else "beyond")
                         + f" {wp.get('window_end', 0):.3g} ms")
        return L


def status_main(argv):
    """main.out --auto_status -f input.json: where an automatic convergence is, without computing (plan R4). Running (a process
    holds the lock), interrupted (the result says running but nobody holds the lock), or how it ended."""
    import argparse
    ap = argparse.ArgumentParser(prog="autoconv status")
    ap.add_argument("input")
    a = ap.parse_args(argv)
    try:
        raw = read_input(Path(a.input))
        work = (Path(a.input).resolve().parent / raw["auto_convergence"]["workdir"]).resolve()
    except (OSError, ValueError, KeyError, TypeError) as e:
        print(f"autoconv status: cannot read the workdir from {a.input}: {e}", file=sys.stderr)
        return 2
    res = work / "result"
    held = False
    if (work / "lock").is_file():
        with open(work / "lock") as f:
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(f, fcntl.LOCK_UN)
            except BlockingIOError:
                held = True
    js = res / "json" if (res / "json" / "summary.json").is_file() else res      # json/ since 2026-10-10
    if not held and not (js / "summary.json").is_file():
        print(f"not started: no result in {work}")
        return 0
    try:
        st = json.loads((js / "summary.json").read_text()).get("status")
    except (OSError, ValueError, AttributeError):
        st = None
    state = "running" if held else ("interrupted (the job was cut; run the same command to continue)" if st == "running"
                                    else (st or "unknown (summary.json missing or unreadable)"))
    print(f"state: {state}")
    def text(name):                       # a fresh run may move result/ away while this reads (82c666-b #11)
        try:
            return (res / name).read_text()
        except OSError:
            return None
    print(text("progress.txt") or "(no progress.txt)\n", end="")
    if not held and st != "running" and text("summary.txt"):
        print("--- summary.txt ---")
        print(text("summary.txt"), end="")
    return 0


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="autoconv run")
    ap.add_argument("input")
    ap.add_argument("--binary", required=True)
    a = ap.parse_args(argv)
    for v in ("PMI_RANK", "PMI_SIZE", "OMPI_COMM_WORLD_RANK", "MPI_LOCALRANKID"):
        if v in os.environ:
            print(f"autoconv run: refused inside an MPI launch ({v} is set); call main.out --auto_conv without mpirun", file=sys.stderr)
            return 2
    try:
        r = Runner(a.input, a.binary)
    except Refused as e:
        try:                                     # an earlier result in the workdir is not left as the present one
            raw = read_input(Path(a.input))
            out = Path(a.input).parent / raw["auto_convergence"]["workdir"] / "result"
            if out.is_dir():                     # written only while holding the workdir's lock (82c666-b #8)
                with open(out.parent / "lock", "a") as f:
                    try:
                        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError:
                        pass                     # another runner holds it: its summary is not ours
                    else:
                        _atomic_write(_json(out) / "summary.json", json.dumps(dict(status="refused", reason=str(e)), indent=1))
                        _atomic_write(out / "summary.txt", f"status: refused ({e})\n")
        except (OSError, ValueError, KeyError, TypeError):
            pass
        print(f"autoconv run: {e}", file=sys.stderr); return 2
    lock = open(r.work / "lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"autoconv run: another runner holds {r.work / 'lock'}", file=sys.stderr); return 2
    if r.resume == "fresh":               # start over: everything of the workdir but the lock goes to <workdir>.old<n> (R1)
        old = [p for p in r.work.iterdir() if p.name != "lock"]
        if old:
            n = 1
            while r.work.with_name(f"{r.work.name}.old{n}").exists():
                n += 1
            keep = r.work.with_name(f"{r.work.name}.old{n}")
            keep.mkdir()
            for p in old:
                shutil.move(str(p), str(keep / p.name))
            print(f"[auto] resume fresh: the earlier work moved to {keep}", file=sys.stderr)
    else:                                  # the lock was free but the result says running: the run before was cut (R2)
        try:
            js = r.work / "result" / "json"
            js = js if (js / "summary.json").is_file() else js.parent                  # json/ since 2026-10-10
            prev = json.loads((js / "summary.json").read_text())
            if prev.get("status") == "running":
                p = json.loads((js / "progress.json").read_text()) if (js / "progress.json").is_file() else {}
                r.interrupted = dict(step=p.get("step"), updated=p.get("updated"))
                print(f"[auto] the run before was interrupted (step {p.get('step')}, last update {p.get('updated')}): continuing it",
                      file=sys.stderr)
        except (OSError, ValueError):
            pass
    def stop(signum, frame):
        raise Stop(f"signal {signum}")
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, stop)
    try:
        return r.execute()
    except Refused as e:
        r.mark("refused", str(e)); print(f"autoconv run: {e}", file=sys.stderr); return 2
    except Stop as e:
        r.mark("stopped", str(e)); print(f"autoconv run: stopped ({e}); call again to continue", file=sys.stderr); return 4
    except (RuntimeError, OSError, ValueError) as e:
        r.mark("error", str(e)); print(f"autoconv run: {e}", file=sys.stderr); return 1
