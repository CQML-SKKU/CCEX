"""Driver configuration of the auto_convergence setup tool ((internal notes), section 0 and the
`auto_convergence` input specification).

The driver reads one JSON file (comments starting with '#' allowed, as in CCEX inputs):

    {
      "base_input": "./ccein.json",
      "auto_convergence": {"tolerance": 0.01, "steps": {"time_window": true, "sampling": true}},
      "bath_source":  {"files": "./baths/config_P1_*"},
      "state_source": {"files": "./states/cfg{k}_state"},
      "server": {"host": "idm", "qos": 1, "workdir": "~/tmp/autoconv/name"}
    }

`load_config` validates the keys and returns a `Settings` with every default filled in. Unknown keys and keys of steps
that are not implemented yet are refused rather than ignored.
"""
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

# Steps the driver can run in this version. rdip and rbath need the streaming executor or the ladder driver (later).
IMPLEMENTED_STEPS = {"time_window", "sampling"}
ALL_STEPS = ("time_window", "rdip", "rbath", "sampling")

# Error budget split (section 0.4 / the input specification): truncation and numerics take 1/3 of the total, sampling 2/3;
# inside the first third rdip : rbath = 50 : 50 (the time step is the user's nstep, not a convergence item; 2026-10-07 user
# decision). Initial values to be checked, not established.
BUDGET_TRUNCATION, BUDGET_SAMPLING = 1 / 3, 2 / 3
SPLIT = {"rdip": 0.5, "rbath": 0.5}

KEYS = {"steps", "tolerance", "rdip_tolerance", "rbath_tolerance", "sampling_tolerance", "error_norm",
        "rdip_start", "rdip_initial_points", "search_nstep", "decay_threshold", "window_margin", "rbath_start",
        "t2_se_target", "truncation_max_tol", "truncation_area_tol", "window_median_points", "t2_fit_points", "window_method"}
# Keys that only make sense once their step exists; giving them while the step is not implemented is refused.
STEP_KEYS = {"rdip": {"rdip_tolerance", "rdip_start", "rdip_initial_points", "search_nstep"},
             "rbath": {"rbath_tolerance", "rbath_start"}}
TOP_KEYS = {"base_input", "auto_convergence", "bath_source", "state_source", "server"}


class ConfigError(ValueError):
    """A driver configuration the tool refuses, with the reason."""


def strip_comments(text):
    """Remove '#' and '!' comments (both accepted by the CCEX input reader) outside JSON strings."""
    out, in_str, esc, skip = [], False, False, False
    for ch in text:
        if skip:
            if ch == "\n":
                skip = False
                out.append(ch)
            continue
        if in_str:
            out.append(ch)
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
            out.append(ch)
        elif ch in "#!":
            skip = True
        else:
            out.append(ch)
    return "".join(out)


def read_json(path):
    """JSON with '#' or '!' comments (CCEX input style)."""
    return json.loads(strip_comments(Path(path).read_text()))


# CCEX reads input keys without regard to case (docs/user/InputFile.md) and the runner looks them up by their written
# names, so a CCEX input's keys are put in those names when it is read (review 20261013-180049 R4): lower case, but for
# the few names written with capitals. Values (paths, names, strings) stay as written.
MIXED_CASE_KEYS = {"ccce": "cCCE", "gcce": "gCCE", "mecce": "meCCE", "delta_a": "delta_A", "d": "D", "e": "E"}
TOP_MIXED_CASE_KEYS = {"qubit": "Qubit", "defect": "Defect"}      # Qubit.qubit stays lower case


class KeyCaseError(ValueError):
    """Two keys of one block that differ only in case."""


def _single_keys(pairs):
    seen = {}
    for k, _ in pairs:
        if k.lower() in seen:
            raise KeyCaseError(f'input keys "{seen[k.lower()]}" and "{k}" are one key to CCEX (case is not read); keep one')
        seen[k.lower()] = k
    return dict(pairs)


def canonical_keys(value, top=True):
    if isinstance(value, list):
        return [canonical_keys(v, False) for v in value]
    if not isinstance(value, dict):
        return value
    return {(TOP_MIXED_CASE_KEYS.get(k.lower()) if top else None) or MIXED_CASE_KEYS.get(k.lower(), k.lower()):
            canonical_keys(v, False) for k, v in value.items()}


def read_input(path):
    """A CCEX input (read_json) with its keys in the names the tools look up; keys that differ only in case are refused."""
    return canonical_keys(json.loads(strip_comments(Path(path).read_text()), object_pairs_hook=_single_keys))


def _positive(name, value):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
        raise ConfigError(f"{name} must be a finite positive number, got {value!r}")
    return float(value)


@dataclass
class Settings:
    steps: dict
    tolerance: float
    rdip_tolerance: float
    rbath_tolerance: float
    sampling_tolerance: float
    error_norm: str
    decay_threshold: float
    window_margin: float
    rdip_start: float | None = None
    rdip_initial_points: int = 5
    search_nstep: int = 50
    rbath_start: float | None = None
    # 2026-10-07 user decisions (auto-truncation-plan.md): sampling is judged on the relative SE of T2, truncation on the
    # fit-line difference (max and area on [0, T_d]); the curve tolerances above stay as diagnostics.
    t2_se_target: float = 0.06
    truncation_max_tol: float = 0.01
    truncation_area_tol: float = 0.005
    window_median_points: int = 9          # rule 2 (2026-10-07): running median of the window envelope, 1 = off
    t2_fit_points: str = "envelope"        # T2 fits on the upper envelope (2026-10-07, modulated curves) or "all" points
    window_method: str = "fit"             # T_d where the T2 fit line reaches the threshold (2026-10-07) or "envelope"
    base_input: dict = field(default_factory=dict)
    bath_source: dict = field(default_factory=dict)
    state_source: dict = field(default_factory=dict)
    server: dict = field(default_factory=dict)
    root: Path = Path(".")


def parse_auto_convergence(block, base_input=None):
    """Validate an auto_convergence block and fill the defaults. base_input (the CCEX input as a dict) is used for the
    step/base-value consistency checks; pass None to skip them."""
    if not isinstance(block, dict):
        raise ConfigError("auto_convergence must be an object")
    unknown = set(block) - KEYS
    if unknown:
        raise ConfigError(f"unknown auto_convergence keys: {sorted(unknown)}")

    raw_steps = block.get("steps", {})
    if not isinstance(raw_steps, dict) or set(raw_steps) - set(ALL_STEPS):
        raise ConfigError(f"steps must be an object with keys from {ALL_STEPS}, got {raw_steps!r}")
    for k, v in raw_steps.items():
        if not isinstance(v, bool):
            raise ConfigError(f"steps.{k} must be true or false, got {v!r}")
    steps = {k: raw_steps.get(k, k in IMPLEMENTED_STEPS) for k in ALL_STEPS}
    for k, on in steps.items():
        if on and k not in IMPLEMENTED_STEPS:
            raise ConfigError(f"step '{k}' is not implemented yet; set steps.{k} to false and give its value in base_input")
    for step, keys in STEP_KEYS.items():
        if step not in IMPLEMENTED_STEPS and keys & set(block):
            raise ConfigError(f"{sorted(keys & set(block))} belong to the step '{step}', which is not implemented yet")

    tol = _positive("tolerance", block.get("tolerance", 0.01))
    trunc = tol * BUDGET_TRUNCATION
    tolerances = {f"{k}_tolerance": _positive(f"{k}_tolerance", block.get(f"{k}_tolerance", trunc * SPLIT[k])) for k in SPLIT}
    tolerances["sampling_tolerance"] = _positive("sampling_tolerance", block.get("sampling_tolerance", tol * BUDGET_SAMPLING))
    parts = sum(tolerances.values())
    if parts > tol * (1 + 1e-12):
        raise ConfigError(f"rdip + rbath + sampling tolerances ({parts:.6g}) exceed the total tolerance ({tol:.6g})")

    norm = block.get("error_norm", "area")
    if norm not in ("area", "max"):
        raise ConfigError(f"error_norm must be 'area' or 'max', got {norm!r}")
    eps = _positive("decay_threshold", block.get("decay_threshold", 0.01))
    if eps >= 1:
        raise ConfigError("decay_threshold must be below 1")
    margin = _positive("window_margin", block.get("window_margin", 1.5))      # 1.5 T_d (user decision 2026-10-07; was 1.2)
    if margin < 1:
        raise ConfigError("window_margin must be at least 1")

    extra = {k: _positive(k, block[k]) for k in ("t2_se_target", "truncation_max_tol", "truncation_area_tol") if k in block}
    if "window_median_points" in block:
        m = block["window_median_points"]
        if not isinstance(m, int) or isinstance(m, bool) or m < 1 or m % 2 == 0:
            raise ConfigError("window_median_points must be a positive odd integer")
        extra["window_median_points"] = m
    if "t2_fit_points" in block:
        if block["t2_fit_points"] not in ("envelope", "all"):
            raise ConfigError("t2_fit_points must be 'envelope' or 'all'")
        extra["t2_fit_points"] = block["t2_fit_points"]
    if "window_method" in block:
        if block["window_method"] not in ("fit", "envelope"):
            raise ConfigError("window_method must be 'fit' or 'envelope'")
        extra["window_method"] = block["window_method"]
    s = Settings(steps=steps, tolerance=tol, error_norm=norm, decay_threshold=eps, window_margin=margin, **tolerances, **extra)
    if base_input is not None:
        check_base_input(s, base_input)
    return s


def check_base_input(s, base):
    """A step that is on must not also have its value in base_input (ambiguous); a step that is off must have it."""
    fixed = {"time_window": ("deltat", "nstep"), "rdip": ("rdip",), "rbath": ("rbath",)}
    for step, keys in fixed.items():
        present = [k for k in keys if k in base]
        if s.steps[step] and step != "time_window" and present:
            raise ConfigError(f"step '{step}' is on but base_input also gives {present}; remove them or turn the step off")
        if not s.steps[step] and len(present) != len(keys):
            raise ConfigError(f"step '{step}' is off, so base_input must give {list(keys)}")
    # The time-window step still needs a trial window (deltat, nstep) to start from; it is replaced by the chosen one.
    if s.steps["time_window"] and not all(k in base for k in ("deltat", "nstep")):
        raise ConfigError("the time_window step needs a trial deltat and nstep in base_input")


def source_capacity(s):
    """Number of distinct existing bath files a fixed file source can give (None for a generator or no source). Absolute
    and ~ patterns are allowed; directories and missing list entries are not counted."""
    import glob
    import os
    files = s.bath_source.get("files")
    if files is None:
        return None
    def resolve(p):
        p = os.path.expanduser(p)
        return p if os.path.isabs(p) else str(Path(s.root) / p)
    paths = [resolve(f) for f in files] if isinstance(files, list) else glob.glob(resolve(files))
    return len({os.path.realpath(p) for p in paths if os.path.isfile(p)})


def read_manifest(path):
    """Pilot/production manifest: {"configurations": [{"id": ..., "run": ..., "states": N or [ids]}, ...],
    "available_configurations": <unused configurations for production>, optional "available_states": <reproducible states per
    configuration>, optional "grid": {"deltat": .., "nstep": ..}}. Run paths are relative to the manifest.
    Returns (runs, ids, expected state counts, available, grid, available_states)."""
    m = read_json(path)
    rows = m.get("configurations")
    if not isinstance(rows, list) or not rows:
        raise ConfigError("manifest: 'configurations' must be a non-empty list")
    runs, ids, counts = [], [], []
    for r in rows:
        if not {"id", "run", "states"} <= set(r):
            raise ConfigError(f"manifest: every configuration needs id, run and states, got {r!r}")
        runs.append(str((Path(path).parent / r["run"]).resolve()))
        ids.append(str(r["id"]))
        counts.append(len(r["states"]) if isinstance(r["states"], list) else int(r["states"]))
    if len(set(ids)) != len(ids):
        raise ConfigError("manifest: configuration IDs repeat")
    return runs, ids, counts, m.get("available_configurations"), m.get("grid"), m.get("available_states")


def expected_run_shape(s):
    """(states, time points) every run must have, from base_input (nstate > 0, nstep); None where not fixed."""
    nstate, nstep = s.base_input.get("nstate"), s.base_input.get("nstep")
    return (nstate if isinstance(nstate, int) and nstate > 0 else None,
            nstep if isinstance(nstep, int) and nstep > 0 else None)


def load_config(path):
    path = Path(path)
    cfg = read_json(path)
    unknown = set(cfg) - TOP_KEYS
    if unknown:
        raise ConfigError(f"unknown driver keys: {sorted(unknown)}")
    if "base_input" not in cfg:
        raise ConfigError("base_input (the CCEX input with the physical settings) is required")
    base_path = (path.parent / cfg["base_input"]).resolve()
    base = read_json(base_path)
    s = parse_auto_convergence(cfg.get("auto_convergence", {}), base)
    bath, state = cfg.get("bath_source", {}), cfg.get("state_source", {})
    if s.steps["sampling"]:
        if not (("files" in bath) ^ ("generator" in bath)):
            raise ConfigError("bath_source needs exactly one of 'files' or 'generator'")
        if not state:
            raise ConfigError("state_source is required: states must be fixed and reproducible (state files or a seed rule)")
    s.base_input, s.bath_source, s.state_source, s.server, s.root = base, bath, state, cfg.get("server", {}), path.parent
    return s
