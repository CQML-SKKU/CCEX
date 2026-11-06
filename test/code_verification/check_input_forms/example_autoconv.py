#!/usr/bin/env python3
"""example_autoconv.py EXAMPLE REPO: the auto_convergence lines of the template `main.out -e` writes are accepted by the
runner of the automatic convergence (A94, review 20261016 R6: the template's lines, uncommented, were refused - keys that
cannot go together on one line).

The commented auto_convergence block of the template is uncommented (the leading // of each of its lines; the trailing
comments stay comments). A key written more than once in one object is a set of alternatives (the line that finds a value
and the line that fixes it); every combination of the alternatives is put into the template's active input (its
configuration.bathfile given a {k} when it has none, as the runner needs) and checked by the runner's own reading of the
input (autoconv/runner.py Runner: the block's keys and the checks between them, before anything runs). A combination with
every section fixed is refused by the runner ("nothing to converge") and is not counted; at least one combination has to
be accepted. Exit 0 when all hold, 1 with the refusals otherwise.
"""
import itertools
import json
import re
import sys
import tempfile
from pathlib import Path

example_path, repo = sys.argv[1:3]
sys.path.insert(0, f"{repo}/example/tools/AutoConvergence")
from autoconv.config import strip_comments            # noqa: E402  (the CCEX comment rules)
from autoconv.runner import Refused, Runner           # noqa: E402

lines = Path(example_path).read_text().splitlines()
start = next((i for i, l in enumerate(lines) if re.match(r'\s*//\s*"auto_convergence"\s*:', l)), None)
problems = []
if start is None:
    print("  the template has no commented auto_convergence block")
    sys.exit(1)
block, depth = [], 0
for l in lines[start:]:
    # a commented line of the input (// then a key, a value or a bracket); a comment line of its own (the continued
    # explanation lines above each key) stays a comment
    text = re.sub(r'^(\s*)//', r'\1', l, count=1) if re.match(r'\s*//\s*["{}\[\]0-9-]', l) else l
    block.append(text)
    code = strip_comments(text)
    depth += code.count("{") - code.count("}")
    if depth == 0 and "{" in "".join(strip_comments(b) for b in block):
        break
src = strip_comments("\n".join(block)).strip().rstrip(",")
src = src[src.index(":") + 1:]


class Multi(list):
    """An object read with every member kept, a key written twice included."""


def expand(v):
    """Every combination of the alternatives (a key written more than once in one object) inside v."""
    if isinstance(v, Multi):
        keys = list(dict.fromkeys(k for k, _ in v))
        choices = [[x for k2, val in v if k2 == k for x in expand(val)] for k in keys]
        return [dict(zip(keys, combo)) for combo in itertools.product(*choices)]
    if isinstance(v, list):
        return [list(c) for c in itertools.product(*[expand(x) for x in v])] if v else [[]]
    return [v]


block_value = json.loads(src, object_pairs_hook=Multi)
combos = expand(block_value)
base = json.loads(strip_comments(Path(example_path).read_text()))
conf = base.setdefault("configuration", {})
if not any("{k}" in str(b) for b in conf.get("bathfile", [])):
    conf["bathfile"] = ["bath_{k}"]
if "states" in src and "s_range" not in conf:      # the state keys of sampling are read with sampled states (s_range)
    conf["s_range"] = [1, 1]
accepted, refused = 0, []
with tempfile.TemporaryDirectory() as tmp:
    for i, combo in enumerate(combos):
        d = Path(tmp) / f"c{i}"
        d.mkdir()
        (d / "input.json").write_text(json.dumps(dict(base, auto_convergence=combo)))
        try:
            Runner(d / "input.json", d / "no_main.out")
        except Refused as e:
            if "nothing to converge" not in str(e):
                refused.append(f"{json.dumps(combo)}: {e}")
            continue
        except OSError:
            pass                       # the block was accepted; the runner then looks for the binary, which is not there
        accepted += 1
print(f"auto_convergence lines of the template: {len(combos)} combinations of the alternatives, {accepted} accepted by the runner")
for r in refused[:5]:
    print("  refused: " + r)
if len(refused) > 5:
    print(f"  ... {len(refused) - 5} more refused")
if not accepted:
    print("  no combination of the lines was accepted")
sys.exit(1 if refused or not accepted else 0)
