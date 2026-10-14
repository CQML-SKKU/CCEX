#!/usr/bin/env python3
"""A Q file from a VASP OUTCAR without a defect spin (review 20261013-180049 T1).

The spin S of the A file's provenance line is read from the OUTCAR's total magnetization. The check ran for every
request, so a Q file alone stopped on an EFG OUTCAR with magnetization 0 or none (an EFG run need not be spin
polarized). On examples/02.hBN_VB_VASP, with the magnetization lines set to 0 and removed: --Q alone runs and writes the
Q file the OUTCAR as it is gives; --A still stops on them, and on a --spin the OUTCAR contradicts.

Usage:  check_qfile_spin.py <work dir>      (exit 0 = pass)
"""
import os, re, shutil, subprocess, sys

WORK = os.path.abspath(sys.argv[1])
TOOL = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(TOOL, "examples/02.hBN_VB_VASP")
MAG = re.compile(r"^( number of electron .*magnetization\s+)\S+", re.M)
fails = 0


def run(name, outcar_text, *options):
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    with open(os.path.join(d, "OUTCAR"), "w") as f:
        f.write(outcar_text)
    args = [a if a != "OUTCAR" else os.path.join(d, "OUTCAR") for a in options]
    p = subprocess.run([sys.executable, os.path.join(TOOL, "main.py"), "VASP", os.path.join(EXAMPLE, "VASP_CONTCAR_flat"),
                        "VB", *args], cwd=d, capture_output=True, text=True)
    with open(os.path.join(d, "log"), "w") as f:
        f.write(p.stdout + p.stderr)
    return p.returncode, p.stdout + p.stderr, d


def check(ok, label, d):
    global fails
    print(f"  {'PASS' if ok else 'FAIL'}  {label}" + ("" if ok else f" (see {d}/log)"))
    fails += not ok


outcar = open(os.path.join(EXAMPLE, "VASP_OUTCAR_flat")).read()
assert MAG.search(outcar), "the example OUTCAR has no magnetization line"
zero = MAG.sub(lambda m: m.group(1) + "0.0000000", outcar)
none = "\n".join(l for l in outcar.split("\n") if not MAG.match(l + "\n"))

rc, _, ref = run("as_is", outcar, "--Q", "OUTCAR")
check(rc == 0 and os.path.isfile(os.path.join(ref, "Qfile")), "Q file from the OUTCAR as it is", ref)
for name, text in (("zero", zero), ("none", none)):
    rc, _, d = run(f"q_{name}", text, "--Q", "OUTCAR")
    same = rc == 0 and os.path.isfile(os.path.join(d, "Qfile")) and \
        open(os.path.join(d, "Qfile")).read() == open(os.path.join(ref, "Qfile")).read()
    check(same, f"Q file alone, magnetization {name}: written, and the same as from the OUTCAR as it is", d)
    rc, out, d = run(f"a_{name}", text, "--A", "OUTCAR")
    check(rc != 0 and ("total magnetization" in out or "no 'magnetization' line" in out) and
          not os.path.isfile(os.path.join(d, "Afile")), f"A file, magnetization {name}: refused, no A file", d)
rc, out, d = run("a_spin", outcar, "--A", "OUTCAR", "--spin", "0.5")
check(rc != 0 and "contradicts the OUTCAR" in out and not os.path.isfile(os.path.join(d, "Afile")),
      "A file with a --spin the OUTCAR contradicts: refused", d)
sys.exit(1 if fails else 0)
