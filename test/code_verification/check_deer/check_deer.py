#!/usr/bin/env python3
"""DEER bath pulses (pulse.bath): checks on the runs check_deer/run.sh made.

Ported from adding_MECCE example/z.Tests/DEER/check_deer.py. Every check is an identity the code must satisfy or a
comparison with the closed form, which at cluster order 1 is exact: the bath spins do not interact, so the echo is a
product over them,

    V(t) = prod_j cos(a_j t min(f, 1-f)),   a_j = 2 pi 52.04 MHz nm^3 (1 - 3 cos^2 theta_j) / r_j^3

over the spins the pump drives, f the fraction of the total time t at which the bath pi lands (f = 0.5: a_j t / 2).

Usage: check_deer.py <work dir> <inputs dir>    (files <work>/<case>_noDiv, <case>_state1_noDiv for nstate 1)
Exit status 0 only when every check passes.
"""
import os
import sys

import numpy as np

work, inputs = sys.argv[1], sys.argv[2]
BATH = os.path.join(inputs, "bath_DEERtest_12spin.txt")
AVAAX = os.path.join(inputs, "avaax_DEERtest")
D_MHZ = 52.04                       # (mu0/4pi) gamma_e^2 hbar / h, two electrons [MHz nm^3]
PUMPED_AXIS = 1                     # the configuration whose detuning (79.975539 MHz) the resonant tone hits
TOL_EXACT = 0.0                     # identities of runs on the same code path: bit for bit
TOL_MERGED = 1e-12                  # a merged sequence (Manual, fractions as edge differences) against the generated
                                    # one: the fractions may differ by an ulp (controller ruling A3 ulp-drift)
TOL_CLOSED = 2e-4                   # closed form vs code (branch check_deer.py)
TOL_REFOCUS = 1e-5                  # refocusing is exact in theory; the propagator's precision with a 2.88 GHz ZFS
TOL_PATH = 1e-9                     # recursive vs exact bound stated in CHANGELOG 2.0.0 (evolution_update)

ENSEMBLE = ("nopulse", "onres", "offres", "multi_hit", "multi_miss", "nearpulse")
PURE = ("path_matrix", "path_matrix_recursive", "path_vector", "path_recursive", "path_nopump")


def out_file(case, kind="noDiv"):
    # nstate 1 runs write <o>_state1_<kind>, the ensemble ones <o>_<kind>
    return os.path.join(work, "{}{}_{}".format(case, "_state1" if case in PURE else "", kind))


def read_curve(case, kind="noDiv"):
    t, v = [], []
    for line in open(out_file(case, kind)):
        w = line.split()
        if len(w) == 2:
            t.append(float(w[0]) * 1e3)                       # ms -> us
            v.append(complex(w[1].replace("+-", "-")))       # the whole complex value (pure states have Im != 0)
    return np.array(t), np.array(v)


def closed_form(t, frac):
    rows = [l.split() for l in open(BATH).read().splitlines()[1:]]
    rows = [r for r in rows if len(r) == 4]
    xyz = np.array([[float(x) for x in r[:3]] for r in rows])
    ax = np.array([int(round(float(x))) for x in open(AVAAX).read().split()[1:]])
    r = np.linalg.norm(xyz, axis=1)
    sel = ax == PUMPED_AXIS
    f = D_MHZ * (1 - 3 * (xyz[sel, 2] / r[sel]) ** 2) / (r[sel] / 10.0) ** 3   # MHz (r in A -> nm)
    eff = min(frac, 1.0 - frac)
    return np.prod(np.cos(2 * np.pi * np.outer(t, f) * eff), axis=1)


def report(name, ok, detail):
    print("  [{}] {:<58} {}".format("PASS" if ok else "FAIL", name, detail))
    return ok


def main():
    missing = [c for c in ENSEMBLE + PURE if not os.path.exists(out_file(c))]
    if missing:
        print("  FAIL  check_deer.py: no output for {}".format(missing))
        return 1

    t, none = read_curve("nopulse")
    _, on = read_curve("onres")
    _, off = read_curve("offres")
    _, mhit = read_curve("multi_hit")
    _, mmiss = read_curve("multi_miss")
    _, near = read_curve("nearpulse")
    ok = []

    print("\n DEER bath pulses, 12 P1 spins, gCCE order 1: ensemble (nstate 0) and one pure state (nstate 1)")

    d = np.abs(none - 1.0).max()
    ok.append(report("order 1 without a pump refocuses exactly", d <= TOL_REFOCUS, "max |V - 1| = {:.2e}".format(d)))
    d = np.abs(off - none).max()
    ok.append(report("off-resonant pump == no pump (merged sequence, 1e-12)", d <= TOL_MERGED,
                     "max |diff| = {:.2e}".format(d)))
    d = np.abs(mhit - on).max()
    ok.append(report("multi-tone [resonant, far] == single resonant", d <= TOL_EXACT, "max |diff| = {:.2e}".format(d)))
    d = np.abs(mmiss - none).max()
    ok.append(report("multi-tone [far, far] == no pump (merged sequence, 1e-12)", d <= TOL_MERGED,
                     "max |diff| = {:.2e}".format(d)))

    # the closed form is real; compare the real part (the ensemble's imaginary part is ~0)
    e = np.abs(on.real - closed_form(t, 0.5)).max()
    ok.append(report("pump at the echo centre matches the closed form", e <= TOL_CLOSED, "max |diff| = {:.2e}".format(e)))
    # a bath pi 0.005 of the sequence away from the qubit pi: the branch once applied it at both edges (two pi cancel)
    e = np.abs(near.real - closed_form(t, 0.505)).max()
    ok.append(report("pump next to the qubit pulse is applied once, not twice", e <= TOL_CLOSED,
                     "max |diff| vs closed form = {:.2e}".format(e)))
    d = np.abs(near - none).max()
    ok.append(report("  ... and is therefore not equal to no pump", d > 1e-3, "max |diff| vs nopulse = {:.2e}".format(d)))

    bad = [c for c in ENSEMBLE if open(out_file(c, "wiDiv")).read() != open(out_file(c, "noDiv")).read()]
    ok.append(report("order 1 has no inclusion-exclusion division (wiDiv == noDiv)", not bad,
                     "differing: {}".format(bad if bad else "none")))

    # the evolution paths on a pure state: the pulse acts on every one of them, and they agree
    _, pm = read_curve("path_matrix")
    _, pv = read_curve("path_vector")
    _, pmr = read_curve("path_matrix_recursive")
    _, pr = read_curve("path_recursive")
    _, pn = read_curve("path_nopump")
    ok.append(report("bath pulse: evolution vector == matrix (exact)", np.abs(pv - pm).max() <= TOL_PATH,
                     "max |diff| = {:.2e}".format(np.abs(pv - pm).max())))
    ok.append(report("bath pulse: recursive == exact (vector)", np.abs(pr - pv).max() <= TOL_PATH,
                     "max |diff| = {:.2e}".format(np.abs(pr - pv).max())))
    ok.append(report("bath pulse: recursive == exact (matrix)", np.abs(pmr - pm).max() <= TOL_PATH,
                     "max |diff| = {:.2e}".format(np.abs(pmr - pm).max())))
    ok.append(report("bath pulse acts on the vector path", np.abs(pv - pn).max() > 1e-3,
                     "max |diff| vs no pump = {:.2e}".format(np.abs(pv - pn).max())))

    n = sum(ok)
    print(" {} / {} checks passed".format(n, len(ok)))
    return 0 if n == len(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
