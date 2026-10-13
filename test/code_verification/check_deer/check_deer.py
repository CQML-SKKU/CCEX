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
TOL_ME = 1e-6                       # megCCE at rate 0 vs gCCE: the bound of check_master_equation (Task B4)
TOL_PCCE = 1e-9                     # cluster_mode pcce, partition_size 1 vs standard (controller ruling R9)
TOL_PERM = 1e-8                     # the spins of a cluster in another order: rounding only, the 2.88 GHz ZFS phase
                                    # (1.8e5 rad at 10 us) through 79 cluster factors at order 2 (295 at order 3);
                                    # measured 6.4e-10 / 1.4e-10 (gCCE / megCCE order 2), 5.1e-10 (order 3) (WSL, MKL).
                                    # The factor 15 above that is arbitrary; a pulse on a wrong spin moves them by ~1
TOL_EXPM = 1e-9                     # gCCE propagator expm vs eigen with the pump: two routes to one propagator, rounding
                                    # only (final review ruling); measured 1.0e-10, the last written decimal (WSL, MKL)

ENSEMBLE = ("nopulse", "onres", "offres", "multi_hit", "multi_miss", "nearpulse")
PURE = ("path_matrix", "path_matrix_recursive", "path_vector", "path_recursive", "path_nopump")
C2_ENSEMBLE = ("meg_onres", "meg_nearpulse", "multi_hit_rev", "pg_onres", "pmeg_onres",   # C2: megCCE, pcce, tones
               "onres_expm")                                                             # final review: expm
C2_PURE = ("o2", "o2_rev", "meg_o2", "meg_o2_rev", "o3", "o3_rot",      # order 2 reversed, order 3 rotated
           "pg_o2")                                                    # final review: pcce at order 2


def out_file(case, kind="noDiv"):
    # nstate 1 runs write <o>_state1_<kind>, the ensemble ones <o>_<kind>
    return os.path.join(work, "{}{}_{}".format(case, "_state1" if case in PURE + C2_PURE else "", kind))


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

    # Task C2: megCCE with the bath pulse (U = U_q x U_b in its pulse superoperator), the tones in the other order,
    # cluster_mode pcce with single-spin partitions (inside pgCCE, pmegCCE), order 2 with the bath file reversed
    print("\n Task C2: megCCE and cluster_mode pcce with bath pulses; the bath file reordered at order 2 and 3")
    curves = {("onres", "noDiv"): on, ("nearpulse", "noDiv"): near}
    for c in C2_ENSEMBLE + C2_PURE + ("path_vector",):
        for kind in ("noDiv", "wiDiv"):
            if os.path.exists(out_file(c, kind)):
                curves[c, kind] = read_curve(c, kind)[1]

    def check(name, a, b, test, kinds=("noDiv",), what="max |diff| = {:.2e}"):   # a missing output fails the check
        keys = [(c, k) for k in kinds for c in (a, b)]
        if any(key not in curves for key in keys):
            ok.append(report(name, False, "no output for {}".format(sorted({c for c, k in keys if (c, k) not in curves}))))
        else:
            d = max(np.abs(curves[a, k] - curves[b, k]).max() for k in kinds)
            ok.append(report(name, test(d), what.format(d)))

    check("megCCE rate 0 with bath pulse == gCCE", "meg_onres", "onres", lambda d: d <= TOL_ME)
    check("megCCE: pump next to the qubit pulse == gCCE (applied once)", "meg_nearpulse", "nearpulse", lambda d: d <= TOL_ME)
    check("multi-tone [far, resonant] == single resonant", "multi_hit_rev", "onres", lambda d: d <= TOL_EXACT)
    check("pcce partition_size 1 + bath pulse == standard (gCCE)", "pg_onres", "onres", lambda d: d <= TOL_PCCE)
    check("pcce partition_size 1 + bath pulse == standard (megCCE)", "pmeg_onres", "meg_onres", lambda d: d <= TOL_PCCE)
    # order 2 (pure state): the pairs act (wiDiv differs from order 1's), and reversing the bath file reverses the two
    # spins of every pair; noDiv and wiDiv must not move beyond rounding (a pulse on the wrong spin of a pair moves them)
    check("order 2 differs from order 1 (pairs act, wiDiv)", "o2", "path_vector", lambda d: d > 1e-3, ("wiDiv",),
          "max |diff| vs path_vector = {:.2e}")
    check("order 2, bath file reversed: same noDiv, wiDiv (gCCE)", "o2", "o2_rev", lambda d: d <= TOL_PERM,
          ("noDiv", "wiDiv"))
    check("order 2, bath file reversed: same noDiv, wiDiv (megCCE)", "meg_o2", "meg_o2_rev", lambda d: d <= TOL_PERM,
          ("noDiv", "wiDiv"))
    check("order 2: megCCE rate 0 with bath pulse == gCCE (noDiv, wiDiv)", "meg_o2", "o2", lambda d: d <= TOL_ME,
          ("noDiv", "wiDiv"))
    # 3-spin clusters that wrap around the rotation keep their spins in rotated order: a U_b in reverse kron order
    # (which the order-2 reversal cannot see) gives another curve
    check("order 3, bath file rotated by 5: same noDiv, wiDiv (gCCE)", "o3", "o3_rot", lambda d: d <= TOL_PERM,
          ("noDiv", "wiDiv"))
    # final review: the expm propagator with the pump, and multi-spin pcce clusters (partition_size 1 at order 2: the
    # standard pairs; the spin order inside a pair may differ, rounding only, as the reorderings above)
    check("gCCE propagator expm + bath pulse == eigen", "onres_expm", "onres", lambda d: d <= TOL_EXPM)
    check("pcce partition_size 1 + bath pulse, order 2 == standard (gCCE)", "pg_o2", "o2", lambda d: d <= TOL_PERM,
          ("noDiv", "wiDiv"))

    n = sum(ok)
    print(" {} / {} checks passed".format(n, len(ok)))
    return 0 if n == len(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
