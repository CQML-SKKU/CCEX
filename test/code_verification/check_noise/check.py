"""method noise (src/calculation/noise.cpp) against references that do not share its code.

model  OU and power-law spectra: chi from the OU closed forms (FID, Hahn; semi-cce classical-bath-noise atten_*) and
       from scipy quadrature of 4 pi^2 int S(f) |Y(2 pi f)|^2 df (CPMG, power law), Y written here anew.
semi   a P1 bath against semi-cce's quantum-bath-noise mode_sum (Hahn echo, order 2), the same bath, state and axis
       files: with the bath mean field, without it, a rotated frame and a second state. semi-cce's hbar*mu0/4pi
       (1.05457172) is set to CCEX's H_BAR (1.054571729) for the comparison; the self-test keeps semi-cce's own
       value and has to fail (the two constants alone differ by up to 6e-8 in W).

Usage: python3 check.py model|semi WORKDIR BIN MPIRUN [SEMI_CCE]; exit 0 pass, 1 fail.
CCEX writes W with 10 decimals, so W is compared to TOL = 2e-10.
"""
import json, os, re, subprocess, sys
import numpy as np
from scipy.integrate import quad

TOL = 2e-10
mode, work, BIN, MPIRUN = sys.argv[1:5]
os.makedirs(work, exist_ok=True)


def ccex(name, d, np_=1):
    json.dump(d, open(f"{work}/{name}.json", "w"), indent=1)
    r = subprocess.run([MPIRUN, "-n", str(np_), BIN, "-f", f"{name}.json", "-o", "out"], cwd=work,
                       capture_output=True, text=True, stdin=subprocess.DEVNULL)
    open(f"{work}/proc_{name}", "w").write(r.stdout + r.stderr)
    if r.returncode:
        raise SystemExit(f"  FAIL  {name}: CCEX exit {r.returncode}")


def coherence(path):
    rows = [l.split() for l in open(path) if l.strip()]
    return np.array([float(r[0]) for r in rows]), np.array([complex(r[1]).real for r in rows])


def report(name, dW, extra=""):
    ok = dW <= TOL
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: max|dW| = {dW:.2e}{extra}")
    return ok


def Y2(w, T, N):   # |int_0^T f e^{iwt} dt|^2 of an ideal CPMG-N (N = 0: FID)
    edges = [0.0] + [T * (2 * k - 1) / (2 * N) for k in range(1, N + 1)] + [T]
    y = 0j
    for j in range(len(edges) - 1):
        a, b = edges[j], edges[j + 1]
        y += (-1) ** j * ((b - a) if w == 0 else (np.exp(1j * w * b) - np.exp(1j * w * a)) / (1j * w))
    return abs(y) ** 2


def chi_quad(S, T, N, fmin, fmax):
    if T == 0:
        return 0.0
    edges = np.r_[fmin, np.arange(fmin, fmax, 0.5 / T)[1:], fmax]
    return 4 * np.pi ** 2 * sum(quad(lambda f: S(f) * Y2(2 * np.pi * f, T, N), a, b, epsabs=0, epsrel=1e-12,
                                     limit=200)[0] for a, b in zip(edges[:-1], edges[1:]) if b > a)


def model():
    base = dict(order=1, bfield=0, rbath=1, rdip=1, Qubit={"qubit": [{"xyz": [0, 0, 0]}]}, method="noise",
                deltat=0.002, nstep=101)
    ok = True
    sig, tc = 10.0, 0.05
    b = 2 * np.pi * sig
    for N in (0, 1, 4):
        name = f"ou_n{N}"
        ccex(name, dict(base, hamiltonian={"noise": {"spectrum": {"model": "OU", "sigma": sig, "tau_c": tc}}},
                        pulse={"count": N}))
        t, W = coherence(f"{work}/out/{name}_wiDiv")
        if N == 0:
            ref = (b * tc) ** 2 * (t / tc + np.exp(-t / tc) - 1)
        elif N == 1:
            x = t / tc
            ref = (b * tc) ** 2 * (x - 3 + 4 * np.exp(-x / 2) - np.exp(-x))
        else:
            t, W = t[::10], W[::10]
            ref = np.array([chi_quad(lambda f: 2 * sig ** 2 * tc / (1 + (2 * np.pi * f * tc) ** 2), T, N, 0, 2e4)
                            for T in t])
        ok &= report(f"OU, {N} pulses", np.abs(W - np.exp(-ref)).max(), f" (W down to {W.min():.2g})")
    A, beta, fmin, fmax = 50.0, 1.0, 0.01, 200.0
    for N in (0, 1, 2):
        name = f"power_n{N}"
        ccex(name, dict(base, hamiltonian={"noise": {"spectrum": {"model": "power", "amplitude": A, "beta": beta,
                                                                  "f_min": fmin, "f_max": fmax}}},
                        pulse={"count": N}))
        t, W = coherence(f"{work}/out/{name}_wiDiv")
        t, W = t[::10], W[::10]
        ref = np.array([chi_quad(lambda f: A * f ** -beta, T, N, fmin, fmax) for T in t])
        ok &= report(f"1/f, {N} pulses", np.abs(W - np.exp(-ref)).max(), f" (W down to {W.min():.2g})")
    return ok


def semi(SEMI):
    data = f"{SEMI}/quantum-bath-noise/examples/data"
    n = 18325   # the bath file's spins
    rng = np.random.default_rng(42)
    np.savetxt(f"{work}/avaax.dat", np.r_[n, rng.integers(1, 13, size=n)], fmt="%g")
    for s in (1, 2):
        np.savetxt(f"{work}/state{s}.dat", np.r_[n, rng.choice([-0.5, 0.5], size=n)], fmt="%g")
    driver = f"{work}/semi_driver.py"   # semi-cce's CLI, with hbar*mu0/4pi set (or not) before it runs
    open(driver, "w").write(f"""import sys, importlib
sys.path.insert(0, {SEMI + '/quantum-bath-noise/src'!r})
hbar = float(sys.argv.pop(1))
if hbar:
    for m in ['constants', 'hamiltonian', 'bath.p1', 'methods.mode_sum']:
        mod = importlib.import_module('quantum_bath_noise.' + m)
        if hasattr(mod, 'HBAR_MU0_O4PI'): mod.HBAR_MU0_O4PI = hbar
from quantum_bath_noise.cli import main
main()
""")
    template = open(f"{SEMI}/quantum-bath-noise/examples/config_p1_mode_sum.txt").read().replace("examples/data", data)

    def semi_run(name, s, mf=True, rot=False, hbar=1.054571729):
        c = re.sub(r"^approach\s*=.*$", "approach = single_sample", template, flags=re.M)
        c = re.sub(r"^avaax_file\s*=.*$", f"avaax_file = {work}/avaax.dat", c, flags=re.M)
        c = re.sub(r"^state_file\s*=.*$", f"state_file = {work}/state{s}.dat", c, flags=re.M)
        c = re.sub(r"^use_mean_field\s*=.*$", f"use_mean_field = {'true' if mf else 'false'}", c, flags=re.M)
        if rot:
            c += "\nrotation = true\nrotation_bath_axis = 0, 0, 1\nrotation_qubit_axis = -1, -1, 1\n"
        open(f"{work}/{name}.txt", "w").write(c)
        r = subprocess.run([sys.executable, driver, str(hbar), "--input", f"{name}.txt", "--outdir", name], cwd=work,
                           capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"  FAIL  semi-cce {name}: {r.stderr[-400:]}")
        return np.loadtxt(f"{work}/{name}/Coh_mode_sum.dat")

    det = [84.946496] * 3 + [0.0] * 3 + [-84.946496] * 3 + [114.0264, 0.0, -114.0264]
    base = {"method": "noise",
            "Qubit": {"qubit": [{"file": f"{data}/p1_defect.dat"}], "alpha_state": {"ms": [-1]}},
            "configuration": {"gyrofile": f"{data}/p1_gyro.dat", "bathfile": [f"{data}/p1_bath.dat"],
                              "s_range": [1, 2], "bath_state_file": [f"{work}/state{{s}}.dat"]},
            "order": 2, "bfield": 500, "rbath": 200, "rdip": 180, "deltat": 0.00012, "nstep": 101,
            "pulse": {"count": 1},
            "Defect": {"species": [{"name": "P1", "n_axis_states": 12,
                                    "detuning": [[i + 1, "e", v] for i, v in enumerate(det)]}]}}
    ok = True
    for name, mf, rot in (("p1", True, False), ("p1_nomf", False, False), ("p1_rot", True, True)):
        d = json.loads(json.dumps(base))
        d["configuration"]["defect_axis_file"] = [f"{work}/avaax.dat"]   # the same axes for both states
        if not mf:
            d["hamiltonian"] = {"common": {"bath_bath_mean_field": "off"}}
        if rot:
            d["frame"] = {"rotation": {"qubit_axis": [-1, -1, 1]}}
        ccex(name, d, 4)
        for s in ((1, 2) if name == "p1" else (1,)):
            _, W = coherence(f"{work}/out/{name}_state{s}_wiDiv")
            ref = semi_run(f"semi_{name}_s{s}", s, mf, rot)
            ok &= report(f"P1 mode_sum, {name}, state {s}", np.abs(W - ref[:, 1]).max(), f" (W down to {ref[:, 1].min():.2g})")
    _, W = coherence(f"{work}/out/p1_state1_wiDiv")
    ref = semi_run("semi_p1_own_hbar", 1, hbar=0)
    dW = np.abs(W - ref[:, 1]).max()
    if dW > TOL:
        print(f"  PASS  self-test: semi-cce's own hbar*mu0/4pi is caught (max|dW| = {dW:.2e})")
    else:
        print(f"  FAIL  self-test: semi-cce's own hbar*mu0/4pi passed (max|dW| = {dW:.2e})")
        ok = False
    return ok


if mode == "model":
    sys.exit(0 if model() else 1)
sys.exit(0 if semi(sys.argv[5]) else 1)
