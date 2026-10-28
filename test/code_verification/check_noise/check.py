"""method noise (src/calculation/noise.cpp) against references that do not share its code.

model  OU and power-law spectra: chi from the OU closed forms (FID, Hahn; semi-cce classical-bath-noise atten_*) and
       from scipy quadrature of 4 pi^2 int S(f) |Y(2 pi f)|^2 df (CPMG, power law), Y written here anew.
semi   against semi-cce's quantum-bath-noise (order 2): a P1 bath by mode_sum (Hahn echo) with the bath mean field,
       without it, a rotated frame and a second state, and by cumulant (exp(-kappa2/2)) for a FID; semi-cce's 13C
       example (DFT hyperfine near the NV) by mode_sum and by cumulant for a FID and CPMG-2, -4. semi-cce's
       hbar*mu0/4pi (1.05457172) is set to CCEX's H_BAR (1.054571729); the self-test keeps semi-cce's own value and has
       to fail (the two constants alone differ by up to 6e-8 in W). hamiltonian.noise.pseudo_secular against mode_sum
       with its noise operator A_zz I_z made A_z . I in the run. Needs semi-cce 2851449 or later: before it, cumulant's
       CPMG with N >= 2 was wrong (|h|^2 off by up to 0.4) and npulse 0 was not a FID.

compute the three values of hamiltonian.noise.compute on one 13C bath agree (C(t) of correlation and spectrum; the
       coherence of spectrum, by FFT, and of coherence).
pycce  compute correlation against PyCCE's noise (CCENoise): the ensemble and a sampled state.
semi   also compute spectrum against semi-cce's psd (C(t), S(f), the coherence; two sets of its options).

Usage: python3 check.py model|compute|pycce|semi WORKDIR BIN MPIRUN [SEMI_CCE]; exit 0 pass, 1 fail.
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
    base = dict(method="noise", deltat=0.002, nstep=101)   # no order, bfield, rbath, rdip, Qubit: a model reads no bath
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
    if 'seq = "FID" if npulse == 0' not in open(f"{SEMI}/quantum-bath-noise/src/quantum_bath_noise/runner.py").read():
        print(f"  FAIL  semi-cce at {SEMI} is older than 2851449 (wrong CPMG filter, no FID): update it")
        return False
    data = f"{SEMI}/quantum-bath-noise/examples/data"
    n = 18325   # the bath file's spins
    rng = np.random.default_rng(42)
    np.savetxt(f"{work}/avaax.dat", np.r_[n, rng.integers(1, 13, size=n)], fmt="%g")
    for s in (1, 2):
        np.savetxt(f"{work}/state{s}.dat", np.r_[n, rng.choice([-0.5, 0.5], size=n)], fmt="%g")
    driver = f"{work}/semi_driver.py"   # semi-cce's CLI, with hbar*mu0/4pi set (or not) before it runs
    open(driver, "w").write("""import sys, importlib
sys.path.insert(0, SEMI_SRC)
hbar = float(sys.argv.pop(1))
if hbar:
    for m in ['constants', 'hamiltonian', 'bath.p1', 'bath.nuclear', 'methods.mode_sum']:
        mod = importlib.import_module('quantum_bath_noise.' + m)
        if hasattr(mod, 'HBAR_MU0_O4PI'): mod.HBAR_MU0_O4PI = hbar
import os
if os.environ.get("NOISE_PSEUDO_SECULAR"):   # CCEX's hamiltonian.noise.pseudo_secular: the noise is A_z . I, not A_zz I_z
    import quantum_bath_noise.methods.mode_sum as mode_sum
    src = open(mode_sum.__file__).read()
    assert src.count("d_bz += bath_spin.A[2, 2] * vectors[local][2]") == 2
    src = src.replace("d_bz += bath_spin.A[2, 2] * vectors[local][2]",
                      "d_bz += sum(bath_spin.A[2, c] * vectors[local][c] for c in range(3))")
    exec(compile(src, mode_sum.__file__, "exec"), mode_sum.__dict__)
from quantum_bath_noise.cli import main   # imports the runner, which then takes the mode_sum above
main()
""".replace("SEMI_SRC", repr(SEMI + "/quantum-bath-noise/src")))
    template = open(f"{SEMI}/quantum-bath-noise/examples/config_p1_mode_sum.txt").read().replace("examples/data", data)

    def semi_run(name, s, mf=True, rot=False, hbar=1.054571729, base=None, result="Coh_mode_sum.dat", pseudo_secular=False,
                 **keys):
        c = re.sub(r"^approach\s*=.*$", "approach = single_sample", base or template, flags=re.M)
        for k, v in keys.items():
            c = re.sub(rf"^{k}\s*=.*$", f"{k} = {v}", c, flags=re.M)
        c = re.sub(r"^avaax_file\s*=.*$", f"avaax_file = {work}/avaax.dat", c, flags=re.M)
        c = re.sub(r"^state_file\s*=.*$", f"state_file = {work}/state{s}.dat", c, flags=re.M)
        c = re.sub(r"^use_mean_field\s*=.*$", f"use_mean_field = {'true' if mf else 'false'}", c, flags=re.M)
        if rot:
            c += "\nrotation = true\nrotation_bath_axis = 0, 0, 1\nrotation_qubit_axis = -1, -1, 1\n"
        open(f"{work}/{name}.txt", "w").write(c)
        r = subprocess.run([sys.executable, driver, str(hbar), "--input", f"{name}.txt", "--outdir", name], cwd=work,
                           capture_output=True, text=True, env=dict(os.environ, NOISE_PSEUDO_SECULAR="1" if pseudo_secular else ""))
        if r.returncode:
            raise SystemExit(f"  FAIL  semi-cce {name}: {r.stderr[-400:]}")
        return np.loadtxt(f"{work}/{name}/{result}")

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
    # a FID: every single cluster's chi passes 745 within the first microseconds (exp(-chi) underflows there)
    d = json.loads(json.dumps(base))
    d["configuration"]["defect_axis_file"] = [f"{work}/avaax.dat"]
    d["configuration"]["s_range"] = [1, 1]
    d.update(pulse={"count": 0}, deltat=0.0001, nstep=31)
    ccex("p1_fid", d, 4)
    _, W = coherence(f"{work}/out/p1_fid_state1_wiDiv")
    ref = semi_run("semi_p1_fid", 1, result="Coh_gauss.dat", method="cumulant", npulse=0, tau_max=0.003, tau_nstep=30)
    ok &= report("P1 cumulant, FID, state 1", np.abs(W - ref[:, 1]).max(), f" (W down to {ref[:, 1].min():.2g})")

    # semi-cce's 13C example, its DFT hyperfine file (no provenance line: x 0.5/S, as semi-cce reads it)
    open(f"{work}/gyro_13C", "w").write("13C 0.5 6.728284\n")   # semi-cce's isotope file value
    nuclear = open(f"{SEMI}/quantum-bath-noise/examples/config_nuclear.txt").read().replace("examples/data", data)
    c13 = {"method": "noise",
           "Qubit": {"qubit": [{"file": f"{data}/nuclear_defect.dat"}], "alpha_state": {"ms": [-1]}},
           "configuration": {"gyrofile": f"{work}/gyro_13C", "bathfile": [f"{data}/nuclear_bath.dat"]},
           "hamiltonian": {"common": {"hyperfine": {"dipole": "dft", "contact": True, "missing_row": "point_dipole",
                                                    "file": f"{data}/nuclear_Atensor.dat"}}},
           "order": 2, "bfield": 500, "rbath": 50, "rdip": 10}
    for N, tmax, method, result in ((1, 3.0, "mode_sum", "Coh_mode_sum.dat"), (0, 0.05, "cumulant", "Coh_gauss.dat"),
                                    (2, 3.0, "cumulant", "Coh_gauss.dat"), (4, 3.0, "cumulant", "Coh_gauss.dat")):
        name = f"c13_n{N}"
        ccex(name, dict(c13, pulse={"count": N}, deltat=tmax / 100, nstep=101), 4)
        _, W = coherence(f"{work}/out/{name}_wiDiv")
        ref = semi_run(f"semi_{name}", 1, mf=False, base=nuclear, result=result, method=method, npulse=N, tau_max=tmax,
                       tau_nstep=100)
        ok &= report(f"13C {method}, {N} pulses", np.abs(W - ref[:, 1]).max(), f" (W down to {ref[:, 1].min():.2g})")
    # hamiltonian.noise.pseudo_secular (A_zx I_x + A_zy I_y + A_zz I_z): the Hahn echo, semi-cce's mode_sum given the same noise
    d = dict(c13, pulse={"count": 1}, deltat=0.03, nstep=101)
    d["hamiltonian"] = dict(c13["hamiltonian"], noise={"pseudo_secular": True})
    ccex("c13_az", d, 4)
    _, W = coherence(f"{work}/out/c13_az_wiDiv")
    ref = semi_run("semi_c13_az", 1, mf=False, base=nuclear, pseudo_secular=True, method="mode_sum", npulse=1, tau_max=3.0,
                   tau_nstep=100)
    _, Wzz = coherence(f"{work}/out/c13_n1_wiDiv")
    ok &= report("13C mode_sum, pseudo_secular", np.abs(W - ref[:, 1]).max(),
                 f" (it changes W by up to {np.abs(W - Wzz).max():.2g})")

    # hamiltonian.noise.compute "spectrum" against semi-cce's psd (C(t), S(f) and the coherence), two sets of its options
    for name, fft in (("plain", {}), ("options", {"window": "tukey", "pad_factor": 2, "subtract_mean": True,
                                                  "clamp_negative": True})):
        d = dict(c13, pulse={"count": 1}, deltat=0.01, nstep=301)
        d["hamiltonian"] = dict(c13["hamiltonian"], noise={"compute": "spectrum",
                                                           "fft": dict({"time_max": 3.0, "nstep": 5000}, **fft)})
        ccex(f"c13_psd_{name}", d, 4)
        keys = {k: (str(v).lower() if isinstance(v, bool) else v) for k, v in fft.items()}
        semi_run(f"semi_c13_psd_{name}", 1, mf=False, base=nuclear, result="Coh_psd.dat", method="psd", npulse=1,
                 tau_max=3.0, tau_nstep=300, noise_max=3.0, noise_nstep=5000, **keys)
        sc = np.loadtxt(f"{work}/semi_c13_psd_{name}/Ct.dat"); cc = np.loadtxt(f"{work}/out/c13_psd_{name}_noise_correlation")
        sp = np.loadtxt(f"{work}/semi_c13_psd_{name}/PSD.dat"); cp = np.loadtxt(f"{work}/out/c13_psd_{name}_noise_spectrum")
        sw = np.loadtxt(f"{work}/semi_c13_psd_{name}/Coh_psd.dat")[:, 1]
        _, W = coherence(f"{work}/out/c13_psd_{name}_wiDiv")
        rc = np.abs(sc[:, 1] - cc[:, 1]).max() / np.abs(sc[:, 1]).max()
        rs = np.abs(sp[:, 1] - cp[:, 1]).max() / np.abs(sp[:, 1]).max() if len(sp) == len(cp) else np.inf
        good = rc <= 1e-10 and rs <= 1e-10   # FFT rounding: S(f) with subtract_mean is small against C(t)
        print(f"  {'PASS' if good else 'FAIL'}  13C spectrum vs psd ({name}): C(t) {rc:.1e}, S(f) {rs:.1e} relative")
        ok &= good
        ok &= report(f"13C spectrum vs psd ({name}), coherence", np.abs(sw - np.clip(W, 0, 1)).max())

    _, W = coherence(f"{work}/out/p1_state1_wiDiv")
    ref = semi_run("semi_p1_own_hbar", 1, hbar=0)
    dW = np.abs(W - ref[:, 1]).max()
    if dW > TOL:
        print(f"  PASS  self-test: semi-cce's own hbar*mu0/4pi is caught (max|dW| = {dW:.2e})")
    else:
        print(f"  FAIL  self-test: semi-cce's own hbar*mu0/4pi passed (max|dW| = {dW:.2e})")
        ok = False
    return ok


NATAB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "check_numerics", "inputs", "natab")


def natab(gyro=f"{NATAB}/gyro_13C"):
    return {"method": "noise", "Qubit": {"qubit": [{"file": f"{NATAB}/bath_Diamond_NV_natab_defect"}],
                                         "alpha_state": {"ms": [-1]}},
            "configuration": {"gyrofile": gyro, "bathfile": [f"{NATAB}/bath_Diamond_NV_natab_1"]},
            "order": 2, "bfield": 500, "rbath": 30, "rdip": 10, "pulse": {"count": 1}}


def compute():
    """the three values of hamiltonian.noise.compute on one 13C bath (ensemble): correlation and spectrum give the same
    C(t), and spectrum (by FFT) the coherence of coherence. The FFT grid has to be fine: 0.1 us over 3 ms gives 2e-10
    (the printed digits), 0.2 us 3e-10, 0.4 us 1.3e-5 (the 13C Larmor period is 1.9 us)"""
    fine = {"time_max": 3.0, "nstep": 30000}
    ccex("cmp_coherence", dict(natab(), deltat=0.002, nstep=501), 4)
    ccex("cmp_spectrum", dict(natab(), deltat=0.002, nstep=501, hamiltonian={"noise": {"compute": "spectrum", "fft": fine}}), 4)
    ccex("cmp_correlation", dict(natab(), deltat=0.0001, nstep=30001, hamiltonian={"noise": {"compute": "correlation"}}), 4)
    a = np.loadtxt(f"{work}/out/cmp_correlation_noise_correlation"); b = np.loadtxt(f"{work}/out/cmp_spectrum_noise_correlation")
    rc = np.abs(a[:, 1] - b[:, 1]).max() / np.abs(a[:, 1]).max()
    ok = rc <= 1e-12
    print(f"  {'PASS' if ok else 'FAIL'}  compute correlation and spectrum: the same C(t) ({rc:.1e} relative)")
    _, wc = coherence(f"{work}/out/cmp_coherence_wiDiv")
    _, ws = coherence(f"{work}/out/cmp_spectrum_wiDiv")
    dW = np.abs(ws - wc).max()   # held to the FFT's accuracy on this grid, 1e-9, not to the printed digits
    print(f"  {'PASS' if dW <= 1e-9 else 'FAIL'}  compute spectrum (FFT) and coherence: the coherence, max|dW| = {dW:.2e} "
          f"(<= 1e-9; W down to {wc.min():.2g})")
    return dW <= 1e-9 and ok


def pycce():
    """compute correlation against PyCCE's noise (CCENoise), the ensemble and a sampled pure state, order 2, PyCCE's
    hbar*mu0/4pi and 13C gyromagnetic ratio as CCEX's / CCEX given PyCCE's"""
    driver = f"{work}/pycce_driver.py"
    open(driver, "w").write("""import os, sys, numpy as np
import pycce as pc, pycce.constants, pycce.bath.array, pycce.h.functions
for m in (pycce.constants, pycce.bath.array, pycce.h.functions): m.HBAR_MU0_O4PI = 1.054571729
natab, out, statefile = sys.argv[1], sys.argv[2], sys.argv[3]
rows = [l.split() for l in open(natab + '/bath_Diamond_NV_natab_1').readlines()[1:] if l.strip()]
pos = np.array([[float(x) for x in r[:3]] for r in rows]); center = np.loadtxt(natab + '/bath_Diamond_NV_natab_defect', max_rows=1)
keep = np.linalg.norm(pos - center, axis=1) <= 30
b = pc.BathArray(ca=pos[keep] - center, sn=['13C'] * int(keep.sum()))
open(out + '_gyro', 'w').write('13C 0.5 %r\\n' % float(b.gyro[0]))
ts = np.arange(1501) * 0.002
for tag in ('ens', 'pure'):
    # the qubit as Sz vectors (an integer alpha is an eigenstate's index in PyCCE): ms -1 and 0
    calc = pc.Simulator(spin=1, position=[0, 0, 0], alpha=[0, 0, 1], beta=[0, 1, 0], bath=b.copy(), r_bath=30, r_dipole=10,
                        order=2, magnetic_field=[0, 0, 500])
    if tag == 'pure':   # density matrices |m><m|: PyCCE 1.1.1's CCENoise stops on state vectors
        st = np.loadtxt(statefile)[1:][keep]
        lookup = {tuple(np.round(p, 6)): s for p, s in zip(pos[keep] - center, st)}
        m = [lookup[tuple(np.round(p, 6))] for p in calc.bath.xyz]
        calc.bath.state = np.array([np.diag([1.0, 0.0]) if x > 0 else np.diag([0.0, 1.0]) for x in m])
    C = calc.compute(ts, quantity='noise', method='cce')
    np.savetxt(out + '_' + tag, np.c_[ts, C.real])   # PyCCE's couplings are in kHz: kHz^2, as CCEX's
sys.stdout.flush()
os._exit(0)   # PyCCE 1.1.1 crashes at interpreter exit on this machine (double free)
""")
    r = subprocess.run([sys.executable, driver, NATAB, f"{work}/pycce", f"{NATAB}/state_1"], capture_output=True, text=True)
    if r.returncode or not os.path.exists(f"{work}/pycce_pure"):
        print(f"  FAIL  PyCCE run: {r.stderr[-400:]}")
        return False
    ok = True
    for tag, extra in (("ens", {}), ("pure", {"s_range": [1, 1], "bath_state_file": [f"{NATAB}/state_{{s}}"]})):
        d = natab(gyro=f"{work}/pycce_gyro")
        d["configuration"].update(extra)
        d.update(deltat=0.002, nstep=1501, hamiltonian={"noise": {"compute": "correlation"}})
        ccex(f"pycce_{tag}", d, 4)
        c = np.loadtxt(f"{work}/out/pycce_{tag}{'_state1' if tag == 'pure' else ''}_noise_correlation")
        p = np.loadtxt(f"{work}/pycce_{tag}")
        rel = np.abs(c[:, 1] - p[:, 1]).max() / np.abs(p[:, 1]).max()
        good = rel <= 1e-9
        print(f"  {'PASS' if good else 'FAIL'}  compute correlation vs PyCCE noise ({'ensemble' if tag == 'ens' else 'sampled state'}): "
              f"{rel:.1e} relative")
        ok &= good
    return ok


if mode == "model":
    sys.exit(0 if model() else 1)
if mode == "compute":
    sys.exit(0 if compute() else 1)
if mode == "pycce":
    sys.exit(0 if pycce() else 1)
sys.exit(0 if semi(sys.argv[5]) else 1)
