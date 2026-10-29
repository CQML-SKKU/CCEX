"""check_pseudospin: the outputs of run.sh against cCCE order 2 and against the channel file's closed form."""
import sys, os
import numpy as np

W = sys.argv[1]
TOL = 1e-9          # the outputs carry 10 decimals
lines = []

def load(name):
    rows = [l.split() for l in open(os.path.join(W, name)) if l.strip() and not l.startswith("#")]
    return np.array([float(r[0]) for r in rows]), np.array([complex(r[1]) for r in rows])

def verdict(label, diff):
    ok = diff <= TOL
    lines.append(f"  {'PASS' if ok else 'FAIL'}  {label} : max |diff| = {diff:.2e}")

for I in ("0.5", "1.5"):
    for n in "012":
        _, c = load(f"out_c_{I}_{n}_wiDiv")
        _, b = load(f"out_b_{I}_{n}_wiDiv")
        verdict(f"spin {I} npulse {n}: block = cCCE2", abs(b - c).max())
        if I == "0.5":
            _, h = load(f"out_h_{I}_{n}_wiDiv")
            verdict(f"spin {I} npulse {n}: channel = cCCE2", abs(h - c).max())

# spin 3/2 Hahn echo: channel = L_1 L_2 (1 - (4/d) sum_c |n_a x n_b|^2 sin^2(w_a t/4) sin^2(w_b t/4))
t, h = load("out_h_1.5_1_wiDiv")
_, o = load("out_o1_wiDiv")
P = np.ones_like(t)
for Om, Da, Db in 2*np.pi*np.loadtxt(os.path.join(W, "out_h_1.5_1_pseudospin_channels"), usecols=(10, 11, 12)):
    wa, wb = np.hypot(Om, Da), np.hypot(Om, Db)
    P -= 4/16*(Om*(Db - Da))**2/(wa*wb)**2*np.sin(wa*t/4)**2*np.sin(wb*t/4)**2
verdict("spin 1.5 npulse 1: channel = L1 L2 P(channel file)", abs(h - o*P).max())

open(os.path.join(W, "check.log"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if all("PASS" in l for l in lines) else 1)
