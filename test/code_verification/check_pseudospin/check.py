"""check_pseudospin: the outputs of run.sh against the pair's cCCE term, Seo's eq. 23 and the channel file's closed form."""
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

# the pair's cCCE term is cCCE order 2 / order 1; compared as h o = c, scaled by 1 + |h| for o's last decimal (under free
# induction o falls to 5e-4 where h is 2e3)
def ratio(h, o, c):
    return (abs(h*o - c)/(1 + abs(h))).max()

# spins 1/2 and 3/2 on the axis (channels grouped, the default)
for I in ("0.5", "1.5"):
    for n in "012":
        _, c = load(f"out_c_{I}_{n}_wiDiv"); _, o = load(f"out_o_{I}_{n}_wiDiv"); _, h = load(f"out_h_{I}_{n}_wiDiv")
        verdict(f"spin {I} npulse {n}: pseudospin x cCCE1 = cCCE2", ratio(h, o, c))

# Seo et al. 2016 SI eq. 23 from the inputs alone (rad/ms, ms, angstrom; H_BAR and the electron's gamma of CCEX):
# L = 1 - K sin^2(sqrt(dA^2 + D^2) t/4) sin^2(D t/4), K = dA^2/(dA^2 + D^2), dA = A_zz,1 - A_zz,2 (alpha m_s = 1, beta 0),
# D the flip-flop splitting; spins on the z axis at 3 and 5.5 A from the qubit, 2.5 A apart
HBAR, GE, GN = 1.054571729, -17608.597050, 6.728284
dA = -2*HBAR*GE*GN*(1/3.0**3 - 1/5.5**3)
D = HBAR*GN*GN/2.5**3
t, h = load("out_h_0.5_1_wiDiv")
seo = 1 - dA**2/(dA**2 + D**2)*np.sin(np.hypot(dA, D)*t/4)**2*np.sin(D*t/4)**2
verdict("spin 0.5 npulse 1: pseudospin = Seo eq. 23", abs(h - seo).max())

# spin 3/2 Hahn echo, channels pairwise: P = 1 - (4/d) sum_c |n_a x n_b|^2 sin^2(w_a t/4) sin^2(w_b t/4)
t, h = load("out_hp_1.5_1_wiDiv")
P = np.ones_like(t)
for Om, Da, Db in 2*np.pi*np.loadtxt(os.path.join(W, "out_hp_1.5_1_pseudospin_channels"), usecols=(12, 14, 15)):
    wa, wb = np.hypot(Om, Da), np.hypot(Om, Db)
    P -= 4/16*(Om*(Db - Da))**2/(wa*wb)**2*np.sin(wa*t/4)**2*np.sin(wb*t/4)**2
verdict("spin 1.5 npulse 1, pairwise: pseudospin = P(channel file)", abs(h - P).max())

# spins 1/2 in the qubit's plane: the flip-flip is a channel of its own
for n in "012":
    _, c = load(f"out_pc_{n}_wiDiv"); _, o = load(f"out_po_{n}_wiDiv"); _, h = load(f"out_ph_{n}_wiDiv")
    verdict(f"plane, spin 0.5 npulse {n}: pseudospin x cCCE1 = cCCE2", ratio(h, o, c))

def channels(name):   # (m1, m2, m1', m2') of every channel with Omega_alpha above 1e-6 kHz
    rows = np.loadtxt(os.path.join(W, name + "_pseudospin_channels"), usecols=(8, 9, 10, 11, 12), ndmin=2)
    return [tuple(r[:4]) for r in rows if r[4] > 1e-6]

def listed(label, ok):
    lines.append(f"  {'PASS' if ok else 'FAIL'}  {label}")

listed("plane, spin 0.5: the flip-flip |1/2,1/2> <-> |-1/2,-1/2> is a channel", (0.5, 0.5, -0.5, -0.5) in channels("out_ph_1"))
double = lambda ch: [x for x in ch if abs(x[0] - x[2]) == 2 and abs(x[1] - x[3]) == 2]
listed("spin 1, axial quadrupole: no double-quantum channel", not double(channels("out_qh_Qaxial")))
listed("spin 1, eta quadrupole: double-quantum channels (|m1-m1'| = |m2-m2'| = 2)", bool(double(channels("out_qh_Qeta"))))

open(os.path.join(W, "check.log"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if all("PASS" in l for l in lines) else 1)
