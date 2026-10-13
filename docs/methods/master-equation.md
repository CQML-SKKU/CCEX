# Master-equation CCE (meCCE, megCCE)

What `method: "meCCE"` and `method: "megCCE"` do: the cluster correlation expansion with Lindblad dissipation, the bath
spins' (both methods) and the qubit's (megCCE only), on the CPU. ME-CCE is the conventional CCE with the qubit projected
out, after Onizhuk et al., *Phys. Rev. Lett.* **132**, 250401 (2024) (DOI 10.1103/PhysRevLett.132.250401, "PRL2024"
below); ME-gCCE is the gCCE with the qubit in the cluster space. Ported from the adding_MECCE branch
(`origin/adding_MECCE`: `src/simulator_mecce.cpp` as of 74a27780b, `src/simulator_megcce.cpp` as of 6022d7998). The key
is [hamiltonian.\<method\>.jump_operators](../user/InputFile_General.md#tag_jump_operators); the test is
`test/code_verification/check_master_equation/`.

## The equations

**Lindblad equation.** With the jump operators $C_k=\sqrt{\gamma_k}\,S^{(k)}_{op}$ of the list,

$$\dot\rho=-i[H,\rho]+\sum_k\Big(C_k\rho C_k^\dagger-\tfrac12\{C_k^\dagger C_k,\rho\}\Big),$$

where $S_{op}$ is $S_-$, $S_+$, $S_z$, $S_x=(S_++S_-)/2$ or $S_y=(S_+-S_-)/2i$ (`"-"`, `"+"`, `"z"`, `"x"`, `"y"`) of
the named spin, of dimension $2S+1$, and $\gamma_k$ is the `rate` in MHz times $2\pi\cdot10^3$: rad·kHz, the unit of the
Hamiltonian, so with $t$ in ms no further $2\pi$ enters. A jump on a bath name acts on every spin of a cluster with that
name, each with its own operator (no collective jumps).

**Liouville space.** $\rho$ is vectorized row by row, so that

$$\mathrm{vec}(A\rho B)=(A\otimes B^T)\,\mathrm{vec}(\rho),$$

$$\mathcal L=-i\big(H\otimes\mathbb 1-\mathbb 1\otimes H^T\big)
+\sum_k\Big(C_k\otimes C_k^*-\tfrac12\,C_k^\dagger C_k\otimes\mathbb 1-\tfrac12\,\mathbb 1\otimes(C_k^\dagger C_k)^T\Big),$$

and an interval $\tau$ of free evolution is $e^{\mathcal L\tau}$ (Eigen's matrix exponential: $\mathcal L$ is not
Hermitian, so the eigendecomposition the cCCE and gCCE use does not apply). A pulse $U$ is $U\otimes U^*$.

**ME-CCE** (PRL2024 Eq. 2). The bath operator $\tilde\rho=\langle\alpha\vert\rho\vert\beta\rangle$ evolves under the two
conditional Hamiltonians of the cCCE,

$$\dot{\tilde\rho}=-iH^{(\alpha)}\tilde\rho+i\tilde\rho H^{(\beta)}+\sum_i\gamma_i\,\mathcal D[L_i](\tilde\rho),\qquad
\mathcal L_{\alpha\beta}=-i\big(H^{(\alpha)}\otimes\mathbb 1-\mathbb 1\otimes H^{(\beta)T}\big)+\mathcal L_{\text{bath}},$$

with jumps on bath spins only. After an odd number of $\pi$ pulses $\alpha$ and $\beta$ are exchanged
($\mathcal L_{\beta\alpha}$, the dissipator unchanged). With $n$ pulses and $\tau=t/2n$, the Hahn echo is
$V=e^{\mathcal L_{\beta\alpha}\tau}e^{\mathcal L_{\alpha\beta}\tau}$, CPMG is $W^{n/2}$ ($n$ even) or
$V\,W^{(n-1)/2}$ ($n$ odd) with $W=e^{\mathcal L_{\alpha\beta}\tau}e^{\mathcal L_{\beta\alpha}\tau}\,V$, and the free
induction ($n=0$) is $e^{\mathcal L_{\alpha\beta}t}$. The cluster's coherence is $L_C(t)=\mathrm{Tr}\,\tilde\rho_C(t)$ from
$\tilde\rho_C(0)=\mathbb 1/d$ (ensemble) or $\vert\psi\rangle\langle\psi\vert$ (one bath state), and
$L=\prod_C\tilde L_C$ as in the cCCE, with $L_\emptyset=1$: there is no 0th cluster.

**ME-gCCE.** The cluster space is the qubit's times the bath's, $d=d_qd_B$, with the gCCE's full Hamiltonian $H$ (no
secular approximation). A bath jump is $\mathbb 1_q\otimes C_B$, a qubit jump (`bath_name` `"qubit"` or the qubit's name)
$C_q\otimes\mathbb 1_B$. The pulses are $U\otimes U^*$ with $U=U_q\otimes\mathbb 1_B$ at the edges of the pulse sequence,
whose fractions are followed (a Manual sequence too):

$$\mathrm{vec}\,\rho(t)=e^{\mathcal L\tau_N}\prod_{i=N-1}^{0}\big(U_i\otimes U_i^*\big)e^{\mathcal L\tau_i}\;
\mathrm{vec}\big(\rho_q(0)\otimes\rho_B(0)\big),\qquad
L_C(t)=\langle\alpha\vert\,\mathrm{Tr}_B\,\rho(t)\,\vert\beta\rangle .$$

The 0th cluster (the qubit alone, with its own jumps) is computed, every cluster is divided by it,
$L=L_\emptyset\prod_C(L_C/L_\emptyset)^{p_C}$, and the result is normalized by $L(0)$, as in the gCCE.

| | ME-CCE (`meCCE`) | ME-gCCE (`megCCE`) |
|---|---|---|
| cluster space | bath, $d_B$; Liouville $d_B^2$ | qubit × bath, $d_qd_B$; Liouville $(d_qd_B)^2$ |
| Hamiltonian | the cCCE's $H^{(\alpha)}$, $H^{(\beta)}$, with `hamiltonian.cCCE` (knight, hfmedi, delta_A) | the gCCE's (`hamiltonian.cCCE` refused) |
| jump operators | bath spins | bath spins and the qubit (one qubit) |
| 0th cluster | none ($L_\emptyset=1$) | computed, with the qubit jumps; divided out |
| pulses | ideal $\pi$ flips (X or Y, 180) from `pulse.count` (FID, Hahn, CPMG); Manual refused | the sequence as given, any rotation (Manual allowed) |
| qubits | one | one or more (a qubit jump needs one) |
| quantity | coherence | coherence |

**Memory.** One complex matrix of the Liouville space is $16N^2$ bytes. For a cluster of $n$ spin-½ bath spins,
$N=4^n$ (ME-CCE) and, with an S = 1 qubit, $N=9\cdot4^n$ (ME-gCCE):

| $n$ | ME-CCE, $N$ | one matrix | ME-gCCE (S = 1 qubit), $N$ | one matrix |
|---|---|---|---|---|
| 4 | 256 | 1 MB | 2304 | 85 MB |
| 5 | 1024 | 16 MB | 9216 | 1.36 GB |

How many such matrices one cluster holds at once (each MPI rank computes its own clusters; counted in
`simulator_mecce.cpp`, `simulator_megcce.cpp` and `superoperator.cpp`, the Hilbert-space matrices of side $\sqrt N$ left
out):

| | held through the time loop | added at a time point | peak |
|---|---|---|---|
| ME-CCE | $\mathcal L$ and the Liouvillian with $\alpha$, $\beta$ exchanged (2) | the exponentials before and after a $\pi$ pulse, the CPMG products | 10 (FID), 11 (Hahn, CPMG) |
| ME-gCCE | $\mathcal L$ (1) | one $e^{\mathcal L\tau}$ per distinct interval of the sequence ($k$: 1 for FID and Hahn, 2 for CPMG, more for a merged DEER sequence with unequal intervals), the pulse superoperator, the running product, two product temporaries | $k+10$ at most: 10 (FID, Hahn), 12 (CPMG) |

The peak is reached while an exponential is computed: Eigen's Padé-13 scaling and squaring (`matrix_exp_pade13` in
`external/eigen/unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h`, Eigen 3.4.1) holds 8 matrices of that size
(the argument $\mathcal L\tau$, the scaled argument, its 2nd, 4th and 6th powers, $U$, $V$ and one work matrix) on top
of what the cluster holds. An ME-gCCE CPMG cluster of 4 bath spins with an S = 1 qubit thus needs about
12 × 85 MB ≈ 1 GB, of 5 spins 12 × 1.36 GB ≈ 16 GB. The time per step is that of a dense exponential of an
$N\times N$ matrix: these methods are meant for small clusters.

## Changes from the branch

- **meCCE refuses qubit jump operators** (the message names megCCE). The branch added a qubit jump to every cluster's
  Liouvillian as a secular scalar $c\,\mathbb 1$ and had no 0th cluster to divide it out, so the factor $e^{ct}$ entered
  once per cluster power: at CCE2 with $N$ spins and $P$ pairs, $e^{(N-P)ct}$ instead of $e^{ct}$ (the 2026-10-13 rebuild
  of the branch gave 0.25 against the right factor on its own test system, 5.5e-9 against $e^{(N-P)ct}$). In the
  projected picture only `z` and, for S = ½, `+`/`-` reduce to such a scalar at all; `x` and `y` flip the qubit
  ($H^{(\alpha)}\leftrightarrow H^{(\beta)}$) and cannot be written there. PyCCE's secular ME-CCE does not take qubit
  jumps either. Pure dephasing of the qubit ($C=\sqrt\gamma\,S_z$) is a factor that does commute with the projected
  evolution: multiply the meCCE result by $e^{-\gamma(m_\alpha-m_\beta)^2t/2}$ ($e^{-\gamma t/2}$ for S = ½), which PRL2024
  p. 1 calls "trivially factored out from the coherence". Any other qubit jump: use megCCE.
- What the kernels stopped on at run time, or ignored, is refused when the input is read (`validateMasterEquation`,
  `src/main.cpp`): a Manual sequence or a rotation other than an ideal $\pi$ flip for meCCE (the branch used `npulse`
  alone and dropped the fractions and angles), more than one qubit for meCCE or for a qubit jump, quantity `dm`,
  savemode `cluster_density_matrix` (the branch's full density-matrix output is not ported), `pulse.per_qubit`,
  `pulse.detuning_factor` ≠ 1, `propagator` other than `eigen`, `cluster_mode` other than `standard` and `pcce`,
  `evolution_update: "recursive"` (the default is `exact` for both), `rdip_ladder`; the branch's top-level
  `jump_operators` and its old operator names `+-`/`-+`, each with the new form named. Once the bath is read
  (`validateJumpOperatorNames`), a `bath_name` that is neither the qubit nor the name of a bath spin of the run (bath
  file and Defect sub-spins, case ignored) is refused with the names there are: the branch matched it to nothing and
  ran that jump at rate 0. A `bath_name` that is both the qubit's name and a bath spin's is refused too (the kernel would
  take it as a qubit jump only).
- The branch's float qubit gyromagnetic ratio and float `deltat` are double in this code (fixed before the port: the
  qubit's gyro in cc2d2c8f2, `Config::deltat` a double).
- DEER bath pulses (`pulse.bath`) run with megCCE as on the branch, $U_q\otimes U_b$ in the pulse superoperator, with
  the gCCE's $U_b$ (`bathPulseOperator`, [deer.md](deer.md)) instead of the branch's own copy; meCCE refuses them.

## Verification

`check_master_equation` (measured 2026-10-13, WSL, g++ with MKL, the CPMG, pure-state and two-qubit rows 2026-10-14; the
suite's README has the cases):

| check | tol | max \|diff\| |
|---|---|---|
| meCCE rate 0 == cCCE, Hahn `_noDiv` / `_wiDiv` | 1e-6 | 7.0e-10 / 7.0e-10 |
| meCCE rate 0 == cCCE, Ramsey (complex) | 1e-8 | 8.1e-10 |
| meCCE rate 0 == cCCE, CPMG with 2 / 3 pulses (complex) | 1e-6 | 7.0e-10 / 1.1e-9 |
| meCCE rate 0 == cCCE, a pure bath state (`check_numerics` p1_cce2_single, complex) | 1e-6 | 4.5e-8 |
| meCCE rate 0 == cCCE with `hamiltonian.cCCE.knight` true (natab, Hahn) | 1e-6 | 0 |
| meCCE bath jumps == the branch's meCCE output | 1e-6 | 1.9e-8 |
| meCCE bath jumps == PyCCE (S = ½ test system) | 1e-6 | 3.6e-9 |
| meCCE, one bath spin, FID == PRL2024 Eq. 6, $2\gamma<a$ / $2\gamma>a$ | 1e-8 | 2.4e-10 / 2.7e-10 |
| megCCE rate 0 == gCCE | 1e-6 | 4.4e-9 |
| megCCE rate 0 == gCCE, two NV qubits (S = 1, two ¹³C, one state) | 1e-6 | 1.0e-10 |
| gCCE == conj(PyCCE) (S = ½) | 1e-6 | 1.9e-8 |
| megCCE bath jumps == conj(PyCCE) (S = ½) | 1e-6 | 5.4e-9 |
| megCCE bath jumps and a qubit `z` jump == conj(PyCCE) (S = ½) | 1e-6 | 4.2e-9 |
| megCCE order 0, qubit `z` / `-` == $e^{-\gamma t/2}$ | 1e-9 | 2.7e-10 / 2.7e-10 |

PRL2024 Eq. 6 is the free induction of one bath spin on the qubit axis ($A_{zx}=A_{zy}=0$, $a=\vert A_{zz}\vert$) at
infinite temperature with $L_\pm=I_\pm$ at the same rate $\gamma$, order 1:

$$L(t)=e^{-\gamma t}\Big(\cosh\frac{\tilde\omega t}{2}+\frac{2\gamma}{\tilde\omega}\sinh\frac{\tilde\omega t}{2}\Big),
\qquad \tilde\omega=\sqrt{4\gamma^2-a^2}.$$

The order-0 checks are the qubit alone (S = ½, B ∥ z, Hahn echo): `z` and `-` (amplitude damping) both damp
$\rho_{\alpha\beta}$ at $\gamma/2$ before and after the $\pi$ pulse, so $L=e^{-\gamma t/2}$; a 0th cluster without the
qubit jump would give 1.

## Relation to PyCCE

PyCCE 1.1.1 has the same three calculations (the conventional CCE with Lindblad jumps on the bath spins, gCCE, and
`LindbladgCCE`); the CCEX steps (Zeeman, dipolar couplings, the initial state, $u(t/2)\,R\,u(t/2)$, the cluster product
and power, the normalization, the row-major dissipator, the jump rate) are the same formulas. Two conventions differ:

1. **The measured element.** CCEX measures $\langle\alpha\vert\rho_q\vert\beta\rangle$ and PyCCE
   $\langle R\alpha\vert\rho_q\vert R\beta\rangle$, $\alpha,\beta$ turned by the pulses. After an odd number of $\pi$
   pulses the two are complex conjugates: gCCE and ME-gCCE give conj(PyCCE), with the same $\vert L\vert$. (The PyCCE paper,
   *Adv. Theory Simul.* **4**, 2100254 (2021), itself writes both, Eq. 5 and Eq. 20.) The ME-CCE values of the test
   system (Hahn echo) have no imaginary part, so there the two cannot be told apart.
2. **The $\pi$ pulse of an S ≥ 1 qubit.** CCEX applies the unitary pulse of the two-level subspace,
   $e^{-i\pi\sigma_x/2}=(\mathbb 1-P)-i\sigma_x$ ($P$ the projector on $\{\alpha,\beta\}$), PyCCE $-i\sigma_x$, which
   removes what is outside the subspace. They agree for S = ½; with leakage to the third level they differ (the branch's
   S = 1 test system: 0.124 for gCCE, 0.081 for ME-gCCE). The external comparisons here therefore use S = ½.
3. The disagreements reported in the branch's time (0.28 for gCCE, 0.052 for ME-gCCE on the S = ½ system) came from its
   float qubit gyromagnetic ratio and float `deltat` on a resonant system (qubit and bath electrons at one Larmor
   frequency); in double precision they are gone (the table above). For the same reason the branch's own gCCE and ME-gCCE
   outputs are not used as references.

## Basis and scope

- ME-CCE is PRL2024 Eqs. 1–4 (the Lindblad equation of the bath, the projected equation of motion Eq. 2 and the CCE
  expansion Eqs. 3–4), with the coherence $\mathrm{Tr}\,\tilde\rho(t)/\mathrm{Tr}\,\tilde\rho(0)$; Eq. 6 is checked
  above.
- Qubit jumps are an extension beyond PRL2024 (which omits "any direct dephasing of the central spin") and are taken by
  megCCE only, where the 0th cluster carries them and every cluster is divided by it, as PyCCE's `LindbladgCCE` does.
- ME-gCCE is not in the paper; it has the structure of PyCCE 1.1.1 `LindbladgCCE`.
- Not supported: collective jumps (PRL2024's $L_{i,\pm}=I^i_\mp\otimes I^{i+1}_\pm$) and the non-Markovian master equation
  of its Supplemental Material. The branch had neither.
- meCCE reads `hamiltonian.cCCE` (`knight`, `hfmedi`, `delta_A`) as the cCCE does and builds the same conditional
  Hamiltonians (`cceClusterHamiltonians`), so at rate 0 it is the cCCE also with those terms. `hfmedi` was made the
  second-order term of PyCCE2021 Eq. 12 (all ordered pairs and $i=j$, every $\beta\ne\alpha$) upstream in 9b0e70b87, and
  meCCE inherits it. megCCE refuses `hamiltonian.cCCE`, as gCCE does.

## Supported and refused

| | meCCE | megCCE |
|---|---|---|
| `cluster_mode` | `standard`, `pcce` (inside pmeCCE; `pcce.cache` refused) | `standard`, `pcce` (inside pmegCCE; same) |
| backend | CPU (`auto` runs on the CPU; `cuda` refused) | CPU (same) |
| `propagator` | `eigen` (the name; the ME step is a matrix exponential); `sparse`, `expm` refused | same |
| `evolution`, `evolution_update` | not used; `vector` and `recursive` refused | same |
| pulses | `pulse.count` 0, 1, n of ideal $\pi$ flips (other rotations and Manual refused); `pulse.bath` refused | any sequence and rotation, Manual included (`per_qubit` false, `detuning_factor` 1); DEER bath pulses `pulse.bath` ([deer.md](deer.md)) |
| jump operators | bath names | bath names, `"qubit"` or the qubit's name (one qubit) |
| `hamiltonian` | `common`, `cCCE`, `meCCE` | `common`, `megCCE` |
| output | coherence, savemode `coherence`, `cluster_coherence` | same |
| refused | qubit jumps, a `bath_name` no bath spin has, two qubits, quantity `dm`, savemode `cluster_density_matrix`, `rdip_ladder`, `-P` | a `bath_name` that is neither the qubit nor a bath spin, two qubits with a qubit jump, quantity `dm`, savemode `cluster_density_matrix`, `rdip_ladder`, `-P`, `hamiltonian.cCCE` |

A `bath_name` is checked against the whole bath once it is read: one that is neither the qubit nor the name of any bath
spin of the run (bath file and Defect sub-spins) is refused, with the names the bath has, and so is one that is both the
qubit's name and a bath spin's. A cluster without that species
gets no jump from it, which is the normal case. No jump operators, or an empty list, is rate 0: meCCE then gives the
cCCE and megCCE the gCCE.

## Where it is in the code

| piece | file |
|---|---|
| reading `hamiltonian.meCCE` / `.megCCE.jump_operators`, the old key | `src/input/option.cpp` (`readJumpOperators`, `refuseOldInputKeys`) |
| the refusals | `src/main.cpp` (`validateMasterEquation`, `validateHamiltonianForMethod`, `validateJumpOperatorNames`), `src/input/option.cpp` (`internalMethod`) |
| superoperators, jump operators | `src/calculation/superoperator.cpp` |
| the ME-CCE kernel | `src/calculation/simulator_mecce.cpp` (`calCoherenceMecce`) |
| the ME-gCCE kernel | `src/calculation/simulator_megcce.cpp` (`calCoherenceMegcce`; `ba == NULL` is the 0th cluster) |
| dispatch, product | `src/calculation/simulator.cpp` (`isMECCE`, `isMEGCCE`, `qubitInMatrix`) |
