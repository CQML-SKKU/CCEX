# Master-equation CCE (meCCE, megCCE)

`method` meCCE and megCCE run the cluster-correlation expansion with Lindblad dissipation, on the CPU. The jump
operators are given in [hamiltonian.\<method\>.jump_operators](../user/InputFile_General.md#tag_jump_operators), and
`test/code_verification/check_master_equation/` compares both methods with the cCCE and gCCE at rate 0, with PyCCE and
with closed forms.

1. [The two methods](#me_methods): what each one keeps of the qubit.
2. [The equations](#me_equations): the Lindblad equation, the Liouville space and the two propagations.
3. [Qubit jumps](#me_qubit): why only megCCE takes them.
4. [Memory](#me_memory): the size of a cluster in the Liouville space.
5. [What runs and what is refused](#me_support): the options each method takes.
6. [Conventions against PyCCE](#me_pycce): the two differences.

<a id="me_methods"></a>
## 1. The two methods

- **ME-CCE** (`meCCE`) is the conventional CCE with the qubit projected out, after Onizhuk et al., PRL 132, 250401
  (2024) ("PRL2024" below). Jumps act on the bath spins.
- **ME-gCCE** (`megCCE`) is the gCCE with the qubit in the cluster space. Jumps act on the bath spins and on the qubit.

| | ME-CCE (`meCCE`) | ME-gCCE (`megCCE`) |
|---|---|---|
| Cluster space | The bath, $d_B$, and $d_B^2$ in Liouville space. | The qubit times the bath, $d_qd_B$, and $(d_qd_B)^2$ in Liouville space. |
| Hamiltonian | The cCCE's $H^{(\alpha)}$ and $H^{(\beta)}$, with `hamiltonian.cCCE` (knight, hfmedi, delta_A). | The gCCE's. `hamiltonian.cCCE` is refused. |
| Jump operators | On bath spins. | On bath spins and on the qubit (one qubit). |
| 0th cluster | None ($L_\emptyset=1$). | Computed with the qubit jumps and divided out. |
| Pulses | Ideal π flips (X or Y, 180°) from `pulse.count` (FID, Hahn, CPMG). | The sequence as given, any rotation, Manual included. |
| Qubits | One. | One or more. A qubit jump needs one. |
| Quantity | Coherence. | Coherence. |

With no jump operators, or an empty list, the rate is 0. meCCE then gives the cCCE, and megCCE the gCCE.

<a id="me_equations"></a>
## 2. The equations

**Lindblad equation.** With the jump operators $C_k=\sqrt{\gamma_k}\,S^{(k)}_{op}$ of the list,

$$\dot\rho=-i[H,\rho]+\sum_k\Big(C_k\rho C_k^\dagger-\tfrac12\{C_k^\dagger C_k,\rho\}\Big).$$

- $S_{op}$ is $S_-$, $S_+$, $S_z$, $S_x=(S_++S_-)/2$ or $S_y=(S_+-S_-)/2i$ (`"-"`, `"+"`, `"z"`, `"x"`, `"y"`) of the
  named spin, of dimension $2S+1$.
- $\gamma_k$ is the `rate` in MHz times $2\pi\cdot10^3$. That is rad·kHz, the unit of the Hamiltonian, so with $t$ in ms
  no further $2\pi$ enters.
- A jump on a bath name acts on every spin of a cluster with that name, each with its own operator. There are no
  collective jumps.

**Liouville space.** $\rho$ is vectorized row by row, so that

$$\mathrm{vec}(A\rho B)=(A\otimes B^T)\,\mathrm{vec}(\rho),$$

$$\mathcal L=-i\big(H\otimes\mathbb 1-\mathbb 1\otimes H^T\big)
+\sum_k\Big(C_k\otimes C_k^*-\tfrac12\,C_k^\dagger C_k\otimes\mathbb 1-\tfrac12\,\mathbb 1\otimes(C_k^\dagger C_k)^T\Big).$$

An interval $\tau$ of free evolution is $e^{\mathcal L\tau}$, computed with Eigen's matrix exponential. $\mathcal L$ is
not Hermitian, so the eigendecomposition used by the cCCE and gCCE does not apply. A pulse $U$ is $U\otimes U^*$.

**ME-CCE** (PRL2024 Eq. 2). The bath operator $\tilde\rho=\langle\alpha\vert\rho\vert\beta\rangle$ evolves under the two
conditional Hamiltonians of the cCCE,

$$\dot{\tilde\rho}=-iH^{(\alpha)}\tilde\rho+i\tilde\rho H^{(\beta)}+\sum_i\gamma_i\,\mathcal D[L_i](\tilde\rho),\qquad
\mathcal L_{\alpha\beta}=-i\big(H^{(\alpha)}\otimes\mathbb 1-\mathbb 1\otimes H^{(\beta)T}\big)+\mathcal L_{\text{bath}},$$

with jumps on bath spins only. After an odd number of $\pi$ pulses, $\alpha$ and $\beta$ are exchanged
($\mathcal L_{\beta\alpha}$, the dissipator unchanged). With $n$ pulses and $\tau=t/2n$:

- Hahn echo: $V=e^{\mathcal L_{\beta\alpha}\tau}e^{\mathcal L_{\alpha\beta}\tau}$.
- CPMG: $W^{n/2}$ ($n$ even) or $V\,W^{(n-1)/2}$ ($n$ odd), with
  $W=e^{\mathcal L_{\alpha\beta}\tau}e^{\mathcal L_{\beta\alpha}\tau}\,V$.
- Free induction ($n=0$): $e^{\mathcal L_{\alpha\beta}t}$.

The cluster's coherence is $L_C(t)=\mathrm{Tr}\,\tilde\rho_C(t)$, starting from $\tilde\rho_C(0)=\mathbb 1/d$ (ensemble)
or $\vert\psi\rangle\langle\psi\vert$ (one bath state). The product is $L=\prod_C\tilde L_C$ as in the cCCE, with
$L_\emptyset=1$, since there is no 0th cluster.

**ME-gCCE.** The cluster space is the qubit's times the bath's, $d=d_qd_B$, with the gCCE's full Hamiltonian $H$ (no
secular approximation).

- A bath jump is $\mathbb 1_q\otimes C_B$.
- A qubit jump (`bath_name` `"qubit"` or the qubit's name) is $C_q\otimes\mathbb 1_B$.
- The pulses are $U\otimes U^*$ with $U=U_q\otimes\mathbb 1_B$, at the edges of the pulse sequence. The fractions of the
  sequence are followed, a Manual sequence too. DEER bath pulses (`pulse.bath`) add $U_b$ to the pulse,
  $U_q\otimes U_b$ ([DEER bath pulses](deer.md)).

$$\mathrm{vec}\,\rho(t)=e^{\mathcal L\tau_N}\prod_{i=N-1}^{0}\big(U_i\otimes U_i^*\big)e^{\mathcal L\tau_i}\;
\mathrm{vec}\big(\rho_q(0)\otimes\rho_B(0)\big),\qquad
L_C(t)=\langle\alpha\vert\,\mathrm{Tr}_B\,\rho(t)\,\vert\beta\rangle .$$

The 0th cluster (the qubit alone, with its own jumps) is computed and every cluster is divided by it,
$L=L_\emptyset\prod_C(L_C/L_\emptyset)^{p_C}$. The result is normalized by $L(0)$, as in the gCCE.

<a id="me_qubit"></a>
## 3. Qubit jumps

meCCE refuses jumps on the qubit, and the message names megCCE. In the projected picture of meCCE only `z` and, for
S = ½, `+` and `-` would reduce to a scalar factor, and `x` and `y` flip the qubit ($H^{(\alpha)}\leftrightarrow
H^{(\beta)}$), which the projection cannot hold. PyCCE's secular ME-CCE does not take qubit jumps either.

Pure dephasing of the qubit ($C=\sqrt\gamma\,S_z$) commutes with the projected evolution. Multiply the meCCE result by
$e^{-\gamma(m_\alpha-m_\beta)^2t/2}$, which is $e^{-\gamma t/2}$ for S = ½. PRL2024 calls this factor "trivially
factored out from the coherence". For any other qubit jump, use megCCE, whose 0th cluster carries the qubit jumps, as
PyCCE's `LindbladgCCE` does.

<a id="me_memory"></a>
## 4. Memory

One complex matrix of the Liouville space is $16N^2$ bytes. For a cluster of $n$ spin-½ bath spins, $N=4^n$ (ME-CCE),
and $N=9\cdot4^n$ with an S = 1 qubit (ME-gCCE):

| $n$ | ME-CCE, $N$ | One matrix | ME-gCCE (S = 1 qubit), $N$ | One matrix |
|---|---|---|---|---|
| 4 | 256 | 1 MB | 2304 | 85 MB |
| 5 | 1024 | 16 MB | 9216 | 1.36 GB |

A cluster holds about 10 to 12 such matrices at its peak, while a matrix exponential is computed. An ME-gCCE CPMG
cluster with an S = 1 qubit thus needs about 1 GB for 4 bath spins and 16 GB for 5. Each step costs a dense exponential
of an $N\times N$ matrix, so these methods are meant for small clusters. Each MPI rank computes its own clusters.

<a id="me_support"></a>
## 5. What runs and what is refused

| | meCCE | megCCE |
|---|---|---|
| `cluster.mode` | `standard` or `pcce` (`cluster.pcce.cache` refused). | The same. |
| `execution.backend` | The CPU. `auto` runs on the CPU and `cuda` is refused. | The same. |
| `execution.propagator` | `eigen` only, as a name (the step is a matrix exponential). | The same. |
| `execution.evolution`, `evolution_update` | Not used. `vector` and `recursive` are refused. | The same. |
| Pulses | `pulse.count` 0, 1 or n ideal π flips. Other rotations, Manual and `pulse.bath` are refused. | Any sequence and rotation, Manual included, and `pulse.bath`. `per_qubit` must be false and `detuning_factor` 1. |
| Jump operators | Bath names. | Bath names, `"qubit"` or the qubit's name (one qubit). |
| `hamiltonian` | `common`, `cCCE` and `meCCE`. | `common` and `megCCE`. |
| Output | Coherence, with `output.savemode` `coherence` or `cluster_coherence`. | The same. |
| Also refused | Qubit jumps and two qubits. | Two qubits with a qubit jump, and `hamiltonian.cCCE`. |

Both also refuse quantity `dm`, `output.savemode` `cluster_density_matrix`, `rdip_ladder` and `output.cluster_only`.

A `bath_name` is checked against the whole bath once it is read. One that is neither the qubit nor the name of any bath
spin of the run (bath file and Defect sub-spins, case ignored) is refused with the names the bath has, and so is one
that is both the qubit's name and a bath spin's. A cluster without that species gets no jump from it, which is the
normal case.

meCCE builds the same conditional Hamiltonians as the cCCE, with `knight`, `hfmedi` and `delta_A`, so at rate 0 it is
the cCCE with those terms as well.

<a id="me_pycce"></a>
## 6. Conventions against PyCCE

PyCCE 1.1.1 has the same three calculations, the conventional CCE with Lindblad jumps on the bath spins, gCCE, and
`LindbladgCCE`, with the same formulas. Two conventions differ.

1. **The measured element.** CCEX measures $\langle\alpha\vert\rho_q\vert\beta\rangle$ and PyCCE
   $\langle R\alpha\vert\rho_q\vert R\beta\rangle$, with $\alpha,\beta$ turned by the pulses. After an odd number of
   $\pi$ pulses the two are complex conjugates, so gCCE and ME-gCCE give conj(PyCCE) with the same $\vert L\vert$. The
   PyCCE paper writes both (Eq. 5 and Eq. 20).
2. **The π pulse of an S ≥ 1 qubit.** CCEX applies the unitary pulse of the two-level subspace,
   $e^{-i\pi\sigma_x/2}=(\mathbb 1-P)-i\sigma_x$, with $P$ the projector on $\{\alpha,\beta\}$. PyCCE applies
   $-i\sigma_x$, which removes what is outside the subspace. They agree for S = ½ and differ with leakage to the third
   level.

Not supported are the collective jumps of PRL2024 ($L_{i,\pm}=I^i_\mp\otimes I^{i+1}_\pm$) and the non-Markovian master
equation of its Supplemental Material. ME-gCCE is not in PRL2024. It has the structure of PyCCE's `LindbladgCCE`.

<a id="me_references"></a>
## References

- M. Onizhuk, Y.-X. Wang, J. Nagura, A. A. Clerk and G. Galli, Understanding central spin decoherence due to
  interacting dissipative spin baths, Phys. Rev. Lett. 132, 250401 (2024), doi:10.1103/PhysRevLett.132.250401.
- M. Onizhuk and G. Galli, PyCCE: A Python package for cluster correlation expansion simulations of spin qubit
  dynamics, Adv. Theory Simul. 4, 2100254 (2021), doi:10.1002/adts.202100254.
