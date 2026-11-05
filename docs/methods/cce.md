# Cluster-correlation expansion: cCCE and gCCE

This page describes what CCEX computes with `method` cCCE and gCCE and the standard clustering (`cluster.mode`
standard, dsj, itb, dsjitb). Partitions and bundles are in [pCCE](pcce.md) and
[Bundle clustering and adaptive refinement](bundle-adaptive.md), the master-equation forms in
[Master-equation CCE](master-equation.md), and the keys in the [input reference](../user/InputFile.md).

1. [The expansion](#cce_expansion): the product over clusters, wiDiv and noDiv, and a diverging product.
2. [Ensemble and sampled bath states](#cce_states): the two starting states of the bath.
3. [cCCE and gCCE](#cce_ccce): the two ways of putting the qubit into a cluster.
4. [The Hamiltonian](#cce_hamiltonian): its terms and units, the mean fields and the corrections of cCCE.
5. [The defect bath](#cce_defect): P1 and other electron spins.
6. [dsj, itb and dsjitb](#cce_dsj): the variants of the standard clustering.

<a id="cce_expansion"></a>
## 1. The expansion

The qubit starts in a superposition of two of its levels, $|\alpha\rangle$ and $|\beta\rangle$ (`Qubit.alpha_state`,
`beta_state`). The coherence function $L(t)$ follows the off-diagonal element $\langle\alpha|\rho_q(t)|\beta\rangle$ of
the qubit's reduced density matrix $\rho_q$ under the pulse sequence. $|L| = 1$ when the bath does not act.

The bath is too large to propagate as a whole. CCE therefore writes $L$ as a product of irreducible contributions over
clusters $C$ of bath spins:

$$L(t) = \prod_C \tilde L_C(t), \qquad \tilde L_C = \frac{L_C}{\prod_{C' \subset C} \tilde L_{C'}}.$$

Here $L_C$ is the coherence computed exactly with the spins of $C$ alone, by building and propagating the cluster's
Hamiltonian. The product in the denominator runs over the proper sub-clusters of $C$. Truncating at clusters of at most
`order` spins gives CCE-`order`.

Written out, the product is $L = \prod_C L_C^{p_C}$, where each power $p_C$ is an integer: +1 for the largest clusters,
and negative or positive for their sub-clusters. The run computes every $L_C$ once, raises it to its power and
multiplies.

The bath is the spins within [rbath](../user/InputFile_General.md#tag_rbath) of the qubit. A cluster is a set of up to
`order` of them that is connected in a graph whose edges are the spin pairs with a distance between the inner and the
outer [rdip](../user/InputFile_General.md#tag_rdip). Clusters of order n are grown from those of order n − 1 by adding
a spin connected to any member. For order 3 and up, [cluster.selection.max_per_order](../user/InputFile_Cluster.md#tag_max_clusters_per_order) keeps
  the strongest clusters of each order, and [add_subclusters](../user/InputFile_Cluster.md#tag_add_subclusters) lists
  their sub-clusters. The strength of a cluster is the sum of $|J_{zz}|$ over its spin pairs (below).

A result has converged when raising `order`, `rdip` and `rbath` no longer changes the curve.
`example/tutorial/advanced/05_P1/cCCE4_Hahn_states/` compares order 2 with order 4, and
[auto_convergence](../user/InputFile_AutoConvergence.md) finds `rdip`, `rbath` and the time window by itself.

<a id="cce_wiDiv"></a>
<a id="cce_divergence"></a>
### 1.1 wiDiv, noDiv and a diverging product

A run writes two curves per bath configuration and state (`output.savemode` `nodiv`, on by default):

- `_wiDiv` ("with division") is $L(t)$, the product above with every power. This is the coherence.
- `_noDiv` ("no division") is the product of the clusters with positive power only. A small denominator cannot make it
  unbounded, but it is not a coherence.

A coherence cannot exceed 1 in magnitude. When `_wiDiv` does, or turns into `nan`/`inf`, the truncated expansion has
not converged at those times. A sub-cluster coherence close to zero sits in a denominator. The analyzer's error
correction (`example/analysis/Analyzer`, [Analysis](../user/Analysis.md)) compares the two curves to find these points.

It happens with strongly coupled spins, for example close pairs of electron spins, a spin-1 qubit without its
zero-field splitting in gCCE at zero field, or cCCE in a tilted field without
[knight](../user/InputFile_General.md#tag_knight) (tutorial `advanced/10_cCCE_corrections/cCCE2_Hahn_knightShift`).
What helps:

- Adding the missing physics, the zero-field splitting or the Knight term.
- A higher `order`.
- Keeping the strongly coupled spins in one cluster ([pCCE](pcce.md), [bundles](bundle-adaptive.md)).
- Another [mean-field](#cce_meanfield) choice.

<a id="cce_states"></a>
## 2. Ensemble and sampled bath states

Without [configuration.s_range](../user/InputFile_General.md#tag_s_range) the bath starts unpolarised. The density
matrix of each cluster's bath part is the identity divided by its dimension (the ensemble).

With `s_range` the bath is taken in product states, one per state s of `s_range`:

- Each spin is in one $I_z$ sublevel, drawn uniformly from its $2I+1$ sublevels with
  [seed](../user/InputFile_General.md#tag_seed) or read from
  [bath_state_file](../user/InputFile_General.md#tag_bath_state_file).
- Each state is expanded on its own and written as `_state<s>`.
- The average over states is left to the analysis.

Sampled states are the usual choice for electron-spin baths, where one bath state differs much from another. The
[mean fields](#cce_meanfield) need them.

<a id="cce_ccce"></a>
## 3. cCCE and gCCE

**cCCE** (conventional CCE, the default `method`). The qubit is not in the cluster. It is projected onto its two
states, and each spin of the cluster sees the field the qubit sets in that state (secular approximation, pure
dephasing). For branch $k \in \{\alpha, \beta\}$, with $m_k$ the $S_z$ sublevel of the branch state, the cluster
Hamiltonian is

$$H_C^{(k)} = H_{\rm bath} + m_k \sum_{n \in C} \mathbf{A}_{z,n} \cdot \mathbf{I}_n \;+\; \text{(corrections below)},$$

with $\mathbf{A}_{z,n}$ the z row of the hyperfine tensor of spin $n$. Then
$L_C(t) = \mathrm{Tr}[U^{(\alpha)}(t)\,\rho_C\,U^{(\beta)\dagger}(t)]$, with the two branch propagators under the pulse
sequence. Each π pulse exchanges the branches. Qubit-bath flip-flops are left out. Every pulse is a π flip of the
qubit, and other pulses are refused with cCCE.

**gCCE** (generalized CCE). The qubit is inside every cluster. The cluster Hamiltonian is the full
$H_q \otimes 1 + 1 \otimes H_{\rm bath} + \sum_n \mathbf{S}\cdot\mathsf{A}_n\cdot\mathbf{I}_n$. The qubit and bath
density matrix is propagated with the pulses as written (any axis and angle), and $L_C$ is the $\alpha\beta$ element of
the cluster's reduced qubit density matrix. The qubit alone is the cluster of order 0, $L_0$, and every cluster enters
normalised by it: $L = L_0 \prod_C (L_C / L_0)^{p_C}$.

gCCE keeps what cCCE drops:

- The qubit's zero-field splitting.
- Its mixing in a tilted or low field.
- Several qubits (`Qubit.qubit`, `intmap`).
- The qubit density matrix (`quantity` dm).

It costs a factor of the qubit dimension in every cluster. `example/tutorial/basic/05_cCCE_vs_gCCE/` compares the two.

<a id="cce_hamiltonian"></a>
## 4. The Hamiltonian

Units:

- Energies are angular frequencies in rad/ms (= rad·kHz), times in ms, positions in Å and fields in G.
- Gyromagnetic ratios $\gamma$ are in rad/ms/G, and the Zeeman term is $-\gamma\,\mathbf{B}\cdot\mathbf{I}$, so the
  electron's $\gamma$ is negative ($-17608.597050$).
- Values read in MHz (the tensor files, the `Defect` block) are multiplied by $2\pi\cdot 10^3$.

The fixed unit of every input number is in the [units table](../user/InputFile.md#Units).

| Term | Form | Where it comes from |
|---|---|---|
| Bath Zeeman | $-\gamma_n \mathbf{B}\cdot\mathbf{I}_n$ | `bfield`, the gyro file or the built-in isotope table. |
| Qubit Zeeman, ZFS, qubit-qubit | $-\gamma_q \mathbf{B}\cdot\mathbf{S} + \mathbf{S}\cdot\mathsf{D}\cdot\mathbf{S} + \mathbf{S}_i\cdot\mathsf{J}_{ij}\cdot\mathbf{S}_j$ | The `Qubit` block (`zfs`, `intmap`). It is in the cluster only with gCCE. |
| Qubit-bath hyperfine | $\mathbf{S}\cdot\mathsf{A}_n\cdot\mathbf{I}_n$ | The point dipole from the positions, or DFT dipolar and contact parts ([hyperfine](../user/InputFile_General.md#tag_hyperfine)). |
| Bath-bath | $\mathbf{I}_i\cdot\mathsf{J}_{ij}\cdot\mathbf{I}_j$ | The point dipole between every pair of the cluster. |
| Quadrupole | $\mathbf{I}_n\cdot\mathsf{Q}_n\cdot\mathbf{I}_n$ for $I > 1/2$ | The [quadrupole](../user/InputFile_General.md#tag_quadrupole) files, or `Defect` `efg` and `quadrupole_moments`. |
| Defect on-site | The detuning $\delta\,I_z$, the zero-field splitting, and the on-site nuclear spins with their hyperfine and quadrupole tensors | The `Defect` block ([below](#cce_defect)). |
| Mean fields | $\sum_j J^{zz}_{nj} m_j I_{z,n}$ and $\sum_j A_{zz,j} m_j S_z$ | [Below](#cce_meanfield). |

The point-dipole tensor between two spins at distance $r$ along $\hat n$ is, as an angular frequency,

$$\mathsf{J} = \frac{\mu_0}{4\pi}\,\frac{\hbar\,\gamma_1\gamma_2}{r^3}\,\left(1 - 3\,\hat n\,\hat n^{T}\right).$$

The strength used to rank clusters and to pair spins in dsj is $|J_{zz}|$.

<a id="cce_meanfield"></a>
### 4.1 Mean fields (sampled states only)

With `s_range` a cluster is computed in the static field of the spins outside it. This is a mean-field approximation of
the rest of the bath in its sampled state $m_j$.

- [bath_bath_mean_field](../user/InputFile_General.md#tag_bath_bath_mean_field) acts on each bath spin of the cluster.
  `"z"` (default) adds $\sum_{j \notin C} J^{zz}_{nj} m_j I_{z,n}$, `"xyz"` adds the whole vector
  $\mathsf{J}_{nj}(0,0,m_j)$, and `"off"` adds nothing. The cluster's own members act through their full coupling
  instead.
- [qubit_bath_mean_field](../user/InputFile_General.md#tag_qubit_bath_mean_field) acts on the qubit (gCCE, megCCE). It
  adds the Overhauser field $\sum_{j \notin C} A_{zz,j} m_j S_z$ of the spins outside the cluster. It is off by
  default. cCCE has no qubit Hamiltonian, so it has no effect there.

In the ensemble neither is used, since the unpolarised bath has no static field.

<a id="cce_corrections"></a>
### 4.2 What cCCE leaves out, and the three corrections

The projection onto the two branches keeps only the $S_z$ part of the hyperfine coupling. The transverse part
$S_x \mathbf{A}_x\cdot\mathbf{I} + S_y \mathbf{A}_y\cdot\mathbf{I}$ changes $m_s$ and acts at second order through
virtual qubit transitions. For a spin-1 qubit with zero-field splitting $\Delta = \tfrac32 D_{zz}$, three pieces of it
can be added to cCCE (`hamiltonian.cCCE`):

| Key | Adds to the branch Hamiltonian | Matters for |
|---|---|---|
| [knight](../user/InputFile_General.md#tag_knight) | The Knight shift of each bath spin by the field across the qubit axis, $-(3\lvert m_s\rvert - 2)(\gamma_q/\Delta)\,\mathbf{B}\cdot\tilde{\mathsf{A}}\cdot\mathbf{I}$, with $\tilde{\mathsf{A}}$ the hyperfine tensor without its z row. | A tilted field. |
| [delta_A](../user/InputFile_General.md#tag_delta_A) | The same-nucleus second-order term $\delta\mathbf{A}_n(m_s)\cdot\mathbf{I}_n$, $\delta\mathbf{A}_n = -m_s(\mathbf{A}_x\times\mathbf{A}_y)/(2\Delta)$. | Nuclei close to the defect at low field. |
| [hfmedi](../user/InputFile_General.md#tag_hfmedi) | The complete second order in the hyperfine coupling, $\sum_{\beta\neq\alpha} E_{\alpha\beta}E_{\alpha\beta}^\dagger/(\omega_\alpha-\omega_\beta)$, which holds the bath-bath coupling mediated by the qubit and the same-nucleus terms. | Close nuclei at low field. |

- knight and delta_A divide by $\Delta$ and are refused without a zero-field splitting.
- delta_A is refused with hfmedi, which already contains it.
- The exact reference for all three is gCCE, which keeps the transverse coupling in the cluster Hamiltonian.

<a id="cce_delta_A"></a>
delta_A is the part of the second-order transverse hyperfine term that is linear in $\mathbf{I}_n$. Its coefficient
$-m_s$ agrees with the coefficient $(2 - 3|m_s|)$ of Maze et al., Eq. (A4), only for $m_s = +1$. It holds with the field
along the qubit axis and $\gamma_e B \ll \Delta$. The qubit Zeeman energy is left out of the denominators, as in knight,
so its relative error is about $\gamma_e B_z/\Delta$ (0.5 % at 5 G, 10 % at 100 G).

<a id="cce_defect"></a>
## 5. The defect bath (P1 and other electron spins)

A bath spin whose name is a `Defect` species is a paramagnetic defect. It is an electron spin with indexed *axis
states*. Each axis state stands for one environment of the electron, such as a Jahn-Teller axis of a P1 centre combined
with a frozen state of its nucleus.

Each defect gets one axis state and the tensors of that state. The axis state comes from
[defect_axis_file](../user/InputFile_General.md#tag_defect_axis_file) or is drawn at random from `seed`. Its nucleus
enters in one of two ways:

- As a static shift of the electron, the [detuning](../user/InputFile_Defect.md#tag_df_detuning) $\delta = A\,m_I$ of
  the axis state, with `n_onsite_spins` 0. A P1 with ¹⁴N has twelve detunings, four axes times three $m_I$.
- As an on-site spin of its own (`n_onsite_spins` > 0), with its hyperfine and electric-field-gradient tensors per axis
  state. It then joins its electron in every cluster that contains the defect. The defect is one site whose dimension
  is the product of its spins'.

Electron baths are usually run with sampled states (`s_range`):

- The electrons' states come from `seed` or `bath_state_file`.
- The on-site spins' states come from `seed` or `defect_state_file`.
- The other bath spins act through the [mean field](#cce_meanfield).

`example/tutorial/basic/06_defect_bath/` runs the ¹⁴N as a detuning, `advanced/05_P1/cCCE2_Hahn_onsite14N/` as an
on-site spin, and [Defect tags](../user/InputFile_Defect.md)
describes the keys.

<a id="cce_dsj"></a>
## 6. dsj, itb and dsjitb (cCCE)

These are three variants of the standard clustering, for cCCE only. gCCE, meCCE and megCCE refuse them.

All three start from the same graph of spin pairs within `rdip` and split the bath into disjoint *groups* of at most
`order` spins. The pairs are taken in order of decreasing $|J_{zz}|$. Two spins (or their groups) are joined when the
joined group still has at most `order` spins. A connected pair whose spins end in different groups, because joining
them would exceed `order`, is an *inter-bath* pair.

| `cluster.mode` | Computes |
|---|---|
| `dsj` | $\prod_G L_G$. Every group is computed exactly as one cluster, and no other cluster is. |
| `itb` | $\prod_{(a,b)} \tilde L_{ab} = \prod_{(a,b)} L_{ab}/(L_a L_b)$ over the inter-bath pairs, the factor that the pairs between groups add. |
| `dsjitb` | Both, $\prod_G L_G \prod_{(a,b)} \tilde L_{ab}$. |

- At `order` 2 the groups are pairs and single spins, every other connected pair is an inter-bath pair, and `dsjitb`
  equals standard CCE2.
- At higher order a group is computed whole, but between groups only pairs enter, so `dsjitb` is not CCE-`order`.
- dsj keeps the strongest couplings inside clusters and needs far fewer clusters than standard CCE. Compare it with
  `standard` at the same `order` and `rdip`, as `example/tutorial/advanced/11_cluster_shapes/cCCE2_Hahn_dsj/` does.

<a id="cce_references"></a>
## References

- W. Yang and R.-B. Liu, Quantum many-body theory of qubit decoherence in a finite-size spin bath, Phys. Rev. B 78,
  085315 (2008), doi:10.1103/PhysRevB.78.085315 (the CCE).
- W. Yang and R.-B. Liu, Quantum many-body theory of qubit decoherence in a finite-size spin bath. II. Ensemble
  dynamics, Phys. Rev. B 79, 115320 (2009), doi:10.1103/PhysRevB.79.115320 (the ensemble CCE).
- M. Onizhuk et al., Probing the coherence of solid-state qubits at avoided crossings, PRX Quantum 2, 010311 (2021),
  doi:10.1103/PRXQuantum.2.010311, and M. Onizhuk and G. Galli, PyCCE: A Python package for cluster correlation
  expansion simulations of spin qubit dynamics, Adv. Theory Simul. 4, 2100254 (2021), doi:10.1002/adts.202100254
  (the generalized CCE).
- J. R. Maze, J. M. Taylor and M. D. Lukin, Electron spin decoherence of single nitrogen-vacancy defects in diamond,
  Phys. Rev. B 78, 094303 (2008), doi:10.1103/PhysRevB.78.094303 (the second-order hyperfine terms).
