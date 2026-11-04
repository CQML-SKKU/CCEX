# Cluster-correlation expansion: cCCE and gCCE

This page describes what CCEX computes with `method` cCCE and gCCE and the standard clustering (`cluster.mode`
standard, dsj, itb, dsjitb). It covers the expansion, the two ways of putting the qubit into a cluster, the
Hamiltonian, the mean fields and the corrections of cCCE.

Other pages:

- partitions and bundles: [pCCE](pcce.md) and [bundle clustering and adaptive refinement](bundle-adaptive.md);
- the master-equation forms meCCE and megCCE: [master-equation CCE](master-equation.md);
- the keys: the [input reference](../user/InputFile.md).

<a id="cce_expansion"></a>
## The expansion

The qubit starts in a superposition of two of its levels, $|\alpha\rangle$ and $|\beta\rangle$ (`Qubit.alpha_state`,
`beta_state`). The coherence function $L(t)$ follows the off-diagonal element $\langle\alpha|\rho_q(t)|\beta\rangle$ of
the qubit's reduced density matrix $\rho_q$ under the pulse sequence. $|L| = 1$ when the bath does not act.

The bath is too large to propagate as a whole. CCE therefore writes $L$ as a product of irreducible contributions over
clusters $C$ of bath spins:

$$L(t) = \prod_C \tilde L_C(t), \qquad \tilde L_C = \frac{L_C}{\prod_{C' \subset C} \tilde L_{C'}},$$

Here $L_C$ is the coherence computed exactly with the spins of $C$ alone, by building and propagating the cluster's
Hamiltonian. The product in the denominator runs over the proper sub-clusters of $C$. Truncating at clusters of at most
`order` spins gives CCE-`order`.

Written out, each computed $L_C$ enters the product with an integer power: +1 for the largest clusters, negative or
positive for their sub-clusters. The run computes the clusters and multiplies these powers.

- A cluster is a set of up to `order` bath spins that is connected in a graph. The edges of the graph are the spin
  pairs whose distance lies between the inner and the outer [rdip](../user/InputFile_General.md#tag_rdip).
- The bath is the spins within [rbath](../user/InputFile_General.md#tag_rbath) of the qubit.
- Clusters of order n are grown from those of order n − 1 by adding a spin connected to any member.
- Convergence is checked by raising `order`, `rdip` and `rbath` until the curve no longer changes.
  `example/tutorial/advanced/01_P1/cCCE4_Hahn_states/` compares order 2 with order 4.
  [auto_convergence](../user/InputFile_AutoConvergence.md) does the `rdip`, `rbath` and time-window part automatically.
- For order 3 and up, [cluster.selection.max_per_order](../user/InputFile_Cluster.md#tag_max_clusters_per_order) keeps
  the strongest clusters of each order, and [add_subclusters](../user/InputFile_Cluster.md#tag_add_subclusters) lists
  their sub-clusters. The strength of a cluster is the sum of $|J_{zz}|$ over its spin pairs (below).

<a id="cce_wiDiv"></a>
### wiDiv and noDiv

A run writes two curves per bath configuration and state (`output.savemode` `nodiv`, on by default):

- `_wiDiv` ("with division") is $L(t)$, the product above with every power. This is the coherence.
- `_noDiv` ("no division") is the product of the clusters with positive power only. A small denominator cannot make it
  unbounded, but it is not a coherence. The analyzer's error correction (`example/analysis/Analyzer`) compares the two
  curves to find the points where `_wiDiv` diverged ([Analysis](../user/Analysis.md)).

<a id="cce_divergence"></a>
### |L| > 1

A coherence cannot exceed 1 in magnitude. When `_wiDiv` does, or turns into `nan`/`inf`, the truncated expansion has
not converged at those times. A sub-cluster coherence close to zero sits in a denominator, and the division amplifies
what the larger cluster does not cancel.

This happens with strongly coupled spins:

- close pairs of electron spins;
- a spin-1 qubit without its zero-field splitting in gCCE;
- cCCE in a tilted field without [knight](../user/InputFile_General.md#tag_knight).

The tutorials `advanced/01_P1/gCCE2_Hahn_zfs` and `basic/04_cCCE_vs_gCCE/cCCE2_Hahn_knightShift` show both cases.

What helps:

- adding the missing physics (the zero-field splitting, the Knight term);
- a higher `order`;
- keeping the strongly coupled spins in one cluster ([pCCE](pcce.md), [bundles](bundle-adaptive.md));
- the [mean-field](#cce_meanfield) choice;
- the analyzer's error correction of the points where `_wiDiv` and `_noDiv` part.

The products carry a separate exponent through the accumulation and the MPI reduction. A large intermediate value
therefore does not overflow when the final one is representable. This does not make the expansion converge.

<a id="cce_states"></a>
### Ensemble and sampled bath states

Without [configuration.s_range](../user/InputFile_General.md#tag_s_range) the bath starts unpolarised: the density
matrix of each cluster's bath part is the identity divided by its dimension (the ensemble).

With `s_range` the bath is taken in product states, one per state s of `s_range`:

- each spin is in one $I_z$ sublevel, drawn uniformly from its $2I+1$ sublevels with
  [seed](../user/InputFile_General.md#tag_seed) or read from
  [bath_state_file](../user/InputFile_General.md#tag_bath_state_file);
- each state is expanded on its own and written as `_state<s>`;
- the average over states is left to the analysis.

Sampled states are the usual choice for electron-spin baths, where one bath state differs much from another. The
[mean fields](#cce_meanfield) need them.

<a id="cce_ccce"></a>
## cCCE and gCCE

**cCCE** (conventional CCE, the default `method`). The qubit is not in the cluster. It is projected onto its two
states, and each spin of the cluster sees the field the qubit sets in that state (secular approximation, pure
dephasing). For branch $k \in \{\alpha, \beta\}$, with $m_k$ the $S_z$ sublevel of the branch state, the cluster
Hamiltonian is

$$H_C^{(k)} = H_{\rm bath} + m_k \sum_{n \in C} \mathbf{A}_{z,n} \cdot \mathbf{I}_n \;+\; \text{(corrections below)},$$

with $\mathbf{A}_{z,n}$ the z row of the hyperfine tensor of spin $n$. Then
$L_C(t) = \mathrm{Tr}[U^{(\alpha)}(t)\,\rho_C\,U^{(\beta)\dagger}(t)]$, with the two branch propagators under the pulse
sequence; each π pulse exchanges the branches. Qubit-bath flip-flops are left out. Every pulse is a π flip of the
qubit; other pulses are refused with cCCE.

**gCCE** (generalized CCE). The qubit is inside every cluster. The cluster Hamiltonian is the full
$H_q \otimes 1 + 1 \otimes H_{\rm bath} + \sum_n \mathbf{S}\cdot\mathsf{A}_n\cdot\mathbf{I}_n$. The qubit and bath
density matrix is propagated with the pulses as written (any axis and angle), and $L_C$ is the $\alpha\beta$ element of
the cluster's reduced qubit density matrix. The qubit alone is the cluster of order 0, $L_0$, and every cluster enters
normalised by it: $L = L_0 \prod_C (L_C / L_0)^{p_C}$.

gCCE keeps what cCCE drops:

- the qubit's zero-field splitting;
- its mixing in a tilted or low field;
- several qubits (`Qubit.qubit`, `intmap`);
- the qubit density matrix (`quantity` dm).

It costs a factor of the qubit dimension in every cluster. `example/tutorial/basic/04_cCCE_vs_gCCE/` compares the two.

<a id="cce_hamiltonian"></a>
## The Hamiltonian

Units:

- energies are angular frequencies in rad/ms (= rad·kHz), times in ms, positions in Å, fields in G;
- gyromagnetic ratios $\gamma$ are in rad/ms/G, and the Zeeman term is $-\gamma\,\mathbf{B}\cdot\mathbf{I}$, so the
  electron's $\gamma$ is negative ($-17608.597050$);
- values read in MHz (the tensor files, the `Defect` block) are multiplied by $2\pi\cdot 10^3$.

The fixed unit of every input number is in the [units table](../user/InputFile.md#Units).

| Term | Form | Where it comes from |
|---|---|---|
| bath Zeeman | $-\gamma_n \mathbf{B}\cdot\mathbf{I}_n$ | `bfield`, the gyro file or the built-in isotope table |
| qubit Zeeman, ZFS, qubit-qubit | $-\gamma_q \mathbf{B}\cdot\mathbf{S} + \mathbf{S}\cdot\mathsf{D}\cdot\mathbf{S} + \mathbf{S}_i\cdot\mathsf{J}_{ij}\cdot\mathbf{S}_j$ | `Qubit` block (`zfs`, `intmap`); in the cluster only with gCCE |
| qubit-bath hyperfine | $\mathbf{S}\cdot\mathsf{A}_n\cdot\mathbf{I}_n$ | point dipole from the positions, or DFT dipolar and contact parts ([hyperfine](../user/InputFile_General.md#tag_hyperfine)) |
| bath-bath | $\mathbf{I}_i\cdot\mathsf{J}_{ij}\cdot\mathbf{I}_j$ | point dipole between every pair of the cluster |
| quadrupole | $\mathbf{I}_n\cdot\mathsf{Q}_n\cdot\mathbf{I}_n$ for $I > 1/2$ | [quadrupole](../user/InputFile_General.md#tag_quadrupole) files, `Defect` `efg` and `quadrupole_moments` |
| defect on-site | detuning $\delta\,I_z$, zero-field splitting, the on-site nuclear spins with their hyperfine and quadrupole tensors | `Defect` block (below) |
| mean fields | $\sum_j J^{zz}_{nj} m_j I_{z,n}$ and $\sum_j A_{zz,j} m_j S_z$ | [below](#cce_meanfield) |

The point-dipole tensor between two spins at distance $r$ along $\hat n$ is, as an angular frequency,

$$\mathsf{J} = \frac{\mu_0}{4\pi}\,\frac{\hbar\,\gamma_1\gamma_2}{r^3}\,\left(1 - 3\,\hat n\,\hat n^{T}\right).$$

The strength used to rank clusters and to pair spins in dsj is $|J_{zz}|$.

<a id="cce_meanfield"></a>
### Mean fields (sampled states only)

With `s_range` a cluster is computed in the static field of the spins outside it. This is a mean-field approximation of
the rest of the bath in its sampled state $m_j$.

- [bath_bath_mean_field](../user/InputFile_General.md#tag_bath_bath_mean_field) acts on each bath spin of the cluster.
  `"z"` (default) adds $\sum_{j \notin C} J^{zz}_{nj} m_j I_{z,n}$; `"xyz"` adds the whole vector
  $\mathsf{J}_{nj}(0,0,m_j)$; `"off"` adds nothing. The cluster's own members act through their full coupling instead.
- [qubit_bath_mean_field](../user/InputFile_General.md#tag_qubit_bath_mean_field) acts on the qubit (gCCE, megCCE). It
  adds the Overhauser field $\sum_{j \notin C} A_{zz,j} m_j S_z$ of the spins outside the cluster. It is off by
  default. cCCE has no qubit Hamiltonian, so it has no effect there.

In the ensemble neither is used, since the unpolarised bath has no static field.

<a id="cce_corrections"></a>
### What cCCE leaves out, and the three corrections

The projection onto the two branches keeps only the $S_z$ part of the hyperfine coupling. The transverse part
$S_x \mathbf{A}_x\cdot\mathbf{I} + S_y \mathbf{A}_y\cdot\mathbf{I}$ changes $m_s$ and acts at second order through
virtual qubit transitions. For a spin-1 qubit with zero-field splitting $\Delta = \tfrac32 D_{zz}$, three pieces of it
can be added to cCCE (`hamiltonian.cCCE`):

| Key | Adds to the branch Hamiltonian | Matters for |
|---|---|---|
| [knight](../user/InputFile_General.md#tag_knight) | the Knight shift of each bath spin by the field across the qubit axis, $-(3\lvert m_s\rvert - 2)(\gamma_q/\Delta)\,\mathbf{B}\cdot\tilde{\mathsf{A}}\cdot\mathbf{I}$ ($\tilde{\mathsf{A}}$: the hyperfine tensor without its z row) | a tilted field |
| [delta_A](../user/InputFile_General.md#tag_delta_A) | the same-nucleus second-order term $\delta\mathbf{A}_n(m_s)\cdot\mathbf{I}_n$, $\delta\mathbf{A}_n = -m_s(\mathbf{A}_x\times\mathbf{A}_y)/(2\Delta)$ | nuclei close to the defect at low field |
| [hfmedi](../user/InputFile_General.md#tag_hfmedi) | the complete second order in the hyperfine coupling, $\sum_{\beta\neq\alpha} E_{\alpha\beta}E_{\alpha\beta}^\dagger/(\omega_\alpha-\omega_\beta)$, which holds the bath-bath coupling mediated by the qubit and the same-nucleus terms | close nuclei, low field |

- knight and delta_A divide by $\Delta$ and are refused without a zero-field splitting.
- delta_A is refused with hfmedi, which already contains it.
- The exact reference for all three is gCCE, which keeps the transverse coupling in the cluster Hamiltonian.

<a id="cce_defect"></a>
## The defect bath (P1 and other electron spins)

A bath spin whose name is a `Defect` species is a paramagnetic defect. It is an electron spin with indexed *axis
states*. Each axis state stands for one environment of the electron, such as a Jahn-Teller axis of a P1 centre combined
with a frozen state of its nucleus.

Each defect gets one axis state and the tensors of that state. The axis state comes from
[defect_axis_file](../user/InputFile_General.md#tag_defect_axis_file) or is drawn at random from `seed`. Its nucleus
enters in one of two ways:

- as a static shift of the electron, the [detuning](../user/InputFile_Defect.md#tag_df_detuning) $\delta = A\,m_I$ of
  the axis state, with `n_onsite_spins` 0. A P1 with ¹⁴N has twelve detunings: four axes × three $m_I$.
- as an on-site spin of its own (`n_onsite_spins` > 0), with its hyperfine and electric-field-gradient tensors per axis
  state. It then joins its electron in every cluster that contains the defect. The defect is one site whose dimension
  is the product of its spins'.

Electron baths are usually run with sampled states (`s_range`):

- the electrons' states come from `seed` or `bath_state_file`;
- the on-site spins' states come from `seed` or `defect_state_file`;
- the other bath spins act through the [mean field](#cce_meanfield).

`example/tutorial/basic/05_defect_bath/` runs both forms of the ¹⁴N.

<a id="cce_dsj"></a>
## dsj, itb and dsjitb (cCCE)

These are three variants of the standard clustering, with cCCE only; gCCE, meCCE and megCCE with them are refused.

All three start from the same graph of spin pairs within `rdip` and split the bath into disjoint *groups* of at most
`order` spins. The pairs are taken in order of decreasing $|J_{zz}|$. Two spins (or their groups) are joined when the
joined group still has at most `order` spins. A connected pair whose spins end in different groups, because joining
them would exceed `order`, is an *inter-bath* pair.

| `cluster.mode` | Computes |
|---|---|
| `dsj` | $\prod_G L_G$: every group computed exactly as one cluster, no other cluster |
| `itb` | $\prod_{(a,b)} \tilde L_{ab} = \prod_{(a,b)} L_{ab}/(L_a L_b)$ over the inter-bath pairs: the factor that the pairs between groups add |
| `dsjitb` | both: $\prod_G L_G \prod_{(a,b)} \tilde L_{ab}$ |

- At `order` 2 the groups are pairs and single spins, every other connected pair is an inter-bath pair, and `dsjitb`
  equals standard CCE2.
- At higher order a group is computed whole, but between groups only pairs enter, so `dsjitb` is not CCE-`order`.
- dsj keeps the strongest couplings inside clusters with far fewer clusters than standard CCE. Compare it with
  `standard` at the same `order` and `rdip` (`example/tutorial/advanced/06_terms_and_models/cCCE2_Hahn_dsj/`).

References:

- W. Yang and R.-B. Liu, Quantum many-body theory of qubit decoherence in a finite-size spin bath, Phys. Rev. B 78,
  085315 (2008), doi:10.1103/PhysRevB.78.085315 (the CCE).
- W. Yang and R.-B. Liu, Quantum many-body theory of qubit decoherence in a finite-size spin bath. II. Ensemble
  dynamics, Phys. Rev. B 79, 115320 (2009), doi:10.1103/PhysRevB.79.115320 (the ensemble CCE).
- M. Onizhuk et al., Probing the coherence of solid-state qubits at avoided crossings, PRX Quantum 2, 010311 (2021),
  doi:10.1103/PRXQuantum.2.010311, and M. Onizhuk and G. Galli, PyCCE: A Python package for cluster correlation
  expansion simulations of spin qubit dynamics, Adv. Theory Simul. 4, 2100254 (2021), doi:10.1002/adts.202100254
  (the generalized CCE).
