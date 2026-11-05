# The pseudospin model of nuclear pairs (method pseudospin)

`"method": "pseudospin"` computes the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops, for
any nuclear spin. It is meant for reading which pairs drive a decay, and it does not give the most accurate curve.
The model for two spins 1/2 is that of Seo et al. (2016), which Lee, Park and Seo (2022) apply to h-BN. The extension
to larger spins and several channels described below is CCEX's own.

The key is [channels](../user/InputFile_General.md#tag_pseudospin_channels), the tutorial is
`example/tutorial/advanced/13_pseudospin/`, and `test/code_verification/check_pseudospin/` holds the exact cases and
the refusals.

1. [The cut Hamiltonian](#ps_cut): which elements of a pair's Hamiltonian are kept.
2. [The pair's term](#ps_term): the channels grouped or pairwise, and the output.
3. [How accurate it is](#ps_accuracy): when the model agrees with cCCE order 2.
4. [The channel table](#ps_table): the file of every channel's Ω and Δ.
5. [Where it runs](#ps_support): what is supported and refused.

<a id="ps_cut"></a>
## 1. The cut Hamiltonian

For each pair and each qubit state (α, β), CCEX writes the pair's Hamiltonian in the product basis of the two spins'
own eigenstates in that state (Zeeman, hyperfine and the whole quadrupole). It keeps the diagonal and every nonzero
element between two product states that differ in both spins. Such a pair of states is a channel:

- A flip-flop, $|m_1, m_2\rangle \leftrightarrow |m_1-1, m_2+1\rangle$.
- A flip-flip, $|m_1, m_2\rangle \leftrightarrow |m_1+1, m_2+1\rangle$.
- With a quadrupole that mixes the $|m\rangle$, also a double-quantum channel such as
  $|m_1, m_2\rangle \leftrightarrow |m_1-2, m_2+2\rangle$.

The elements that change one spin only are dropped.

The α and β product states of one index are taken to be the same state. This leaves out their overlap, which is the
pair's share of the ESEEM.

<a id="ps_term"></a>
## 2. The pair's term

With the nuclear state $I/d$, $d = (2I_1+1)(2I_2+1)$, each product state has weight $1/d$. Starting from the diagonal
coherence $L_\text{diag}$, the channels replace the diagonal values of their states:

- `"grouped"` (the default). The channels join the product states into blocks $B$. A state that sits in two channels
  joins them into a chain, for example $|1,-1\rangle \leftrightarrow |0,0\rangle \leftrightarrow |-1,1\rangle$ for two
  spins 1. Each block of $k_B$ states evolves as a whole,
  $$L_\text{pair} = L_\text{diag} + \sum_B \frac{k_B}{d}\,\bigl(l_B - l_{B,\text{diag}}\bigr).$$
  Within the cut Hamiltonian this is exact.
- `"pairwise"`. Each channel is a two-level system with its own $\Omega$ and $\Delta$, and
  $$L_\text{pair} = L_\text{diag} + \frac{2}{d} \sum_c \bigl(l_c - l_{c,\text{diag}}\bigr).$$
  This misses the coupling inside a chain.

The two agree when no two channels share a state.

The cluster's value is $P = L_\text{pair}/(L_1 L_2)$. The single spins $L_1$ and $L_2$ come from their eigenvalues
alone (1 under a Hahn echo or CPMG), and a single spin's cluster is 1. The output $\prod_\text{pairs} P$ is therefore
the pairs' part of the decay, without the single spins' ESEEM. Multiplied by cCCE order 1 it can be compared with
cCCE order 2.

For two spins 1/2 with the secular hyperfine and the flip-flop alone, both choices give Seo's eq. 23,
$P = 1 - K \sin^2(\sqrt{\Delta A^2 + D^2}\, t/4) \sin^2(D t/4)$ with $K = \Delta A^2/(\Delta A^2 + D^2)$.

<a id="ps_accuracy"></a>
## 3. How accurate it is

The comparison is P times cCCE order 1 against cCCE order 2.

- When the cut drops nothing (two spins on the qubit axis, no quadrupole), `"grouped"` agrees with cCCE order 2 and
  `"pairwise"` does not, since it misses the chains.
- A quadrupole that mixes the $|m\rangle$ can put both choices far from cCCE order 2.
- On the h-BN case the two choices agree with each other, and both move further from cCCE order 2 as the field is
  lowered.

<a id="ps_table"></a>
## 4. The channel table

Each run writes `{output.file}_pseudospin_channels`, one line per channel, with:

- The two spins and their positions.
- The $|m\rangle$ labels of the two product states (by the largest weight in α).
- $\Omega_\alpha$ and $\Omega_\beta$, the splitting $2|\langle m_1,m_2|H|m_1',m_2'\rangle|$, in kHz.
- $\Delta_\alpha$ and $\Delta_\beta$, the level gap, in kHz.

The file can be large, about 100 MB for the 459 nuclei of the h-BN tutorial.

<a id="ps_support"></a>
## 5. Where it runs

- One qubit, `quantity` coherence and the ensemble (no `s_range`).
- `order` 2, its default. Other orders are refused.
- `cluster.mode` standard.
- The field along the qubit axis.
- Ideal π pulses from `pulse.count`, with `pulse.per_qubit` false.
- `execution.propagator` eigen, on the CPU.

It refuses `hamiltonian.cCCE` (knight, hfmedi, delta_A) and Defect species with on-site spins.

<a id="ps_references"></a>
## References

- H. Seo, A. L. Falk, P. V. Klimov, K. C. Miao, G. Galli and D. D. Awschalom, Quantum decoherence dynamics of
  divacancy spins in silicon carbide, Nat. Commun. 7, 12935 (2016), doi:10.1038/ncomms12935, Supplementary eqs. 20–24
  (the pseudospin model of two spins 1/2).
- J. Lee, H. Park and H. Seo, First-principles theory of extending the spin qubit coherence time in hexagonal boron
  nitride, npj 2D Mater. Appl. 6, 60 (2022), doi:10.1038/s41699-022-00336-2 (the model applied to h-BN).
