# The pseudospin model of nuclear pairs (method pseudospin)

`"method": "pseudospin"` computes the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops, for
any nuclear spin. It is meant for reading which pairs drive a decay; it does not give the most accurate curve.

The model for two spins 1/2 is H. Seo et al., Nat. Commun. 7, 12935 (2016), Supplementary eqs. 20–24; J. Lee, H. Park and H. Seo, npj 2D
Mater. Appl. 6, 60 (2022) apply it to h-BN. The extension to larger spins and several channels described below is
CCEX's own.

The key is [channels](../user/InputFile_General.md#tag_pseudospin_channels); the tutorial is
`example/tutorial/advanced/10_pseudospin/`.

## The cut Hamiltonian

For each pair and each qubit state (α, β), CCEX writes the pair's Hamiltonian in the product basis of the two spins'
own eigenstates in that state (Zeeman, hyperfine and the whole quadrupole). It keeps the diagonal and every nonzero
element between two product states that differ in both spins. Such a pair of states is a channel:

- a flip-flop, $|m_1, m_2\rangle \leftrightarrow |m_1-1, m_2+1\rangle$;
- a flip-flip, $|m_1, m_2\rangle \leftrightarrow |m_1+1, m_2+1\rangle$;
- with a quadrupole that mixes the $|m\rangle$, also a double-quantum channel such as $|m_1, m_2\rangle \leftrightarrow |m_1-2, m_2+2\rangle$.

The elements that change one spin only are dropped.

The α and β product states of one index are taken to be the same state. This leaves out their overlap, which is the
pair's share of the ESEEM.

## The pair's term

With the nuclear state $I/d$, $d = (2I_1+1)(2I_2+1)$, each product state has weight $1/d$. Starting from the diagonal
coherence $L_\text{diag}$, the channels replace the diagonal values of their states:

- `"grouped"` (the default): the channels join the product states into blocks $B$. A state that sits in two channels
  joins them into a chain, for example $|1,-1\rangle \leftrightarrow |0,0\rangle \leftrightarrow |-1,1\rangle$ for two
  spins 1. Each block of $k_B$ states evolves as a whole:
  $$L_\text{pair} = L_\text{diag} + \sum_B \frac{k_B}{d}\,\bigl(l_B - l_{B,\text{diag}}\bigr).$$
  Within the cut Hamiltonian this is exact.
- `"pairwise"`: each channel is a two-level system with its own $\Omega$ and $\Delta$, and
  $$L_\text{pair} = L_\text{diag} + \frac{2}{d} \sum_c \bigl(l_c - l_{c,\text{diag}}\bigr).$$
  This misses the coupling inside a chain.

The two agree when no two channels share a state.

The cluster's value is $P = L_\text{pair}/(L_1 L_2)$. The single spins $L_1$, $L_2$ come from their eigenvalues alone
(1 under a Hahn echo or CPMG), and a single spin's cluster is 1. The output is therefore $\prod_\text{pairs} P$: the
pairs' part of the decay, without the single spins' ESEEM. Multiplied by cCCE order 1 it can be compared with cCCE
order 2.

For two spins 1/2 with the secular hyperfine and the flip-flop alone, both choices give Seo's eq. 23,
$P = 1 - K \sin^2(\sqrt{\Delta A^2 + D^2}\, t/4) \sin^2(D t/4)$ with $K = \Delta A^2/(\Delta A^2 + D^2)$.

## How accurate it is

The table compares P times cCCE order 1 with cCCE order 2.

| Case | `"grouped"` | `"pairwise"` |
|---|---|---|
| two spins on the qubit axis, no quadrupole (the cut drops nothing), spins 1, 3/2, 3 | within 1e-10 | 0.43, 0.71, 0.65 off |
| spins 1, quadrupole that mixes the $\lvert m\rangle$ (η = 0.5) | 0.74 off | 0.74 off |
| h-BN case of `check_numerics`, at 3, 0.3 and 0.03 T | within 2.0e-3, 2.2e-2, 4.2e-2 | within 2.0e-3, 2.2e-2, 4.2e-2 |

A quadrupole that splits a chain's levels brings the two close (spins 1, axial quadrupole: pairwise 2.9e-6 off). In
the mixing case the likely cause is that the α and β bases differ; this has not been tested separately. Apart from
the spin-3/2 case, which `check_pseudospin` holds, these comparisons were run once each.

## The channel table

Each run writes `{output.file}_pseudospin_channels`, one line per channel with:

- the two spins and their positions;
- the $|m\rangle$ labels of the two product states (by the largest weight in α);
- $\Omega_\alpha$ and $\Omega_\beta$, the splitting $2|\langle m_1,m_2|H|m_1',m_2'\rangle|$, in kHz;
- $\Delta_\alpha$ and $\Delta_\beta$, the level gap, in kHz.

The file can be large: about 100 MB for the 459 nuclei of the h-BN tutorial.

## Where it runs

- one qubit;
- the ensemble (no `s_range`);
- `order` 2 (its default; other orders are refused);
- `cluster.mode` standard;
- the field along the qubit axis;
- ideal π pulses from `pulse.count`;
- on the CPU.

## Checks

`test/code_verification/check_pseudospin/` holds these exact cases within 1e-9:

- two spins on the qubit axis (spins 1/2, and 3/2 grouped, under free induction, Hahn and CPMG-2) against
  cCCE order 2 / order 1;
- Seo's eq. 23 from the positions;
- the closed form of the channel file for spins 3/2 pairwise;
- two spins in the qubit's plane (a flip-flip);
- the double-quantum channels of an asymmetric quadrupole.

It also checks the refusals.
