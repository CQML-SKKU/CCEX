# Advanced 13: Pseudospin pairs

Folder: `example/tutorial/advanced/13_pseudospin/`

Method `pseudospin` computes the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops (Seo et al.,
Nat. Commun. 7, 12935 (2016), and Lee, Park, Seo, npj 2D Mater. Appl. 6, 60 (2022)), for any nuclear spin. Its
output, the pairs' part of the decay without the single spins' ESEEM, times cCCE order 1 is comparable with cCCE
order 2. It is meant for reading which pairs drive a decay, not for the best curve. The model, its channels and its
accuracy are in [Pseudospin model of nuclear pairs](../../../../docs/methods/pseudospin.md), and the key `channels`
in [General tags](../../../../docs/user/InputFile_General.md#tag_pseudospin_channels).

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Times
are for a desktop CPU.

This page has one case:

1. [The pairs of an h-BN bath against cCCE2](#ps_vscce) (`pseudospin2_Hahn_vsCCE/`).

<a id="ps_vscce"></a>
## 1. The pairs of an h-BN bath against cCCE2

Folder: `pseudospin2_Hahn_vsCCE/`

<a id="a13_vscce_input"></a>
### 1.1 The input

The system is the V_B⁻ of `advanced/06_hBN`, with one configuration of the flat layer, DFT hyperfine and quadrupole,
3 T, a Hahn echo, the ensemble and 459 nuclei within 20 Å. The rest is as in
`advanced/06_hBN/cCCE2_Hahn_tensors` (`ccein_hf3_qd2.jsonc`). Four inputs differ in the method:

```jsonc
"order"  : 2                                                                       // ccein_cce2.jsonc, cCCE2
"order"  : 1                                                                       // ccein_cce1.jsonc, single nuclei
"method" : "pseudospin", "hamiltonian" : {"pseudospin": {"channels": "grouped"}}   // ccein_grouped.jsonc
"method" : "pseudospin", "hamiltonian" : {"pseudospin": {"channels": "pairwise"}}  // ccein_pairwise.jsonc
```

For `pseudospin`, `order` is 2 by default and the only one accepted. The pseudospin inputs keep the
`hamiltonian.common` block of the others beside `hamiltonian.pseudospin`.

<a id="a13_vscce_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

It runs the four inputs, multiplies each pseudospin curve by cCCE1 and plots it with cCCE2. It also plots the
pseudospin curve alone, the pairs' part of the decay. The four runs take 12 s and 0.06 GB per rank. Each pseudospin
run also writes `results/ccein_<name>_pseudospin_channels` (about 100 MB here), one line per channel with the two
nuclei, their |m⟩ labels, Ω (the splitting 2|⟨m₁,m₂|H|m₁′,m₂′⟩|) and Δ (the level gap) in kHz, in the qubit's α and β
states.

<a id="a13_vscce_result"></a>
### 1.3 The result

![pseudospin2_Hahn_vsCCE: expected result](pseudospin2_Hahn_vsCCE/expected.png)

T2 is 58.0 µs by cCCE2 and 58.1 µs by pseudospin × cCCE1, with both `channels` within 2.3e-3 of cCCE2. The pairs
alone (58.2 µs) are within 3.9e-2 of cCCE2, so at 3 T the decay is the pairs'. Here `grouped` and `pairwise` agree,
because the quadrupole splits the levels of a chain. Without a quadrupole, `pairwise` misses the coupling inside a
chain and `grouped` is exact within the cut Hamiltonian
([General tags](../../../../docs/user/InputFile_General.md#tag_pseudospin_channels)).

Next, try `bfield` 3000 (0.3 T). The single nuclei's ESEEM grows, and pseudospin × cCCE1 is farther from cCCE2.

Inputs: [`ccein_cce2.jsonc`](pseudospin2_Hahn_vsCCE/ccein_cce2.jsonc), [`ccein_cce1.jsonc`](pseudospin2_Hahn_vsCCE/ccein_cce1.jsonc), [`ccein_grouped.jsonc`](pseudospin2_Hahn_vsCCE/ccein_grouped.jsonc), [`ccein_pairwise.jsonc`](pseudospin2_Hahn_vsCCE/ccein_pairwise.jsonc)
