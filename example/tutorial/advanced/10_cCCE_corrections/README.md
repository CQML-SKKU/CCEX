# Advanced 10: Corrections to cCCE

Folder: `example/tutorial/advanced/10_cCCE_corrections/`

cCCE projects the qubit onto its two states and keeps only the $S_z$ part of the hyperfine coupling. The keys of
`hamiltonian.cCCE` add back second-order pieces of the rest, for a spin-1 qubit with a zero-field splitting: the
hyperfine-mediated coupling of two nuclei through the defect (`hfmedi`, on the h-BN bath of
[Advanced 06: V_B⁻ in hexagonal BN](../06_hBN/README.md)) and the Knight shift in a tilted field (`knight`). How they
are made is in [cCCE and gCCE](../../../../docs/methods/cce.md#cce_corrections), and gCCE is their exact reference
([Basic 05: cCCE against gCCE](../../basic/05_cCCE_vs_gCCE/README.md)).
Each folder: `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every key
not written in an input has its default (`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

## cCCE2_Hahn_hfMediated

- **What it shows:** the hyperfine-mediated bath-bath coupling (`hamiltonian.cCCE.hfmedi`): at 1000 G with the
  V_B⁻ zero-field splitting, two nuclei coupled through the defect's transverse hyperfine; one configuration, ensemble.
- **Keys that matter:** `hamiltonian.cCCE.hfmedi` true / absent, the ZFS as a `Qubit.intmap` self entry (kHz), `bfield` 1000, `rbath` 14 Å,
  `rdip` 8 Å, `deltat` 5 ns × `nstep` 50.
- **How to run:** `NP=4 bash run.sh`: 1.2 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 86 ns with the term; without it |L| is still 0.66 at 245 ns.
- **What to change next:** a higher `bfield` (the term falls off with the field).

![cCCE2_Hahn_hfMediated: expected result](cCCE2_Hahn_hfMediated/expected.png)

Inputs: [`ccein_wi_hfmedi.jsonc`](cCCE2_Hahn_hfMediated/ccein_wi_hfmedi.jsonc), [`ccein_wo_hfmedi.jsonc`](cCCE2_Hahn_hfMediated/ccein_wo_hfmedi.jsonc)


## cCCE2_Hahn_knightShift

- **What it shows:** in a field tilted from the NV axis ([50, 0, 100] G) the qubit's transverse field shifts the ¹³C
  precession; `hamiltonian.cCCE.knight` adds that shift to cCCE (`ccein_wi_knight.jsonc`) against without
  (`ccein_wo_knight.jsonc`).
- **Keys that matter:** `hamiltonian.cCCE.knight` true, `bfield` [50, 0, 100], the ZFS (`Qubit.intmap` self entry).
- **How to run:** `NP=4 bash run.sh`: 1.2 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.082 ms with the Knight shift; without it the divided curve (`_wiDiv`)
  diverges (|L| up to 16, cut at 1.2 in the plot) - the analyzer's error correction
  (`example/analysis/Analyzer`, `-err`) repairs such curves against `_noDiv`.

![cCCE2_Hahn_knightShift: expected result](cCCE2_Hahn_knightShift/expected.png)

Inputs: [`ccein_wi_knight.jsonc`](cCCE2_Hahn_knightShift/ccein_wi_knight.jsonc), [`ccein_wo_knight.jsonc`](cCCE2_Hahn_knightShift/ccein_wo_knight.jsonc)
