# 04_hBN

The boron vacancy V_B⁻ (spin 1) in hexagonal BN with natural-abundance ¹⁰B, ¹¹B and ¹⁴N, every nuclear spin with
its DFT hyperfine tensor and, for spin > 1/2, its electric-field gradient (quadrupole). Each folder:
`cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every key not written in
an input has its default (`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

`bath/flat/` holds 50 bath configurations of a flat layer, the defect position (`bath_hBN_VB_flat_Defect`) and the
tensor files the inputs read: `hBN_VB_flat_ATensor_HSE_CellPBE_v3` (hyperfine of the S = 1 defect, tensor format v3;
nearest ¹⁴N $A_{zz}$ 48.6 MHz) and `hBN_VB_flat_QTensor_HSE_CellPBE` (EFG), made with `tools/SpinTensorReader`; the
gyromagnetic ratios are `bath/hBN_natab_gyro`, beside the `do_*` scripts that show how the tensors were prepared.

## cCCE2_Hahn_ensemble

- **What it shows:** what the DFT tensors change: cCCE2 Hahn echo at 3 T over five configurations (`k_range`), with
  the hyperfine and quadrupole tensors (`ccein_hf3_qd2.jsonc`), the quadrupole only (the hyperfine a point dipole,
  `hf0_qd2`) and the hyperfine only (`hf3_qd0`).
- **Keys that matter:** `hamiltonian.common.hyperfine` {`dipole` dft, `contact` true, `file`},
  `hamiltonian.common.quadrupole` {`structure` plain, `file`}, `bfield` 30000 G, `rbath` 20 Å, `rdip` 10 Å, `deltat`
  3 µs × `nstep` 50, `k_range` [1, 5].
- **How to run:** `NP=4 bash run.sh`: 8 s, 0.07 GB per rank.
- **Expected result:** `expected.png`; T2 56.5 µs with both tensors, 56.6 µs quadrupole only, 34.0 µs hyperfine only.
- **What to change next:** `k_range` [1, 50]; `rbath` 30 to check convergence.

![cCCE2_Hahn_ensemble: expected result](cCCE2_Hahn_ensemble/expected.png)

Inputs: [`ccein_hf0_qd2.jsonc`](cCCE2_Hahn_ensemble/ccein_hf0_qd2.jsonc), [`ccein_hf3_qd0.jsonc`](cCCE2_Hahn_ensemble/ccein_hf3_qd0.jsonc), [`ccein_hf3_qd2.jsonc`](cCCE2_Hahn_ensemble/ccein_hf3_qd2.jsonc)

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

## gCCE2_Hahn_states

- **What it shows:** gCCE2 of one sampled bath state, the qubit density matrix (`quantity` dm) over 0.5 µs in a small
  bath (`rbath` 5 Å); the plot shows the element between ms +1 and 0, divided by its value at t = 0.
- **Keys that matter:** `method` gCCE, `quantity` dm, `s_range` [1, 1], `seed` 1, the two tensor files, the
  zero-field splitting as a `Qubit.intmap` self term (D = 3.48 GHz, E = 50 MHz, as in `gCCE2_Hahn_Overhauser`).
- **How to run:** `NP=4 bash run.sh`: 0.6 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 0.25 µs, with the oscillation of the nearest ¹⁴N's hyperfine resolved by
  the 1 ns step (before 2.1.0 the input had no zero-field splitting and 10 points of 0.1 µs, and gave 0.44 µs). Over
  a longer window this one state diverges (|L| > 1 from about 0.85 µs).
- **What to change next:** `rbath` 10 and `s_range` [1, 10].

![gCCE2_Hahn_states: expected result](gCCE2_Hahn_states/expected.png)

Input: [`ccein.jsonc`](gCCE2_Hahn_states/ccein.jsonc)

## gCCE2_Hahn_Overhauser

- **What it shows:** the Overhauser field of the bath on the qubit (`hamiltonian.common.qubit_bath_mean_field`) at
  zero field, between the ms ±1 superposition (`alpha_state` vector [1, 0, 1]) and ms 0.
- **Keys that matter:** `qubit_bath_mean_field` true / false, `Qubit.alpha_state`, `beta_state`, the `intmap` ZFS
  (`"unit": "MHz"`, D ≈ 3.48 GHz, E ≈ 50 MHz), `bfield` 0, `deltat` 10 ps × `nstep` 50.
- **How to run:** `NP=4 bash run.sh`: 1.7 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; in both a shallow dip to |L| 0.966 near 0.28 ns (about 1/D) and back to 0.998
  at 0.49 ns; they differ by 9e-5 in this short window.
- **What to change next:** a longer window (`nstep` 1000) to see where the two separate.

![gCCE2_Hahn_Overhauser: expected result](gCCE2_Hahn_Overhauser/expected.png)

Inputs: [`ccein_with_Overhauser.jsonc`](gCCE2_Hahn_Overhauser/ccein_with_Overhauser.jsonc), [`ccein_without_Overhauser.jsonc`](gCCE2_Hahn_Overhauser/ccein_without_Overhauser.jsonc)
