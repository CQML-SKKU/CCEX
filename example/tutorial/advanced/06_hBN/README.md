# Advanced 06: V_B⁻ in hexagonal BN

Folder: `example/tutorial/advanced/06_hBN/`

The qubit here is the boron vacancy V_B⁻ (spin 1) in hexagonal BN with natural-abundance ¹⁰B, ¹¹B and ¹⁴N. Every
nuclear spin has its DFT hyperfine tensor and, for spin > 1/2, its electric-field gradient (quadrupole).

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

`bath/flat/` holds 50 bath configurations of a flat layer, the defect position (`bath_hBN_VB_flat_Defect`) and the
tensor files the inputs read, both made with `example/tools/SpinTensorReader`:

- `hBN_VB_flat_ATensor_HSE_CellPBE_v3`, the hyperfine of the S = 1 defect in tensor format v3 (the nearest ¹⁴N has
  $A_{zz}$ 48.6 MHz).
- `hBN_VB_flat_QTensor_HSE_CellPBE`, the electric-field gradients.

The gyromagnetic ratios are in `bath/hBN_natab_gyro`, beside the `do_*` scripts that show how the tensors were
prepared. [Advanced 07: DFT hyperfine and quadrupole tensors](../07_DFT_tensors/README.md) makes such files and
compares the coherence with and without them, and the hyperfine-mediated term on this bath is in
[Advanced 10: Corrections to cCCE](../10_cCCE_corrections/README.md).

This page has two cases:

1. [The qubit density matrix of one sampled bath state](#hbn_states) (`gCCE2_Hahn_states/`).
2. [The Overhauser field of the bath on the qubit](#hbn_overhauser) (`gCCE2_Hahn_Overhauser/`).

<a id="hbn_states"></a>
## 1. One sampled bath state

Folder: `gCCE2_Hahn_states/`

<a id="a06_states_input"></a>
### 1.1 The input

gCCE2 of one sampled bath state writes the qubit density matrix over 0.5 µs in a small bath (`rbath` 5 Å). The
zero-field splitting is a `Qubit.intmap` self term (D = 3.48 GHz, E = 50 MHz, as in `gCCE2_Hahn_Overhauser`), and the
two tensor files give the hyperfine and the quadrupole:

```jsonc
"method"        : "gCCE",
"quantity"      : "dm",
"Qubit"         : {"qubit": [{"file": "../bath/flat/bath_hBN_VB_flat_Defect"}],
                   "intmap": [{"between": ["q0", "q0"], "tensor": [[-1080, 0, 0], [0, -1180, 0], [0, 0, 2320]],
                               "unit": "MHz"}]},
"configuration" : {"gyrofile": "../bath/hBN_natab_gyro", "bathfile": ["../bath/flat/bath_hBN_VB_flat_1"],
                   "s_range": [1, 1], "seed": 1},
"rbath"         : 5,
"hamiltonian"   : {"common": {"hyperfine": {"dipole": "dft", "contact": true,
                                            "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"},
                              "quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}}}
```

<a id="a06_states_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

It takes 0.6 s and 0.05 GB per rank.

<a id="a06_states_result"></a>
### 1.3 The result

![gCCE2_Hahn_states: expected result](gCCE2_Hahn_states/expected.png)

The plot shows the element between ms +1 and 0, divided by its value at t = 0. The T2 is 0.25 µs, with the
oscillation of the nearest ¹⁴N's hyperfine resolved by the 1 ns step. Over a longer window this one state diverges,
with |L| > 1 from about 0.85 µs.

Next, try `rbath` 10 and `s_range` [1, 10].

Input: [`ccein.jsonc`](gCCE2_Hahn_states/ccein.jsonc)

<a id="hbn_overhauser"></a>
## 2. The Overhauser field

Folder: `gCCE2_Hahn_Overhauser/`

<a id="a06_overhauser_input"></a>
### 2.1 The input

The Overhauser field of the bath on the qubit is the key `hamiltonian.common.qubit_bath_mean_field`, true in
`ccein_with_Overhauser.jsonc` and false in `ccein_without_Overhauser.jsonc`. Both run at zero field, between the
ms ±1 superposition and ms 0, with `deltat` 10 ps × `nstep` 50:

```jsonc
"Qubit"       : {"intmap": [{"between": ["VB", "VB"], "tensor": [[-1080, 0, 0], [0, -1180, 0], [0, 0, 2320]],
                             "unit": "MHz"}],
                 "alpha_state": {"vector": [1, 0, 1]}, "beta_state": {"vector": [0, 1, 0]}, ...},
"bfield"      : 0,
"hamiltonian" : {"common": {"qubit_bath_mean_field": true, ...}}
```

The `intmap` tensor is the ZFS with D ≈ 3.48 GHz and E ≈ 50 MHz.

<a id="a06_overhauser_run"></a>
### 2.2 Running it

```sh
NP=4 bash run.sh
```

It takes 1.7 s and 0.06 GB per rank.

<a id="a06_overhauser_result"></a>
### 2.3 The result

![gCCE2_Hahn_Overhauser: expected result](gCCE2_Hahn_Overhauser/expected.png)

Both show a shallow dip to |L| 0.966 near 0.28 ns (about 1/D) and come back to 0.998 at 0.49 ns. They differ by 9e-5
in this short window.

Next, a longer window (`nstep` 1000) shows where the two separate.

Inputs: [`ccein_with_Overhauser.jsonc`](gCCE2_Hahn_Overhauser/ccein_with_Overhauser.jsonc), [`ccein_without_Overhauser.jsonc`](gCCE2_Hahn_Overhauser/ccein_without_Overhauser.jsonc)
