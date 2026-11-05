# Examples

The tutorials are in `example/tutorial/`, in two levels. `basic/` covers the input itself and runs in seconds to a
minute, so read it in order. `advanced/` builds on it, so pick the topic you need.

1. [Basic](#ex_basic): the input, six topics in order.
2. [Advanced](#ex_advanced): fourteen topics in four groups, each on its own.
3. [What should I use for my case?](#WhatToUse): the method and the tutorial for a situation.

<a id="ex_basic"></a>
## 1. Basic

| Page | What it covers |
|---|---|
| [Basic 01: NV in ¹³C, from a bath to T2](../../example/tutorial/basic/01_first_run/README.md) | An NV in natural-abundance ¹³C, a cCCE2 Hahn echo in seconds on a laptop. |
| [Basic 02: Bath configurations, sampled states and isotopes](../../example/tutorial/basic/02_bath_and_states/README.md) | Averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`), and the built-in isotope table. |
| [Basic 03: Checking convergence](../../example/tutorial/basic/03_convergence/README.md) | `rbath`, `rdip`, `order`, the number of configurations and the time axis, each changed until the curve no longer moves. |
| [Basic 04: Pulse sequences](../../example/tutorial/basic/04_pulses/README.md) | The Hahn echo, CPMG, and a sequence written out (`Manual`). |
| [Basic 05: cCCE against gCCE](../../example/tutorial/basic/05_cCCE_vs_gCCE/README.md) | The two methods, the same at a strong field and different near a level crossing. |
| [Basic 06: An electron-spin bath (P1)](../../example/tutorial/basic/06_defect_bath/README.md) | The `Defect` block, sampled states and the mean field of the bath. |

<a id="ex_advanced"></a>
## 2. Advanced

<a id="ex_adv_run"></a>
### 2.1 Running and handling a calculation

| Page | What it covers |
|---|---|
| [Advanced 01: Automatic convergence](../../example/tutorial/advanced/01_auto_convergence/README.md) | The time window, rdip, rbath and the configurations found by the runner. |
| [Advanced 02: Analysis of the output](../../example/tutorial/advanced/02_analysis/README.md) | The analyzer (error correction, average, T2 fit) and the cluster reports. |
| [Advanced 03: Runs on a cluster (Slurm)](../../example/tutorial/advanced/03_server_runs/README.md) | A production run as a Slurm job, its size, one job and a job array. |
| [Advanced 04: The GPU backend](../../example/tutorial/advanced/04_GPU/README.md) | The CUDA backend against the CPU, and the automatic GPU options. |

<a id="ex_adv_systems"></a>
### 2.2 Other systems and inputs

| Page | What it covers |
|---|---|
| [Advanced 05: Electron-spin baths further](../../example/tutorial/advanced/05_P1/README.md) | Making a P1 bath, the ¹⁴N as a detuning or a spin (2e, 2e2n), order 2 checked at order 4 and a DEER pump. |
| [Advanced 06: V_B⁻ in h-BN with DFT tensors](../../example/tutorial/advanced/06_hBN/README.md) | Making an h-BN bath, the tensor files from a VASP run (`example/tools/SpinTensorReader`), and the coherence without and with them. |
| [Advanced 07: Multiple qubits](../../example/tutorial/advanced/07_multi_qubits/README.md) | Two and three NVs as one register. |
| [Advanced 08: Coordinate rotation](../../example/tutorial/advanced/08_coordinate_rotation/README.md) | A bath written in another frame (`frame.rotation`), and a `Defect` species in that frame. |

<a id="ex_adv_cce"></a>
### 2.3 Hamiltonian terms and cluster shapes

| Page | What it covers |
|---|---|
| [Advanced 09: The terms of the cluster Hamiltonian](../../example/tutorial/advanced/09_hamiltonian_terms/README.md) | The corrections to cCCE (the hyperfine-mediated term, the Knight shift) against gCCE, and the mean fields. |
| [Advanced 10: The cluster modes (dsj, pCCE, bundles) and adaptive refinement](../../example/tutorial/advanced/10_cluster_modes/README.md) | Disjoint clusters, partitions (pCCE), bundles and the adaptive refinement. |

<a id="ex_adv_methods"></a>
### 2.4 Other models of the bath

| Page | What it covers |
|---|---|
| [Advanced 11: Dissipation (meCCE and megCCE)](../../example/tutorial/advanced/11_dissipation/README.md) | Lindblad relaxation of the bath spins (meCCE) and of the qubit (megCCE) in every cluster. |
| [Advanced 12: Semiclassical noise (method noise)](../../example/tutorial/advanced/12_noise/README.md) | The Gaussian coherence of the phase noise, from the bath or from a model spectrum, and its correlation C(t) and spectrum S(f). |
| [Advanced 13: Pseudospin pairs](../../example/tutorial/advanced/13_pseudospin/README.md) | The pairs of a nuclear bath in the pseudospin model, against cCCE2. |

<a id="WhatToUse"></a>
## 3. What should I use for my case?

| Situation | Use | Example |
|---|---|---|
| A defect in a nuclear-spin bath, with the field along the defect axis. | cCCE (default). | [Basic 01: NV in ¹³C, from a bath to T2](../../example/tutorial/basic/01_first_run/README.md) |
| An electron-spin bath (P1, ...). | The `Defect` block and sampled states. | [Basic 06: An electron-spin bath (P1)](../../example/tutorial/basic/06_defect_bath/README.md) |
| Hyperfine and quadrupole from DFT. | The tensor files. | [Advanced 06: V_B⁻ in h-BN with DFT tensors](../../example/tutorial/advanced/06_hBN/README.md) |
| Near a level crossing, a low or tilted field, or the density matrix. | gCCE. | [Basic 05: cCCE against gCCE](../../example/tutorial/basic/05_cCCE_vs_gCCE/README.md) |
| Several qubits. | gCCE with `intmap`. | [Advanced 07: Multiple qubits](../../example/tutorial/advanced/07_multi_qubits/README.md) |
| Strongly coupled bath spins. | pCCE or bundles. | [Advanced 10: The cluster modes (dsj, pCCE, bundles) and adaptive refinement](../../example/tutorial/advanced/10_cluster_modes/README.md) |
| Large clusters. | The GPU backend. | [Advanced 04: The GPU backend](../../example/tutorial/advanced/04_GPU/README.md) |
| A bath written in the crystal frame. | `frame.rotation`. | [Advanced 08: Coordinate rotation](../../example/tutorial/advanced/08_coordinate_rotation/README.md) |
| A noise spectrum, or a model noise without a bath. | Method noise. | [Advanced 12: Semiclassical noise (method noise)](../../example/tutorial/advanced/12_noise/README.md) |
| Which nuclear pairs drive the decay. | Method pseudospin. | [Advanced 13: Pseudospin pairs](../../example/tutorial/advanced/13_pseudospin/README.md) |
| Relaxation of the bath spins or of the qubit. | meCCE or megCCE. | [Advanced 11: Dissipation (meCCE and megCCE)](../../example/tutorial/advanced/11_dissipation/README.md) |
