# Examples

The tutorials are in `example/tutorial/`, in two levels. `basic/` covers the input itself and runs in seconds to a minute; read
it in order. `advanced/` builds on it; pick the topic you need.

**Basic**

| Page | Folder | Topic |
|---|---|---|
| [Basic 01: NV in ¹³C, from a bath to T2](../../example/tutorial/basic/01_first_run/README.md) | `basic/01_first_run` | NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds on a laptop |
| [Basic 02: Bath configurations, sampled states and isotopes](../../example/tutorial/basic/02_bath_and_states/README.md) | `basic/02_bath_and_states` | averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`); the built-in isotope table |
| [Basic 03: Checking convergence](../../example/tutorial/basic/03_convergence/README.md) | `basic/03_convergence` | `rbath`, `rdip`, `order`, the number of configurations and the time axis: change one, see whether the curve still moves |
| [Basic 04: Pulse sequences](../../example/tutorial/basic/04_pulses/README.md) | `basic/04_pulses` | Hahn echo, CPMG, a sequence written out (`Manual`) |
| [Basic 05: cCCE against gCCE](../../example/tutorial/basic/05_cCCE_vs_gCCE/README.md) | `basic/05_cCCE_vs_gCCE` | the two methods: the same at a strong field, different near a level crossing |
| [Basic 06: An electron-spin bath (P1)](../../example/tutorial/basic/06_defect_bath/README.md) | `basic/06_defect_bath` | an electron-spin bath (P1): the `Defect` block, sampled states, the mean field of the bath |

**Advanced**

| Page | Folder | Topic |
|---|---|---|
| [Advanced 01: The P1 bath further](../../example/tutorial/advanced/01_P1/README.md) | `advanced/01_P1` | the P1 bath further: the ¹⁴N as a spin, order 4, gCCE at zero field, the density matrix, Vp, a DEER pump |
| [Advanced 02: V_B⁻ in hexagonal BN](../../example/tutorial/advanced/02_hBN/README.md) | `advanced/02_hBN` | V_B⁻ in h-BN with DFT tensors: sampled states, the Overhauser field |
| [Advanced 03: Several NV centres](../../example/tutorial/advanced/03_multiNV/README.md) | `advanced/03_multiNV` | several NV centres as one register |
| [Advanced 04: DFT hyperfine and quadrupole tensors](../../example/tutorial/advanced/04_DFT_tensors/README.md) | `advanced/04_DFT_tensors` | the hyperfine and quadrupole files from a VASP run (`example/tools/SpinTensorReader`), and what they change |
| [Advanced 05: Coordinate rotation](../../example/tutorial/advanced/05_coordinate_rotation/README.md) | `advanced/05_coordinate_rotation` | a bath written in another frame (`frame.rotation`) |
| [Advanced 06: More terms and models](../../example/tutorial/advanced/06_terms_and_models/README.md) | `advanced/06_terms_and_models` | disjoint clusters, the hyperfine-mediated term, meCCE, the Knight shift |
| [Advanced 07: pCCE and bundles](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) | `advanced/07_pCCE_and_bundle` | clusters of groups: partitions (pCCE), bundles, adaptive refinement |
| [Advanced 08: Analysis of the output](../../example/tutorial/advanced/08_analysis/README.md) | `advanced/08_analysis` | the analyzer (error correction, average, T2 fit) and the cluster reports |
| [Advanced 09: Phase noise (method noise)](../../example/tutorial/advanced/09_noise/README.md) | `advanced/09_noise` | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum; its correlation C(t) and spectrum S(f) |
| [Advanced 10: Pseudospin pairs](../../example/tutorial/advanced/10_pseudospin/README.md) | `advanced/10_pseudospin` | method `pseudospin`: the pairs of a nuclear bath in the pseudospin model, against cCCE2 |
| [Advanced 11: The GPU backend](../../example/tutorial/advanced/11_GPU/README.md) | `advanced/11_GPU` | the CUDA backend against the CPU, and the automatic GPU options |
| [Advanced 12: Automatic convergence](../../example/tutorial/advanced/12_auto_convergence/README.md) | `advanced/12_auto_convergence` | the automatic convergence of the time window, rdip, rbath and the configurations |
| [Advanced 13: Runs on a cluster (Slurm)](../../example/tutorial/advanced/13_server_runs/README.md) | `advanced/13_server_runs` | a production run as a Slurm job: sizing it, one job, a job array |

<a id="WhatToUse"></a>
## What should I use for my case?

| Situation | Use | Example |
|---|---|---|
| a defect in a nuclear-spin bath, field along the defect axis | cCCE (default) | [Basic 01: NV in ¹³C, from a bath to T2](../../example/tutorial/basic/01_first_run/README.md) |
| an electron-spin bath (P1, ...) | the `Defect` block, sampled states | [Basic 06: An electron-spin bath (P1)](../../example/tutorial/basic/06_defect_bath/README.md) |
| hyperfine and quadrupole from DFT | the tensor files | [Advanced 04: DFT hyperfine and quadrupole tensors](../../example/tutorial/advanced/04_DFT_tensors/README.md) |
| near a level crossing, a low or tilted field, density matrix | gCCE | [Basic 05: cCCE against gCCE](../../example/tutorial/basic/05_cCCE_vs_gCCE/README.md) |
| several qubits | gCCE with `intmap` | [Advanced 03: Several NV centres](../../example/tutorial/advanced/03_multiNV/README.md) |
| strongly coupled bath spins | pCCE or bundles | [Advanced 07: pCCE and bundles](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) |
| large clusters | the GPU backend | [Advanced 11: The GPU backend](../../example/tutorial/advanced/11_GPU/README.md) |
| a bath written in the crystal frame | `frame.rotation` | [Advanced 05: Coordinate rotation](../../example/tutorial/advanced/05_coordinate_rotation/README.md) |
| a noise spectrum, or a model noise without a bath | method noise | [Phase noise](../../example/tutorial/advanced/09_noise/README.md) |
| which nuclear pairs drive the decay | method pseudospin | [Advanced 10: Pseudospin pairs](../../example/tutorial/advanced/10_pseudospin/README.md) |
