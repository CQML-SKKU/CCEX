# Examples

Twelve worked examples live in `example/tutorial/`, in the order to learn them; each has a page here with its text,
its input files and its expected figures. Each topic folder holds one folder per calculation, named
`<method><order>_<pulse>_<approach>` (`cCCE2_Hahn_ensemble`: cCCE to order 2, Hahn echo, the bath ensemble; the
approach is `ensemble`, `states` or what the calculation compares). A calculation folder has its inputs
(`ccein*.jsonc`, JSON with `//` comments, only the keys that matter; `main.out -e` writes every key), a `run.sh` and the
`expected.png` a rerun is compared with. Each page says, for every calculation, what it shows, the keys that matter,
how to run it (time, ranks, memory), the expected result and what to change next.

| Page | Topic |
|---|---|
| [00_quickstart](../../example/tutorial/00_quickstart/README.md) | NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds on a laptop |
| [01_ensemble_and_states](../../example/tutorial/01_ensemble_and_states/README.md) | averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`); the built-in isotope table |
| [02_pulses](../../example/tutorial/02_pulses/README.md) | Hahn echo, CPMG, a sequence written out (`Manual`) |
| [03_P1](../../example/tutorial/03_P1/README.md) | an electron-spin bath: P1 at 1 and 10 ppm, Vp, a DEER pump, the mean field of the bath |
| [04_hBN](../../example/tutorial/04_hBN/README.md) | V_B⁻ in h-BN with DFT hyperfine and quadrupole tensors |
| [05_multiNV](../../example/tutorial/05_multiNV/README.md) | several NV centres as one register |
| [06_cCCE_vs_gCCE](../../example/tutorial/06_cCCE_vs_gCCE/README.md) | the two methods, the options of cCCE (dsj, Knight shift), meCCE |
| [07_pCCE_and_bundle](../../example/tutorial/07_pCCE_and_bundle/README.md) | clusters of groups: partitions (pCCE), bundles, adaptive refinement, the cluster reports |
| [08_GPU](../../example/tutorial/08_GPU/README.md) | the CUDA backend against the CPU |
| [09_auto_convergence](../../example/tutorial/09_auto_convergence/README.md) | the automatic convergence of the time window, rdip, rbath and the configurations |
| [10_coordinate_rotation](../../example/tutorial/10_coordinate_rotation/README.md) | a bath written in another frame (`frame.rotation`) |
| [11_noise](../../example/tutorial/11_noise/README.md) | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum |

<a id="RunningOne"></a>
## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/00_quickstart
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` finds the binary through `git rev-parse --show-toplevel` (`CCEX_BIN` names another one, e.g. a CMake build's
`build/bin/main.out`; `CCEX_GPU_BIN` the GPU binary of `08_GPU`), runs every input with `mpirun -n $NP ... -o results`
(the outputs are named after the input: `ccein_wiDiv`, `ccein_k1_state3_wiDiv`, ...; see [Output files](OutputFiles.md))
and then `example/tutorial/plot_coherence.py`, which writes `results/coherence.png` and prints one number per curve
(T2 = the last time |L| ≥ 1/e). A few `run.sh` also run a variant of the input written by
`example/tools/InputKeys/set_input_keys.py` into `results/`. Inputs name their files relative to the calculation
folder (`../bath/`, `../../00_quickstart/bath/`), so `run.sh` is run from there. Git ignores `results/`.

Shared by the tutorials: `example/tutorial/plot_coherence.py` (the figures and numbers) and
`example/tutorial/make_natab_baths.sh` (natural-abundance ¹³C configurations made with `example/tools/BathGenerator`,
for `01` and `09`). Each `expected.png` was made by that folder's `run.sh` on a workstation (i7-14700K, WSL, 4 ranks);
the three calculations of `08_GPU` and `05_multiNV/gCCE2_HahnRamsey_NV6` need a server and have no `expected.png` yet.
Times and memory on the pages are from the same workstation.

<a id="WhatToUse"></a>
## What should I use for my case?

| Situation | Method and keys | Example |
|---|---|---|
| a defect in a nuclear-spin bath (¹³C, ²⁹Si, ...), high field along the defect axis | cCCE (the default), `order` 2, ensemble | [00_quickstart](../../example/tutorial/00_quickstart/README.md), [01_ensemble_and_states](../../example/tutorial/01_ensemble_and_states/README.md) (`cCCE2_Hahn_ensemble`) |
| several bath configurations from one input | `configuration.bathfile` with `{k}` and `configuration.k_range` | [01_ensemble_and_states](../../example/tutorial/01_ensemble_and_states/README.md) (`cCCE2_Hahn_ensemble`) |
| the bath in sampled pure states (electron baths, polarised or single-shot) | `configuration.s_range` and `seed` (or `bath_state_file`, `example/tools/StateFiles`) | [01_ensemble_and_states](../../example/tutorial/01_ensemble_and_states/README.md) (`cCCE2_Hahn_states`), [03_P1](../../example/tutorial/03_P1/README.md) |
| dynamical decoupling | `pulse.count` (1 Hahn, N CPMG-N), `pulse.name` `Manual` with `pulse.sequence` for any other π sequence | [02_pulses](../../example/tutorial/02_pulses/README.md) |
| an electron-spin bath (P1, other paramagnetic defects) | the `Defect` block (axes, detunings, on-site nuclei), `s_range`, check `order` 4 | [03_P1](../../example/tutorial/03_P1/README.md) |
| DFT hyperfine and quadrupole tensors (2D materials, defects with strong contact terms) | `hamiltonian.common.hyperfine` / `quadrupole` with the files of `example/tools/SpinTensorReader` ([Tools](Tools.md)) | [04_hBN](../../example/tutorial/04_hBN/README.md) |
| low or tilted field, a spin-1 qubit whose zero-field splitting matters, the qubit density matrix | gCCE (`method` gCCE, `Qubit.intmap` self entry), `quantity` dm | [06_cCCE_vs_gCCE](../../example/tutorial/06_cCCE_vs_gCCE/README.md) (`gCCE2_Hahn_vsCCE`), [03_P1](../../example/tutorial/03_P1/README.md) (`gCCE2_Hahn_zfs`) |
| cCCE in a tilted field | `hamiltonian.cCCE.knight` | [06_cCCE_vs_gCCE](../../example/tutorial/06_cCCE_vs_gCCE/README.md) (`cCCE2_Hahn_knightShift`) |
| several qubits | `Qubit.qubit` list with `intmap`, gCCE | [05_multiNV](../../example/tutorial/05_multiNV/README.md) |
| strongly coupled groups of bath spins (dense electron baths) | pCCE (`cluster.mode` pcce, `partition_size`) or bundles (`cluster.mode` bundle, `bundle.max_dim`, optionally `cluster.adaptive`) | [07_pCCE_and_bundle](../../example/tutorial/07_pCCE_and_bundle/README.md) |
| large clusters (Hilbert dimension of hundreds or more) and a GPU | `execution.backend` auto with `bin/main_gpu.out` | [08_GPU](../../example/tutorial/08_GPU/README.md) |
| not sure which rbath, rdip, time window or number of configurations | the `auto_convergence` block | [09_auto_convergence](../../example/tutorial/09_auto_convergence/README.md) |
| a bath file written in the crystal frame (defect axis not along z) | the `frame.rotation` block, `qubit_axis` | [10_coordinate_rotation](../../example/tutorial/10_coordinate_rotation/README.md) |

Which combinations of method, cluster mode and backend run at all is the table
[What runs where](Overview.md#WhatRunsWhere).
