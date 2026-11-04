# Examples

`example/tutorial/` holds the worked examples in two levels, each with a page here with its text, its input files and
its expected figures.

- `basic/` is the input itself: the qubit, the bath, the pulses, the two main methods and an electron-spin bath. Each
  runs in seconds on a laptop; read them in order.
- `advanced/` builds on them: other defects and baths, DFT tensors, more methods, the analysis of the output and runs on
  a cluster. Pick the topic you need.
- Each topic folder holds one folder per calculation, named `<method><order>_<pulse>_<approach>`. For example,
  `cCCE2_Hahn_ensemble` is cCCE to order 2, Hahn echo, the bath ensemble. The approach is `ensemble`, `states` or what
  the calculation compares.
- A calculation folder has its inputs, a `run.sh` and the `expected.png` a rerun is compared with. The inputs
  (`ccein*.jsonc`) are JSON with `//` comments and hold only the keys that matter; `main.out -e` writes every key.
- For every calculation, the page says what it shows, the keys that matter, how to run it (time, ranks, memory), the
  expected result and what to change next.

**Basic**

| Page | Topic |
|---|---|
| [basic/01_first_run](../../example/tutorial/basic/01_first_run/README.md) | NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds on a laptop |
| [basic/02_bath_and_states](../../example/tutorial/basic/02_bath_and_states/README.md) | averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`); the built-in isotope table |
| [basic/03_pulses](../../example/tutorial/basic/03_pulses/README.md) | Hahn echo, CPMG, a sequence written out (`Manual`) |
| [basic/04_cCCE_vs_gCCE](../../example/tutorial/basic/04_cCCE_vs_gCCE/README.md) | the two methods; the Knight shift of cCCE |
| [basic/05_defect_bath](../../example/tutorial/basic/05_defect_bath/README.md) | an electron-spin bath: P1 with the `Defect` block, sampled states, the mean field of the bath |

**Advanced**

| Page | Topic |
|---|---|
| [advanced/01_P1](../../example/tutorial/advanced/01_P1/README.md) | the P1 bath further: order 4, gCCE at zero field, the density matrix, Vp, a DEER pump |
| [advanced/02_hBN](../../example/tutorial/advanced/02_hBN/README.md) | V_B⁻ in h-BN with DFT tensors: sampled states, the Overhauser field |
| [advanced/03_multiNV](../../example/tutorial/advanced/03_multiNV/README.md) | several NV centres as one register |
| [advanced/04_DFT_tensors](../../example/tutorial/advanced/04_DFT_tensors/README.md) | the hyperfine and quadrupole files from a VASP run (`example/tools/SpinTensorReader`), and what they change |
| [advanced/05_coordinate_rotation](../../example/tutorial/advanced/05_coordinate_rotation/README.md) | a bath written in another frame (`frame.rotation`) |
| [advanced/06_terms_and_models](../../example/tutorial/advanced/06_terms_and_models/README.md) | disjoint clusters, the hyperfine-mediated term, meCCE |
| [advanced/07_pCCE_and_bundle](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) | clusters of groups: partitions (pCCE), bundles, adaptive refinement |
| [advanced/08_analysis](../../example/tutorial/advanced/08_analysis/README.md) | the analyzer (error correction, average, T2 fit) and the cluster reports |
| [advanced/09_noise](../../example/tutorial/advanced/09_noise/README.md) | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum; its correlation C(t) and spectrum S(f) |
| [advanced/10_pseudospin](../../example/tutorial/advanced/10_pseudospin/README.md) | method `pseudospin`: the pairs of a nuclear bath in the pseudospin model, against cCCE2 |
| [advanced/11_GPU](../../example/tutorial/advanced/11_GPU/README.md) | the CUDA backend against the CPU, and the automatic GPU options |
| [advanced/12_auto_convergence](../../example/tutorial/advanced/12_auto_convergence/README.md) | the automatic convergence of the time window, rdip, rbath and the configurations |
| [advanced/13_server_runs](../../example/tutorial/advanced/13_server_runs/README.md) | a production run as a Slurm job: sizing it, one job, a job array |

<a id="RunningOne"></a>
## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/basic/01_first_run
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` does three things:

1. It finds the binary through `git rev-parse --show-toplevel`. `CCEX_BIN` names another one, e.g. a CMake build's
   `build/bin/main.out`; `CCEX_GPU_BIN` names the GPU binary of `advanced/11_GPU`.
2. It runs every input with `mpirun -n $NP ... -o results`. The outputs are named after the input: `ccein_wiDiv`,
   `ccein_k1_state3_wiDiv`, ... (see [Output files](OutputFiles.md)).
3. It runs `example/tutorial/plot_coherence.py`, which writes `results/coherence.png` and prints one number per curve:
   T2, the last time $|L| \ge 1/e$.

A few `run.sh` also run a variant of the input, written into `results/` by
`example/tools/InputKeys/set_input_keys.py`. Inputs name their files relative to the calculation folder (`../bath/`,
`../../../basic/01_first_run/bath/`), so run `run.sh` from there. Git ignores `results/`.

Two files are shared by the tutorials:

- `example/tutorial/plot_coherence.py` makes the figures and numbers.
- `example/tutorial/make_natab_baths.sh` makes the natural-abundance ¹³C configurations of `basic/02_bath_and_states` and
  `advanced/12_auto_convergence` with
  `example/tools/BathGenerator`.

Each `expected.png` was made by that folder's `run.sh` on a workstation (i7-14700K, WSL, 4 ranks). Those of `advanced/11_GPU`
were made on a GPU server node (idm4-gpu). Times and memory on the pages are from the same workstation.

<a id="WhatToUse"></a>
## What should I use for my case?

| Situation | Method and keys | Example |
|---|---|---|
| a defect in a nuclear-spin bath (¹³C, ²⁹Si, ...), high field along the defect axis | cCCE (the default), `order` 2, ensemble | [basic/01_first_run](../../example/tutorial/basic/01_first_run/README.md), [basic/02_bath_and_states](../../example/tutorial/basic/02_bath_and_states/README.md) (`cCCE2_Hahn_ensemble`) |
| several bath configurations from one input | `configuration.bathfile` with `{k}` and `configuration.k_range` | [basic/02_bath_and_states](../../example/tutorial/basic/02_bath_and_states/README.md) (`cCCE2_Hahn_ensemble`) |
| the bath in sampled pure states (electron baths, polarised or single-shot) | `configuration.s_range` and `seed` (or `bath_state_file`, `example/tools/StateFiles`) | [basic/02_bath_and_states](../../example/tutorial/basic/02_bath_and_states/README.md) (`cCCE2_Hahn_states`), [basic/05_defect_bath](../../example/tutorial/basic/05_defect_bath/README.md) |
| dynamical decoupling | `pulse.count` (1 Hahn, N CPMG-N), `pulse.name` `Manual` with `pulse.sequence` for any other π sequence | [basic/03_pulses](../../example/tutorial/basic/03_pulses/README.md) |
| an electron-spin bath (P1, other paramagnetic defects) | the `Defect` block (axes, detunings, on-site nuclei), `s_range`, check `order` 4 | [basic/05_defect_bath](../../example/tutorial/basic/05_defect_bath/README.md) |
| DFT hyperfine and quadrupole tensors (2D materials, defects with strong contact terms) | `hamiltonian.common.hyperfine` / `quadrupole` with the files of `example/tools/SpinTensorReader` ([Tools](Tools.md)) | [advanced/04_DFT_tensors](../../example/tutorial/advanced/04_DFT_tensors/README.md), [advanced/02_hBN](../../example/tutorial/advanced/02_hBN/README.md) |
| low or tilted field, a spin-1 qubit whose zero-field splitting matters, the qubit density matrix | gCCE (`method` gCCE, `Qubit.intmap` self entry), `quantity` dm | [basic/04_cCCE_vs_gCCE](../../example/tutorial/basic/04_cCCE_vs_gCCE/README.md) (`gCCE2_Hahn_vsCCE`), [advanced/01_P1](../../example/tutorial/advanced/01_P1/README.md) (`gCCE2_Hahn_zfs`) |
| cCCE in a tilted field | `hamiltonian.cCCE.knight` | [basic/04_cCCE_vs_gCCE](../../example/tutorial/basic/04_cCCE_vs_gCCE/README.md) (`cCCE2_Hahn_knightShift`) |
| several qubits | `Qubit.qubit` list with `intmap`, gCCE | [advanced/03_multiNV](../../example/tutorial/advanced/03_multiNV/README.md) |
| strongly coupled groups of bath spins (dense electron baths) | pCCE (`cluster.mode` pcce, `partition_size`) or bundles (`cluster.mode` bundle, `bundle.max_dim`, optionally `cluster.adaptive`) | [advanced/07_pCCE_and_bundle](../../example/tutorial/advanced/07_pCCE_and_bundle/README.md) |
| large clusters (Hilbert dimension of hundreds or more) and a GPU | `execution.backend` auto with `bin/main_gpu.out` | [advanced/11_GPU](../../example/tutorial/advanced/11_GPU/README.md) |
| not sure which rbath, rdip, time window or number of configurations | the `auto_convergence` block | [advanced/12_auto_convergence](../../example/tutorial/advanced/12_auto_convergence/README.md) |
| averaging, the T2 fit, the division errors of `_wiDiv` | the analyzer of `example/analysis` (`-err`, `-ea`, `-fit_w`) | [advanced/08_analysis](../../example/tutorial/advanced/08_analysis/README.md) |
| a run that needs a cluster | a Slurm job or a job array; `configuration.seed` one per configuration to split a run | [advanced/13_server_runs](../../example/tutorial/advanced/13_server_runs/README.md) |
| a bath file written in the crystal frame (defect axis not along z) | the `frame.rotation` block, `qubit_axis` | [advanced/05_coordinate_rotation](../../example/tutorial/advanced/05_coordinate_rotation/README.md) |
| a decay from many weak noise sources, a noise spectrum S(f), or a model noise without a bath | `method` noise, `hamiltonian.noise` (`spectrum`, `compute`) | [advanced/09_noise](../../example/tutorial/advanced/09_noise/README.md) |
| which nuclear pairs drive a decay (any nuclear spin, with quadrupoles) | `method` pseudospin, its channel table | [advanced/10_pseudospin](../../example/tutorial/advanced/10_pseudospin/README.md) |

The table [What runs where](Overview.md#WhatRunsWhere) shows which combinations of method, cluster mode and backend
run at all.
