# example/tutorial/

Worked examples in two levels. `basic/` is the input itself: the qubit, the bath, the pulses, the two main methods and
an electron-spin bath; each runs in seconds on a laptop, and they are meant to be read in order. `advanced/` builds on
them: other defects and baths, DFT tensors, more methods, the analysis of the output and runs on a cluster; pick the
topic you need. Each numbered folder is one topic and holds one folder per calculation, named
`<method><order>_<pulse>_<approach>`. For example, `cCCE2_Hahn_ensemble` is cCCE to order 2 with a Hahn echo on the bath
ensemble; the approach is `ensemble`, `states` or what the calculation compares.

A calculation folder has its inputs (`ccein*.jsonc`, JSON with `//` comments that explain every key), a `run.sh` and an
`expected.png` to compare a rerun with. The topic's `README.md` says, for each calculation, what it shows, which keys
matter, how to run it (time, ranks, memory), the expected result and what to try next.

**basic/**

| Folder | Topic | Calculations |
|---|---|---|
| `01_first_run/` | NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds on a laptop | the folder itself |
| `02_bath_and_states/` | averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`) | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_states`, `cCCE2_Hahn_builtinIsotopes` |
| `03_convergence/` | checking `rbath`, `rdip`, `order`, the number of configurations and the time axis | `cCCE2_Hahn_rbath`, `cCCE2_Hahn_rdip`, `cCCE_Hahn_order`, `cCCE2_Hahn_configurations`, `cCCE2_Hahn_timeAxis` |
| `04_pulses/` | Hahn echo, CPMG, a sequence written out (`Manual`) | `cCCE2_CPMG_ensemble`, `cCCE2_Manual_ensemble` |
| `05_cCCE_vs_gCCE/` | the two methods: the same at a high field, different near the level crossing | `gCCE2_Hahn_vsCCE`, `gCCE2_Hahn_levelCrossing` |
| `06_defect_bath/` | an electron-spin bath: P1 with the `Defect` block, sampled states, the mean field | `cCCE2_Hahn_states`, `cCCE2_Hahn_meanField` |

**advanced/**, in four groups:

*1. Running and handling a calculation*

| Folder | Topic | Calculations |
|---|---|---|
| `01_auto_convergence/` | the automatic convergence of the time window, rdip, rbath, the configurations and the states: found, given as candidates, or fixed | `cCCE2_Hahn_ensemble_noCandidates`, `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_ensemble_fixed`, `cCCE2_Hahn_states` |
| `02_analysis/` | the analyzer (error correction, average, T2 fit) and the cluster reports | `cCCE2_Hahn_analyzer`, `cCCE2_Hahn_clusterReport` |
| `03_server_runs/` | a production run as a Slurm job: sizing it, one job, a job array | `cCCE2_Hahn_slurm` |
| `04_GPU/` | the CUDA backend against the CPU, and the automatic GPU options | `cCCE5_Hahn_backends`, `cCCE5_Hahn_autoOptions`, `gCCE4_Hahn_backends` |

*2. Other systems and inputs*

| Folder | Topic | Calculations |
|---|---|---|
| `05_P1/` | electron-spin baths further: making a P1 bath, the ¹⁴N as a detuning or a spin (2e, 2e2n), order 2 checked at order 4, a DEER pump | `cCCE2_Hahn_makeBath`, `cCCE2_Hahn_onsite14N`, `cCCE4_Hahn_states`, `gCCE2_Hahn_DEER` |
| `06_hBN/` | V_B⁻ in h-BN with DFT tensors: sampled states, the Overhauser field | `gCCE2_Hahn_states`, `gCCE2_Hahn_Overhauser` |
| `07_DFT_tensors/` | the hyperfine and quadrupole files from a VASP run, and what they change | `makeTensors_hBN`, `cCCE2_Hahn_tensors` |
| `08_multi_qubits/` | several NV centres as one register | `gCCE2_Hahn_oneVsTwoNV`, `gCCE2_HahnRamsey_NV2`, `gCCE2_HahnRamsey_NV3` |
| `09_coordinate_rotation/` | a bath written in another frame (`frame.rotation`) | `cCCE2_Hahn_rotation` |

*3. Hamiltonian terms and cluster shapes*

| Folder | Topic | Calculations |
|---|---|---|
| `10_cCCE_corrections/` | corrections to cCCE: the hyperfine-mediated term, the Knight shift | `cCCE2_Hahn_hfMediated`, `cCCE2_Hahn_knightShift` |
| `11_cluster_shapes/` | cluster shapes: disjoint clusters (dsj), partitions (pCCE), bundles, adaptive refinement | `cCCE2_Hahn_dsj`, `pCCE2_Hahn_states`, `cCCE2_Hahn_bundle`, `cCCE2_Hahn_adaptive`, `cCCE2_Hahn_regroup` |

*4. Other models of the bath*

| Folder | Topic | Calculations |
|---|---|---|
| `12_dissipation/` | dissipation: meCCE, Lindblad relaxation of the bath spins | `meCCE2_Hahn_bathRelaxation` |
| `13_noise/` | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum; C(t) and S(f) | `noise2_Hahn_vsCCE`, `noiseModel_OU`, `noiseSpectrum_13C` |
| `14_pseudospin/` | method `pseudospin`: the pairs of a nuclear bath in the pseudospin model, against cCCE2 | `pseudospin2_Hahn_vsCCE` |

## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/basic/01_first_run
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` finds the binary from the repository root (`git rev-parse --show-toplevel`); `CCEX_BIN` names another one,
such as a CMake build's `build/bin/main.out`, and `CCEX_GPU_BIN` the GPU binary for `advanced/04_GPU`. It runs every input with
`mpirun -n $NP ... -o results`, so the outputs are named after the input (`ccein_wiDiv`, `ccein_k1_state3_wiDiv`, ...).
Then `plot_coherence.py` writes `results/coherence.png` and prints one number per curve: T2, the last time |L| ≥ 1/e.
Some `run.sh` also run a variant of the input that `example/tools/InputKeys/set_input_keys.py` writes into `results/`.
Inputs name their files relative to the calculation folder, so run `run.sh` from there. Git ignores `results/`.

Two scripts here are shared: `plot_coherence.py` (the figures and numbers) and `make_natab_baths.sh`
(natural-abundance ¹³C configurations made with `tools/BathGenerator`, for `basic/02_bath_and_states` and `advanced/01_auto_convergence`). Each `expected.png` was
made by its folder's `run.sh` on a desktop CPU with 4 ranks, except those of `advanced/04_GPU`, made on a GPU
server node.

## Folders before 2.1.0

2.1.0 split the tutorials into `basic/` and `advanced/`. The folders of 2.0 (`00_quickstart` ... `12_pseudospin`) are in
`git show v2.0.1:example/tutorial/README.md`, and the table of where the folders before 2.0.0 went in
`git show 2c14eb918:example/tutorial/README.md`.
