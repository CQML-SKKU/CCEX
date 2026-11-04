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
| `03_pulses/` | Hahn echo, CPMG, a sequence written out (`Manual`) | `cCCE2_CPMG_ensemble`, `cCCE2_Manual_ensemble` |
| `04_cCCE_vs_gCCE/` | the two methods; the Knight shift of cCCE | `gCCE2_Hahn_vsCCE`, `cCCE2_Hahn_knightShift`, `gCCE1_Hahn_handwritten` |
| `05_defect_bath/` | an electron-spin bath: P1 with the `Defect` block, sampled states, the mean field | `cCCE2_Hahn_states`, `cCCE2_Hahn_onsite14N`, `cCCE2_Hahn_meanField` |

**advanced/**

| Folder | Topic | Calculations |
|---|---|---|
| `01_P1/` | the P1 bath further: order 4, gCCE at zero field, the density matrix, Vp, a DEER pump | `cCCE4_Hahn_states`, `gCCE2_Hahn_zfs`, `gCCE2_Ramsey_population`, `cCCE2_Hahn_VpBath`, `gCCE2_Hahn_DEER` |
| `02_hBN/` | V_B⁻ in h-BN with DFT tensors: sampled states, the Overhauser field | `gCCE2_Hahn_states`, `gCCE2_Hahn_Overhauser` |
| `03_multiNV/` | several NV centres as one register | `gCCE2_Hahn_oneVsTwoNV`, `gCCE2_HahnRamsey_NV2`, `gCCE2_HahnRamsey_NV3` |
| `04_DFT_tensors/` | the hyperfine and quadrupole files from a VASP run, and what they change | `makeTensors_hBN`, `cCCE2_Hahn_tensors` |
| `05_coordinate_rotation/` | a bath written in another frame (`frame.rotation`) | `cCCE2_Hahn_rotation` |
| `06_terms_and_models/` | disjoint clusters, the hyperfine-mediated term, meCCE | `cCCE2_Hahn_dsj`, `cCCE2_Hahn_hfMediated`, `meCCE2_Hahn_bathRelaxation` |
| `07_pCCE_and_bundle/` | clusters of groups: partitions (pCCE), bundles, adaptive refinement | `pCCE2_Hahn_states`, `cCCE2_Hahn_bundle`, `cCCE2_Hahn_adaptive`, `cCCE2_Hahn_regroup` |
| `08_analysis/` | the analyzer (error correction, average, T2 fit) and the cluster reports | `cCCE2_Hahn_analyzer`, `cCCE2_Hahn_clusterReport` |
| `09_noise/` | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum; C(t) and S(f) | `noise2_Hahn_vsCCE`, `noiseModel_OU`, `noiseSpectrum_13C` |
| `10_pseudospin/` | method `pseudospin`: the pairs of a nuclear bath in the pseudospin model, against cCCE2 | `pseudospin2_Hahn_vsCCE` |
| `11_GPU/` | the CUDA backend against the CPU, and the automatic GPU options | `cCCE5_Hahn_backends`, `cCCE5_Hahn_autoOptions`, `gCCE4_Hahn_backends`, `pCCE2_Hahn_backends` |
| `12_auto_convergence/` | the automatic convergence of the time window, rdip, rbath and the configurations | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_ensemble_noCandidates` |
| `13_server_runs/` | a production run as a Slurm job: sizing it, one job, a job array | `cCCE2_Hahn_slurm` |

## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/basic/01_first_run
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` finds the binary from the repository root (`git rev-parse --show-toplevel`); `CCEX_BIN` names another one,
such as a CMake build's `build/bin/main.out`, and `CCEX_GPU_BIN` the GPU binary for `advanced/11_GPU`. It runs every input with
`mpirun -n $NP ... -o results`, so the outputs are named after the input (`ccein_wiDiv`, `ccein_k1_state3_wiDiv`, ...).
Then `plot_coherence.py` writes `results/coherence.png` and prints one number per curve: T2, the last time |L| ≥ 1/e.
Some `run.sh` also run a variant of the input that `example/tools/InputKeys/set_input_keys.py` writes into `results/`.
Inputs name their files relative to the calculation folder, so run `run.sh` from there. Git ignores `results/`.

Two scripts here are shared: `plot_coherence.py` (the figures and numbers) and `make_natab_baths.sh`
(natural-abundance ¹³C configurations made with `tools/BathGenerator`, for `basic/02_bath_and_states` and `advanced/12_auto_convergence`). Each `expected.png` was
made by its folder's `run.sh` on a workstation (i7-14700K, WSL, 4 ranks), except those of `advanced/11_GPU`, made on a GPU
server node.

## Folders before 2.1.0

2.1.0 split the tutorials into `basic/` and `advanced/`. The folders of 2.0 (`00_quickstart` ... `12_pseudospin`) are in
`git show v2.0.1:example/tutorial/README.md`, and the table of where the folders before 2.0.0 went in
`git show 2c14eb918:example/tutorial/README.md`.
