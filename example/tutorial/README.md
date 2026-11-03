# example/tutorial/

Thirteen worked examples, in the order to learn them. Each numbered folder is one topic and holds one folder per
calculation, named `<method><order>_<pulse>_<approach>`. For example, `cCCE2_Hahn_ensemble` is cCCE to order 2 with a
Hahn echo on the bath ensemble; the approach is `ensemble`, `states` or what the calculation compares.

A calculation folder has its inputs (`ccein*.jsonc`, JSON with `//` comments that explain every key), a `run.sh` and an
`expected.png` to compare a rerun with. The topic's `README.md` says, for each calculation, what it shows, which keys
matter, how to run it (time, ranks, memory), the expected result and what to try next.

| Folder | Topic | Calculations |
|---|---|---|
| `00_quickstart/` | NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds on a laptop | the folder itself |
| `01_ensemble_and_states/` | averaging over bath configurations (`k_range`) and over sampled bath states (`s_range`) | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_states`, `cCCE2_Hahn_builtinIsotopes` |
| `02_pulses/` | Hahn echo, CPMG, a sequence written out (`Manual`) | `cCCE2_CPMG_ensemble`, `cCCE2_Manual_ensemble` |
| `03_P1/` | an electron-spin bath: P1 at 1 and 10 ppm, Vp | `cCCE2_Hahn_states`, `cCCE4_Hahn_states`, `cCCE2_Hahn_onsite14N`, `gCCE2_Hahn_zfs`, `gCCE2_Ramsey_population`, `cCCE2_Hahn_VpBath`, `gCCE2_Hahn_DEER`, `cCCE2_Hahn_meanField` |
| `04_hBN/` | V_B⁻ in h-BN with DFT hyperfine and quadrupole tensors | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_hfMediated`, `gCCE2_Hahn_states`, `gCCE2_Hahn_OH` |
| `05_multiNV/` | several NV centres as one register | `gCCE2_Hahn_oneVsTwoNV`, `gCCE2_HahnRamsey_NV2`, `gCCE2_HahnRamsey_NV6` |
| `06_cCCE_vs_gCCE/` | the two methods, and the options of cCCE | `gCCE2_Hahn_vsCCE`, `cCCE2_Hahn_dsj`, `cCCE2_Hahn_knightShift`, `gCCE1_Hahn_handwritten`, `meCCE2_Hahn_bathRelaxation` |
| `07_pCCE_and_bundle/` | clusters of groups: partitions (pCCE), bundles, adaptive refinement | `pCCE2_Hahn_states`, `cCCE2_Hahn_bundle`, `cCCE2_Hahn_adaptive`, `cCCE2_Hahn_regroup`, `cCCE2_Hahn_clusterReport` |
| `08_GPU/` | the CUDA backend against the CPU, and the automatic GPU options | `cCCE5_Hahn_backends`, `cCCE5_Hahn_autoOptions`, `gCCE4_Hahn_backends`, `pCCE2_Hahn_backends` |
| `09_auto_convergence/` | the automatic convergence of the time window, rdip, rbath and the configurations | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_ensemble_noCandidates` |
| `10_coordinate_rotation/` | a bath written in another frame (`frame.rotation`) | `cCCE2_Hahn_rotation` |
| `11_noise/` | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum; C(t) and S(f) | `noise2_Hahn_vsCCE`, `noiseModel_OU`, `noiseSpectrum_13C` |
| `12_pseudospin/` | method `pseudospin`: the pairs of a nuclear bath in the pseudospin model, against cCCE2 | `pseudospin2_Hahn_vsCCE` |

## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/00_quickstart
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` finds the binary from the repository root (`git rev-parse --show-toplevel`); `CCEX_BIN` names another one,
such as a CMake build's `build/bin/main.out`, and `CCEX_GPU_BIN` the GPU binary for `08_GPU`. It runs every input with
`mpirun -n $NP ... -o results`, so the outputs are named after the input (`ccein_wiDiv`, `ccein_k1_state3_wiDiv`, ...).
Then `plot_coherence.py` writes `results/coherence.png` and prints one number per curve: T2, the last time |L| ≥ 1/e.
Some `run.sh` also run a variant of the input that `example/tools/InputKeys/set_input_keys.py` writes into `results/`.
Inputs name their files relative to the calculation folder, so run `run.sh` from there. Git ignores `results/`.

Two scripts here are shared: `plot_coherence.py` (the figures and numbers) and `make_natab_baths.sh`
(natural-abundance ¹³C configurations made with `tools/BathGenerator`, for `01` and `09`). Each `expected.png` was
made by its folder's `run.sh` on a workstation (i7-14700K, WSL, 4 ranks), except those of `08_GPU`, made on a GPU
server; `05_multiNV/gCCE2_HahnRamsey_NV6` has none yet.

## Folders before 2.0.0

The tutorials were reorganised in learning order for 2.0.0. The table of where each earlier folder went is in the
history: `git show 2c14eb918:example/tutorial/README.md`.
