# example/tutorial/

Twelve worked examples in the order to learn them. Each numbered folder is one topic and holds one folder per
calculation, named `<method><order>_<pulse>_<approach>` (`cCCE2_Hahn_ensemble`: cCCE to order 2, Hahn echo, the bath
ensemble; the approach is `ensemble`, `states` or what the calculation compares). A calculation folder has its inputs
(`ccein*.jsonc`, JSON with `//` comments, only the keys that matter; `main.out -e` writes every key), a `run.sh` and the
`expected.png` a rerun is compared with. The topic's `README.md` says, for every calculation, what it shows, the keys
that matter, how to run it (time, ranks, memory), the expected result and what to change next.

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
| `08_GPU/` | the CUDA backend against the CPU | `cCCE5_Hahn_backends`, `gCCE4_Hahn_backends`, `pCCE2_Hahn_backends` |
| `09_auto_convergence/` | the automatic convergence of the time window, rdip, rbath and the configurations | `cCCE2_Hahn_ensemble`, `cCCE2_Hahn_ensemble_noCandidates` |
| `10_coordinate_rotation/` | a bath written in another frame (`frame.rotation`) | `cCCE2_Hahn_rotation` |
| `11_noise/` | method `noise`: the Gaussian coherence of the phase noise, from the bath or from a model spectrum | `noise2_Hahn_vsCCE`, `noiseModel_OU` |

## Running one

```
bash do_compile.sh wsl-gnu                     # once, from the repository root; makes bin/main.out
cd example/tutorial/00_quickstart
NP=4 bash run.sh                               # results/; needs python3 with numpy and matplotlib for the plot
```

`run.sh` finds the binary through `git rev-parse --show-toplevel` (`CCEX_BIN` names another one, e.g. a CMake build's
`build/bin/main.out`; `CCEX_GPU_BIN` the GPU binary of `08_GPU`), runs every input with `mpirun -n $NP ... -o results`
(the outputs are named after the input: `ccein_wiDiv`, `ccein_k1_state3_wiDiv`, ...) and then
`plot_coherence.py`, which writes `results/coherence.png` and prints one number per curve (T2 = the last time |L| ≥
1/e). A few `run.sh` also run a variant of the input written by `example/tools/InputKeys/set_input_keys.py` into
`results/`. Inputs name their files relative to the calculation folder (`../bath/`, `../../00_quickstart/bath/`), so
`run.sh` is run from there. Git ignores `results/`.

Shared here: `plot_coherence.py` (the figures and numbers) and `make_natab_baths.sh` (natural-abundance ¹³C
configurations with `tools/BathGenerator`, for `01` and `09`). Each `expected.png` was made by that folder's `run.sh`
(i7-14700K, WSL, 4 ranks), except the three of `08_GPU` and `05_multiNV/gCCE2_HahnRamsey_NV6`, which need a server
and will be added.

## Where the earlier folders went

The calculations were reorganised in learning order for 2.0.0; `git mv` kept their history. The trimmed inputs
give byte-identical `_wiDiv` to the old ones (28 inputs compared, same keys set where `run.sh` used to set them). The
idm `job` files were dropped (they ran `run.sh` with `NCONF` and `NSTATE`, which `configuration.k_range` and
`configuration.s_range` in the inputs replace); the old tree is `git show c3c862972:example/tutorial/...`.

| Before | Now | Note |
|---|---|---|
| `01_Arbitrary_system/` | `06_cCCE_vs_gCCE/gCCE1_Hahn_handwritten/` | its `qubit` file dropped: the input uses the `Qubit` block and never read it |
| `02_Diamond_NV_natab/bath/` | `00_quickstart/bath/` | the shared ¹³C bath (00, 01, 02, 05, 06, 07, 10; `genbath/` for 01 and 09) |
| `02_Diamond_NV_natab/ensemble_CCE/` | `01_ensemble_and_states/cCCE2_Hahn_ensemble/` | `k_range` over four generated configurations; its input on the committed configuration is `00_quickstart` |
| `02_Diamond_NV_natab/single-sample_CCE/` | `01_ensemble_and_states/cCCE2_Hahn_states/` | compared with the ensemble (no `s_range`) |
| `02_Diamond_NV_natab/builtin_isotopes/` | `01_ensemble_and_states/cCCE2_Hahn_builtinIsotopes/` | |
| `02_Diamond_NV_natab/ensemble_Dsj/` | `06_cCCE_vs_gCCE/cCCE2_Hahn_dsj/` | `ccein_cce.json` is `ccein_standard.jsonc` |
| `02_Diamond_NV_natab/knight_shift/` | `06_cCCE_vs_gCCE/cCCE2_Hahn_knightShift/` | |
| `02_Diamond_NV_natab/ensemble_gCCE/` | `05_multiNV/gCCE2_Hahn_oneVsTwoNV/` | |
| (new) | `02_pulses/cCCE2_CPMG_ensemble/`, `02_pulses/cCCE2_Manual_ensemble/`, `06_cCCE_vs_gCCE/gCCE2_Hahn_vsCCE/` | Hahn/CPMG-2/CPMG-8; a Manual CPMG-2 against the built-in one; cCCE2 against gCCE2 |
| `03_Diamond_multiNV_natab/NV2_13C2/` | `05_multiNV/gCCE2_HahnRamsey_NV2/` | inputs `ccein_hahn.jsonc`, `ccein_ramsey.jsonc`; the Hahn input order 1 → 2 (as the Ramsey one) |
| `03_Diamond_multiNV_natab/NV6_13C2/` | `05_multiNV/gCCE2_HahnRamsey_NV6/` | reads NV2's bath, gyro and state files (identical copies dropped); `nstep` 300 → 50 (out of memory at 300) |
| `03_Diamond_multiNV_natab/bath/` | dropped | read by nothing (a copy of `genbath/` and one configuration) |
| `04_Diamond_NV_P1Bath/bath/` | `03_P1/bath/` | 1 and 10 ppm × 5 configurations kept (99 MB → 0.9 MB); the rest is in git history |
| `04_Diamond_NV_P1Bath/SpinEcho_CCE2_1ppm/` | `03_P1/cCCE2_Hahn_states/` | `k_range` [1, 2] × `s_range` [1, 10] (recorded: 50 × 50 on idm) |
| `04_Diamond_NV_P1Bath/SpinEcho_CCE4_1ppm/` | `03_P1/cCCE4_Hahn_states/` | as above |
| `04_Diamond_NV_P1Bath/SpinEcho_CCE2e2n_1ppm/` | `03_P1/cCCE2_Hahn_onsite14N/` | as above |
| `04_Diamond_NV_P1Bath/SpinEcho_gCCE2_10ppm_ZFS/` | `03_P1/gCCE2_Hahn_zfs/` | `ccein_with_ZFS_v2.json` (the same system with the qubit's position from a file) and its `qubitfile` dropped |
| `04_Diamond_NV_P1Bath/Population_gCCE2_10ppm/` | `03_P1/gCCE2_Ramsey_population/` | its unread `qubitfile` dropped |
| `04_Diamond_NV_P1Bath/SpinEcho_pCCE2_1ppm/` | `07_pCCE_and_bundle/pCCE2_Hahn_states/` | four inputs (`partition_size` 1-4) → one, `run.sh` sets the size; one configuration × 4 states |
| `05_Diamond_NV_VpBath/` | `03_P1/cCCE2_Hahn_VpBath/`, `03_P1/bath/Vp/` | `bath_Diamond_NV_eleBath_1` and its `_defect` dropped (read by nothing) |
| `06_Pulse_Manual/HahnEcho_natab_{ensemble,single}_{origin,wahuha}/` | dropped | plain Hahn echoes (`pulse.count` 1, no sequence): `00_quickstart` and `01/cCCE2_Hahn_states` on another configuration; the wahuha copies were byte-identical |
| `06_Pulse_Manual/CPMG_P1_1ppm/`, `SpinEcho_P1_1ppm/` | dropped | CPMG-8 and Hahn by `pulse.count` on the P1 bath, three of the five inputs byte-identical copies; `02_pulses` shows the pulses, `03_P1` the bath |
| `06_Pulse_Manual/WAHUHA_P1_1ppm/` | dropped | gCCE with the WAHUHA sequence written out on P1; the `Manual` form is `02_pulses/cCCE2_Manual_ensemble` |
| `06_Pulse_Manual/WAHUHA_natab_single_ccce/` | dropped | refused by the code (cCCE with pulses that are not π flips) |
| `06_Pulse_Manual/WAHUHA_natab_single_gcce/` | dropped | its `_wiDiv` diverges (\|L\| about 6e172 at 2 ms on a rerun; `check_numerics` leaves WAHUHA out for the same reason) |
| `06_Pulse_Manual/bath/` | dropped | a second ¹³C configuration, read only by the dropped runs |
| `07_hBN_VB_natab/bath/` | `04_hBN/bath/` | |
| `07_hBN_VB_natab/ensemble_CCE/` | `04_hBN/cCCE2_Hahn_ensemble/` | `k_range` [1, 5] |
| `07_hBN_VB_natab/hyperfine_mediated_term/` | `04_hBN/cCCE2_Hahn_hfMediated/` | |
| `07_hBN_VB_natab/single-sample_gCCE/` | `04_hBN/gCCE2_Hahn_states/` | `seed` 1 added (the state was drawn from the clock) |
| `07_hBN_VB_natab/single-sample_gCCE_DPS/` | `04_hBN/gCCE2_Hahn_OH/` | |
| `08_Bundle_clustering/CCE2_bundle/` | `07_pCCE_and_bundle/cCCE2_Hahn_bundle/` | |
| `08_Bundle_clustering/CCE2_adaptive/` | `07_pCCE_and_bundle/cCCE2_Hahn_adaptive/` | |
| `08_Bundle_clustering/CCE2_regroup/` | `07_pCCE_and_bundle/cCCE2_Hahn_regroup/` | `bath/bath_p1_6` → `07_pCCE_and_bundle/bath/`; the copies of the ¹³C files dropped |
| `09_GPU_backend/CCE5_P1/`, `gCCE_P1/`, `pCCE_P1/` | `08_GPU/cCCE5_Hahn_backends/`, `gCCE4_Hahn_backends/`, `pCCE2_Hahn_backends/` | three inputs per case (differing in `backend`) → one |
| `10_Coordinate_rotation/HahnEcho_CCE2/` | `10_coordinate_rotation/cCCE2_Hahn_rotation/` | the copies of the ¹³C files dropped (`make_cubic_frame.py` reads `00_quickstart/bath/`) |
| `11_Auto_convergence/ensemble_CCE2/` | `09_auto_convergence/cCCE2_Hahn_ensemble/` | |
