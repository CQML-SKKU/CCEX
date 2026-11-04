# basic/02_bath_and_states

The NV of [basic/01_first_run](../01_first_run/README.md) averaged two ways: over bath configurations (positions of the ¹³C), with
`configuration.k_range`, and over sampled bath states (`configuration.s_range`) against the ensemble state. Each folder:
`cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every key not written in
an input has its default (`main.out -e` writes every key). Times: i7-14700K, WSL, 4 ranks.

## cCCE2_Hahn_ensemble

- **What it shows:** the coherence depends on where the ¹³C sit; a measured ensemble averages over configurations.
  `k_range` [1, 4] runs the input once per k with `{k}` replaced in the file names (`ccein_k1_wiDiv` ... `_k4_`).
- **Keys that matter:** `configuration.bathfile` `["results/bath/config_{k}"]` with `configuration.k_range` [1, 4];
  the rest as `basic/01_first_run`.
- **How to run:** `NP=4 bash run.sh` makes the four configurations first (`../../../make_natab_baths.sh`, 110 Å cells,
  seed k, about 5 s with NP=4) and then runs: 7 s in all, 0.06 GB per rank.
- **Expected result:** `expected.png` (the four curves and their mean); T2 of the mean 0.92 ms (per configuration
  0.88 to 1.06 ms).
- **What to change next:** more configurations (`k_range` [1, 20] and `make_natab_baths.sh 20`); how many are enough is
  what [advanced/12_auto_convergence](../../advanced/12_auto_convergence/README.md) decides.

![cCCE2_Hahn_ensemble: expected result](cCCE2_Hahn_ensemble/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_ensemble/ccein.jsonc)

## cCCE2_Hahn_states

- **What it shows:** `s_range` [1, 20] replaces the unpolarised bath by 20 product states drawn from `seed`
  (`ccein_state1_wiDiv` ...); their mean approaches the ensemble curve (no `s_range`), while one state alone does not
  decay the same way.
- **Keys that matter:** `configuration.s_range` [1, 20], `configuration.seed` 1 (the states are drawn from it); `run.sh` runs
  the same input without `s_range` (the ensemble) as `results/ensemble.json` (`tools/InputKeys/set_input_keys.py`), and
  `ccein_files.jsonc`: the 20 states read from files (`bath_state_file` `["results/states/state{s}"]`) that
  [`tools/StateFiles/make_state_files.py`](../../../../docs/user/Tools.md#tag_tool_statefiles) writes first
  (`results/coherence_files.png`).
- **How to run:** `NP=4 bash run.sh`: 5 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 ensemble 1.068 ms, mean of 20 states 1.09 ms (largest |L| difference 0.09),
  state 1 alone 1.50 ms.
- **What to change next:** `s_range` [1, 100] (the difference from the ensemble shrinks); `bath_state_file` to give the states
  from files instead of `seed`.

![cCCE2_Hahn_states: expected result](cCCE2_Hahn_states/expected.png)

Inputs: [`ccein.jsonc`](cCCE2_Hahn_states/ccein.jsonc), [`ccein_files.jsonc`](cCCE2_Hahn_states/ccein_files.jsonc)

## cCCE2_Hahn_builtinIsotopes

- **What it shows:** without `configuration.gyrofile` the bath species are taken from the built-in isotope table
  (EasySpin, 154 isotopes): ¹³C γ 6.72828544 against the file's 6.728284 rad/ms/G. The log says where each species came
  from.
- **Keys that matter:** `configuration.gyrofile` present (`ccein_gyrofile.jsonc`) or absent (`ccein_builtin.jsonc`).
- **How to run:** `NP=4 bash run.sh`: 1.3 s, 0.05 GB per rank.
- **Expected result:** `expected.png` (the two curves on top of each other); largest difference 1.1e-4 (the 2e-7
  relative difference of γ over 2 ms at 500 G); with the table's value written in a gyro file the outputs are
  byte-identical.
- **What to change next:** a bath of another isotope (²⁹Si, ¹⁵N) without a gyro file; species that are not isotopes
  (P1, Vp) still need one ([basic/05_defect_bath](../05_defect_bath/README.md)).

![cCCE2_Hahn_builtinIsotopes: expected result](cCCE2_Hahn_builtinIsotopes/expected.png)

Inputs: [`ccein_builtin.jsonc`](cCCE2_Hahn_builtinIsotopes/ccein_builtin.jsonc), [`ccein_gyrofile.jsonc`](cCCE2_Hahn_builtinIsotopes/ccein_gyrofile.jsonc)
