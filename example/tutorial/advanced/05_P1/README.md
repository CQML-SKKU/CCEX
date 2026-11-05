# Advanced 05: The P1 bath further

Folder: `example/tutorial/advanced/05_P1/`

More on the P1 bath of [Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md), whose `bath/` every input here
reads: the ¹⁴N of each P1 as a nuclear spin, the expansion to order 4, gCCE at zero field (where the NV zero-field splitting is needed), the qubit density
matrix, a Vp bath at three densities and a DEER pump on the bath. Each folder: `cd <folder> && NP=4 bash run.sh`, then
compare `results/coherence.png` with `expected.png`. Every key not written in an input has its default (`main.out -e`).
Times: i7-14700K, WSL, 4 ranks.

## cCCE2_Hahn_onsite14N

- **What it shows:** each P1 written with its ¹⁴N as an on-site nuclear spin (`n_onsite_spins` 1, four Jahn-Teller
  axes with the ¹⁴N hyperfine and quadrupole tensors), so a cluster of one P1 is its electron and its ¹⁴N.
- **Keys that matter:** the `Defect` species: `n_onsite_spins`, `types`, `spins`, `gyros`, `quadrupole_moments`,
  `n_axis_states` 4, `hyperfine`, `efg`.
- **How to run:** `NP=4 bash run.sh`: 65 s, 0.07 GB per rank.
- **Expected result:** `expected.png`; T2 of the mean 266 µs (against 307 µs with the detunings of `basic/06_defect_bath/cCCE2_Hahn_states`).
- **What to change next:** the on-site states from a file (`configuration.defect_state_file`).

![cCCE2_Hahn_onsite14N: expected result](cCCE2_Hahn_onsite14N/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_onsite14N/ccein.jsonc)

## cCCE4_Hahn_states

- **What it shows:** the same to order 4, the convergence check of the expansion; `cluster.selection.max_per_order`
  keeps the 5000 strongest three-spin and 10000 strongest four-spin clusters.
- **Keys that matter:** `order` 4, `cluster.selection.max_per_order` [0, 0, 0, 5000, 10000].
- **How to run:** `NP=2 bash run.sh` to stay within 2 GB: 0.84 GB per rank (the cluster list); 34 s on 2 ranks.
- **Expected result:** `expected.png`; T2 of the mean 307 µs, the same as order 2 (converged at order 2 here).
- **What to change next:** `max_per_order` 0 (all clusters; much more memory).

![cCCE4_Hahn_states: expected result](cCCE4_Hahn_states/expected.png)

Input: [`ccein.jsonc`](cCCE4_Hahn_states/ccein.jsonc)

## gCCE2_Hahn_zfs

- **What it shows:** at zero field the NV must be a spin 1 with its zero-field splitting (D = 2.88 GHz, the `intmap`
  tensor in kHz) in gCCE; without it the ms ±1 levels are degenerate and the divided curve (`_wiDiv`) diverges.
- **Keys that matter:** `method` gCCE, `bfield` 0, the `Qubit.qubit` entry and its `intmap` self tensor (present in
  `ccein_with_ZFS.jsonc`, absent in `ccein_without_ZFS.jsonc`), `s_range` [1, 1], `seed` 8810, 10 ppm, `deltat` 50 ns.
- **How to run:** `NP=4 bash run.sh`: 1.5 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 with the ZFS 3.59 µs; without it |L| reaches 3.8e15 (cut at 1.2 in the
  plot).
- **What to change next:** a small field (`bfield` 10) without the ZFS; more states.

![gCCE2_Hahn_zfs: expected result](gCCE2_Hahn_zfs/expected.png)

Inputs: [`ccein_with_ZFS.jsonc`](gCCE2_Hahn_zfs/ccein_with_ZFS.jsonc), [`ccein_without_ZFS.jsonc`](gCCE2_Hahn_zfs/ccein_without_ZFS.jsonc)

## gCCE2_Ramsey_population

- **What it shows:** the NV prepared in ms = 0 (`alpha_state` = `beta_state` = `{"ms": [0]}`), no pulse, the qubit density matrix
  (`quantity` dm) of gCCE2 over 50 ms in 10 ppm P1 at 300 G: the ms = 0 population, which flip-flops with the bath
  would lower.
- **Keys that matter:** `quantity` dm, `alpha_state` / `beta_state` `{"ms": [0]}`, `pulse.count` 0, `bfield` 300, `deltat` 10 µs ×
  `nstep` 5000.
- **How to run:** `NP=4 bash run.sh`: 76 s, 0.06 GB per rank.
- **Expected result:** `expected.png` (the y range from the data); the population stays within 0.997 to 1.003.
- **What to change next:** a scan of `bfield`: where an NV transition comes into resonance with the bath, the
  population drops.

![gCCE2_Ramsey_population: expected result](gCCE2_Ramsey_population/expected.png)

Input: [`ccein.jsonc`](gCCE2_Ramsey_population/ccein.jsonc)

## cCCE2_Hahn_VpBath

- **What it shows:** the Hahn echo in a bath of Vp defects (spin 1/2 without on-site nuclei: no `Defect` block) at
  10, 12 and 50 ppm, one configuration and one sampled state each.
- **Keys that matter:** `configuration.gyrofile` (`Vp/DiaVp_gyro`), `bathfile`, `s_range` [1, 1], `seed`, `rbath`, `rdip`
  and `deltat` per concentration.
- **How to run:** `NP=4 bash run.sh`: 1.7 s, 0.06 GB per rank.
- **Expected result:** `expected.png`, the time axis on a log scale (`plot_coherence.py --log-time`) since the T2
  differ tenfold: 5.2 µs (10 ppm), 13.5 µs (12 ppm, a longer window), 2.3 µs (50 ppm, a 10 ns step) for this one state.
- **What to change next:** more states (`s_range` [1, 20]) before comparing concentrations.

![cCCE2_Hahn_VpBath: expected result](cCCE2_Hahn_VpBath/expected.png)

Inputs: [`ccein_10ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_10ppm.jsonc), [`ccein_12ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_12ppm.jsonc), [`ccein_50ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_50ppm.jsonc)

## gCCE2_Hahn_DEER

- **What it shows:** a DEER pump (`pulse.bath`): a π pulse on the P1 spins at the time of the NV's π pulse
  (`ccein_pump.jsonc`) against the plain Hahn echo (`ccein_nopump.jsonc`), gCCE2, NV ms +1 / 0, 1 ppm, ten sampled
  states. The tone `energy_shift` 0 MHz drives only the P1 whose `detuning` is within `bandwidth` of it (the central
  ¹⁴N line, axis states 4, 5, 6 and 11); their coupling to the NV is no longer refocused, so the echo decays faster.
  (At 500 G the NV's ms −1 / 0 transition is close to the P1 lines and its gCCE2 curve diverges, hence ms +1.)
- **Keys that matter:** `method` gCCE, `pulse.bath` {`count` 1, `defect` "P1", `energy_shift` [0.0], `sequence`
  [["1/2", "X", 180], ["1/2", "I", 0]]}, the `Defect` detunings, `s_range` [1, 10], `seed`.
- **How to run:** `NP=4 bash run.sh`: 9 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 of the mean of the 10 states 150 µs (Hahn echo) and 36 µs (with the pump).
- **What to change next:** another tone (`energy_shift` [84.946496] drives the axis states 1-3), or the pump at another
  time (`pulse.bath.sequence` [["0.3", "X", 180], ["0.7", "I", 0]]); see [DEER bath pulses](../../../../docs/methods/deer.md).

![gCCE2_Hahn_DEER: expected result](gCCE2_Hahn_DEER/expected.png)

Inputs: [`ccein_nopump.jsonc`](gCCE2_Hahn_DEER/ccein_nopump.jsonc), [`ccein_pump.jsonc`](gCCE2_Hahn_DEER/ccein_pump.jsonc)

