# Advanced 05: The P1 bath further

Folder: `example/tutorial/advanced/05_P1/`

This topic goes further with the P1 bath of
[Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md), and every input here reads its `bath/`.

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

This page has six cases:

1. [The ¹⁴N of each P1 as a nuclear spin](#p1_onsite) (`cCCE2_Hahn_onsite14N/`).
2. [The expansion to order 4](#p1_order4) (`cCCE4_Hahn_states/`).
3. [gCCE at zero field, where the NV zero-field splitting is needed](#p1_zfs) (`gCCE2_Hahn_zfs/`).
4. [The qubit density matrix](#p1_population) (`gCCE2_Ramsey_population/`).
5. [A Vp bath at three densities](#p1_vp) (`cCCE2_Hahn_VpBath/`).
6. [A DEER pump on the bath](#p1_deer) (`gCCE2_Hahn_DEER/`).

<a id="p1_onsite"></a>
## 1. The ¹⁴N as a spin

Folder: `cCCE2_Hahn_onsite14N/`

<a id="a05_onsite_input"></a>
### 1.1 The input

Each P1 is written with its ¹⁴N as an on-site nuclear spin, so a cluster of one P1 is its electron and its ¹⁴N. The
`Defect` species gives the ¹⁴N and, for each of the four Jahn-Teller axes, its hyperfine and quadrupole tensors:

```jsonc
"Defect" : {"species": [{"name": "P1", "n_onsite_spins": 1, "types": ["14N"], "spins": [1],
                         "quadrupole_moments": [2.044], "gyros": [1.9337792], "n_axis_states": 4,
                         "hyperfine": [[1, "14N", [81.312, 0.0, 0.0, 0.0, 81.312, 0.0, 0.0, 0.0, 114.0264]], ...],
                         "efg":       [[1, "14N", [0.52191, 0.0, 0.0, 0.0, 0.522771, -8.6e-05, 0.0, -8.6e-05, -1.044681]], ...]}]}
```

The bath is two configurations (`k_range` [1, 2]) of ten sampled states each (`s_range` [1, 10]).

<a id="a05_onsite_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

It takes 65 s and 0.07 GB per rank.

<a id="a05_onsite_result"></a>
### 1.3 The result

![cCCE2_Hahn_onsite14N: expected result](cCCE2_Hahn_onsite14N/expected.png)

The T2 of the mean is 266 µs, against 307 µs with the detunings of `basic/06_defect_bath/cCCE2_Hahn_states`.

Next, the on-site states can be read from a file (`configuration.defect_state_file`).

Input: [`ccein.jsonc`](cCCE2_Hahn_onsite14N/ccein.jsonc)

<a id="p1_order4"></a>
## 2. Order 4

Folder: `cCCE4_Hahn_states/`

<a id="a05_order4_input"></a>
### 2.1 The input

The bath of Basic 06's `cCCE2_Hahn_states` is computed to order 4, to check the convergence of the expansion.
`cluster.selection.max_per_order` keeps the 5000 strongest three-spin and 10000 strongest four-spin clusters:

```jsonc
"order"   : 4,
"cluster" : {"selection": {"max_per_order": [0, 0, 0, 5000, 10000]}}
```

<a id="a05_order4_run"></a>
### 2.2 Running it

```sh
NP=2 bash run.sh
```

Two ranks keep the run within 2 GB. The cluster list takes 0.84 GB per rank, and the run 34 s on 2 ranks.

<a id="a05_order4_result"></a>
### 2.3 The result

![cCCE4_Hahn_states: expected result](cCCE4_Hahn_states/expected.png)

The T2 of the mean is 307 µs, the same as at order 2, so the expansion has converged at order 2 here.

Next, `max_per_order` 0 keeps all clusters and needs much more memory.

Input: [`ccein.jsonc`](cCCE4_Hahn_states/ccein.jsonc)

<a id="p1_zfs"></a>
## 3. Zero field needs the zero-field splitting

Folder: `gCCE2_Hahn_zfs/`

<a id="a05_zfs_input"></a>
### 3.1 The input

At zero field the NV must be a spin 1 with its zero-field splitting in gCCE, D = 2.88 GHz, given as the `intmap`
self tensor in kHz. `ccein_with_ZFS.jsonc` has it, and `ccein_without_ZFS.jsonc` is the same input without it:

```jsonc
"method" : "gCCE",
"Qubit"  : {"qubit": [{"name": "NV", "xyz": [523.336, 524.936, 521.003]}],
            "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]},
            "intmap": [{"between": ["NV", "NV"], "tensor": [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]}]},
"bfield" : 0
```

Both read one sampled state (`s_range` [1, 1], `seed` 8810) of the 10 ppm bath, with `deltat` 50 ns.

<a id="a05_zfs_run"></a>
### 3.2 Running it

```sh
NP=4 bash run.sh
```

It takes 1.5 s and 0.06 GB per rank.

<a id="a05_zfs_result"></a>
### 3.3 The result

![gCCE2_Hahn_zfs: expected result](gCCE2_Hahn_zfs/expected.png)

The T2 with the ZFS is 3.59 µs. Without it the ms ±1 levels are degenerate and the divided curve (`_wiDiv`) diverges,
with |L| up to 3.8e15 (cut at 1.2 in the plot).

Next, try a small field (`bfield` 10) without the ZFS, or more states.

Inputs: [`ccein_with_ZFS.jsonc`](gCCE2_Hahn_zfs/ccein_with_ZFS.jsonc), [`ccein_without_ZFS.jsonc`](gCCE2_Hahn_zfs/ccein_without_ZFS.jsonc)

<a id="p1_population"></a>
## 4. The density matrix

Folder: `gCCE2_Ramsey_population/`

<a id="a05_population_input"></a>
### 4.1 The input

The NV is prepared in ms = 0 and has no pulse, and gCCE2 writes the qubit density matrix over 50 ms in 10 ppm P1 at
300 G. Its ms = 0 population is what flip-flops with the bath would lower:

```jsonc
"method"   : "gCCE",
"quantity" : "dm",
"Qubit"    : {"alpha_state": {"ms": [0]}, "beta_state": {"ms": [0]}, ...},
"bfield"   : 300,
"deltat"   : 0.01,
"nstep"    : 5000,
"pulse"    : {"count": 0}
```

<a id="a05_population_run"></a>
### 4.2 Running it

```sh
NP=4 bash run.sh
```

It takes 76 s and 0.06 GB per rank.

<a id="a05_population_result"></a>
### 4.3 The result

![gCCE2_Ramsey_population: expected result](gCCE2_Ramsey_population/expected.png)

The y axis spans only the range of the data. The population stays within 0.997 to 1.003.

Next, a scan of `bfield` shows where an NV transition comes into resonance with the bath, since the population drops
there.

Input: [`ccein.jsonc`](gCCE2_Ramsey_population/ccein.jsonc)

<a id="p1_vp"></a>
## 5. A Vp bath

Folder: `cCCE2_Hahn_VpBath/`

<a id="a05_vp_input"></a>
### 5.1 The input

This case computes the Hahn echo in a bath of Vp defects at 10, 12 and 50 ppm, one configuration and one sampled
state each. A Vp is a spin 1/2 without on-site nuclei, so the input has no `Defect` block, and the gyro file gives it:

```jsonc
"configuration" : {"gyrofile": "../../../basic/06_defect_bath/bath/Vp/DiaVp_gyro",
                   "bathfile": ["../../../basic/06_defect_bath/bath/Vp/10ppm/bath_DiaVp_10ppm_1"],
                   "s_range": [1, 1], "seed": 1}
```

`rbath`, `rdip` and `deltat` differ per concentration:

| Input | `rbath` | `rdip` | `deltat` |
|---|---|---|---|
| `ccein_10ppm.jsonc` | 370 Å | 330 Å | 30 ns |
| `ccein_12ppm.jsonc` | 370 Å | 330 Å | 0.5 µs |
| `ccein_50ppm.jsonc` | 240 Å | 220 Å | 10 ns |

<a id="a05_vp_run"></a>
### 5.2 Running it

```sh
NP=4 bash run.sh
```

It takes 1.7 s and 0.06 GB per rank.

<a id="a05_vp_result"></a>
### 5.3 The result

![cCCE2_Hahn_VpBath: expected result](cCCE2_Hahn_VpBath/expected.png)

The time axis is on a log scale (`plot_coherence.py --log-time`), since the T2 differ tenfold. For this one state
they are 5.2 µs at 10 ppm, 13.5 µs at 12 ppm (a longer window) and 2.3 µs at 50 ppm (a 10 ns step).

Next, compute more states (`s_range` [1, 20]) before comparing concentrations.

Inputs: [`ccein_10ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_10ppm.jsonc), [`ccein_12ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_12ppm.jsonc), [`ccein_50ppm.jsonc`](cCCE2_Hahn_VpBath/ccein_50ppm.jsonc)

<a id="p1_deer"></a>
## 6. A DEER pump

Folder: `gCCE2_Hahn_DEER/`

<a id="a05_deer_input"></a>
### 6.1 The input

A DEER pump (`pulse.bath`) applies a π pulse to the P1 spins at the time of the NV's π pulse. `ccein_pump.jsonc` has
it, and `ccein_nopump.jsonc` is the plain Hahn echo. Both are gCCE2 with the NV in ms +1 / 0, 1 ppm and ten sampled
states (`s_range` [1, 10]):

```jsonc
"pulse" : {"count": 1, "bath": {"count": 1, "defect": "P1", "energy_shift": [0.0],
                                "sequence": [["1/2", "X", 180.0], ["1/2", "I", 0.0]]}}
```

The tone `energy_shift` 0 MHz drives only the P1 whose `detuning` is within `bandwidth` of it, the central ¹⁴N line
(axis states 4, 5, 6 and 11 of the `Defect` detunings). The NV is in ms +1 because at 500 G its ms −1 / 0 transition is
close to the P1 lines and its gCCE2 curve diverges.

<a id="a05_deer_run"></a>
### 6.2 Running it

```sh
NP=4 bash run.sh
```

It takes 9 s and 0.06 GB per rank.

<a id="a05_deer_result"></a>
### 6.3 The result

![gCCE2_Hahn_DEER: expected result](gCCE2_Hahn_DEER/expected.png)

The pumped spins' coupling to the NV is no longer refocused, so the echo decays faster. The T2 of the mean of the 10
states is 150 µs for the Hahn echo and 36 µs with the pump.

Next, try another tone (`energy_shift` [84.946496] drives the axis states 1–3), or the pump at another time
(`pulse.bath.sequence` [["0.3", "X", 180], ["0.7", "I", 0]]). [DEER bath pulses](../../../../docs/methods/deer.md)
explains the pulses.

Inputs: [`ccein_nopump.jsonc`](gCCE2_Hahn_DEER/ccein_nopump.jsonc), [`ccein_pump.jsonc`](gCCE2_Hahn_DEER/ccein_pump.jsonc)
