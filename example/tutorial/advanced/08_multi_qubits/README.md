# Advanced 08: Multiple qubits

Folder: `example/tutorial/advanced/08_multi_qubits/`

A run can hold several qubits as one register, written as a list in `Qubit.qubit`, with their zero-field splittings
and couplings in `Qubit.intmap`. The examples here are NV centres in a shared ¹³C bath, computed with gCCE. Every
qubit is a spin 1 with its zero-field splitting (the `intmap` self tensors, in kHz), and the qubits couple through
their dipolar interaction from the positions.

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

1. [One NV against two in the same bath](#mq_one_two) (`gCCE2_Hahn_oneVsTwoNV/`).
2. [The density matrix of two NVs and a ¹³C pair](#mq_nv2) (`gCCE2_HahnRamsey_NV2/`).
3. [The same with a third NV](#mq_nv3) (`gCCE2_HahnRamsey_NV3/`).

<a id="mq_one_two"></a>
## 1. One NV against two

Folder: `gCCE2_Hahn_oneVsTwoNV/`

- **What it shows:** the Hahn echo of one NV against two NVs 25 Å apart along z, in the bath of Basic 01. The single
  NV (`ccein_single.jsonc`) takes its position from a file (`Qubit.qubit[].file`) and has its ZFS. The two NVs
  (`ccein_multi2.jsonc`) are the `Qubit.qubit` list, with the coherence between |+1,+1⟩ and |0,0⟩.
- **Keys that matter:** `method` gCCE, `Qubit.qubit` (two entries with `xyz`) and `Qubit.intmap` (one ZFS tensor
  each), or a `Qubit.qubit[].file` entry with its ZFS as a `Qubit.intmap` self entry.
- **How to run:** `NP=4 bash run.sh`. It takes 2 min (the two-NV run) and 0.07 GB per rank.
- **Expected result:** `expected.png`. T2 is 1.084 ms for one NV and 0.956 ms for two.
- **What to change next:** move the second NV further away (`xyz`), or choose another pair of register states with
  `alpha_state` and `beta_state`.

![gCCE2_Hahn_oneVsTwoNV: expected result](gCCE2_Hahn_oneVsTwoNV/expected.png)

Inputs: [`ccein_multi2.jsonc`](gCCE2_Hahn_oneVsTwoNV/ccein_multi2.jsonc), [`ccein_single.jsonc`](gCCE2_Hahn_oneVsTwoNV/ccein_single.jsonc)

<a id="mq_nv2"></a>
## 2. The density matrix of two NVs

Folder: `gCCE2_HahnRamsey_NV2/`

- **What it shows:** two NVs 50 Å apart and a pair of ¹³C near the first. The run writes the 2-qubit density matrix
  (`quantity` dm) for the bath state in the file `state1`, for the Hahn echo (`ccein_hahn.jsonc`) and for Ramsey
  (`ccein_ramsey.jsonc`, `pulse.count` 0), over 30 ns. The plot shows ρ between |+1,+1⟩ and |0,0⟩ divided by its value
  at t = 0. The bath, gyro and state files of this folder are written by hand (`bath` holds a count line, then x y z
  name).
- **Keys that matter:** `quantity` dm, `Qubit.qubit` and `intmap`, `configuration.bath_state_file` ["state"] with
  `s_range` [1, 1], `rbath` and `rdip` large enough to take every spin, and `deltat` 0.1 ns × `nstep` 300.
- **How to run:** `NP=4 bash run.sh`. It takes under 1 s and 0.05 GB per rank.
- **Expected result:** `expected.png`. |L| is 1.000 (Hahn) and 0.996 (Ramsey) at 30 ns.
- **What to change next:** a longer window (`deltat`), or order 1 to see what the pair cluster adds.

![gCCE2_HahnRamsey_NV2: expected result](gCCE2_HahnRamsey_NV2/expected.png)

Inputs: [`ccein_hahn.jsonc`](gCCE2_HahnRamsey_NV2/ccein_hahn.jsonc), [`ccein_ramsey.jsonc`](gCCE2_HahnRamsey_NV2/ccein_ramsey.jsonc)

<a id="mq_nv3"></a>
## 3. Three NVs

Folder: `gCCE2_HahnRamsey_NV3/`

- **What it shows:** the calculation of `gCCE2_HahnRamsey_NV2` with a third NV, on a window 100 times longer. The
  register has dimension 3³ = 27, and with the two ¹³C the clusters have dimension 108. The qubits' dipolar coupling
  (50 Å apart) makes L oscillate. The run reads the bath, gyro and state files of `../gCCE2_HahnRamsey_NV2/`.
- **Keys that matter:** three `Qubit.qubit` entries, three `intmap` self tensors (the zero-field splittings), and
  `deltat` 10 ns × `nstep` 300.
- **How to run:** `NP=4 bash run.sh`. It takes 1.2 s and writes 30 MB of output, since each file holds nstep × 27²
  complex numbers.
- **Expected result:** `expected.png`. The Ramsey curve falls to 0.41 at 0.46 µs and comes back to 1.0 at 0.92 µs,
  about three periods in the window. The Hahn echo removes the first oscillation and then follows it, ending at 0.85.
- **What to change next:** a fourth qubit, with one more `qubit` and `intmap` entry. The qubit space grows threefold
  with every NV, and the output ninefold.

![gCCE2_HahnRamsey_NV3: expected result](gCCE2_HahnRamsey_NV3/expected.png)

Inputs: [`ccein_hahn.jsonc`](gCCE2_HahnRamsey_NV3/ccein_hahn.jsonc), [`ccein_ramsey.jsonc`](gCCE2_HahnRamsey_NV3/ccein_ramsey.jsonc)
