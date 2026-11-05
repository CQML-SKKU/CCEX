# Advanced 12: Dissipation (meCCE and megCCE)

Folder: `example/tutorial/advanced/12_dissipation/`

The master-equation CCE adds Lindblad dissipation to every cluster. meCCE is the cCCE with jump operators on the bath
spins, and megCCE the gCCE with jump operators on the bath spins and on the qubit. How they work is in
[Master-equation CCE](../../../../docs/methods/master-equation.md).

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

1. [meCCE2_Hahn_bathRelaxation](#ds_bath): the relaxation of the bath spins against cCCE.

<a id="ds_bath"></a>
## 1. meCCE2_Hahn_bathRelaxation: the bath spins relax

- **What it shows:** `method` meCCE, the cCCE with a Lindblad master equation in each cluster, on the system of
  Basic 01 with 500 points of 4 µs. Jump operators S₋ and S₊ on every ¹³C, at γ = 0.01 rad/ms each (`rate`
  1.5915e-6 MHz × 2π·10³), relax the bath spins (`ccein_meCCE.jsonc`), compared with cCCE (`ccein_cCCE.jsonc`). The
  rates of all the ¹³C add up, so even this small rate halves T2. Without jump operators meCCE gives the cCCE curve.
- **Keys that matter:** `method` meCCE and `hamiltonian.meCCE.jump_operators` [{`bath_name` "13C", `operator` "-" or
  "+", `rate` (MHz)}].
- **How to run:** `NP=4 bash run.sh`. It takes 18–20 s, nearly all of it meCCE, and 0.06 GB per rank.
- **Expected result:** `expected.png`. T2 is 1.08 ms with cCCE and 0.57 ms with meCCE.
- **What to change next:** `operator` "z" (pure dephasing of the ¹³C), or `method` megCCE with a qubit jump
  (`bath_name` "qubit", in `hamiltonian.megCCE`).

![meCCE2_Hahn_bathRelaxation: expected result](meCCE2_Hahn_bathRelaxation/expected.png)

Inputs: [`ccein_cCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_cCCE.jsonc), [`ccein_meCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_meCCE.jsonc)
