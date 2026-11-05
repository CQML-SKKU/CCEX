# Advanced 11: Dissipation (meCCE and megCCE)

Folder: `example/tutorial/advanced/11_dissipation/`

The master-equation CCE adds Lindblad dissipation to every cluster. meCCE is the cCCE with jump operators on the bath
spins, and megCCE the gCCE with jump operators on the bath spins and on the qubit. How they work is in
[Master-equation CCE](../../../../docs/methods/master-equation.md).

In each folder run `cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every
key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

This page has one case:

1. [The relaxation of the bath spins against cCCE](#ds_bath) (`meCCE2_Hahn_bathRelaxation/`).

<a id="ds_bath"></a>
## 1. The bath spins relax

Folder: `meCCE2_Hahn_bathRelaxation/`

<a id="a11_bath_input"></a>
### 1.1 The input

The system is that of Basic 01, with 500 points of 4 µs. `ccein_cCCE.jsonc` is the plain cCCE run, and
`ccein_meCCE.jsonc` adds the method and the jump operators:

```jsonc
"method"      : "meCCE",
"hamiltonian" : {"meCCE": {"jump_operators": [{"bath_name": "13C", "operator": "-", "rate": 1.5915e-06},
                                              {"bath_name": "13C", "operator": "+", "rate": 1.5915e-06}]}}
```

meCCE is the cCCE with a Lindblad master equation in each cluster. The jump operators S₋ and S₊ act on every ¹³C at
γ = 0.01 rad/ms each (`rate` 1.5915e-6 MHz × 2π·10³) and relax the bath spins. Without jump operators meCCE gives the
cCCE curve.

<a id="a11_bath_run"></a>
### 1.2 Running it

```sh
NP=4 bash run.sh
```

It runs both inputs and plots them. It takes 18–20 s, nearly all of it meCCE, and 0.06 GB per rank.

<a id="a11_bath_result"></a>
### 1.3 The result

![meCCE2_Hahn_bathRelaxation: expected result](meCCE2_Hahn_bathRelaxation/expected.png)

T2 is 1.08 ms with cCCE and 0.57 ms with meCCE. The rates of all the ¹³C add up, so even this small rate halves T2.

Next, try `operator` "z" (pure dephasing of the ¹³C), or `method` megCCE with a qubit jump (`bath_name` "qubit", in
`hamiltonian.megCCE`).

Inputs: [`ccein_cCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_cCCE.jsonc), [`ccein_meCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_meCCE.jsonc)
