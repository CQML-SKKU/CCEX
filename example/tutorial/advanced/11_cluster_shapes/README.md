# Advanced 11: Cluster shapes (dsj, pCCE and bundles)

Folder: `example/tutorial/advanced/11_cluster_shapes/`

Other ways of cutting the bath into clusters than the standard one. **dsj** (`cluster.mode` dsj, cCCE only) splits
the bath into disjoint groups of the coupling graph and computes each group as one cluster. Clusters can also be made
of groups instead of single spins. **pCCE** (`cluster.mode` pcce) groups the bath into partitions of
`partition_size` spins by k-means on the positions and runs the expansion over partitions (order 2 = pairs of
partitions). **Bundle clustering** (`cluster.mode` bundle) bundles strongly coupled spins,
from the coupling strengths and under the dimension cap `cluster.bundle.max_dim`, and every bundle enters a cluster whole;
`order` counts bundles. **Adaptive refinement** (the `cluster.adaptive` block,
[Cluster tags](../../../../docs/user/InputFile_Cluster.md#tag_adaptive)) repeats the bundled calculation, measures the error at the bundle boundaries
and re-bundles. The cluster reports (the contribution of each cluster, the plan) are in [Advanced 02: Analysis of the output](../02_analysis/README.md). Each folder: `cd <folder> && NP=4 bash
run.sh`, then compare `results/coherence.png` with `expected.png`. Every key not written in an input has its default
(`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

## cCCE2_Hahn_dsj

- **What it shows:** `cluster.mode` dsj (the clusters are the disjoint groups of the coupling graph within `rdip`)
  against the standard pairs (`ccein_standard.jsonc`).
- **Keys that matter:** `cluster.mode` dsj.
- **How to run:** `NP=4 bash run.sh`: 1.1 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.158 ms (dsj) against 1.084 ms; largest difference 0.065.
- **What to change next:** `cluster.mode` itb or dsjitb (cCCE only).

![cCCE2_Hahn_dsj: expected result](cCCE2_Hahn_dsj/expected.png)

Inputs: [`ccein.jsonc`](cCCE2_Hahn_dsj/ccein.jsonc), [`ccein_standard.jsonc`](cCCE2_Hahn_dsj/ccein_standard.jsonc)

## pCCE2_Hahn_states

- **What it shows:** pCCE2 of the NV in 1 ppm P1 ([Basic 06: An electron-spin bath (P1)](../../basic/06_defect_bath/README.md)'s bath, one configuration, four sampled states) with
  partitions of 1, 2 and 3 P1 (`partition_size` 1 is ordinary CCE2 over single spins).
- **Keys that matter:** `cluster.mode` pcce, `cluster.pcce.partition_size` (2 in the input; `run.sh` sets 1, 2, 3 or
  `SIZES`), `trials`, `max_iterations`.
- **How to run:** `NP=4 bash run.sh`: 8 s, 0.06 GB per rank. `SIZES="1 2 3 4" bash run.sh` adds 4 (much longer).
- **Expected result:** `expected.png`; T2 260 µs (size 1), 183 µs (2), 202 µs (3) for these four states.
- **What to change next:** more states (`s_range`) before reading the differences; `cluster.pcce.connections` spins.

![pCCE2_Hahn_states: expected result](pCCE2_Hahn_states/expected.png)

Input: [`ccein.jsonc`](pCCE2_Hahn_states/ccein.jsonc)

## cCCE2_Hahn_bundle

- **What it shows:** basic/01_first_run's NV (`rbath` 30) with ordinary clusters (`ccein_standard.jsonc`) and with bundles of
  up to three ¹³C (`bundle.max_dim` 8, `ccein_bundle.jsonc`); the `[adaptive]` log line says how the 199 spins were
  bundled (84 bundles).
- **Keys that matter:** `cluster.mode` bundle, `cluster.bundle.max_dim` 8.
- **How to run:** `NP=4 bash run.sh`: 2.7 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.165 ms (standard) and 1.169 ms (bundles); largest difference 0.015.
- **What to change next:** `max_dim` 4 (pairs) or 16.

![cCCE2_Hahn_bundle: expected result](cCCE2_Hahn_bundle/expected.png)

Inputs: [`ccein_bundle.jsonc`](cCCE2_Hahn_bundle/ccein_bundle.jsonc), [`ccein_standard.jsonc`](cCCE2_Hahn_bundle/ccein_standard.jsonc)

## cCCE2_Hahn_adaptive

- **What it shows:** the bundled calculation (`rbath` 20, 2 ms) with adaptive refinement, merge policy, `atol` 1e-3:
  two rounds (27 → 26 bundles), then "no feasible improving merge" (exit status 2); `results/adaptive/` holds one
  folder per step and `summary.json`. `adaptive` needs an explicit `seed`: every round recomputes the same system.
- **Keys that matter:** `cluster.adaptive` {`workdir`, `atol`, `rtol`, `max_rounds`}, `configuration.seed`.
- **How to run:** `NP=4 bash run.sh`: 48 s, 0.12 GB per process.
- **Expected result:** `expected.png` (first and final run); |L| 0.38 at 2 ms in both, largest difference 1.6e-5 (the
  merges hardly move this curve).
- **What to change next:** `atol` 1e-4; `max_rounds` 8.

![cCCE2_Hahn_adaptive: expected result](cCCE2_Hahn_adaptive/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_adaptive/ccein.jsonc)

## cCCE2_Hahn_regroup

- **What it shows:** the NV in six P1 electron spins (`bath/bath_p1_6`, five sampled states, bundles up to dimension 4)
  with `"regroup": true`: the bundles [0,3] [1,5] [2,4] give three states with |L| > 1.01; the culprit pairs are found
  and one change, splitting [2,4], is accepted; the rerun has no unphysical state (status `physical`, exit 0). The NV
  position and the P1 gyromagnetic ratio are [Advanced 04: The GPU backend](../04_GPU/README.md)'s.
- **Keys that matter:** `cluster.adaptive.regroup` true, `cluster.bundle.max_dim` 4, `s_range` [1, 5], `seed`.
- **How to run:** `NP=4 bash run.sh`: 1 s, 0.06 GB per rank.
- **Expected result:** `expected.png` (state 4 before and after); the printed rounds: culprit pairs [0,3]-[1,5]
  (g 1.63), [0,3]-[2,4] (0.10), [1,5]-[2,4] (5e-4), final bundles [[0,3], [1,5], [2], [4]].
- **What to change next:** `max_dim` 8; more states.

![cCCE2_Hahn_regroup: expected result](cCCE2_Hahn_regroup/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_regroup/ccein.jsonc)

