# Advanced 01: Automatic convergence

Folder: `example/tutorial/advanced/01_auto_convergence/`

With an `auto_convergence` block, `mpirun -n N main.out -f ccein.jsonc` starts the runner instead of one calculation.
The runner (`example/tools/AutoConvergence/`, Python with numpy, scipy and matplotlib) runs CCEX itself, chooses the
time window, rdip, rbath and the number of bath configurations, and then computes the production. These are the
checks of [Basic 03: Checking convergence](../../basic/03_convergence/README.md), done by the runner. How it works is
in [Automatic convergence](../../../../docs/methods/auto-convergence.md), and its keys are in
[Automatic convergence tags](../../../../docs/user/InputFile_AutoConvergence.md). Every key not written in an input
has its default (`main.out -e`).

1. [cCCE2_Hahn_ensemble](#ac_given): the rdip and rbath candidates given.
2. [cCCE2_Hahn_ensemble_noCandidates](#ac_found): the candidates found by the runner.
3. [The two cases compared](#ac_compared): the values found, the times and the combination check.

<a id="ac_given"></a>
## 1. cCCE2_Hahn_ensemble: the candidates given

- **What it shows:** the NV of Basic 01, a cCCE2 Hahn echo of the ensemble. The rdip candidates 4–12 Å and the rbath
  candidates 30–45 Å are given, and the runner finds everything else. It finds the window (1.5 × the decay time),
  rdip, rbath and how many configurations the mean needs (SE(A)/A ≤ 5 %).
- **Keys that matter:** `auto_convergence` {`workdir`, `rdip.candidates`, `rbath.candidates`, `sampling`},
  `configuration.bathfile` ["results/bath/config_{k}"] (every existing k is used), `nstep` 101 (the points of every
  run) and `deltat` (only the first window).
- **How to run:** `NP=4 bash run.sh`. It first makes 40 configurations (`../../../make_natab_baths.sh`, 110 Å cells,
  enough for rbath up to 53 Å, and `NBATH` sets how many), then runs the runner. It takes about 5 min in all, most of it
  making the baths, and under 0.2 GB per process. Calling it again continues, keeping the baths and reusing the
  finished runs and decisions. The runner exits 0 when `converged` and 3 when `unresolved`, and the summary then says
  what is missing. `run.sh` prints the status and `results/auto/result/summary.txt`.
- **Expected result:** `expected.png` (the runner's `6_converged_coherence.png`). The status is converged, with the
  window 2.58 ms, rdip 8 Å, rbath 35 Å and 25 production configurations (σA/A 4.0 %), and T2 0.82 ms
  (stretched-exponential fit, n 1.37). `results/auto/result/` also has `1_window.png`, `2_rdip.png`, `3_rbath.png`,
  `5_configurations.png` and `coherence_mean.txt`.
- **What to change next:** raise `nstep`. The 101 points over 2.58 ms are 26 µs apart, coarser than the ¹³C
  precession at 500 G (about 2 µs), so the ESEEM ripple is aliased. The runner keeps `nstep` and does not check it.

![cCCE2_Hahn_ensemble: expected result](cCCE2_Hahn_ensemble/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_ensemble/ccein.jsonc)

<a id="ac_found"></a>
## 2. cCCE2_Hahn_ensemble_noCandidates: the candidates found

- **What it shows:** the same calculation with the `rdip` and `rbath` blocks left out. The runner counts the ¹³C of
  the test configurations and starts from 4 rbath candidates 1.25 apart, the largest where 300 spins are inside
  (33.4 Å), and from 7 rdip candidates 0.8 apart below 1.25² × the pair-phase estimate (8.2 Å). It then raises the
  largest of each until the comparisons agree (`example/tools/AutoConvergence/autoconv/autostart.py`).
- **Keys that matter:** `auto_convergence` {`workdir`, `sampling`} only. The rest is as in `cCCE2_Hahn_ensemble`.
- **How to run:** `NP=4 bash run.sh`. The baths are made in 140 Å cells (rbath up to 68 Å). Candidates the runner
  found itself need two agreeing comparisons in a row, so rbath converges only against 65.2 Å. With the 110 Å cells
  of `cCCE2_Hahn_ensemble` (53 Å) the run stops as `truncation-unresolved` and asks for baths covering 65.5 Å.
- **Expected result:** `expected.png` (the runner's `6_converged_coherence.png`). The status is converged.
- **What to change next:** compare the values it found with those of `cCCE2_Hahn_ensemble` ([section 3](#ac_compared)).

![cCCE2_Hahn_ensemble_noCandidates: expected result](cCCE2_Hahn_ensemble_noCandidates/expected.png)

Input: [`ccein.jsonc`](cCCE2_Hahn_ensemble_noCandidates/ccein.jsonc)

<a id="ac_compared"></a>
## 3. The two cases compared

Times are for 4 ranks of a desktop CPU. Most of the time is making the 40 baths, which is done once.

| | rdip | rbath | Window | Configurations | T2 (n) | Time in all (baths + runner) | Peak memory |
|---|---|---|---|---|---|---|---|
| `cCCE2_Hahn_ensemble` (candidates given) | 8 Å | 35 Å | 2.58 ms | 25 | 0.82 ms (1.37) | 4.8 min (3.9 + 0.8) | 0.18 GB |
| `cCCE2_Hahn_ensemble_noCandidates` | 8.24 Å | 41.7 Å | 2.98 ms | 30 | 0.81 ms (1.81) | 6.9 min (5.6 + 1.2) | 0.19 GB |

The T2 agree within 2 %, but the exponents do not (1.37 and 1.81). The bath configurations are not the same, since the
cells have another size.

The final combination check compares the chosen rbath and rdip with the largest computed on three production
configurations. It differs by 9 % in the area for `cCCE2_Hahn_ensemble` (rbath 35, rdip 8 against rbath 45, rdip 15)
and by 5.7 % for `cCCE2_Hahn_ensemble_noCandidates`, more than the 2 % threshold of each step. The runner reports this
as a warning in the summary and does not count it against convergence, so check it before using the result.
