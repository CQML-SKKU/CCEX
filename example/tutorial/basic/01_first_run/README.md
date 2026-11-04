# basic/01_first_run

    bash do_compile.sh wsl-gnu                 # 1. build (from the repository root; or: cmake -B build && cmake --build build)
    cd example/tutorial/basic/01_first_run && NP=4 bash run.sh    # 2. run: results/ccein_wiDiv (about 1 s)
    #                                            3. plot: run.sh writes results/coherence.png and prints T2
    #                                            4. compare results/coherence.png with expected.png (T2 1.068 ms)

**What it shows.** The Hahn-echo coherence |L(t)| of an NV centre in diamond with natural-abundance ¹³C (1.1 %):
conventional CCE to order 2 (cCCE2), the bath in its unpolarised (ensemble) state. The fast ripple is the ¹³C
precession at 500 G (ESEEM); the slow decay is the ¹³C flip-flops.

**Keys that matter** (`ccein.jsonc`; everything else is the default). `main.out -e` writes an input template,
`ccex_example.jsonc`, with these settings and every other key a commented line; it runs once its `bathfile` names a
bath file (`bath/bath_Diamond_NV_natab_1` here). Keys: `Qubit.qubit[].file` (the
NV position; spin 1, ms +1 / 0 by default), `configuration.gyrofile` and `bathfile` (`bath/`: one 150 Å configuration
made with `bath/genbath/`), `order` 2, `bfield` 500 G, `rbath` 50 Å, `rdip` 10 Å, `deltat` 1 µs × `nstep` 2000,
`pulse.count` 1.

**How to run.** `NP=4 bash run.sh`: 0.8 s on 4 ranks (i7-14700K, WSL), 0.05 GB per rank. Writes `results/ccein_wiDiv`
(the coherence), `ccein_noDiv` (without the division by the sub-clusters, for the analyzer's error correction),
`ccein.log` and `coherence.png`.

**Expected result.** `expected.png`; T2 = 1.068 ms (the last time |L| ≥ 1/e, from `plot_coherence.py`).

![basic/01_first_run: expected result](expected.png)

Input: [`ccein.jsonc`](ccein.jsonc)

**What to change next.** `pulse.count` 2 or 8 (CPMG, [basic/03_pulses](../03_pulses/README.md)); `configuration.s_range` [1, 20] for sampled bath states
([basic/02_bath_and_states](../02_bath_and_states/README.md)); `rbath` and `rdip` larger to check convergence (or let [advanced/12_auto_convergence](../../advanced/12_auto_convergence/README.md) choose them);
`method` gCCE with the NV zero-field splitting ([basic/04_cCCE_vs_gCCE](../04_cCCE_vs_gCCE/README.md)).
