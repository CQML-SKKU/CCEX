# check_numerics/: does this build still give the numbers an earlier one gave

Twenty-seven calculations, cut from the tutorials so that each runs in seconds
(`cases/<name>/ccein.json`, the bath, gyro, tensor, state and axis files they
read in `inputs/`), run with `bin/main.out` and compared against the outputs
recorded from earlier builds on the same machine (`reference/<version>/`).

    bash test/code_verification/check_numerics/run.sh [NP]     # needs bin/main.out (bash do_compile.sh wsl-gnu), mpirun, numpy

`run.sh` writes to `_work/`, which is ignored.

## How it decides

One table per reference. The verdict of a case is the largest |difference| over
its `_noDiv` and density-matrix files. The `_wiDiv` coherence is printed but not
judged: it divides by the free evolution and amplifies rounding wherever that is small.

- The exit status is 0 only when the newest `reference/head-*` (by name) shows no difference above 1e-9. That
  reference is this line of development's own, so a change that moves any number is seen. A case it does not have
  fails.
- Against the tagged versions the tolerance is 1e-4. `expected.json` names the cases that differ for a known reason,
  and the table reports them as "expected" instead of FAIL.
- The earlier `head-*` folders are kept, are printed with 1e-9 and do not decide. Against every reference except the
  deciding one, a case added after it is "new", not a failure.
- A case with no output, or with a reference file the run did not write, fails whatever `expected.json` says.

A binary built with MKL (`do_compile.sh` and CMake turn it on when they find it; the run log's header says
`MKL for large clusters : on`) is judged against the deciding `head-*` reference relative to it instead:
|run − reference| ≤ 1e-12 |reference| element by element, 1e-12 absolute where the reference is 0. The references
were recorded without MKL, and MKL changes the last digits of the clusters it takes (zgemm in the cCCE ensemble
Hahn echo's recursive path from bath dimension 21, zheevd from dimension 63: the bath for cCCE, qubit × bath for gCCE).
Held relative, a diverged element (|L| of order 1e5, where 1e-9 absolute is a 1e-14 change) is measured against its
own size, and every value below 1e3 is held tighter than by 1e-9. `compare.py` reads the run logs to know which rule applies.

## The cases

| name | system | method | approach | pulse | quantity |
|---|---|---|---|---|---|
| natab_cce2_hahn / _ramsey / _cpmg8 | NV + natural-abundance ¹³C, rbath 30 | CCE2 | ensemble (nstate 0) | Hahn / Ramsey (`pulse.count` 0) / CPMG-8 | coherence |
| natab_cce2_single | same | CCE2 | single (nstate 2, state files) | Hahn | coherence |
| natab_gcce2_hahn, natab_gcce2_bathadjust | same, one NV from the `Qubit.file` form; the second with the bath file given twice, the copy shifted 50 Å along z (`configuration.bath_shift`). `make_cases.py` cuts both from tutorial 02's `ensemble_CCE` input | gCCE2 | ensemble | Hahn | coherence |
| natab_dsj2 | same | dsj | ensemble | Hahn | coherence |
| natab_knight_shift | same, tilted field | CCE2 | ensemble | Hahn | coherence |
| natab_manual_hahn | same, the Hahn echo written as a `sequence` | CCE2 | ensemble | manual | coherence |
| p1_cce2_ensemble / _single | NV + P1 1 ppm, rbath 500 (103 spins), JT axes from a file | CCE2 | nstate 0 / nstate 2 | Hahn | coherence |
| p1_cce2e2n, p1_cce4, p1_pcce2_sk2 | same (rbath 400 for CCE2e2n) | CCE2e2n / CCE4 / pCCE2 `pcce.partition_size` 2 | single | Hahn | coherence |
| p1_gcce2_zfs | NV + P1 10 ppm, B = 0, zero-field splitting | gCCE2 | single | Hahn | coherence |
| p1_gcce2_population | NV + P1 10 ppm | gCCE2 | single | Ramsey | dm |
| p1_cpmg8_manual, p1_manual_hahn | NV + P1 1 ppm, `sequence` | gCCE2 | single | CPMG-8 / Hahn | coherence |
| vp_cce2 | NV + Vp 10 ppm | CCE2 | single | Hahn | coherence |
| hbn_cce2_hf3_qd2 / _hf3_qd0, hbn_hfmedi | V_B⁻ + h-BN, DFT hyperfine and quadrupole tensors (`inputs/hbn/hBN_VB_flat_ATensor_HSE_CellPBE_v3`) | CCE2 | ensemble | Hahn | coherence |
| hbn_gcce2_dm, hbn_gcce2_oh | same | gCCE2 | single | Hahn | dm / coherence (OH term) |
| arbitrary_gcce1_dm | the hand-written three-file system of tutorial 01 | gCCE1 | single | Hahn | dm |
| p1_noise2_single, natab_noise2_cpmg8 | p1_cce2_single and natab_cce2_cpmg8 with `method` `noise` (the second with `hamiltonian.noise.pseudo_secular` true); written by hand from those two, not by `make_cases.py`. No tagged version has the method, so only the `head-*` references from 2026-10-28 on hold them | noise | single / ensemble | Hahn / CPMG-8 | coherence |

Everything a case draws at random in a tutorial is a file here: the bath states
(`state_1`, `state_2`), the P1 Jahn-Teller axes (`avaax_nav12`, `avaax_nav4`),
the ¹⁴N sublevels (`exstate_*`). Two binaries that read the same files can be
compared. With `seed` alone they cannot, because the random sequence changed on
this branch. `make_cases.py` wrote the cases and these files from the tutorials
and can write them again.

The full density-matrix files of the arbitrary case (`savemode` cluster_density_matrix) are named after the cluster:
`out_0_full_state1` (the bath spin 0) and `out_full_state1` (the zeroth cluster). Older references call them
`out_0_1` and `out_0` (head) or `out_0_0` and `out_0` (tags); `expected.json` maps the names (`renames`) for every
reference, and `compare.py` compares the values (1.6e-8 against both tags on 2026-10-02). The per-cluster files (`out_state1`, `out_0_state1`) have no reference.

Not in the pool: the two-NV density matrices of tutorial 03 (`05_multiNV/gCCE2_HahnRamsey_NV2`). Every version up to
64f0a1a41, and the tutorial's own recorded result, held `nan` in them (elements the free evolution keeps at zero,
divided by zero). From 64f0a1a41 on they are finite, and the case can join the pool once the tags' `nan` is named in
`expected.json`. Also out: the WAHUHA sequences of tutorial 06, whose gCCE coherence grows past 1e4 under every
version. Both are calculations to look at, not references.

## The fixed-size path, bit for bit

`test_fixed_size.cpp` (built and run by `make cpu-test`, CMake `--target cpu-test`) calls the single-state cCCE
evolution of `simulator_cce.cpp` with the fixed-size types it uses for bath dimension 2 and 4 and with the dynamic
ones, `evolution_update: exact`. It requires every output to have the same bits, the sign of a zero included (the
fixed size once gave (1, +0) where the dynamic path gives (1, -0)). It runs that case and 2,000 random Hermitian
pairs, Hahn and CPMG-2. `recursive` is not held to this: its fixed size keeps Eigen's unrolled products and may
differ from the dynamic path in the last bits.

The bits are held only on the builds where they were seen to agree, chosen by compiler macros as the batched
kernel's table chooses: icpx 2023.2 with AVX-512 (idm4, gpu), icpx 2026.1 with AVX2 (wsl-icpx), g++ 13 with AVX2
tuned for the WSL machine (`__tune_alderlake__`: wsl-gnu, generic there). Any other build (the CI's generic preset
on its runner gave 1-ulp differences, CI run 38487968260; the portable CMake build) holds each output within 1e-12
relative to the larger magnitude or 1e-15 absolute, and prints how many outputs have other bits and the largest
difference. The outputs are overlaps of unit vectors, so rounding leaves a few ulp of 1 on any value, near 0 included.

## The batched small clusters, bit for bit

`test_cce_batch.cpp` (also in `make cpu-test` and CMake `--target cpu-test`) holds the batched CPU evolution of
bath dimension 2 and 4 (`include/calculation/cce_batch.h`) to the scalar path's bits. It uses the dumped clusters of
`test/gpu_verification/check_gpu/data/small_cluster/*.bin` (128 production P1 CR60 sigma4 clusters and 51 pCCE2
clusters of check_rdip_ladder's p1_pcce2_spins). Each file runs as dumped and with the other `evolution_update`, and
a one-state file also as the ensemble Hahn echo. Both sides start from the dumped Hamiltonians (`cceSetup` +
`cceEvolve` against `cceSetup` + the kernel).

The kernel reproduces the arithmetic of the builds it was checked on (g++ 13 AVX2, icpx 2026.1 AVX2, icpx 2023.2
AVX-512). In any other build it is off and the check reports every output as "scalar only". Run it in every build the
batch is to be used in: a compiler change can move where the scalar path fuses its products. `CCEX_CPU_BATCH=0` turns
the batch off in a run (for A/B comparisons).

## The column blocks

`test_cce_columns.cpp` (also in `make cpu-test` and CMake `--target cpu-test`) holds the one-state cCCE evolution
in column blocks (`cceEvolveColumns`, input block `column_evolution`) to the one-time-point path it replaces. Inputs:
random Hermitian Hamiltonians of bath dimension 300 and 600 (the block's `min_dim` set to the dimension), `exact` and
`recursive`, Ramsey, Hahn echo, CPMG-8 and two `Manual` sequences (four different intervals; repeated ones), block
sizes 32 (not a divisor of nstep 77 and 41), 1 and more than nstep.

- Every output within 1e-12 relative to the larger magnitude. The two sum the same products in another order
  (1.9e-15 with MKL, 2.1e-15 without on WSL).
- The phase rule is held separately, since the relative test cannot see it (`exact` and `recursive` differ by about
  1e-16 there). After a run both paths leave the last step's phases of every distinct interval in the workspace,
  compared for last steps around the recursive anchors (16, 32, 64) and the block edges. `exact` and the recursive
  step factors must be the bits of the path before; `recursive` at an anchor step the bits of `exact`; between anchors
  carried (some bit differs from `exact`, all within 1e-11) and within 1e-14 of the path before.
- Not bit for bit between anchors: the carried product is a complex multiplication whose fused multiply-adds the
  compiler places per translation unit (wsl-gnu: the same bits with MKL, 280 of 480 intervals 1 ulp apart, at most
  8.1e-16, with `CCEX_BLAS=off`).
- An anchor of 15 instead of 16, or sin and cos recomputed at each block's first column, fails it.
- With the block on but the dimension below `min_dim`, the outputs are the bits of the path before.

## The references

| folder | what | how it was made |
|---|---|---|
| `head-<date><time>-<commit>/` (older ones `head-<date>-<commit>/`) | this branch at that commit, on WSL, 2 ranks, MKL off (`CCEX_BLAS=off`) | `bash make_reference.sh current` |
| `v1.1.0/`, `v1.2.2/`, `v2.0.1/` | the tagged releases, built from the tag on the same WSL | `bash make_reference.sh v2.0.1` (a worktree of the tag under /tmp or `CCEX_WORKTREE_DIR`, `do_compile.sh wsl-gnu`, removed afterwards; a case the tag cannot run, such as method noise before 2.1.0, is left out) |

The name carries the time so that two references of one day sort in the order they were made (only the newest by
name decides). v1.0.0 does not build with the current compiler (`src/bath.cpp: taking address of temporary array`)
and has no reference. v2.0.1 gives this branch's numbers exactly in all 25 cases it runs. v1.1.0 and v1.2.2 give the
same numbers as each other. Against them this branch differs on purpose in the cases `expected.json` names: ensemble
cCCE returned the complex conjugate for an even pulse count until 8185aa31e (Ramsey and CPMG change the sign of the
imaginary part), deltat was a float until 39d03aa5c (the fast P1 phases and the population at 5 ms move), the seeded
pCCE partition search draws a different partition, and the later fixes listed below.

When a change is meant to move a number, make the new head reference after it (`bash make_reference.sh current`)
and say why in the commit. `expected.json` records, per case and per earlier reference, why the case differs and by
how much. The head references made for a deliberate change:

| reference | change | cases that moved (noDiv/dm, wiDiv, against the previous head) |
|---|---|---|
| the three reference folders of 2026-10-02 (re-recorded) | inputs that wrote `"alphams": -1` at the top of the qubitfile form ran with the NV branch +1; the inputs now give the branch as `"alpha_state": {"ms": [-1]}`, and CCEX refuses the old key. All three were run again with the same command (NP 2; the tags built from the tag) | `p1_cce2_ensemble`, `p1_cce2_single`, `p1_cce2e2n`, `p1_cce4`, `p1_cpmg8_manual`, `p1_manual_hahn`, `p1_pcce2_sk2`, `vp_cce2`: up to 1.26 (`p1_cce2_single`), 8e-6 for the ensemble case |
| `head-20261013-8a775921f` | the h-BN hyperfine file `hBN_VB_flat_ATensor_HSE_CellPBE_tag` (v1 data, the OUTCAR hyperfine x 2, with a `VASP` first line, read x 1: every tensor 2x) replaced by `inputs/hbn/hBN_VB_flat_ATensor_HSE_CellPBE_v3` (tensor format v3, the v1 numbers / 2, read x 1) | `hbn_cce2_hf3_qd0` 5.37e-2 / 5.37e-2, `hbn_cce2_hf3_qd2` 0.694 / 3.46e-2, `hbn_gcce2_dm` 0.544 / 0.619, `hbn_gcce2_oh` 0.614 / 0.117, `hbn_hfmedi` 0.880 / 0.775; the other twenty byte-identical |
| `head-202610132130-2b2cc3e16` | the dimension 2 and 4 cCCE eigensolves by a fixed four-sweep complex Jacobi (`CCEX_SMALL_EIG=0` gives the earlier bits) | `p1_cce2_single` 9.2e-10 / 9.9e-10; `natab_cce2_hahn`, `natab_manual_hahn`, `p1_cce2_ensemble`, `vp_cce2` 1e-10; `natab_cce2_cpmg8`, `natab_knight_shift` wiDiv 1e-10; `p1_cce4` wiDiv 4.8e-4 and `p1_pcce2_sk2` wiDiv 0.38 on diverging curves (about 2e-7 of the value); the rest byte-identical |
| `head-202610152038-5381794ce` | the qubit Hamiltonian subtracts the static Overhauser field of the cluster's own spins only with `hamiltonian.common.qubit_bath_mean_field` true | the six gCCE cases with nstate > 0 and the mean field off: `arbitrary_gcce1_dm` 4.1e-8 / 4.1e-8, `hbn_gcce2_dm` 0.103 / 0.246, `p1_cpmg8_manual` 9.47e-2 (wiDiv diverges), `p1_gcce2_population` 1.71e-3 / 1.89e-3, `p1_gcce2_zfs` 1.81e-3 / 2.50e-2, `p1_manual_hahn` 1.76e-4 (wiDiv diverges); the other nineteen byte-identical |
| `head-202610161458-9d22ab633` | `hbn_gcce2_oh`: the V_B⁻ zero-field splitting (`Qubit.intmap` self term `[[-1080,0,0],[0,-1180,0],[0,0,2320]]`) is MHz, as in tutorial `04_hBN/gCCE2_Hahn_OH`, so the input says `"unit": "MHz"` (D 3.48 GHz, E 50 MHz) | `hbn_gcce2_oh` noDiv 0.85, wiDiv 4.69e-2 (\|L\| at 0.5 ns 0.054 → 0.81, the dip near 1/D 0.29 ns); the other twenty-four byte-identical |

Keep `head-20261012-aa883d0d4`: it was recorded with the `_tag` tensor file, and
`../check_input_forms/tensor_format.sh` needs it to show that an old file is read as before.

The tags read neither the current input keys nor `"unit"` (v1.1.0), so `cases/hbn_gcce2_oh/ccein_tags.json` is the
same calculation in their input form (the ZFS in kHz, x 1000; the v3 tensor file, which they read x 1 as a `VASP`
file with `CorrTotSpin` 0). Its tag outputs were recorded from the case folder with each tag's binary (built from the
tag with its own `do_compile.sh wsl`; the tags have no `wsl-gnu` preset), 2 ranks, `-o <folder>/out`. The same procedure with the 2026-09-24 input and its `_tag` tensor file gives the
recorded bytes of both tags. Both tags agree
with the current head within 2.9e-8 (noDiv) / 8.3e-10 (wiDiv).
