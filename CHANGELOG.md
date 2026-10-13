# Changelog

## 2.0.0 — unreleased (changes since v1.2.2, 2026-09-04)

### Added
- **DEER bath pulses** (the block `pulse.bath`: `count`, `sequence`, `defect`, `energy_shift`, `spin`, `alphams`,
  `betams`) for gCCE (matrix and vector, exact and recursive) and megCCE ($U_q\otimes U_b$ in its pulse superoperator)
  on the CPU, with `cluster_mode` standard or pcce. Ported from the adding_DEER/adding_MECCE branches, whose flat keys
  (`bnpulse`, ...) are refused with the new names. See `docs/methods/deer.md`.
- **ME-CCE / ME-gCCE** (`method`: `meCCE`, `megCCE`; `hamiltonian.meCCE.jump_operators`,
  `hamiltonian.megCCE.jump_operators`) — Lindblad dissipation of bath spins (meCCE, megCCE) and of the qubit (megCCE
  only), CPU. Ported from adding_MECCE; meCCE refuses qubit jump operators, which the branch multiplied into every
  cluster. Method: Onizhuk et al., Phys. Rev. Lett. 132, 250401 (2024). See `docs/methods/master-equation.md`.
- **`cluster_mode: "pcce"` with gCCE, meCCE and megCCE** (inside pgCCE, pmeCCE, pmegCCE): the pCCE partitions with the
  gCCE and master-equation kernels, on the CPU; `pcce.cache` stays with cCCE. The adding_MECCE branch's method names
  `pgCCE`, `pmeCCE` and `pmegCCE` are refused with the new form (method + `"cluster_mode": "pcce"`).
- **Reserved method name `noise`**: `method` (and `-m`) accepts the name and refuses it as not implemented yet;
  `hamiltonian.noise` is an empty section like `hamiltonian.gCCE`.
- **`hamiltonian.common.bath_bath_mean_field: "off"`** (nstate > 0): the static field of the bath spins outside a cluster
  left out (`"z"` stays the default, `"xyz"` is the former `vector_mean_field: true`).
- **`hamiltonian.cCCE.delta_A: false`**: the hyperfine tensors of a cluster's spins replaced by their mean, so the pair
  flip-flop of different Knight fields is gone (a comparison calculation; clusters of one spin unchanged).
- **`example/tools/InputKeys/convert_input_keys.py`**: rewrites inputs written with the old keys (comments kept, the
  changed keys listed, any number of files; `--check`, `--stdout`).
- **GPU backends** (`backend`: `cpu`, `cuda`, `auto`; settings in the `gpu_options` block) for cCCE, gCCE and pCCE, with
  `propagator` (`eigen`, `expm`), `evolution` (`vector`, `matrix`) and `gcce_refinement` for gCCE. See `docs/methods/gpu.md`.
- **Adaptive clustering** (`cluster_mode: adaptive` with an `Adaptive` input block): bundles of strongly coupled spins with
  the merge policy, and `"regroup": true` in that block (command line `--regroup`) to regroup where a bundle expansion
  gives |L| > 1.
  See `docs/methods/adaptive-cluster-mode.md` and `docs/methods/adaptive-input-block.md`.
- **pCCE**: a deterministic partition search, a partition cache (`pcce_cache_file`, `pcce_cache_mode`) and
  `pcce_connections` (`centers`, `spins`).
- **Automatic convergence**: an `auto_convergence` block in the input (`main.out -f input.json`; `"enabled": false` runs
  the ordinary calculation; command line `--auto_conv`) chooses rbath, rdip, the time window and the number of bath
  configurations and states (`example/tools/AutoConvergence/README.md`). Neighbouring candidates are compared on the
  area under a T2 fit of |L|; where the fit fails, on the area under |L| itself (same 2 %), and the summary says so.
- **`rdip_ladder`**: one order-2 run (cCCE, gCCE or pCCE) also writes the curves of smaller rdip values.
- **`evolution_update`**: `exact` or `recursive` (the phases by a recurrence; 1.0–1.6x faster on the tested systems,
  where the curves stay within 1e-9 of `exact`; not bit-identical). The default is `recursive` where it exists (see Changed).
- **`propagator: "sparse"`** (opt-in; also `{"method": "sparse", "min_dim": 4096, "tol": 1e-12}`): one sampled state of
  cCCE and gCCE `vector` without diagonalizing the cluster Hamiltonian (Krylov action on the state, H sparse), for
  clusters of dimension ≥ `min_dim`. It pays only where the dense eigendecomposition does not fit in memory (dim 16384:
  about 36 GB dense against about 0.2 GB sparse). It is much slower at high field or with a wide spectrum and long times:
  the cost grows with spectral width × time and with the square of the number of time points (P1 at 500 G, dim 4096,
  151 points: about 18 h against 5 min dense, extrapolated from the first 11 points). `propagator` is accepted by cCCE now (`eigen`, `sparse`; `expm` stays gCCE).
- **Built-in isotope table**: `gyrofile` is optional; species it does not name come from the table.
- **`savemode.nodiv`**: `false` writes `_wiDiv` alone (see the `savemode` block under Changed).
- **Packed GPU kernel for small cCCE clusters** (`gpu_options.packed_max_dim`: 4 by default in `bin/main_gpu.out` with
  `backend` `cuda` or `auto` and method CCE, 0 with `backend` `cpu` and in `bin/main.out`; 2 or 4 may be given): every cCCE
  cluster of bath dimension 2 and 4 is computed by one GPU thread, thousands per kernel, for the whole node through one CUDA
  context (node rank 0 feeds the GPU from node shared memory), the rest unchanged; single state and ensemble, Hahn, CPMG,
  Ramsey, `exact` and `recursive`; `bin/main_gpu.out` and `method` CCE. Within 1e-8 + 2e-6 |CPU| of the CPU, not
  bit-identical. On idm4-gpu (32 ranks, one GPU) the production P1 cCCE2 runs took 0.59–0.84 of the CPU's time; a run of a
  few seconds gains nothing (CR0 ensemble, 1.8 s: 0.95). See `docs/user/InputFile_Backend.md#tag_gpu_options` and
  (internal notes).
- **Warp GPU kernel for gCCE and cCCE clusters up to dimension 27** (`gpu_options.warp_max_dim`: 27 by default in
  `bin/main_gpu.out` with `backend` `cuda` or `auto` and method CCE or gCCE, 0 with `backend` `cpu` and in `bin/main.out`;
  0–32 may be given): one cluster per GPU warp (shared-memory Jacobi, the lanes over the time points, exact phases), through
  the packed kernel's node feeder; gCCE with a qubit of dimension 3 (total dimension 6, 9, 12, 18, 21, 24, 27), single state
  and ensemble, coherence and dm; cCCE bath dimension 8 and 16. Runs it does not take (`gcce_refinement`, `pulseiter`,
  savemode `cluster_density_matrix`, ...) switch it off by themselves. On idm4-gpu (32 ranks + MPS, one GPU) P1 gCCE2 took 0.60–0.69 and
  h-BN gCCE2 0.86 of the earlier GPU build's time. See `docs/user/InputFile_Backend.md#tag_gpu_options` and
  (internal notes).
- **Builds**: `do_compile.sh` presets (`wsl-gnu`, `wsl-icpx`, `idm4`, `gpu`, `nurion`, `mac`, `generic`) and a CMake build.
- **MKL for large clusters, when the build finds it** (`do_compile.sh` every preset and CMake, `CCEX_BLAS` auto by
  default; `off` without it, `mkl` to insist; the build and the run log's header say on or off): the cCCE ensemble Hahn step's D³ product from bath dimension 21 by zgemm and the cCCE/gCCE cluster
  eigensolve from dimension 63 by zheevd (MKL lp64, sequential). idm, 64 ranks: h-BN cCCE2 ensemble 0.80–0.88, gCCE3
  0.72–0.80, gCCE2 0.98 of the time; P1 production inputs unchanged. The curves differ in the last digits (cCCE ≤ 2e-10,
  h-BN gCCE late times up to 9e-8). (internal notes) A46.

### Changed
- **Input keys renamed and grouped** (A47, 2026-10-13; the table of old and new names is in `docs/user/InputFile.md`):
  `method` is `cCCE` or `gCCE` and the cluster shape is `cluster_mode` (`standard`, `pcce`, `bundle` - the former
  `adaptive` -, `dsj`, `itb`, `dsjitb`); the Hamiltonian settings are one `hamiltonian` block (`common`:
  `qubit_bath_mean_field`, `bath_bath_mean_field`, `hyperfine {dipole, contact, file, missing_row}`, `quadrupole`; `cCCE`:
  `knight`, `hfmedi`, `delta_A`); `pcce`, `bundle`, `adaptive`, `pulse`, `gcce_refinement` and `cluster_contributions`
  are blocks; `Qubit.file` replaces `qubitfile` (with `spin`, `zfs`, `alphams`, `betams`), `Qubit.alpha_state` /
  `beta_state` (`{"index"}` or `{"vector"}`, each side on its own) replace `alphaidx`/`betaidx`/`psia`/`psib`, the length of
  `Qubit.qubit` replaces `nqubit`; `rbath`/`rdip` take `[inner, outer]` (`rbathcut`/`rdipcut`); `bath_shift`,
  `bath_state_file`, `add_subclusters`, `max_clusters_per_order`; `Defect` is `{state_file, axis_file, species: [...]}`
  with `name`, `n_nuclear_spins`, `n_configurations`, `quadrupole_moments`, `nuclear_positions`, `hyperfine`. Numeric
  codes became words and 0/1 switches true/false. `-m` takes `cCCE` or `gCCE` and keeps the input's `cluster_mode`.
  Results are unchanged: every input of the repository, converted, gives the bytes it gave before.
- **Faster CPU paths with unchanged results** (byte-identical output): the cCCE and gCCE product assembly, the state preparation,
  fixed-size paths for small clusters, the order-2 cluster list built without duplicated work, the small-cluster
  Hamiltonians, and the pCCE partition search (up to 9x faster on large baths). See Performance below.
- **Small cCCE clusters batched on the CPU** (byte-identical output): clusters of bath dimension 2 and 4 are evolved
  4 (AVX2) or 8 (AVX-512) at a time, in the builds whose arithmetic the batched kernel reproduces (g++ 13 AVX2,
  icpx 2026.1 AVX2, icpx 2023.2 AVX-512; any other build keeps the one-by-one path). One state with any pulse plan and the
  ensemble Hahn echo, `exact` and `recursive`. Production P1 cCCE2 on idm (64 ranks): 0.60–0.64 of the time with
  `recursive`, 0.74–0.77 with `exact`. See (internal notes).
- **Repository layout**: `src/` and `include/` in five parts (input, physical, clustering, calculation, common);
  `example/{tutorial,analysis,tools}`; `external/`; `test/{code_verification,cluster_verification,gpu_verification}`.
  The tutorials were rewritten as numbered folders (twelve; 11 and 12 added on 2026-10-12) with `run.sh`.
- In input files, `#` and `!` start a comment only outside JSON strings.
- **GPU settings moved into one block, `gpu_options`**: `min_dim` (was `gpu_min_dim`), `evolution_batch_size`
  (`gpu_batch_clusters`), `evolution_batch_steps` (`gpu_batch_steps`), `cluster_batch_mb` (`gpu_cluster_batch_mb`),
  `cache_mb` (`gpu_cache_mb`), and the new `packed_max_dim`. The old flat keys are refused with a message naming the new
  form (e.g. `gpu_min_dim is now gpu_options.min_dim`); `backend` keeps `cpu`, `cuda`, `auto`.
- **`savemode` is a string or a block `{"mode": ..., "nodiv": true|false}`, with new names**: `coherence` (was `normal`,
  the default), `coherence_tagged` (was `info`), `cluster_coherence` (was `all`), `cluster_density_matrix` (was
  `allfull`); `nodiv` (default `true`) replaces the top-level `save_nodiv`. The files written are the same as before for
  the same setting. See `docs/user/InputFile_Output.md#tag_savemode`.
- **`bin/main_gpu.out` with `backend` `cuda` or `auto` runs the cCCE clusters of bath dimension 2 and 4 in the packed
  kernel by default** (`gpu_options.packed_max_dim` 4; `"packed_max_dim": 0` gives the earlier path). Its results agree with
  the CPU within 1e-8 + 2e-6 |CPU| and are not bit-identical to those of the earlier GPU path. `backend` `cpu` and
  `bin/main.out` are unchanged.
- **`bin/main_gpu.out` with `backend` `cuda` or `auto` runs the gCCE and cCCE clusters of dimension up to 27 in the warp
  kernel by default** (`gpu_options.warp_max_dim` 27; `"warp_max_dim": 0` gives the earlier path). Its results agree with
  the CPU within 1e-8 + 2e-6 |CPU| outside the stretch where a product diverges (|L| > 1 and between) and are not
  bit-identical to those of the earlier GPU path; `backend` `cpu` and `bin/main.out` are unchanged.
- **Short runs leave the packed and warp kernels off** (A40 follow-up, 2026-10-13): with neither `packed_max_dim` nor
  `warp_max_dim` given, a run whose cluster evaluations for them (clusters x states) are fewer than 20000 per rank over all
  ranks does not use the card (its early start is closed again)
  and computes everything on the CPU path (byte-identical to both set to 0), with one `off: ... fewer than the threshold`
  log line; and the warp kernel's shared-memory region is allocated only when some rank has a warp cluster (P1 CR0 cCCE2
  had been 4.7 % slower from it). Reasoning: (internal notes).
- **The default of `evolution_update` is `recursive` where it exists** (the CPU cCCE, single state and ensemble, the packed
  kernel, and the CPU gCCE with `propagator` = `eigen`, `evolution` `vector` or `matrix`), and `exact` elsewhere (the GPU
  large-cluster path, gCCE `expm`), without a message, per path: in `bin/main_gpu.out` with `backend` `auto` the clusters
  on the CPU and the packed kernel run `recursive` and the clusters on the GPU `exact`, with `cuda` the packed kernel
  `recursive` and the rest `exact`; `backend` `auto` in `bin/main.out` computes everything on the CPU and so runs
  `recursive`. The run log says which and why. The outputs are therefore not bit-identical to earlier runs (within
  1e-9 on the tested curves at |L| ≤ 1). `"evolution_update": "exact"` gives the earlier bits.
- **Clusters of an order are split over the ranks by cost when their sizes differ** (dimension D, cost D²(D + nstep); each
  rank takes a contiguous run of the list): mixed baths such as h-BN run 1.6–3x faster on 64 ranks (cCCE2 1.6x, gCCE2 2.3x,
  cCCE3 3.1x). Orders whose clusters all have one size (P1, ¹³C) keep the equal-count split, byte for byte; in a mixed bath
  the product of each rank's clusters is grouped differently, so the last printed digit can change (h-BN gCCE2 dm: 5 of
  36,000 numbers, 1e-10). Data (internal notes).
- **`gpu_options.min_dim` defaults to 768 (was 256)**: on idm4-gpu (RTX PRO 6000, 32 ranks; A37, job 699) the GPU was
  faster than the CPU from a cluster dimension of about 670–750 on, on every path (cCCE and gCCE, single state and
  ensemble), and slower at 648 and below. So in `bin/main_gpu.out` with `backend` `auto` and no `min_dim` given, clusters
  of dimension 256–767 now run on the CPU: their last digits can change (within the GPU tolerance 1e-8 + 2e-6 |CPU|),
  and the A35 cCCE4 ensemble case ran 22.3 → 9.6 s (measured at `min_dim` 512, which routes it as 768 does). An explicit `min_dim` keeps its value. Data
  (internal notes).
- **`backend` `auto` no longer stops on a GPU that cannot be opened or that runs out of memory**: the packed kernel is
  turned off when the GPU cannot be opened, and a cluster that runs out of GPU memory is computed on the CPU, with every
  later cluster of its dimension or larger on that rank (`[gpu-memory-fallback]` in the log). `backend` `cuda` still
  stops in both cases.
- **Each rank keeps only its own order-2 cluster rows** (cCCE/gCCE order 2; `-P`, `cluster_contributions` and adaptive
  work keep the full list): peak memory per rank on the P1 production CR60 σ4 bath 334 → 181 MiB, output byte-identical,
  same speed. Data (internal notes).
- **Hyperfine tensor files, format v3**: SpinTensorReader writes line 1 `VASP` / `QUANTUMESPRESSO` and line 2
  `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`,
  and the numbers are the hyperfine of the spin-S defect for both codes (QE's S = 1/2 coupling x 1/(2S), `--spin S`
  required; VASP's OUTCAR as it is, S from its magnetization). CCEX reads such a file x 1: `DefectTotSpin` is only
  checked against S and `CorrTotSpin` is not used. A file without the line is read by the earlier rule (first line `VASP`
  -> x 1, otherwise x 0.5/`DefectTotSpin`), with a warning in the log; its results are unchanged.
  `--fidentify` is removed (the lines are always written). See `docs/user/InputFile_General.md#tag_hyperfine_file`.

### Compared with v1.2.2 (same inputs and node, idm Xeon Gold 6530, 64 ranks, median of 3)

| Case | v1.2.2 | 2.0.0 (`exact`) | 2.0.0 (`recursive`) |
|---|---|---|---|
| P1 production cCCE2, CR0 | 76.0 s | 15.2 s (5.0x) | 9.3 s (8.1x) |
| P1 production cCCE2, CR60 sigma 4 | 463.5 s | 79.0 s (5.9x) | 45.4 s (10.2x) |
| P1 production gCCE2, CR0 | 90.0 s | 39.6 s (2.3x) | 30.0 s (3.0x) |
| P1 production gCCE2, CR60 sigma 4 | 544.9 s | 216.4 s (2.5x) | 161.6 s (3.4x) |

The coherence differs from v1.2.2 only where a fix in this release changes it (deltat read as a float in v1.2.2,
the sampled states and `seed`, the ensemble half interval in double, the density-matrix `nan`, the pCCE partition).

### Performance (measured during the October 2026 work)

Each row compares the code before and after that change, on the same inputs and server (idm, 64 ranks unless noted);
the results are byte-identical, or within the stated bound for `recursive`. For v1.2.2 see the table above.

| What | Case | Before -> after | Speed-up | Baseline |
|---|---|---|---|---|
| P1 production cCCE2, 30 states | CR60 sigma 4 | 229 s -> 85.9 s | 2.7x | code of 2026-10-09 |
| | CR0 | 21.3 s -> 16.5 s | 1.3x | code of 2026-10-11 |
| P1 production gCCE2, 30 states | CR60 sigma 4 | 313.3 s -> 224.7 s | 1.4x | code of 2026-10-11 |
| | CR0 | 56.1 s -> 41.2 s | 1.4x | code of 2026-10-11 |
| pCCE partition (P1, 30,000 trials) | rbath 483 A | 1,616 s -> 173 s | 9.3x | code of 2026-10-11 |
| Automatic convergence, h-BN gCCE2 ensemble | whole run | 14,406 s -> 4,281 s (CCEX runs 596 -> 103) | 3.4x | code of 2026-10-11 |
| Automatic convergence, pCCE runs differing only in rdip | NV-13C, 4 ranks (local) | 8.2 s -> 0.3 s | 27x | code of 2026-10-11 |
| Product assembly, P1 production cCCE2, 30 states, `recursive` | CR60 sigma 4 | 55.7 s -> 45.7 s (assembly per rank 12.9 s -> 4.9 s) | 1.2x | code of 2026-10-12 |
| | CR0 | 10.8 s -> 9.4 s (assembly per rank 2.3 s -> 0.9 s) | 1.15x | code of 2026-10-12 |
| Product assembly, P1 production gCCE2, 30 states, `recursive` | CR60 sigma 4 | 167.9 s -> 162.2 s (assembly per rank 17.3 s -> 10.0 s) | 1.04x | code of 2026-10-12 |
| | CR0 | 31.0 s -> 29.9 s (assembly per rank 3.1 s -> 1.8 s) | 1.04x | code of 2026-10-12 |
| `evolution_update: recursive` (default where it exists) | 16 inputs (NV-13C, h-BN, P1) | | 1.0-1.6x | `exact` |

### Fixed
- Ensemble cCCE: the conjugation of the coherence and the half interval computed in double precision.
- Density-matrix output: a component that divides noise by zero is set to 0 and reported, instead of `nan`.
- `seed` had no effect on sampled states when `nstate` > 0.
- `savemode: "allfull"` (now `cluster_density_matrix`) wrote no result of the expansion (the save step knew only `normal`, `all` and `info`), named
  each cluster's full density matrix by positions 0, 1, … so that clusters of the same size overwrote one another (and
  every rank wrote the zeroth cluster), and with cCCE wrote nothing. It now writes what `all` writes and, with gCCE,
  each cluster's full density matrix to `<outfile>_<bath indices>_full[_state<n>]_wiDiv/noDiv`. Other savemodes unchanged.
- Builds with Intel `icpx`: the adaptive code is compiled without fast floating point, which had folded `isnan`/`isinf`.
- The analyzer's plots of several curves (`-pl_m`) failed with matplotlib 3.11 (`cm.get_cmap` removed); the numbers it writes are unchanged.
- **`hfmedi` (hyperfine-mediated interaction)**: the term is now the complete second-order correction of each branch,
  H_med^(α) = Σ_{β≠α} E_αβ E_αβ† / (ω_α − ω_β) with E_αβ = Σ_i ⟨α|S|β⟩·A_i·I_i: every intermediate eigenstate β of the
  qubit Hamiltonian (it used only the other qubit level, missing for example the |+1⟩ path of the NV |0⟩ branch), every
  ordered pair (i, j) (it summed i < j only, which was not Hermitian for a spin-1 qubit: the CPU used one triangle of it
  and the GPU refused it), and the i = j terms. Hermitian by construction. Results with `hfmedi: true` change.
  An intermediate level degenerate with the branch (within 1e-12 of the spread of the levels, so an energy shift does not
  change the decision) is skipped when the coupling E_αβ vanishes (spin-1 at zero field: +1 and -1, which gave `nan`;
  or a zero hyperfine) and refused when it does not.
  The h-BN inputs that use it (tutorial 07 `hyperfine_mediated_term`) had an asymmetric `qzfs` typo, now
  `[[-1080000,0,0],[0,-1180000,0],[0,0,2320000]]`; with both, the `check_numerics` case `hbn_hfmedi` moves by up to 0.95.
- **`knight` (Knight shift)**: the term multiplies by the nuclear Zeeman vector −γ_n B instead of the bare field B, and
  the (3|m_s| − 2) factor takes |m_s| of each branch as ⟨S_z²⟩ (it was always 1). Zero for a field along the qubit axis,
  as before; in a tilted field results with `knight: true` change (`natab_knight_shift`: up to 0.98).
- **The h-BN tutorial and test hyperfine tensors were read 2x**: `hBN_VB_flat_ATensor_HSE_CellPBE_tag` (tutorial 07 and
  the `check_numerics` copy) was SpinTensorReader v1 data (the OUTCAR x 2) with a `VASP` first line, which CCEX reads x 1,
  so every hyperfine tensor was twice the DFT value (nearest ¹⁴N A_zz 97.3 MHz against 47 measured) and the
  hyperfine-mediated term 4x. They are replaced by `hBN_VB_flat_ATensor_HSE_CellPBE_v3` (format v3, the numbers / 2,
  A_zz 48.6 MHz); the `check_numerics` h-BN cases change by up to 0.88 and were re-recorded, and tutorial 07 and its
  figures were rerun. SpinTensorReader example 03 (MoS2) had the same mistake. Studies run on the old file
  ((internal notes)…`, `a20-…`, `a27-…`, `a33-…`) used the doubled hyperfine.
  (internal notes).
- SpinTensorReader's QE writer used names it never defined (`sys`, `defectposition`, `fout_Afile`/`fout_Qfile` in its
  summary), so its error paths crashed instead of stopping with their message; and its tag was misspelt
  `QUANTUMESSPRESSO`.

### Refused now (inputs that used to be accepted)
- **The input keys before A47** (2026-10-13): an old name, value (`method` `CCE`, `pCCE`, `dsj`, ...; `cluster_mode`
  `adaptive`) or form (`Defect` as a list, `gcce_refinement: true`) stops the run with the new one named and the
  converter pointed to; `DefectTotSpin` (the qubit spin is used) and `CorrTotSpin` (old-format VASP tensor files only)
  are removed. Also refused: `gCCE` with `cluster_mode` `dsj`/`itb`/`dsjitb` (not implemented yet), a
  `hamiltonian.cCCE` setting in a gCCE run, the `adaptive` block outside `pcce`/`bundle`, the `pcce` block outside `pcce`,
  a qubit state with both `index` and `vector`, 0/1 for the switches.
- `alphams` at the top of a `qubitfile` (it had been ignored since v1.1.0).
- Arrays shorter than needed (Defect types shorter than `naddspin`, `bathadjust` rows with fewer than 3 numbers).
- An `avaax` axis beyond the defect's `navaax`.
- The `savemode` names `normal`, `info`, `all`, `allfull` and the top-level key `save_nodiv`, each with a message naming
  the new form (e.g. `savemode "all" is now "cluster_coherence"`, `save_nodiv is now savemode.nodiv`); also an unknown key
  of the `savemode` block and a value of the wrong type.
- A hyperfine tensor file of format v3 that contradicts itself or the input: a first line that is not the tag of the
  provenance line's source, the provenance line elsewhere than line 2 or twice, an unknown format version or source, an S
  that is missing or not a half-integer, or an S other than `DefectTotSpin`.
- SpinTensorReader `--version v1` (it wrote the OUTCAR x 2, right only untagged and for S = 1), and QE `--A` without
  `--spin S`.
- Combinations a method does not support; `test/code_verification/check_refusals` lists them.
- `pulse.bath` (DEER) where it cannot act: methods cCCE and meCCE, `cluster_mode` bundle, `backend` `cuda`,
  `propagator` `sparse`, `gpu_options.warp_max_dim` > 0, `pulse.detuning_factor` ≠ 1, `pulse.per_qubit` true; an
  unknown key of `pulse.bath`; the branch flat keys (`bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift`,
  `bspin`, `balphams`, `bbetams`); a last bath row that is not identity and two pulses of one sequence on one edge; a
  `pulse.bath.defect` no bath spin has, a `pulse.bath.spin` other than the driven spins' S, and `alphams`/`betams` that
  are not two different levels of it (a pump that reaches no spin of that name only warns).
- meCCE / megCCE where they cannot act: the adding_MECCE top-level `jump_operators` (named with its new place) and the
  old operator names `+-`/`-+`; a `hamiltonian.meCCE`/`.megCCE` section in a run of another method (also after `-m`), a
  key other than `jump_operators` in it, jump operators in both sections, an entry without `bath_name`, `operator` and
  `rate`, an unknown operator, a rate that is negative or not finite, a `bath_name` that is neither the qubit nor a
  bath spin of the run (checked once the bath is read; the branch ran it at rate 0); for meCCE a Manual sequence, a
  rotation other than an ideal π flip and qubit jump operators (megCCE named); more than one qubit for meCCE or with a
  qubit jump; for both, `pulse.per_qubit` true, `pulse.detuning_factor` ≠ 1, `propagator` other than `eigen`,
  `cluster_mode` other than `standard` and `pcce`, `evolution_update: "recursive"` given, quantity `dm`, savemode
  `cluster_density_matrix`, `rdip_ladder`; megCCE with a `hamiltonian.cCCE` setting.

### Removed
- The MKL dependency, `hf_cutoff`, the bundled `zlib/`, the top-level `main.cpp` and `config/` (moved into `src/` and the tools).
