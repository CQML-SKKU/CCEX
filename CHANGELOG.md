# Changelog

## Unreleased (feature/bath-drive)

### Added

- **Continuous bath drive, `pulse.bath_drive`** (cCCE, the ensemble): electron-like spin-1/2 bath spins on the lines of a
  Defect table driven by one or more tones, each spin in the rotating frame of its tone or of its own line, with the
  rotating-wave approximation (García Oliván et al., arXiv 2512.06948). The qubit coupling A_zx, A_zy is kept exactly by
  tilted branch frames (`qubit_transverse`, default true), `shifts` adds the Bloch-Siegert, AC Stark, Van Vleck and
  drive-assisted cross-line flip-flop terms. A port of pyCCEX `drive.py` (b3ce43e), with the same outputs to the printed
  digits on a 20 ppm P1 bath (Fig. 3's four protocols, Ramsey, Hahn, CPMG-2 and -3, order 1 and 2). New kernel
  `src/calculation/simulator_cce_drive.cpp`, page `docs/methods/bath-drive.md`.
- **The mean field of the driven bath** (`pulse.bath_drive` with `configuration.s_range`): every spin's sampled m along
  the axis of its driven-frame Hamiltonian, the Hartree field of the cluster's pair Hamiltonian from the spins outside the
  cluster and from a shell file (`mean_field` pfull, projected, z, none; `mean_field_shell`, `mean_field_species`,
  `mean_field_shell_state_file`), the coherence averaged over the samples inside every cluster (one pass, the ensemble's
  output names). The same outputs as pyCCEX `drive.py` given the same samples.
- **The qubit's transverse terms near resonance** (`pulse.bath_drive`): the run prints the qubit's admixture of flipped
  states, sum_k |A_perp,k|^2 (1/gap_k^2 + 1/w_q^2), and warns above 0.01, where the cCCE that drops these terms is not
  reliable (an NV and P1 about 512 G, the NV level anticrossing about 1024 G). The accuracy of the rotating-wave model
  away from them is in `docs/methods/bath-drive.md` (Limits).

## 2.0.1 — 2026-10-22 (changes since 2.0.0)

### Fixed

- `findProperty` (`src/clustering/cluster_hash.cpp`, the subcluster lookup of `solveTilde`, which sets the power of each
  cluster in the product) returned an uninitialized pointer when there was no cluster of the asked order at all; it
  returns NULL then, as for a cluster that is not in its order's table. The orders are grown from the lower ones, so no
  run is known to have reached it; clang's `-Wsometimes-uninitialized` found it.
- **GPU gCCE with a qubit dimension of 8 or more (A97).** Built with CUDA 12.8 for the RTX PRO 6000 (compute capability
  12.0), the pulse kernel `applyPulse` (`src/calculation/gpu/gpu_backend.cu`) read through an address register that
  `ptxas -O3` had just overwritten, and the gCCE run stopped with `an illegal memory access` (NV6 gCCE2 `dm`, on 1 rank
  as on 16 ranks with MPS; the CPU binary ran it). The faulty instructions are in the loop body unrolled by 8, which only
  a qubit dimension q >= 8 (the product of 2S+1 over the qubits) reaches: `test_gcce_cuda` passed for q = 3 to 7 and
  failed for 8 and 9 on every repeat. The source is unchanged: CUDA 13.2.2 compiles it right (compute-sanitizer finds
  no error, NV6 gives the CPU's outputs exactly, and 08_GPU's cCCE5, gCCE4 and pCCE2 agree with the 12.8 build to 1e-10,
  the last printed digit). The `gpu` preset and `test/gpu_verification/server/environment.sh` load `cuda/cuda-13.2`
  (`/usr/local/cuda-13.2` on idm4-gpu), `do_compile.sh` and CMake take only nvcc 13.2 or newer (an older one found on
  `PATH` by `generic` or by CMake's `CCEX_CUDA=AUTO` is left out with a warning; one named by the `gpu` preset,
  `CUDA_HOME`, `NVCC` or `CCEX_CUDA=ON` stops the build), and `test_gcce_cuda` has a q = 9 case, so `make gpu-test` catches it. Where 12.8
  must stay, compiling `gpu_backend.cu` alone with `-Xptxas -O0`, or embedding only its PTX for the driver to compile,
  also ran NV6 correctly.

### Builds

- **idm4-gpu builds with Intel MPI, as idm4 does.** The `gpu` preset takes `mpiicpc` with icpx and MKL from the module
  `intel/oneapi-2023.2` (oneAPI 2023.2, Intel MPI 2021.10, MKL 2023.2: the paths of idm4's `intel-23.2/icc-23.2`), with
  GCC 11.5 and CUDA 13.2; it took Open MPI's `mpicxx` with icpx from `setvars.sh` before, and stops now with a message
  when the module is not there. `test/gpu_verification/server/environment.sh` loads the same module, so the jobs'
  `mpirun`, the one the checks take (the Intel MPI of `I_MPI_ROOT` first, `test/README.md`) and the binaries' MPI are
  one; with the Open MPI binaries a check run in that environment needed `CCEX_MPI`. The module conflicts with
  `gcc/openmpi-*` and `nvhpc`, so a job has one MPI. Intel MPI reduces gCCE's density matrices faster there (one
  node with 32 cores, one Xeon 6731P; one rank per core): `MPI_Reduce` with `MPI_PROD` of 729 × 729 complex numbers
  on 16 ranks took 1.55 ms against 12.98 ms with Open MPI
  4.1.5 and 14.00 ms with 5.0.10, and NV6 gCCE2 `dm` (16 ranks with MPS) spent 0.41 s in it on rank 0 against
  1.1-1.3 s, with the same outputs (job 17695). Rank 0 then wrote the 1.14 GB of the two `dm` files in 17.1 s instead
  of 9.7 s in 11 of 21 such runs with Intel MPI and in none of 9 with Open MPI: the same instructions at a lower rate
  (2.33 against 4.18 per cycle) through the whole write, with no other task on its core, the same clock and the same
  share of transparent huge pages; `I_MPI_PIN=off` did not change it (jobs 17695-17697). The cause is not known. A
  `coherence` run writes a few numbers per step, so its write is short.
- **The CUDA files compile with `nvcc -split-compile`.** `warp_cluster_gpu.cu` instantiates 30 kernels (2.17 million
  SASS instructions), and nvcc optimized them one after another: with CUDA 13.2 the file took 1217 s on idm4-gpu
  (one node with 32 cores, one Xeon 6731P), nearly all of a 1331 s build at `CCEX_JOBS=32` (job 17690).
  `do_compile.sh` gives the three CUDA files
  `-split-compile=$CCEX_JOBS` (default 5; the make variable `NVCC_SPLIT`, which `NVCC_SPLIT=` on the make line
  empties), CMake `-split-compile=${CCEX_CUDA_SPLIT_COMPILE}` (default 5; 1 is off). At 32 the file takes 284 s and
  the build 301 s (job 17692). The machine code differs (every kernel of the three files), the results do not: NV6
  gCCE2 `dm` on 1 rank and on 16 ranks with MPS gives the CPU binary's outputs exactly, 08_GPU's cCCE5, gCCE4 and
  pCCE2 those of the build without it, and `make gpu-test` passes 203 of 203 (jobs 17692, 17694). With 2 threads,
  what `test/gpu_verification/server/build.sbatch` gives on its two CPUs, that job took 788 s (the build 765 s; at most
  4.9 GB of memory with the page cache), within its 20 minutes and 8 GiB; the default 5 built in 393 s on five CPUs.
  The machine code also depends on the number of threads, and both builds gave the outputs above exactly (job 17708).
- **macOS (Apple Silicon).** The `mac` preset builds at `-O3 -mcpu=native` (it said `-O2`) and links Apple's Accelerate
  for the large clusters in MKL's place: the same Fortran `zgemm_` and `zheevd_` (`src/common/blas_option.cpp`), with the
  system on every Mac. `CCEX_BLAS=auto` (the default) takes it on macOS, `accelerate` insists, `off` leaves it out, and
  `mkl` stops there (there is no MKL for Apple Silicon). CMake does the same (`-DCCEX_BLAS=AUTO|ACCELERATE|OFF`;
  `CCEX_NATIVE` gives `-mcpu=native` on Apple Silicon). Until now CMake's AUTO found Intel's macOS MKL 2023 under
  `/opt/intel`, x86_64 only, and the arm64 link failed. The run log says `MKL for large clusters : on (Apple Accelerate
  in its place; ...)`, which the checks read as an external library, and `Accelerate on` in its `Build` line.
  Measured on an M4 (4 ranks, median of three; h-BN cCCE3 ensemble Hahn with clusters to dimension 343, h-BN gCCE3 to
  1029, P1 pCCE2 `partition_size` 4 to 256): 76.0 / 19.5 / 26.5 s with the old preset, 71.1 / 19.3 / 26.0 s at `-O3`,
  32.8 / 14.6 / 18.3 s with Accelerate (2.3, 1.3 and 1.5 times faster). `-O3` gives the bits of `-O2`; with Accelerate
  the curves agree with the run without it to the last printed digit or so (noDiv within 7.6e-9; a wiDiv element
  diverged to 1.9e5 moved by 1.5e-9 relative), and the 25 cases of `check_numerics` give the same bits but one
  non-finite element. `VECLIB_MAXIMUM_THREADS=1` changed nothing on 8 ranks and was slower on 4, so it is not set.
- **The memory check on macOS** (A67) reads the rank's resident size (`task_info`) and the free and inactive pages
  (`host_statistics64`, what psutil calls available) where Linux reads `/proc/self/status` and `/proc/meminfo`. Without
  `/proc` both read 0: the estimate left out what the rank held and the node's memory was "unknown (not checked)", so
  a run that does not fit was not stopped (`check_refusals`' `rj_memory_not_enough` was accepted).
- The g++ warning `-Wno-format-overflow` is no longer given to clang (the `mac` preset, CMake with a compiler other
  than g++), which said once per object that it does not know it.
- `test_scaled_product` holds a component 2^-1200 below the other only where `long double` reaches that far (x87, quad;
  it is zero where `long double` is `double`, as on Apple Silicon), and the sign of a NaN result only where `long double`
  is x87's (clang -O3 on Apple Silicon gives the two assembly paths' NaNs different signs; every finite value keeps its
  bits).
- `check_numerics/run.sh`, `make_reference.sh` and `check_rdip_ladder/run.sh` run where `timeout` is missing (macOS),
  with perl's `alarm` in its place (exported, for the runs `xargs` starts).
- The checks run with macOS's `/bin/bash` 3.2: `check_input_forms/legacy.sh` names its systems with a `case` instead of
  an associative array (`declare -A`, bash 4), and `check_master_equation` and `check_refusals` keep their jump-operator
  lists in variables, because bash 3.2 brace-expands a list written inside `"$( ... )"` (`[{'a':1,'b':2}]` reached
  Python as `['a':1]`). `check_refusals`' `write_full_*` (result files linked to `/dev/full`) are skipped where there
  is no `/dev/full` (macOS), and `check_input_forms`' autoconv A10 (a tool writable by the owner's private group) where
  there is no `getent` to tell a private group (macOS, whose primary group `staff` is every user's: it ran there and
  failed, rightly refused).

## 2.0.0 — 2026-10-17 (changes since v1.2.2)

### Summary

- **A new input format.** The keys are grouped in blocks (`configuration`, `Qubit`, `Defect`, `hamiltonian`, `pulse`,
  `cluster`, `execution`, `output`, `frame`, `auto_convergence`); inputs may be `.jsonc` files with `//` comments; many
  numbers take a unit (`{"value": 5, "unit": "us"}`); `main.out -e` writes a commented template with every key. The
  command line is only `-f`, `-o DIR`, `-v`, `-e`, `-h`. Old inputs are refused with the new name, and a converter
  rewrites them (see "Updating old inputs").
- **Faster.** The test inputs run up to about 35 times faster than v1.2.2 (P1 cCCE2 production: 442 s → 22 s on 64
  ranks), mostly with unchanged results; MKL is used for large clusters when the build finds it; a CUDA backend runs cCCE, gCCE
  and pCCE on a GPU.
- **New methods and options.** ME-CCE and ME-gCCE (Lindblad dissipation), DEER pulses on bath spins, pCCE with gCCE and
  the master-equation methods, bundle clustering with an adaptive refinement, the automatic convergence of rbath, rdip,
  the time window and the number of configurations, several bath configurations from one input.
- **Physics fixes that change results**: `hfmedi` (complete second-order term), `knight` (nuclear Zeeman vector, |m_s|
  of each branch), the static Overhauser field of gCCE with sampled states, the h-BN tutorial hyperfine tensors (read
  twice too large), and a memory corruption of pCCE with state files.
- **Safer inputs.** A misspelt key, a file that does not open, a combination a method cannot do or a run that does not
  fit in memory stops before the calculation, with one message on one rank; a key that is read but not used is
  reported.
- **Tutorials and documentation.** Eleven tutorials in learning order (`example/tutorial/00_quickstart` …
  `10_coordinate_rotation`), each with `run.sh` and `expected.png`; a rewritten manual (`docs/user/`), method pages
  (`docs/methods/`) and a doxygen site.

### Updating old inputs

`python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites inputs in place, keeping their comments and
listing what it changed (`--check` lists only, `--stdout` prints the new text). A rewritten input gives the same results
as before. The table below lists every renamed or removed key; `docs/user/InputFile.md` (section "Old names") has the
same table with more detail.

| Old | New |
|---|---|
| **Method and cluster shape** | |
| `"method": "CCE"` | `"method": "cCCE"` |
| `"method": "pCCE"`, `"dsj"`, `"itb"`, `"dsjitb"` | `"method": "cCCE"` with `"cluster": {"mode": "pcce"}` (`"dsj"`, `"itb"`, `"dsjitb"`) |
| `"method": "pgCCE"`, `"pmeCCE"`, `"pmegCCE"` (adding_MECCE branch) | `gCCE`, `meCCE`, `megCCE` with `"cluster": {"mode": "pcce"}` |
| `cluster_mode`; its value `adaptive` | `cluster.mode`; `bundle` |
| `addsubclus`, `add_subclusters` | `cluster.selection.add_subclusters` (true/false) |
| `nk`, `max_clusters_per_order` | `cluster.selection.max_per_order` |
| `rbathcut`, `rdipcut` | `rbath: [inner, outer]`, `rdip: [inner, outer]` |
| **Bath configuration** | |
| `gyrofile`, `bathfile`, `seed` (top level) | `configuration.gyrofile`, `configuration.bathfile`, `configuration.seed` |
| `nstate` N (top level or in `configuration`) | `configuration.s_range: [1, N]`; `nstate` 0 (the ensemble): no `s_range` |
| a `seed` number with `k_range`, meaning the same seed for every k | a list of that number, one per configuration (a number is now the first configuration's seed) |
| `bathadjust`, `bath_shift` (top level) | `configuration.bath_shift` |
| `statefile`, `bath_state_file` (top level, a string) | `configuration.bath_state_file: [...]`, one per bath file (`{s}` the state) |
| `avaaxfile`, `Defect.axis_file` | `configuration.defect_axis_file: [...]`, one per bath file |
| `exstatefile`, `Defect.state_file` | `configuration.defect_state_file: [...]`, one per bath file |
| **Qubit** | |
| `qubitfile`, `Qubit.file`, `Qubit.spin`, `Qubit.name`, `qspin` | an entry of `Qubit.qubit`: `[{"file": ..., "spin": ..., "name": ...}]` (`file` in place of `xyz`) |
| `Qubit.nqubit` | removed: the length of `Qubit.qubit` |
| `qalphams`/`qbetams`, `Qubit.alphams`/`betams`, `Qubit.qubit[].alphams`/`betams` | `Qubit.alpha_state: {"ms": [m_1, m_2, ...]}`, `Qubit.beta_state: {"ms": [...]}` (one m_s per qubit) |
| `Qubit.alphaidx`, `Qubit.betaidx` | `Qubit.alpha_state: {"index": n}`, `Qubit.beta_state: {"index": n}` |
| `Qubit.psia`, `Qubit.psib` | `Qubit.alpha_state: {"vector": [...]}`, `Qubit.beta_state: {"vector": [...]}` |
| `qzfs`, `Qubit.zfs` (kHz) and its `Qubit.frame` | a `Qubit.intmap` self entry `{"between": ["q0", "q0"], "tensor": ..., "unit": "kHz", "frame": ...}` |
| `Qubit.qubit[].detuning: [x]` | `x` (a number or `{"value", "unit"}`) |
| `Qubit.overhaus` | `hamiltonian.common.qubit_bath_mean_field` |
| **Hamiltonian** | |
| `vector_mean_field` false / true | `hamiltonian.common.bath_bath_mean_field` `"z"` / `"xyz"` (`"off"` is new) |
| `hf_readmode` 0 / 1 / 2 / 3 | `hamiltonian.common.hyperfine`: `dipole` `"point"` / `"point"` / `"dft"` / `"dft"`, `contact` false / true / false / true |
| `hf_tensorfile`, `hf_ignore_oor` 0 / 1 | `hamiltonian.common.hyperfine.file`, `.missing_row` `"stop"` / `"point_dipole"` |
| `qd_readmode` 0; 1; 2 / 3 / 4 | `hamiltonian.common.quadrupole: "off"`; not available; `.structure` `"plain"` / `"strained"` / `"bubble"` |
| `qd_tensorfile`, `qd_tensorfile_woqubit`, `qd_cellpara` | `hamiltonian.common.quadrupole.file`, `.file_without_defect`, `.cell` |
| `knight`, `hfmedi` | `hamiltonian.cCCE.knight`, `hamiltonian.cCCE.hfmedi` |
| `jump_operators` (top level, adding_MECCE branch) | `hamiltonian.meCCE.jump_operators` or `hamiltonian.megCCE.jump_operators` (not rewritten by the converter) |
| jump `operator` `"+-"` / `"-+"` | `"-"` / `"+"` |
| `DefectTotSpin` | removed: the qubit spin is used |
| `CorrTotSpin` | removed: remake old-format VASP tensor files with `example/tools/SpinTensorReader` |
| `hf_cutoff` | removed |
| `magnetic_field_axis` | removed: write `bfield` as [0, 0, Bz] and turn the bath with `frame.rotation` |
| **Pulse** | |
| `npulse`, `pulsename`, `sequence`, `pulseiter`, `pulse_time`, `detuning_factor` | `pulse.count`, `pulse.name`, `pulse.sequence`, `pulse.per_qubit`, `pulse.duration`, `pulse.detuning_factor` |
| `bnpulse`, `bsequence`, `bpulse_defect`, `bpulse_energy_shift` (adding_DEER branch) | `pulse.bath.count`, `.sequence`, `.defect`, `.energy_shift` (not rewritten by the converter) |
| `bspin`, `pulse.bath.spin` | removed: the driven spins' own S |
| `balphams`/`bbetams`, `pulse.bath.alphams`/`betams` | `pulse.bath.levels: [m_a, m_b]` |
| **Defect** | |
| `"Defect": [ ... ]` | `"Defect": {"species": [ ... ]}` |
| `dfname`, `eqs`, `rxyzs`, `hypf` | `name`, `quadrupole_moments`, `nuclear_positions`, `hyperfine` |
| `naddspin`, `n_nuclear_spins` | `n_onsite_spins` |
| `navaax`, `n_configurations` | `n_axis_states` |
| `apprx` | removed (it changed nothing) |
| `coordinate_frame` | `frame` |
| species `axis` | `zfs.axis`: `"zfs": {"D": ..., "axis": [...], "unit": ...}` |
| `zfs.E` (0 only) | removed: `zfs.tensor` for E ≠ 0 |
| a unit in another case (`"mhz"`) | the exact spelling (`"MHz"`) |
| **pCCE and bundles** | |
| `sK`, `max_trial`, `max_iter`, `kmeans_pp` | `cluster.pcce.partition_size`, `.trials`, `.max_iterations`, `.kmeans_pp` |
| `iter_detail`, `pcce_detail_limit`, `print_trials` true/false + `print_limit` | `cluster.pcce.print_trials`: the number of trials printed |
| `pcce_partition`, `pcce_connections`, `pcce_cache_mode`, `pcce_cache_file` | `cluster.pcce.partition`, `.connections`, `.cache.mode`, `.cache.file` |
| `bundle_max_dim`, `bundle_merge_cutoff`, `bundle_connections`, `bundle_connection_budget` | `cluster.bundle.max_dim`, `.merge_cutoff`, `.connections`, `.connection_budget` |
| `expansion_max_clusters`, `expansion_max_work` | `cluster.bundle.max_clusters`, `.max_work` |
| top-level `pcce`, `bundle` blocks | `cluster.pcce`, `cluster.bundle` |
| `Adaptive`, `Refine` or top-level `adaptive` block | `cluster.adaptive` |
| `adaptive_work_file`, `probe_memory_gb` | `cluster.adaptive.work_file`, `cluster.adaptive.probe_memory_gb` |
| `cluster.adaptive.output` | `cluster.adaptive.workdir` |
| adaptive keys with `-` (`max-rounds`) | `_` (`max_rounds`) |
| `cluster.adaptive.policy` | removed (one policy) |
| **Execution** | |
| `backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `column_evolution`, `profile` (top level) | the same names in `execution` |
| `gpu_min_dim`, `gpu_batch_clusters`, `gpu_batch_steps`, `gpu_cluster_batch_mb`, `gpu_cache_mb` | `execution.gpu_options.min_dim`, `.evolution_batch_size`, `.evolution_batch_steps`, `.cluster_batch_mb`, `.cache_mb` |
| `gcce_refinement` true/false, `gcce_refinement_max_sweeps` | `execution.gcce_refinement: {"enabled": ..., "max_sweeps": ...}` |
| **Output** | |
| `outfile` | `output.file` |
| `savemode` `normal` / `info` / `all` / `allfull` | `output.savemode` `coherence` / `coherence_tagged` / `cluster_coherence` / `cluster_density_matrix` |
| `save_nodiv` | `output.savemode: {"mode": ..., "nodiv": ...}` |
| `cluster_contributions` true/false, `cluster_contribution_max_mb` | `output.cluster_contributions: {"enabled": ..., "max_mb": ...}` |
| `cluster_only` (top level) | `output.cluster_only` |
| **Frame** | |
| `coordinate_frame_rotation`, `bath_coordinate_rotation` | `frame.rotation`, on when the block is given (a block that was off is removed) |
| `hf_tensor_frame`, `qd_tensor_frame` | `hamiltonian.common.hyperfine.frame`, `hamiltonian.common.quadrupole.frame` |
| `Qubit.intmap[].tensor_frame` | `Qubit.intmap[].frame` |
| frame value `"bath"` | `"file"` |
| no frame with the rotation on (a Defect species, a single qubit's intmap entry, the qubit ZFS) | `"frame": "qubit"` (the old hidden default; now required) |
| `bath_axis`, `reference_qubit`, `qubit_position_frame`, `defect_axis_reference` | removed (x is `qubit_axis` × [0, 0, 1], the origin the first qubit) |
| **Automatic convergence** | |
| `auto_convergence.configurations` | `auto_convergence.sampling`; its file names are the `configuration` block's |
| `converge` (each section) | removed: a section whose value is given is fixed |
| `time_window.unit`, `time_window.threshold` | `time_window.value: {"value", "unit"}`, `time_window.decay_level` |
| `sampling.number`, `sampling.test_number` | `sampling.configurations`, `sampling.test_configurations` |
| `launcher`; `output.path`, `output.file` with `auto_convergence` | removed (the runs use `mpirun -n N` of `main.out` and go to `auto_convergence.workdir`) |
| **Command line** | |
| `-I`, `-N`, `-m`, `-q`, `-B` | `configuration.bathfile`, `configuration.s_range`, `method`, `quantity`, `bfield` |
| `-s`, `-S`, `-a` | `configuration.bath_state_file`, `defect_state_file`, `defect_axis_file` |
| `-P` | `output.cluster_only` |
| `--auto_conv`, `--adaptive` | the `auto_convergence` and `cluster.adaptive` blocks |
| `--auto_status` | `python -m autoconv status FILE` in `example/tools/AutoConvergence` |
| `-o FILE` (an output name) | `-o DIR` is the output folder; the name is `output.file` |

### Input and command line

- **Blocks.** `configuration` {`path`, `k_range`, `gyrofile`, `bathfile`, `bath_shift`, `bath_state_file`,
  `defect_axis_file`, `defect_state_file`, `s_range`, `seed`}; `hamiltonian` {`common`: `qubit_bath_mean_field`,
  `bath_bath_mean_field`, `hyperfine` {`dipole`, `contact`, `file`, `missing_row`, `frame`}, `quadrupole`; one section
  per method: `cCCE` {`knight`, `hfmedi`, `delta_A`}, `meCCE` and `megCCE` {`jump_operators`}, `gCCE` and `noise`
  empty}; `cluster` {`mode`, `selection` {`add_subclusters`, `max_per_order`}, `pcce`, `bundle`, `adaptive`};
  `execution` {`backend`, `gpu_options`, `propagator`, `evolution`, `evolution_update`, `gcce_refinement`,
  `column_evolution`, `profile`}; `output` {`path`, `file`, `savemode`, `cluster_contributions`, `cluster_only`};
  `frame` {`rotation`}. `order`, `rbath`, `rdip`, `quantity`, `bfield`, `deltat` and `nstep` stay at the top.
- `method` is `cCCE`, `gCCE`, `meCCE` or `megCCE`; the cluster shape is `cluster.mode` (`standard`, `pcce`, `bundle`,
  `dsj`, `itb`, `dsjitb`). `noise` is a reserved method name, refused as not implemented yet.
- Numeric codes became words and 0/1 switches became true/false.
- **The command line** is `-f FILE`, `-o DIR`, `-v`, `-e` and `-h`/`--help`; every other setting is an input key.
  `-o DIR` names the output folder (made when missing) and takes the place of `output.path`.
- **`main.out -e`** writes the input template `ccex_example.jsonc` (in the working directory, or in `DIR` with
  `-o DIR`; an existing file is not overwritten) and stops. Its active part is the setting of
  `example/tutorial/00_quickstart` with a `bathfile` you fill in (it does not run until it names a bath file) and
  `gyrofile` commented (the built-in table); every other key is a commented line with its allowed values and default,
  and each section points to its manual page and a tutorial that uses it.
- **`-h`** is a short usage: one line per option and where the input keys and tutorials are.
- **`//` comments and `.jsonc` inputs.** `//` outside a JSON string starts a comment to the end of the line, as `#` and
  `!` do (`"http://x"` in a string stays as written), in CCEX, the adaptive refinement, the automatic convergence and the
  converter. An input may be named `.jsonc`; the default output name drops `.jsonc` as it drops `.json`.
- **A unit per key**: `deltat`, `bfield`, `rbath`, `rdip`, `configuration.bath_shift`, `Qubit.qubit[].xyz`,
  `Qubit.qubit[].detuning`, `Qubit.qubit[].gyro`, `pulse.duration`, `pulse.bath.energy_shift`, `pulse.bath.bandwidth`,
  `hamiltonian.<method>.jump_operators[].rate`, `hamiltonian.common.quadrupole.cell`, `cluster.bundle.merge_cutoff` and
  `auto_convergence.time_window.value` take `{"value": ..., "unit": "..."}` (a number or a list) beside the plain
  number, which keeps its unit. Time `s`, `ms`, `us`, `ns`; field `G`, `mT`, `T`; length `A`, `nm`; frequency `Hz`,
  `kHz`, `MHz`, `GHz`; rates `rad/s`, `rad/ms`, `rad/us`. Units are case-sensitive. One exact operation per value (×
  or / a power of ten), so `{"value": 5, "unit": "us"}` is `0.005` ms bit for bit. A plain `Qubit.qubit[].detuning` or
  intmap tensor (kHz) is noted in the log, beside the Defect's MHz. See `docs/user/InputFile.md#Units`.
- **Several bath configurations in one run**: `configuration.k_range: [first, last]` runs every k as a configuration of
  its own, one after another, with `{k}` in the file names (`bathfile`, `bath_state_file`, `defect_axis_file`,
  `defect_state_file`, `cluster.pcce.cache.file`) replaced by k (`{{k}}` is a literal `{k}`). Each writes
  `<output.file>_k<k>_...` (no average), byte for byte what its files give alone. `configuration.path` is the folder of
  the block's file names.
- **Sampled states**: `configuration.s_range: [first, last]` (absent: the ensemble). States drawn from the seed are
  drawn for 1 .. last and calculated for first .. last, so [11, 20] continues [1, 10]. `{s}` in `bath_state_file` and
  `defect_state_file` names the state (without it the state is appended).
- **Seeds**: a `configuration.seed` number is the first configuration's seed (k uses seed + k − first); a list gives one
  per configuration. Without a seed, configuration k uses the generated seed + (k − first), printed per k in the log
  (`[seed] configuration k=...: seed S`), so every configuration draws its own bath states.
- **The state and axis files are lists, one per bath file**, read in the order of the bath files into the merged bath;
  `""` is a bath file without such a file (drawn from the seed). The defect state file takes several on-site spins per
  line. The log names a merge of bath files (`Bath = fileA (n spins) + fileB (m spins)`).
- **The qubit states in one place**: `Qubit.alpha_state` and `Qubit.beta_state` are each exactly one of
  `{"ms": [m_1, m_2, ...]}` (one m_s per qubit in the order of `Qubit.qubit`), `{"index": n}` or `{"vector": [...]}`.
  Vector and `psi0` components may be complex, `[re, im]`. An `{"index": n}` eigenvector is multiplied by the phase that
  makes its largest component real and positive (the first of equal ones), so psi0 = |α⟩ + |β⟩ no longer depends on
  the eigensolver's phase (it matters for gCCE with E ≠ 0).
- **One qubit from a file** is an entry of `Qubit.qubit`, `{"file": ..., "spin": ..., "name": ...}` in place of `xyz`
  (name default `q0`). A qubit's zero-field splitting is a `Qubit.intmap` self entry in both forms.
- **Frames**: `frame.rotation` holds `qubit_axis` and `enabled` only and is on when the block is given (`enabled`
  defaults to true). Each set of directional numbers has a key `frame`, `"file"` (written in the bath file's frame,
  turned with the positions) or `"qubit"` (written in the computational frame, z = `qubit_axis`, used as written):
  `hamiltonian.common.hyperfine.frame`, `hamiltonian.common.quadrupole.frame`, `Qubit.intmap[].frame`,
  `Defect.species[].frame`. With the rotation on every such place needs its `frame`; with it off a written `frame` is
  reported as unused. The computational x is `qubit_axis` × [0, 0, 1] and the origin the first qubit. `bfield` is read in
  the computational frame, so `[0, 0, B]` is along `qubit_axis` in the lab.
- **Defect**: `zfs` is `{"D", "axis", "unit"}` or `{"tensor", "unit"}`; `frame` is needed only for a species with
  `nuclear_positions`, `hyperfine`, `efg` or `zfs`; `{"value", "unit"}` is taken beside `{"values", "unit"}`.
- **pCCE and bundles**: `cluster.pcce.print_trials` is the number of trials printed; `cluster.bundle.connections` takes
  `checked` (the mode of a bundle input without rdip), and the log names the mode; `cluster.adaptive.workdir` holds the
  refinement's steps, and adaptive keys are written with `_` only.
- **A file the input names has to open**: a missing state, axis, defect state, tensor, gyro, qubit, bath or pCCE cache
  file stops the run before the calculation, naming the key and the file (it used to be a warning, and states and axes
  were drawn from the seed). `""` is still no file.
- **A key the run reads and does not use** gets one line: `Warning: input key "..." is not used by this run: ...`
  (`psi0` with cCCE, the mean fields in the ensemble, state files in the ensemble, `gpu_options` in `main.out`, unused
  `cluster.selection.max_per_order` entries 0–2, ...).
- **Error messages** name the keys by their current paths; a refused input stops all MPI ranks together with one message
  on rank 0 and exit code 1.
- **Merged bath files**: spins of different bath files within 0.5 Å of each other give a warning naming the files and
  positions; two at the same position (closer than 1e-6 Å), whose coupling is infinite, are refused before the
  calculation.
- **Memory check before the calculation** (no input key): the peak memory of a rank is estimated from its largest cluster
  on the path the run takes (dense matrices 16 D², meCCE/megCCE 16 D⁴, the column blocks, the dm buffers,
  `cluster_density_matrix`) plus what it holds, and compared with what a rank of the node may use
  (`SLURM_MEM_PER_NODE`, `SLURM_MEM_PER_CPU`, else `MemAvailable`, shared by the node's ranks). One log line; a run that
  does not fit is refused before the first cluster, with the cluster dimension, both sizes and what to change.
- **An output folder that cannot be written** is refused at the start (rank 0 writes and removes a file there).
- `output.path`: the folder of the outputs, made when missing (`-o DIR` > `output.path` > the working directory). The
  output name is `output.file`, else the input file's name without `.json`/`.jsonc`. An absolute `output.file` with
  `output.path` or `-o` is refused, and so is `output.path` with the adaptive refinement.
- In input files, `#`, `!` and `//` start a comment only outside JSON strings. Keys are read without regard to case,
  also in the `savemode`, `propagator` and `adaptive` blocks.

### Calculation and physics

- **ME-CCE and ME-gCCE** (`method` `meCCE`, `megCCE`; `hamiltonian.meCCE.jump_operators`,
  `hamiltonian.megCCE.jump_operators`): Lindblad dissipation of bath spins (both) and of the qubit (megCCE only), on the
  CPU. meCCE refuses qubit jump operators. Method: Onizhuk et al., Phys. Rev. Lett. 132, 250401 (2024). See
  `docs/methods/master-equation.md`.
- **DEER bath pulses** (the block `pulse.bath`: `count`, `sequence`, `defect`, `energy_shift`, `levels`, `bandwidth`)
  for gCCE (matrix and vector, `exact` and `recursive`) and megCCE (U_q ⊗ U_b in its pulse superoperator), on the CPU,
  with `cluster.mode` standard or pcce. `pulse.bath.bandwidth` (MHz) is the half width of the frequency window of the
  bath pulse; its default is the fixed value used before (10 rad·kHz). See `docs/methods/deer.md`.
- **pCCE with gCCE, meCCE and megCCE** (`cluster.mode` pcce): the pCCE partitions with the gCCE and master-equation
  kernels, on the CPU; `cluster.pcce.cache` stays with cCCE.
- **pCCE**: a deterministic partition search, a partition cache (`cluster.pcce.cache` {`file`, `mode`}) and
  `cluster.pcce.connections` (`centers`, `spins`). When the bath is not divisible by `partition_size`, the log says
  which spins are left out (`pCCE spins left out : ...`).
- **Bundle clustering and the adaptive refinement** (`cluster.mode` bundle with a `cluster.adaptive` block): bundles of
  strongly coupled spins with the merge policy, and `"regroup": true` to regroup where a bundle expansion gives
  |L| > 1. See `docs/user/InputFile_Cluster.md`.
- **`hamiltonian.cCCE.delta_A: true`** (default false): adds to each bath spin the second-order hyperfine term
  δA(m_s)·I = −m_s (A_x × A_y)·I / (2Δ) of Maze, Taylor, Lukin, Phys. Rev. B 78, 094303 (2008), Eqs. (A3), (A4), beside
  `knight` (δg) and `hfmedi`. The coefficient −m_s differs from Maze's (2 − 3|m_s|): a direct derivation, time-reversal
  symmetry and exact diagonalization give −m_s. Refused with `hfmedi` (which already contains it), without one spin-1
  qubit with a zero-field splitting, and with `propagator` sparse; a warning when the field across the qubit axis mixes
  the levels. See `docs/user/InputFile_General.md#tag_delta_A`.
- **`hamiltonian.common.bath_bath_mean_field: "off"`** (sampled states): leaves out the static field of the bath spins
  outside a cluster (`"z"` stays the default, `"xyz"` is the former `vector_mean_field: true`).
- **Built-in isotope table**: `configuration.gyrofile` is optional; species it does not name come from the table.
- **Hyperfine tensor files, format v3**: SpinTensorReader writes line 1 `VASP` / `QUANTUMESPRESSO` and line 2
  `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`,
  and the numbers are the hyperfine of the spin-S defect for both codes (QE's S = 1/2 coupling × 1/(2S), `--spin S`
  required; VASP's OUTCAR as it is, S from its magnetization). CCEX reads such a file × 1. A file without the line is
  read by the earlier rule (first line `VASP` → × 1, otherwise × 0.5/S), with a warning; its results are unchanged.
  See `docs/user/InputFile_General.md#tag_hyperfine_file`.
- **The automatic convergence**: an `auto_convergence` block in the input (`mpirun -n N main.out -f input.jsonc`;
  `"enabled": false` runs the ordinary calculation) chooses rbath, rdip, the time window and the number of bath
  configurations and states, then runs the production. Neighbouring candidates are compared on the area under a T2 fit
  of |L|; where the fit fails, on the area under |L| itself (same 2 %), and the summary says so. A section whose value
  is given is fixed, one without it is found. See `docs/user/InputFile_AutoConvergence.md`.
- `rdip_ladder`, an internal key of the automatic convergence: one order-2 run (cCCE, gCCE or pCCE) also writes the
  curves of smaller rdip values.
- **`pulse.name`** is Ramsey, HahnEcho, CPMG, WAHUHA or Manual and must fit `pulse.count`; WAHUHA is the four-pulse
  sequence (`count` 4).
- The cost of a cluster in the split over the ranks counts the qubit dimension for megCCE, pgCCE and pmegCCE as for
  gCCE (results unchanged).

### Speed

Unless noted, the results are byte-identical to those before each change.

- **The CPU paths**: the cCCE and gCCE product assembly, the state preparation, fixed-size paths for small clusters, the
  order-2 cluster list built without duplicated work, the small-cluster Hamiltonians, and the pCCE partition search (up
  to 9× faster on large baths).
- **Small cCCE clusters batched on the CPU**: clusters of bath dimension 2 and 4 are evolved 4 (AVX2) or 8 (AVX-512) at
  a time, in the builds whose arithmetic the batched kernel reproduces (g++ 13 AVX2, icpx 2026.1 AVX2, icpx 2023.2
  AVX-512); one state with any pulse plan and the ensemble Hahn echo. Each process first checks the batch against the
  one-by-one path on a few made-up clusters (about 1 ms) and turns it off, with one `[cpu-batch] off` line, where a bit
  differs (e.g. a build with `-ffp-contract=off`, or another CPU's `-march=native` tuning); results are the one-by-one
  path's everywhere. Its window storage is made only when a window is used, and a time grid whose coherences for one
  group of lanes pass 4 MiB runs one by one. Production P1 cCCE2 on 64 ranks: 0.60–0.64 of the time with `recursive`,
  0.74–0.77 with `exact`.
- **`evolution_update`**: `exact` or `recursive` (the phases by a recurrence; 1.0–1.6× faster on the tested systems,
  within 1e-9 of `exact` at |L| ≤ 1, not bit-identical). The default is `recursive` where it exists (the CPU cCCE,
  single state and ensemble, the packed GPU kernel, and the CPU gCCE with `propagator` `eigen`), `exact` elsewhere (the
  GPU large-cluster path, gCCE `expm`); the run log says which. `"evolution_update": "exact"` gives the earlier bits.
- **Clusters of an order are split over the ranks by cost when their sizes differ** (dimension D, cost D²(D + nstep);
  each rank takes a contiguous run of the list): mixed baths such as h-BN run 1.6–3× faster on 64 ranks (cCCE2 1.6×,
  gCCE2 2.3×, cCCE3 3.1×). Orders whose clusters all have one size (P1, ¹³C) keep the equal-count split, byte for
  byte; in a mixed bath the last printed digit can change (h-BN gCCE2 dm: 5 of 36,000 numbers, 1e-10).
- **Each rank keeps only its own order-2 cluster rows** (cCCE/gCCE; `cluster_only`, `cluster_contributions` and adaptive
  work keep the full list): peak memory per rank on the P1 production CR60 σ4 bath 334 → 181 MiB, same speed.
- **MKL for large clusters, when the build finds it** (`CCEX_BLAS` auto by default, `off` without it, `mkl` to insist;
  the build and the run log say on or off): the cCCE ensemble Hahn step's D³ product from bath dimension 21 by zgemm and
  the cCCE/gCCE cluster eigensolve from dimension 63 by zheevd (MKL lp64, sequential). 64 ranks: h-BN cCCE2 ensemble
  0.80–0.88, gCCE3 0.72–0.80, gCCE2 0.98 of the time; P1 production inputs unchanged. The curves differ in the last
  digits (cCCE ≤ 2e-10, h-BN gCCE late times up to 9e-8).
- **`execution.column_evolution` for large one-state clusters** (on by default, `min_dim` 256, `columns` 64;
  `{"enabled": false}` turns it off): in the single-state cCCE/pCCE evolution of clusters of dimension ≥ `min_dim`, the
  time points are carried as matrix columns, so each pulse interval is one matrix product instead of one per time point.
  Curves agree with the path before to about 2e-15 relative, not bit-identical. 64 ranks, pCCE sK5/sK3o3/sK4 (clusters
  of dimension 1024/512/256, Hahn and CPMG-8, nstep 200 and 1000; 12 inputs): with MKL 0.40 of the time (geometric mean;
  sK5 CPMG-8 nstep 1000 114.3 s → 13.4 s), with Eigen alone 0.61; peak memory per rank −8 to +11 MB.
- **`propagator: "sparse"`** (opt-in; also `{"method": "sparse", "min_dim": 4096, "tol": 1e-12}`): one sampled state of
  cCCE and gCCE `vector` without diagonalizing the cluster Hamiltonian (Krylov action on the state), for clusters of
  dimension ≥ `min_dim`. It pays only where the dense eigendecomposition does not fit in memory (dim 16384: about 36 GB
  dense against about 0.2 GB sparse); it is much slower at high field or with a wide spectrum and long times (P1 at
  500 G, dim 4096, 151 points: about 18 h against 5 min dense, extrapolated from the first 11 points). cCCE accepts
  `propagator` now (`eigen`, `sparse`; `expm` stays gCCE).
- megCCE builds the pulse superoperators of a cluster once per distinct pulse instead of at every time point, and the
  DEER bath pulse computes its rotation once per sequence row instead of once per driven spin.
- **The GPU backend** (`execution.backend` `cpu`, `cuda`, `auto`; `execution.gpu_options`) for cCCE, gCCE and pCCE in
  `bin/main_gpu.out`, with `propagator` (`eigen`, `expm`), `evolution` (`vector`, `matrix`) and `gcce_refinement` for
  gCCE. See `docs/methods/gpu.md`.
  - **Packed kernel for small cCCE clusters** (`gpu_options.packed_max_dim`, 4 by default with `backend` `cuda` or
    `auto`, 2 or 4 may be given): every cCCE cluster of bath dimension 2 and 4 is computed by one GPU thread, thousands
    per kernel, for the whole node through one CUDA context (node rank 0 feeds the GPU from node shared memory); single
    state and ensemble, Hahn, CPMG, Ramsey, `exact` and `recursive`. Within 1e-8 + 2e-6 |CPU| of the CPU, not
    bit-identical. With 32 ranks and one GPU the production P1 cCCE2 runs took 0.59–0.84 of the CPU's time; a run of a
    few seconds gains nothing (CR0 ensemble, 1.8 s: 0.95). With `packed_max_dim` not given and a run the CPU batch takes, it is on for one state
    only from 40000 such clusters per rank per state and off for the ensemble (break-even measured at about 39.6 k, 32
    ranks); the clusters it does not take go to the CPU batch, and runs where it is off write the bytes of
    `bin/main.out`.
  - **Warp kernel for gCCE and cCCE clusters up to dimension 27** (`gpu_options.warp_max_dim`, 27 by default, 0–32 may
    be given): one cluster per GPU warp (shared-memory Jacobi, the lanes over the time points, exact phases); gCCE with a
    qubit of dimension 3 (total dimension 6–27), cCCE bath dimension 8 and 16, single state and ensemble, coherence and
    dm, through the packed kernel's node feeder. Within 1e-8 + 2e-6 |CPU| outside the stretch where a product diverges,
    not bit-identical. Runs it does not take (`gcce_refinement`, `per_qubit`, `cluster_density_matrix`, ...) switch it
    off by themselves. With 32 ranks + MPS and one GPU, P1 gCCE2 took 0.60–0.69 and h-BN gCCE2 0.86 of the time of the
    GPU path without it.
  - **Short runs leave both kernels off**: with neither `packed_max_dim` nor `warp_max_dim` given, a run with fewer than
    20000 cluster evaluations (clusters × states) per rank does not use the card and writes the bytes of the CPU path,
    with one log line. The early CUDA and shared-memory start is skipped where the input alone rules the kernels out,
    and the warp kernel's shared-memory region is allocated only when some rank has a warp cluster (P1 CR0 cCCE2 had
    been 4.7 % slower from it).
  - **`gpu_options.min_dim` defaults to 768**: on an RTX PRO 6000 (32 ranks) the GPU was faster than the CPU from a
    cluster dimension of about 670–750 on, on every path, and slower at 648 and below; clusters of dimension below 768
    run on the CPU with `backend` `auto` (their last digits can change within the GPU tolerance; a cCCE4 ensemble case
    ran 22.3 → 9.6 s). An explicit `min_dim` keeps its value.
  - **`backend` `auto` does not stop on a GPU that cannot be opened or runs out of memory**: the packed kernel is turned
    off when the GPU cannot be opened, and a cluster that runs out of GPU memory is computed on the CPU, with every later
    cluster of its dimension or larger on that rank (`[gpu-memory-fallback]` in the log). `backend` `cuda` still stops.
  - With `configuration.k_range` a GPU run keeps its CUDA context from one configuration to the next.

Measured during development, each row against the development code just before that change (same inputs and server,
64 ranks unless noted; byte-identical results, or within the stated bound for `recursive`):

| What | Case | Before → after | Speed-up |
|---|---|---|---|
| P1 production cCCE2, 30 states | CR60 σ4 | 229 s → 85.9 s | 2.7× |
| | CR0 | 21.3 s → 16.5 s | 1.3× |
| P1 production gCCE2, 30 states | CR60 σ4 | 313.3 s → 224.7 s | 1.4× |
| | CR0 | 56.1 s → 41.2 s | 1.4× |
| pCCE partition (P1, 30,000 trials) | rbath 483 Å | 1,616 s → 173 s | 9.3× |
| Automatic convergence, h-BN gCCE2 ensemble | whole run | 14,406 s → 4,281 s (CCEX runs 596 → 103) | 3.4× |
| Automatic convergence, pCCE runs differing only in rdip | NV-¹³C, 4 ranks (workstation) | 8.2 s → 0.3 s | 27× |
| Product assembly, P1 production cCCE2, 30 states, `recursive` | CR60 σ4 | 55.7 s → 45.7 s (assembly per rank 12.9 s → 4.9 s) | 1.2× |
| | CR0 | 10.8 s → 9.4 s (assembly per rank 2.3 s → 0.9 s) | 1.15× |
| Product assembly, P1 production gCCE2, 30 states, `recursive` | CR60 σ4 | 167.9 s → 162.2 s (assembly per rank 17.3 s → 10.0 s) | 1.04× |
| | CR0 | 31.0 s → 29.9 s (assembly per rank 3.1 s → 1.8 s) | 1.04× |
| `evolution_update: recursive` | 16 inputs (NV-¹³C, h-BN, P1) | against `exact` | 1.0–1.6× |

One later set of changes (commit f9b0ba8f3 against ccc16e507), on the GPU node (Xeon 6731P + RTX PRO 6000, 32 ranks,
MPS; `do_compile.sh gpu`, so the later build has MKL), one job, median of 3 (AB BA AB). CPU = `bin/main.out`, GPU =
`bin/main_gpu.out` with `backend: auto`; times in seconds. Curves: the CPU build byte-identical (P1 gCCE2) or in the last
printed digit (cCCE, 1e-10 to 8e-10), h-BN gCCE dm within 1.5e-8 relative; the GPU P1 gCCE2 difference is the warp kernel
against the CPU path (relative 1.3e-6, 1.6e-5 on the diverging CR60 curve).

| Case | CPU before → after | CPU speed-up | GPU before → after | GPU speed-up |
|---|---|---|---|---|
| P1 production cCCE2, 30 states, CR0 | 11.24 → 6.24 | 1.8× | 8.22 → 6.46 | 1.3× |
| P1 production cCCE2, ensemble, CR0 | 1.63 → 1.20 | 1.4× | 1.49 → 1.24 | 1.2× |
| P1 production cCCE2, 30 states, CR60 σ4 | 59.11 → 30.73 | 1.9× | 27.11 → 27.45 | 0.99× |
| P1 production cCCE2, ensemble, CR60 σ4 | 5.72 → 2.93 | 2.0× | 3.23 → 3.07 | 1.05× |
| P1 production gCCE2, 30 states, CR0 | 39.79 → 39.93 | 1.0× | 39.69 → 22.34 | 1.8× |
| P1 production gCCE2, 30 states, CR60 σ4 | 219.59 → 220.29 | 1.0× | 218.51 → 105.57 | 2.1× |
| h-BN cCCE2 ensemble Hahn | 14.28 → 10.70 | 1.3× | 14.65 → 10.74 | 1.4× |
| h-BN gCCE2 dm, 10 states | 10.61 → 7.30 | 1.5× | 10.71 → 7.62 | 1.4× |
| h-BN gCCE3 dm, 1 state | 17.97 → 10.56 | 1.7× | 17.97 → 10.98 | 1.6× |
| NV-¹³C cCCE2 Hahn (rbath 70 Å, rdip 15 Å) | 1.49 → 1.43 | 1.04× | 1.76 → 1.50 | 1.2× |
| h-BN cCCE2 ensemble Ramsey | 14.87 → 4.42 | 3.4× | 15.15 → 4.60 | 3.3× |
| Geometric mean of the 11 | | 1.53× | | 1.47× |

### Output and log

- **`output.savemode`** is a string or a block `{"mode": ..., "nodiv": true|false}` with new names: `coherence` (the
  default), `coherence_tagged`, `cluster_coherence`, `cluster_density_matrix`; `nodiv` false writes `_wiDiv` alone.
  The files written are the same as before for the same setting. See `docs/user/InputFile_Output.md#tag_savemode`.
- **`output.cluster_only: true`**: writes the cluster plan `<output.file>.clusters.json` and no dynamics.
- **`output.cluster_contributions`**: a table of each cluster's contribution; its header line carries `seed=` only when
  the run draws from a seed.
- **A shorter run log** (results unchanged): the progress of a state is at most one line per tenth of rank 0's clusters
  (`At rank 0, N of M clusters (x%) are computed at t s`), and every state ends with
  `State done : i of n, C clusters, t s` (`ensemble` for the ensemble). The bath states, mean fields, the steps of a
  state and their wall times are printed with `-v` only; the qubit states (psia, psib, psi0) for the first state of the
  run. A configuration of `k_range` after the first prints only the files it read, the spin counts, the tagged lines and
  warnings (everything with `-v`), and every run prints `Clusters per order : 1: n1, 2: n2 (total, bath spins, ranks)`
  before the calculation. Measured: a P1 cCCE2 run, 3 configurations × 30 states on 2 ranks, 321 KB → 121 KB; a 1.2 MB
  production log (10 states, 92 % progress lines) comes to about 42 KB by extrapolation.
- **A `Build :` line in the banner**: `Build : CPU | <compiler version> | <MPI library> | MKL on/off | preset <name> |
  commit <hash>[-dirty]`, or `Build : GPU (CUDA x.y, sm_NN) | ...` for `main_gpu.out`. The compiler from its macros, the
  MPI from `MPI_Get_library_version`, the preset from `do_compile.sh` (`cmake` from CMakeLists.txt), the commit from a
  header both builds write at build time (`-dirty` when the tree differs from it, `unknown` outside git).
- `cluster.selection.max_per_order` entries 1 and 2 (which are not used: only orders 3 and up are cut) give one warning.

### Tools

- **`example/tools/InputKeys/convert_input_keys.py`** rewrites inputs written with the old keys (comments `#`, `!` and
  `//` kept, `.json` and `.jsonc`, CRLF line ends kept, the changed keys listed, any number of files; `--check`,
  `--stdout`). A key written twice with the same value is kept once; different values (or keys that differ only in case)
  stop the file before it is written; it takes JSON booleans only and writes through a temporary file. Tests: `example/tools/InputKeys/tests`.
- **`example/tools/StateFiles/make_state_files.py`** writes the bath state, axis and defect state files of an input
  (the automatic convergence uses the same rule).
- **`example/tools/AutoConvergence/`**: the runner of the automatic convergence; `python -m autoconv status FILE` shows
  where a run is; `analyze` reads the same `auto_convergence` block as the runner. Keys are read without regard to case,
  and two keys that differ only in case are refused.
- **SpinTensorReader** writes tensor format v3 (above); `--fidentify` is removed (the lines are always written).
- **Builds**: `do_compile.sh` presets (`wsl-gnu`, `wsl-icpx`, `idm4`, `gpu`, `nurion`, `mac`, `generic`) and a CMake
  build (`cmake -B build && cmake --build build`, everything in `build/`).
- `example/tutorial/make_natab_baths.sh` makes the natural-abundance ¹³C configurations of the tutorials with MPI
  (the configurations divided over `NP` ranks).

### Documentation and tutorials

- **Tutorials** in learning order, one folder per topic: `00_quickstart` (NV-¹³C cCCE2 Hahn echo, seconds),
  `01_ensemble_and_states`, `02_pulses`, `03_P1` (with the Vp bath, DEER and the mean fields), `04_hBN`, `05_multiNV`,
  `06_cCCE_vs_gCCE` (with meCCE), `07_pCCE_and_bundle`, `08_GPU`, `09_auto_convergence`, `10_coordinate_rotation`. One
  folder per calculation named `<method><order>_<pulse>_<approach>` (e.g. `cCCE2_Hahn_ensemble`), inputs `.jsonc` with
  `//` comments and only the keys that matter, a `run.sh`, an `expected.png` per calculation (`plot_coherence.py`), and
  `configuration.k_range` instead of loops in `run.sh`. `example/README.md` has the five-line start and "what should I
  use for my case"; `example/tutorial/README.md` maps every old folder to its new place.
- **The manual** (`docs/user/`): getting started (building, the first run, the run log, exit codes), one page per input
  block, output files, analysis, tools, the units, and the table of old names. **Method pages** (`docs/methods/`): CCE,
  pCCE, bundles and the adaptive refinement, DEER, the master equation, GPU, the coordinate rotation. The doxygen site
  (`cd docs/doxygen && doxygen Doxyfile`) holds the manual, the method pages, the tutorials with their inputs and figures,
  and the API.
- **Repository layout**: `src/` and `include/` in five parts (input, physical, clustering, calculation, common);
  `example/{tutorial,analysis,tools}`; `external/`; `test/{code_verification,cluster_verification,gpu_verification}`.

### Fixed

- **`hfmedi` (hyperfine-mediated interaction)** is the complete second-order correction of each branch,
  H_med^(α) = Σ_{β≠α} E_αβ E_αβ† / (ω_α − ω_β) with E_αβ = Σ_i ⟨α|S|β⟩·A_i·I_i: every intermediate eigenstate β of the
  qubit Hamiltonian (it used only the other qubit level, missing for example the |+1⟩ path of the NV |0⟩ branch), every
  ordered pair (i, j) (it summed i < j only, which was not Hermitian for a spin-1 qubit: the CPU used one triangle of it
  and the GPU refused it), and the i = j terms. An intermediate level degenerate with the branch (within 1e-12 of the
  spread of the levels) is skipped when the coupling vanishes (spin-1 at zero field: +1 and −1, which gave `nan`; or a
  zero hyperfine) and refused when it does not. Results with `hfmedi: true` change (a test case by up to 0.95); the h-BN tutorial input that uses it
  had an asymmetric qubit ZFS typo, now `[[-1080000,0,0],[0,-1180000,0],[0,0,2320000]]` kHz.
- **`knight` (Knight shift)** multiplies by the nuclear Zeeman vector −γ_n B instead of the bare field B, and the
  (3|m_s| − 2) factor takes |m_s| of each branch (it was always 1). Zero for a field along the qubit axis, as before; in
  a tilted field results with `knight: true` change (a test case by up to 0.98).
- **gCCE and megCCE with sampled bath states**: the qubit Hamiltonian subtracted the static Overhauser field of the
  cluster's own spins, Σ A_zz m S_z, also with `qubit_bath_mean_field` false, where no mean field had been added, so the
  qubit lost the static field of the spins in the cluster (no T2* phase in Ramsey). The subtraction now happens only
  with the mean field on; that case and ensemble runs are unchanged. Hahn and CPMG change slightly (up to 1.8e-3 in the
  coherence of the P1 test cases, 9.5e-2 for CPMG-8 at 20 us, 0.10 in an h-BN density-matrix element).
- **The h-BN tutorial and test hyperfine tensors were read 2× too large**: the file
  `hBN_VB_flat_ATensor_HSE_CellPBE_tag` was SpinTensorReader v1 data (the OUTCAR × 2) with a `VASP` first line, which
  CCEX reads × 1 (nearest ¹⁴N A_zz 97.3 MHz against 47 measured; the hyperfine-mediated term 4×). It is replaced by
  `hBN_VB_flat_ATensor_HSE_CellPBE_v3` (format v3, A_zz 48.6 MHz); the h-BN test cases change by up to 0.88, and the tutorial was rerun. SpinTensorReader example 03
  (MoS2) had the same mistake.
- **The h-BN `gCCE2_Hahn_OH` tutorial** wrote the qubit ZFS without a unit, so it was read in kHz (D 3.48 MHz instead
  of 3.48 GHz); it now says `"unit": "MHz"`, and its results changed.
- **pCCE with a bath not divisible by `partition_size`**: after pCCE dropped the last `nspin mod partition_size` spins,
  reading a bath-state file still wrote the states of the dropped spins past the end of the shortened bath (memory
  corruption; a crash or silently wrong memory). Completed runs compared before and after the fix gave the same bytes.
- **`configuration.seed`** had no effect on sampled states.
- **Ensemble cCCE**: the conjugation of the coherence for an even pulse count, and the half interval computed in double
  precision; `deltat` is read as a double (v1.2.2 rounded it to a float).
- **Density-matrix output**: a component that divides noise by zero is set to 0 and reported, instead of `nan`.
- **`savemode` `cluster_density_matrix`** (formerly `allfull`) wrote no result of the expansion, named each cluster's full
  density matrix by positions so that clusters of the same size overwrote one another, and with cCCE wrote nothing. It
  now writes what `cluster_coherence` writes and, with gCCE, each cluster's full density matrix to
  `<output.file>_<bath indices>_full[_state<n>]_wiDiv/noDiv`.
- **`pulse.name` WAHUHA** with a `pulse.count` other than 4 wrote past the sequence array (a crash for 0–3, a wrong curve
  for 5); it is refused now.
- **Bundle clustering** with `connections` `checked` (no rdip) crashed (segfault) when bundles of several sites formed.
- **`execution.backend` `auto` with `configuration.k_range`**: from the second configuration on, the ranks that had
  opened the GPU skipped a collective node split that the others entered, and the run hung.
- **A refused input under several MPI ranks** gave exit code 1 or 255 and zero, one or two messages at random; it now
  writes the message once on rank 0 and stops all ranks with exit code 1. A refused input while the GPU card was being
  opened in the background ended in `std::terminate` (exit 255); it now closes the card first.
- **The adaptive refinement**: the step inputs set `output.file`, `output.savemode` and `cluster.adaptive.work_file` in
  the blocks as the input spells them (instead of adding a second block that CCEX did not read), and the input manifest
  follows `Qubit.qubit[].file`, `{s}` names and no longer takes `""` entries as paths.
- **The input file given as `-fFILE`** skipped the `adaptive` block: a bundle input with one ran as a plain bundle run
  (its curve off by up to 2.84 on a six-spin case).
- **The automatic convergence** read its keys with their case, unlike CCEX (`"Method": "gCCE"` ran as cCCE), and
  started CCEX's Python side with the working directory first on `sys.path` (now `PYTHONSAFEPATH=1`, Python 3.11 and
  later).
- **SpinTensorReader**: the QE writer used names it never defined, so its error paths crashed instead of stopping with
  their message; its tag was misspelt `QUANTUMESSPRESSO`; the VASP reader read the defect spin for a Q file alone too
  and stopped on an OUTCAR with magnetization 0 or none.
- **The analyzer**: `-fit_w` with `-ea` raised an IndexError (the average now gets one fit); plots of several curves
  (`-pl_m`) failed with matplotlib 3.11.
- Builds with Intel `icpx`: the adaptive code is compiled without fast floating point, which had folded
  `isnan`/`isinf`.
- `-P` (now `output.cluster_only`) with neither `outfile` nor `-o` stopped with "basic_string: construction from null".
- A refusal of a bad `pulse.bath.sequence` fraction named `sequence[i]`; the refusal of a plan with gCCE, meCCE or
  megCCE and pcce said those plans were supported.

### Refused now (inputs that used to be accepted)

- Every old key, value or form in the table above, with its new name and a pointer to the converter.
- `pulse.name` other than Ramsey, HahnEcho, CPMG, WAHUHA or Manual, or one that does not fit `pulse.count`; WAHUHA with
  `count` ≠ 4; `pulse.sequence` without `pulse.name` Manual; a Manual sequence whose row count is not `count + 1`.
- `pulse.duration` or `detuning_factor` ≠ 1 with cCCE (and `duration` with meCCE, megCCE); `pulse.per_qubit` with an
  `{"index"}` or `{"vector"}` qubit state.
- `hamiltonian.cCCE.knight` or `delta_A` without a zero-field splitting of the qubit; `delta_A` with `hfmedi`; a
  `hamiltonian.cCCE` setting in a gCCE or megCCE run.
- A unit in the wrong case; a Defect value object with both `value` and `values`; `Qubit.qubit[].detuning` `[x]`.
- A seed list of the wrong length; `k_range` without `{k}` in the file names, a file of any k that cannot be read, an
  enabled pCCE cache whose file has no `{k}`, `k_range` with an `adaptive` block; a string where `bath_state_file`,
  `defect_axis_file` or `defect_state_file` takes a list, or a list of another length than `bathfile`.
- `cluster.bundle.connections` `checked` with rdip.
- `gCCE` with `cluster.mode` `dsj`/`itb`/`dsjitb` (not implemented); the `adaptive` block outside `pcce`/`bundle`, the
  `pcce` block outside `pcce`; a qubit state with more than one of `ms`, `index`, `vector`; 0/1 for the switches.
- `alphams` at the top of a qubit file (it had been ignored since v1.1.0).
- Arrays shorter than needed (Defect types shorter than `n_onsite_spins`, `bath_shift` rows with fewer than 3 numbers);
  an axis state beyond the defect's `n_axis_states`.
- An unknown key of the `savemode` block and a value of the wrong type.
- A hyperfine tensor file of format v3 that contradicts itself or the input: a first line that is not the tag of the
  provenance line's source, the provenance line elsewhere than line 2 or twice, an unknown format version or source, an
  S that is missing or not a half-integer, or an S other than the qubit's.
- SpinTensorReader `--version v1` (it wrote the OUTCAR × 2, right only untagged and for S = 1), and QE `--A` without
  `--spin S`.
- `pulse.bath` (DEER) where it cannot act: methods cCCE and meCCE, `cluster.mode` bundle, `backend` `cuda`, `propagator`
  `sparse`, `gpu_options.warp_max_dim` > 0, `pulse.detuning_factor` ≠ 1, `pulse.per_qubit` true; an unknown key of
  `pulse.bath`; a last bath row that is not identity and two pulses of one sequence on one edge; a `pulse.bath.defect` no
  bath spin has, and `levels` that are not two different levels of the driven spins (a pump that reaches no spin of
  that name only warns).
- meCCE / megCCE where they cannot act: a `hamiltonian.meCCE`/`.megCCE` section in a run of another method, a key other
  than `jump_operators` in it, jump operators in both sections, an entry without `bath_name`, `operator` and `rate`, an
  unknown operator, a rate that is negative or not finite, a `bath_name` that is neither the qubit nor a bath spin of the
  run or is both; for meCCE a Manual sequence, a rotation other than an ideal π flip and qubit jump operators; more than
  one qubit for meCCE or with a qubit jump; for both, `pulse.per_qubit` true, `pulse.detuning_factor` ≠ 1, `propagator`
  other than `eigen`, `cluster.mode` other than `standard` and `pcce`, `evolution_update: "recursive"` given, `quantity`
  `dm`, savemode `cluster_density_matrix`, `rdip_ladder`.
- `-e` with `-f`. Other combinations a method does not support: `test/code_verification/check_refusals` lists them.

### Removed

- The command-line options other than `-f`, `-o DIR`, `-v`, `-e`, `-h` (see the table above).
- MKL as a requirement (it is optional now), `hf_cutoff`, `DefectTotSpin`, `CorrTotSpin`, the bundled `zlib/`, the
  top-level `main.cpp` and `config/` (moved into `src/` and the tools).
- `example/tools/Runner/` (job folders per variable value) and `example/analysis/ConvTest/`: `configuration.k_range`
  and the automatic convergence took their place (git history has both). The adaptive-refinement study scripts moved to
  (internal notes).
- The P1 bath of the tutorials keeps 1 and 10 ppm × 5 configurations (99 MB → 0.9 MB); the other configurations are in
  git history (`git show c3c862972:example/tutorial/04_Diamond_NV_P1Bath/bath/...`).

### Compared with v1.2.2

Same inputs and node (one node with 64 cores, two 32-core Xeon Gold 6530; 64 ranks, one per core), median of 3 runs.
v1.2.2 built with icpx 2023.2 without MKL; 2.0.0 (a
development build, commit 7b33ca7d5) built with `do_compile.sh idm4`, MKL on, so the times include the effect of MKL.
Wall time of each run in seconds, start-up included.

| Case | v1.2.2 (s) | 2.0.0 (s) |
|---|---|---|
| P1 production cCCE2, CR60 σ4, 30 states | 442.4 | 22.3 |
| P1 production gCCE2, CR60 σ4, 30 states | 544.0 | 161.7 |
| P1 production gCCE2, CR0, 30 states | 96.0 | 31.4 |
| P1 production cCCE2, CR0, 30 states | 73.5 | 5.3 |
| multi-NV (6 NV) gCCE2, density matrix | 55.2 | 21.9 |
| h-BN gCCE3, density matrix, 1 state | 44.1 | 6.4 |
| h-BN gCCE2, density matrix, 10 states | 27.7 | 5.9 |
| P1 pCCE2 (partition 2), 10 states | 5.9 | 1.3 |
| ¹³C ensemble cCCE2 | 2.2 | 1.6 |
| tutorial-size inputs (13 cases) | 1.0–3.2 | 1.0–1.6 |

Where the coherence converges, it agrees with v1.2.2 within 1e-4; with `deltat` rounded to a float as v1.2.2 keeps it,
within 1e-9 to 1e-10 (P1 cCCE2 and gCCE2, ¹³C). Larger differences come only from fixes in this release (`deltat` read
as a float, the ensemble half interval in double, the complex conjugate for an even pulse count, the static Overhauser
field of gCCE, the density-matrix `nan`, the pCCE partition) or from points where the curve already diverges
(|L| > 1), where v1.2.2 often wrote `nan`. The scripts and tables are under (internal notes).

#### On the GPU node

The released code (job 17715), built as the `gpu` preset builds it after 2.0.0 (CUDA 13.2, Intel MPI 2021.10, icpx
2023.2; see 2.0.1), on the GPU node (one node with 32 cores, one Xeon 6731P, and one
RTX PRO 6000; 32 ranks and 16 for the 6 NV, one per core; CUDA MPS): v1.2.2
(CPU only, built with icpx 2023.2 and the same MPI) against 2.0.0 `main.out` (CPU) and `main_gpu.out` (GPU,
`backend auto`) on 23 of A97's 27 inputs: its four P1 production inputs (CR0 and CR60, cCCE2 and gCCE2, 30 states) were
not at hand. Three runs per build in the order AB BA AB, the median in seconds with start-up; the last column as in the
table after it.

| Case | v1.2.2 (s) | 2.0.0 CPU (s) | 2.0.0 GPU (s) | Largest coherence difference to v1.2.2 |
|---|---|---|---|---|
| h-BN gCCE3, density matrix, 1 state | 79.1 | 10.5 (7.5×) | 10.9 (7.3×) | |
| P1 pCCE2, sK 3, 10 states | 63.6 | 2.4 (26.1×) | 2.5 (25.0×) | |
| multi-NV (6 NV) gCCE2, density matrix | 52.4 | 29.4 (1.8×) | 14.8 (3.5×) | 2e-6 |
| h-BN gCCE2, density matrix, 10 states | 24.8 | 7.2 (3.4×) | 7.7 (3.2×) | |
| P1 pCCE2, sK 2, 10 states | 7.4 | 1.2 (6.1×) | 1.3 (5.6×) | |
| P1 cCCE2, 10 states | 3.1 | 1.0 (3.1×) | 1.1 (2.7×) | 1e-6 |
| natural ¹³C cCCE2 Hahn, 10 states | 2.9 | 1.2 (2.4×) | 1.3 (2.3×) | 8e-5 |
| tutorial-size inputs (12 cases) | 0.8–1.9 | 0.8–1.3 | 1.0–1.5 | |
| 08_GPU cCCE5 (one cluster of dimension 1024) | | 10.4 | 1.9 (5.3× the CPU) | |
| 08_GPU gCCE4 | | 34.6 | 3.1 (11.3× the CPU) | |
| 08_GPU pCCE2 | | 1.0 | 1.2 | |

Geometric mean of the speed-up over the 19 inputs that v1.2.2 runs: 2.0× CPU / 1.8× GPU. The GPU gives the CPU's
coherence to the printed digits but in cCCE5, gCCE4 and h-BN gCCE3, which differ by 1e-10 (relative where the h-BN curve
diverges). Against the 2.0.0 build (CUDA 12.8, Open MPI 4.1.5; A97, the same inputs), the outputs are the same bytes but
for cCCE5 on the GPU (1e-10) and the points where h-BN gCCE2 and gCCE3 diverge (2e-11 relative at most; the 6 NV's
1.14 GB not kept), and the times are within 0.2 s or 5% but for the 6 NV: with CUDA 12.8 its GPU run stopped
(2.0.1, Fixed), and its CPU run took 29.4 s against 22.4 s in all three runs (the GPU run in one of three, 22.1 s
against 14.8 s), the cluster done at 11.1 s in both; the rest is mostly rank 0 writing the 1.14 GB of the `dm` files,
the slower write with Intel MPI of 2.0.1, Builds.

An earlier development build of 2.0.0 was measured
there (Xeon 6731P + one RTX PRO 6000, 32 ranks, one run each): v1.2.2 (CPU only) against 2.0.0 `main.out` (CPU) and
`main_gpu.out` (GPU, `backend auto`), both built with icpx 2023.2 and OpenMPI 4.1.5, on the inputs of the 2.0.0 checks
and the v1.2.2 examples. One run per build, so small differences in time are within the noise. Times in seconds; the
last column is the largest difference of the coherence (`wiDiv`), an empty cell is not compared.

| Case | v1.2.2 (s) | 2.0.0 CPU (s) | 2.0.0 GPU (s) | Largest coherence difference to v1.2.2 |
|---|---|---|---|---|
| P1 production cCCE2, CR0 | 90.4 | 6.4 (14.1×) | 6.4 (14.2×) | 3e-7 |
| P1 production cCCE2, CR60 | 526.4 | 31.3 (16.8×) | 27.4 (19.2×) | 9e-6 |
| P1 production ensemble cCCE2, CR0 | 12.7 | 1.2 (10.7×) | 1.2 (10.4×) | 6e-8 |
| P1 production ensemble cCCE2, CR60 | 80.5 | 3.0 (27.0×) | 3.0 (26.5×) | 2e-8 |
| P1 production gCCE2, CR0 | 117.0 | 40.7 (2.9×) | 22.3 (5.3×) | |
| P1 production gCCE2, CR60 | 681.7 | 219.5 (3.1×) | 108.2 (6.3×) | |
| P1 example cCCE2 e2n, 50 states | 480.8 | 23.5 (20.4×) | 23.7 (20.3×) | |
| P1 example pCCE2, sK 3 | 312.7 | 8.8 (35.7×) | 8.7 (35.8×) | |
| h-BN ensemble cCCE2 | 67.7 | 10.8 (6.3×) | 10.8 (6.3×) | 4e-6 |
| h-BN ensemble Ramsey | 31.9 | 4.4 (7.2×) | 4.7 (6.8×) | |
| h-BN gCCE3 | 78.3 | 10.6 (7.4×) | 11.0 (7.1×) | |
| ¹³C ensemble cCCE2 | 3.9 | 2.0 (1.9×) | 2.0 (1.9×) | 8e-6 |

Geometric mean of the speed-up over the 57 inputs (most of the v1.2.2 examples run in about a second, where start-up
dominates): cCCE (18) 3.2× CPU / 3.1× GPU, gCCE (18) 1.9× / 1.7×, pCCE (3) 14.8× / 14.4×, ensemble (18) 2.8× / 2.5×,
all 2.8× / 2.6×.
