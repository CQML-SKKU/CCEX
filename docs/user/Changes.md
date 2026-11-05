# Changes

This page summarizes what changed in each version, in the current key names. The file `CHANGELOG.md` at the root of
the CCEX repository lists every change, with its details and measurements.

## 2.1.0

This section lists what changed in 2.1.0 (released 2026-11-04) since 2.0.1. Version 2.0.1 (2026-10-22) fixed GPU gCCE
and the builds of 2.0.0.

- **Two new methods.**
  - `noise` gives the Gaussian coherence exp(−χ) of the qubit's phase noise. It works from the modes of the bath's cCCE
    clusters or from a model spectrum (Ornstein–Uhlenbeck or a power law). With `hamiltonian.noise.compute` it also
    gives the correlation function C(t) and the spectrum S(f). Its mode sums also run on the GPU
    ([method](InputFile_General.md#tag_method), tutorial
    [Advanced 13: Semiclassical noise (method noise)](../../example/tutorial/advanced/13_noise/README.md)).
  - `pseudospin` gives the bath's pairs alone, in the pseudospin model of nuclear flip-flops, for any nuclear spin. It
    writes a table of every channel's Ω and Δ ([channels](InputFile_General.md#tag_pseudospin_channels), tutorial
    [Advanced 14: Pseudospin pairs](../../example/tutorial/advanced/14_pseudospin/README.md)).
- **GPU runs on a full node** ([Execution tags](InputFile_Backend.md#tag_backend)).
  - `execution.gpu_options.ranks_per_gpu` lets only some ranks of a node use the GPU while the others stay on the CPU.
    The clusters are split between them by measured cost.
  - A GPU batch of cCCE clusters is diagonalized in one call.
  - For cCCE with `backend` `auto` or `cuda`, the `gpu_options` an input leaves out are chosen at the start, with a
    `[gpu-auto]` line in the log.
- **Builds.** An nvcc older than 13.2 is accepted for cards below compute capability 12.0.
- **Fixed.**
  - The MKL builds write the same last digits from run to run.
  - The quadrupole structure `"strained"` gives a far spin the tensor of `file_without_defect`, as `"bubble"` does. It
    took the defect file's tensor again.
  - `"strained"` and `"bubble"` with a vertex-format `file_without_defect` no longer hang.
  - cCCE with `quantity` dm is refused at input. It stopped in the middle of the run.
- **Tutorials and documentation.** The tutorials are split into `basic/` (the input itself, six topics to read in
  order) and `advanced/` (fourteen topics); [Examples](Examples.md) lists both. Every tutorial input explains its keys
  in comments. New: the pseudospin tutorial, and calculations for the noise spectrum, the automatic GPU options, the
  tensor files from a VASP run, the analyzer and a run on a Slurm cluster. The site is reorganized: one page per input
  block with a small table per key, a page for the old names ([Upgrading old inputs](InputFile_OldNames.md)), and the
  method pages under [Theory and implementation](../methods/Methods.md), new for noise, pseudospin and the automatic
  convergence.
- **Inputs.** No key of 2.0.1 was renamed or removed, so a 2.0.1 input runs as it is. In `bin/main_gpu.out` with
  `backend` `auto`, the automatic `gpu_options` change which clusters go to the card. The outputs therefore move
  within rounding.

## 2.0.0

This section lists what changed in 2.0.0 (released 2026-10-17) since v1.2.2, in the current key names.

### Summary

- **A new input format.**
  - The keys are grouped in blocks: `configuration`, `Qubit`, `Defect`, `hamiltonian`, `pulse`, `cluster`,
    `execution`, `output`, `frame`, `auto_convergence`.
  - Inputs may be `.jsonc` files with `//` comments.
  - Many numbers take a unit (`{"value": 5, "unit": "us"}`).
  - `main.out -e` writes a commented template with every key.
  - The command line is only `-f`, `-o DIR`, `-v`, `-e`, `-h`.
  - Old inputs are refused with the new name, and a converter rewrites them (see "Updating old inputs").
- **Faster.** The test inputs run up to about 35 times faster than v1.2.2, mostly with unchanged results. MKL is used for large clusters when the build finds it. A CUDA backend runs
  cCCE, gCCE and pCCE on a GPU.
- **New methods and options.**
  - ME-CCE and ME-gCCE (Lindblad dissipation).
  - DEER pulses on bath spins.
  - pCCE with gCCE and the master-equation methods.
  - Bundle clustering with an adaptive refinement.
  - The automatic convergence of rbath, rdip, the time window and the number of configurations.
  - Several bath configurations from one input.
- **Physics fixes that change results.**
  - `hfmedi`: the complete second-order term.
  - `knight`: the nuclear Zeeman vector, with the $|m_s|$ of each branch.
  - The static Overhauser field of gCCE with sampled states.
  - The h-BN tutorial hyperfine tensors, which were read twice too large.
  - A memory corruption of pCCE with state files.
- **Safer inputs.** A run stops before the calculation, with one message on one rank, for a key that looks like a
  misspelling of one the run reads (two letters or fewer apart), a file that does not open, a combination a method
  cannot do, or a run that does not fit in memory. Any other key the run does not read gets a warning, and a key that
  is read but not used is reported.
- **Tutorials and documentation.** Eleven tutorials in learning order (`example/tutorial/00_quickstart` …
  `10_coordinate_rotation`), each with `run.sh` and `expected.png`. A rewritten manual (`docs/user/`), method pages
  (`docs/methods/`) and a doxygen site.

### Updating old inputs

`python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites inputs in place. It keeps their comments and
lists what it changed. `--check` lists only, and `--stdout` prints the new text. A rewritten input gives the same
results as before.

[Upgrading old inputs](InputFile_OldNames.md) lists every renamed or removed key, by block, with what it is now.
