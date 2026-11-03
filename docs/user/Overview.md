# Overview

CCEX computes the coherence of spin qubits in a bath of nuclear and electron spins with the cluster-correlation
expansion (CCE). It is a C++ and MPI code written by the CQML group at Sungkyunkwan University:
<https://github.com/CQML-SKKU/CCEX>.

## What it computes

The input names one or more qubits, the bath around them and a pulse sequence. From the positions, gyromagnetic ratios
and interaction tensors, CCEX builds the spin Hamiltonian, divides the bath into clusters and solves each cluster
exactly. The cluster results are then combined into one quantity as a function of time:

- The quantity is the qubit coherence `L(t)` under the pulse sequence, or the qubit density matrix with gCCE
  (`quantity`).
- The method is conventional CCE (cCCE), generalized CCE (gCCE) or one of their master-equation versions, meCCE and
  megCCE (`method`). The clusters can also be partitions (pCCE), bundles, or the disjoint-cluster and inter-bath
  variants dsj, itb and dsjitb (`cluster.mode`).
- Two further methods are models for interpreting a decay. `noise` gives the Gaussian coherence exp(−χ) of the qubit's
  phase noise, computed from the modes of the bath's clusters or from a model spectrum, and can also write the noise
  correlation function and spectrum (`hamiltonian.noise`). `pseudospin` keeps only the bath's pairs, in the pseudospin
  model of nuclear flip-flops, and writes a table of every flip-flop channel (`hamiltonian.pseudospin`).
- The sequence can be a Hahn echo, a Ramsey experiment, CPMG (`pulse.count`) or any sequence written out pulse by
  pulse (`pulse.sequence`). DEER pulses can act on the bath (`pulse.bath`).
- A run covers one bath configuration or several. Electron-spin baths are usually computed with sampled bath states
  (`configuration.s_range`).

A bath can hold nuclear spins with point-dipole or DFT hyperfine and quadrupole tensors, electron-spin defects such as
P1 centres with their own on-site nuclei and axes, or both. Besides the standard clustering by order and distance,
CCEX has bundle clustering (`cluster.mode` bundle, `cluster.bundle.max_dim`) and an adaptive refinement of the
clusters (the `cluster.adaptive` block). The CUDA backend runs the cluster propagation on a GPU
(`execution.backend`). For cCCE it chooses the `gpu_options` an input leaves out at the start of the run, and some of a
node's ranks can use the GPU while the rest stay on the CPU. A coordinate frame rotation (`frame.rotation`) lets the
bath file stay in the crystal frame.

<a id="WhatCCEComputes"></a>
## The cluster-correlation expansion in short

The qubit coherence in the whole bath is written as a product over clusters $C$ of bath spins. The clusters are all
those of up to `order` spins that the clustering keeps, and each contributes its irreducible part $\tilde L_C$:

$$L(t) = \prod_{C} \tilde L_C(t), \qquad \tilde L_C(t) = \frac{L_C(t)}{\prod_{C' \subset C} \tilde L_{C'}(t)},$$

Here $L_C$ is the coherence computed exactly with only the spins of $C$ in the bath, and the product in the denominator
runs over the proper sub-clusters of $C$. Order 1 keeps single spins, order 2 adds pairs, and so on. A result can be
trusted once it stops changing when the order, `rbath` or `rdip` is raised (the two radii decide which spins and which
clusters are included).

The `_wiDiv` output is this product, the CCE coherence. Written out, it is a product of the cluster coherences $L_C$
raised to positive and negative integer powers. The `_noDiv` output keeps only the factors with a positive power, so
it has no division. It never diverges, but it is not a coherence; the analyzer uses it to repair points where a division
blew up ([Analysis](Analysis.md), [Output files](OutputFiles.md)).

Without `configuration.s_range`, the bath starts in its unpolarised (ensemble) state. With `s_range` [first, last],
each state s is a product state of the bath spins, drawn from `configuration.seed` or read from a file, and gets its own
output. The mean over the states approaches the ensemble result.

A physical coherence has $|L| \le 1$. If a CCE product gives $|L| > 1$ or nan, the expansion did not converge for that
run: a sub-cluster term in a denominator came close to zero, or a strongly coupled group of spins was split between
clusters. The usual fixes are a higher order, gCCE, pCCE or bundles, or the analyzer's error correction.

Details are on the methods pages: [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundles and the adaptive refinement](../methods/bundle-adaptive.md).

<a id="WhatRunsWhere"></a>
## What runs where

The table crosses `method` with `cluster.mode`; each entry says which `execution.backend` works. "CPU, GPU" means
`cuda` and `auto` may put clusters on the GPU. "CPU" means `cuda` is refused and `auto` runs on the CPU. "refused"
means the run stops with a message.

| method | `standard` | `pcce` | `bundle` | `dsj`, `itb`, `dsjitb` |
|---|---|---|---|---|
| `cCCE` (default) | CPU, GPU | CPU, GPU | CPU, GPU | CPU |
| `gCCE` | CPU, GPU¹ | CPU | CPU, GPU¹ | refused |
| `meCCE` (cCCE + Lindblad jumps on bath spins) | CPU | CPU | refused | refused |
| `megCCE` (gCCE + jumps on bath spins and the qubit) | CPU | CPU | refused | refused |
| `noise` (Gaussian phase noise) | CPU, GPU² | refused | refused | refused |
| `pseudospin` (the pairs in the pseudospin model) | CPU | refused | refused | refused |

¹ With `execution.propagator` `eigen` (the default); `expm` runs on the CPU.
² With `backend` `cuda`, the sums over the bath's modes run on the GPU and the eigensystems on the CPU; `auto` runs on
the CPU.

DEER bath pulses (`pulse.bath`, gCCE and megCCE with `standard` or `pcce`) and a `Manual` sequence on the ensemble (no
`s_range`) also run only on the CPU. `quantity` `dm` is for gCCE; meCCE, megCCE, noise and pseudospin compute
`coherence` only. GPU runs need the CUDA binary `bin/main_gpu.out`. Details:
[method](InputFile_General.md#tag_method), [cluster.mode](InputFile_Cluster.md#tag_cluster_mode),
[backend](InputFile_Backend.md#tag_backend), [GPU execution](../methods/gpu.md),
[master-equation CCE](../methods/master-equation.md).

To find the method and example that fit your system, see "What should I use for my case?" on the
[Examples](Examples.md#WhatToUse) page.

## Building

Build with `bash do_compile.sh <preset>` or `cmake -B build && cmake --build build -j 8`. The presets, the
requirements and the options are in [Getting started](GettingStarted.md#Building).

## Running

```
mpirun -n 4 bin/main.out -f ccein.jsonc
```

The input file is JSON (`.json` or `.jsonc`) and may contain `//`, `#` and `!` comments. Its keys are described in
[Input file](InputFile.md). Every setting of a run is an input key. The command line only names the input (`-f`), the
output folder (`-o`, which replaces `output.path`), a longer log (`-v`), the input template (`-e`) and the usage
(`-h`); see [Getting started](GettingStarted.md#HowToRun). A run writes the coherence or density matrix of each bath
configuration to text files ([Output files](OutputFiles.md)). The analyzer in `example/analysis/` corrects errors in
them, averages over the ensemble and plots the result.

## Where to go next

- [Getting started](GettingStarted.md): building, the first calculation, the command-line options, the run log, the
  environment variables and the order in which a run reads its input.
- [Input file](InputFile.md): every key by group, and the unit of every number.
- [Output files](OutputFiles.md): every file a run writes, with its columns.
- [Examples](Examples.md): thirteen tutorials in learning order, each with its inputs and expected figures. They start
  with a quickstart that runs in seconds and go on to the GPU backend, the automatic convergence, the frame rotation
  and the noise and pseudospin models.
- [Analysis](Analysis.md) and [Tools](Tools.md): the programs that process the output (`example/analysis/`) and that
  prepare the inputs (`example/tools/`).
- Methods ([overview](../methods/Methods.md)): [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
  [bundles and the adaptive refinement](../methods/bundle-adaptive.md),
  [the coordinate frame rotation](../methods/coordinate-frame-rotation.md),
  [DEER bath pulses](../methods/deer.md), [master-equation CCE](../methods/master-equation.md),
  [Gaussian phase noise](../methods/noise.md), [the pseudospin model](../methods/pseudospin.md) and
  [GPU execution](../methods/gpu.md).
- [Changes](Changes.md): what changed in 2.1.0 and in 2.0.0.
- Data structures and Source files: the code documentation, for developers.

<a id="Citing"></a>
## How to cite / references

`CITATION.cff` at the repository root lists the authors and the version of the code.

The methods CCEX implements are described in:

- W. Yang and R.-B. Liu, Phys. Rev. B 78, 085315 (2008) and Phys. Rev. B 79, 115320 (2009): the CCE and its ensemble
  form.
- M. Onizhuk et al., PRX Quantum 2, 010311 (2021), and M. Onizhuk and G. Galli, Adv. Theory Simul. 4, 2100254 (2021):
  the generalized CCE.
- M. Onizhuk et al., Phys. Rev. Lett. 132, 250401 (2024): the master-equation CCE.
- H. Seo et al., Nat. Commun. 7, 12935 (2016), and J. Lee, H. Park, H. Seo, npj 2D Mater. Appl. 6, 60 (2022): the
  pseudospin model of nuclear pairs.

[cCCE and gCCE](../methods/cce.md) lists them with their DOIs. A paper describing CCEX itself is to be added by the
CCEX authors.
