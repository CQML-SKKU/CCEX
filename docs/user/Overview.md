# Overview

CCEX computes the coherence of spin qubits in a bath of nuclear and electron spins with the cluster-correlation
expansion (CCE). It is a C++ and MPI code written by the CQML group at Sungkyunkwan University:
<https://github.com/CQML-SKKU/CCEX>.

1. [What it computes](#WhatItComputes): the input, the methods, the pulse sequences and the baths.
2. [The cluster-correlation expansion in short](#WhatCCEComputes): the formula, the two output files, and when a
   result can be trusted.
3. [What runs where](#WhatRunsWhere): which method runs with which clustering, on the CPU or the GPU.
4. [Where to go next](#WhereNext): the other pages.
5. [How to cite](#Citing): the papers behind the methods.

<a id="WhatItComputes"></a>
## 1. What it computes

The input names one or more qubits, the bath around them and a pulse sequence. From the spin positions, gyromagnetic
ratios and interaction tensors, CCEX builds the spin Hamiltonian, divides the bath into clusters, solves each cluster
exactly and combines the results into the qubit coherence `L(t)`, or with gCCE the qubit density matrix
(`quantity`).

| Choice | Options | Key |
|---|---|---|
| Method | conventional CCE (cCCE), generalized CCE (gCCE), and their master-equation versions meCCE and megCCE | `method` |
| Clusters | by order and distance; partitions (pCCE); bundles of strongly coupled spins; the variants dsj, itb and dsjitb; an adaptive refinement | `cluster.mode`, `cluster.adaptive` |
| Sequence | Hahn echo, Ramsey, CPMG, or any sequence written pulse by pulse; DEER pulses on the bath | `pulse` |
| Bath | nuclear spins with point-dipole or DFT tensors, electron-spin defects such as P1 centres with their own nuclei and axes, or both | `configuration`, `Defect` |
| Bath state | the unpolarised ensemble, or sampled product states (usual for electron-spin baths) | `configuration.s_range` |
| Hardware | CPU with MPI; the cluster propagation also on a GPU | `execution.backend` |

Two more methods interpret a decay rather than compute it cluster by cluster. `noise` gives the Gaussian coherence
exp(−χ) of the qubit's phase noise, and can write the noise spectrum. `pseudospin` keeps only the pairs of bath
spins, in the pseudospin model of nuclear flip-flops, and lists every flip-flop channel.

A coordinate frame rotation (`frame.rotation`) lets the bath file stay in the crystal frame.

A run reads the input and the files it names, finds the clusters and spreads them over the MPI ranks. It then
computes and saves each bath state in turn:

![The steps of a CCEX run](../images/CCEX_workflow.png)

<a id="WhatCCEComputes"></a>
## 2. The cluster-correlation expansion in short

CCE writes the qubit coherence in the whole bath as a product over clusters $C$ of bath spins. The clusters are those
of up to `order` spins that the clustering keeps, and each contributes its irreducible part $\tilde L_C$:

$$L(t) = \prod_{C} \tilde L_C(t), \qquad \tilde L_C(t) = \frac{L_C(t)}{\prod_{C' \subset C} \tilde L_{C'}(t)},$$

Here $L_C$ is the coherence computed exactly with only the spins of $C$ in the bath, and the product in the denominator
runs over the proper sub-clusters of $C$. Order 1 keeps single spins, order 2 adds pairs, and so on. A result can be
trusted once it stops changing when the order, `rbath` or `rdip` is raised; the two radii decide which spins and which
pairs are included ([Basic 03](../../example/tutorial/basic/03_convergence/README.md)).

A run writes two files. `_wiDiv` is this product, the CCE coherence. Written out, it is a product of the cluster
coherences $L_C$ raised to positive and negative integer powers. `_noDiv` keeps only the factors with a positive power,
so it has no division and never diverges, but it is not a coherence. The analyzer uses it to repair points where a
division blew up ([Analysis](Analysis.md), [Output files](OutputFiles.md)).

Without `configuration.s_range`, the bath starts in its unpolarised (ensemble) state. With `s_range` [first, last],
each state s is a product state of the bath spins, drawn from `configuration.seed` or read from a file, with its own
output. The mean over many states approaches the ensemble result.

A physical coherence has $|L| \le 1$. A product with $|L| > 1$ or nan means the expansion did not converge in that
run: a sub-cluster term in a denominator came close to zero, or a strongly coupled group of spins was split between
clusters. A higher order, gCCE, pCCE or bundles usually help, and the analyzer's error correction repairs single
points. The methods pages have the details: [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundles and the adaptive refinement](../methods/bundle-adaptive.md).

<a id="WhatRunsWhere"></a>
## 3. What runs where

Each entry says where a `method` runs with a `cluster.mode`:

- **CPU, GPU**: `execution.backend` `cuda` and `auto` may put clusters on the GPU.
- **CPU**: `cuda` stops the run with an error, and `auto` runs on the CPU.
- **refused**: the run stops with a message.

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

Further limits:

- DEER bath pulses (`pulse.bath`, gCCE and megCCE with `standard` or `pcce`) run only on the CPU, and so does a
  `Manual` sequence on the ensemble (no `s_range`) with cCCE or meCCE.
- `quantity` `dm` is for gCCE with `cluster.mode` `standard` or `pcce`. meCCE, megCCE, noise and pseudospin compute
  the coherence only.
- GPU runs need the CUDA binary `bin/main_gpu.out` ([GPU execution](../methods/gpu.md)).

The keys are described under [method](InputFile_General.md#tag_method),
[cluster.mode](InputFile_Cluster.md#tag_cluster_mode) and [backend](InputFile_Backend.md#tag_backend). To choose a
method and an example for your system, see [What should I use for my case?](Examples.md#WhatToUse).

<a id="WhereNext"></a>
## 4. Where to go next

- [Getting started](GettingStarted.md): what to install, how to build, and a first calculation. Start here.
- [Building and running](BuildingAndRunning.md): the build options, the command-line options, the exit codes, the run
  log, the environment variables and troubleshooting.
- [Theory and implementation](../methods/Methods.md): what each method computes and how CCEX does it.
- [Input file](InputFile.md): every key by group, and the unit of every number.
- [Output files](OutputFiles.md): every file a run writes, with its columns.
- [Examples](Examples.md): the tutorials, with their inputs and expected figures. The six basic ones run in seconds;
  the thirteen advanced ones cover other baths, DFT tensors, the analysis, the GPU, the automatic convergence, the
  noise and pseudospin models and runs on a cluster.
- [Analysis](Analysis.md): the scripts that average, fit and plot the output.
- [Tools](Tools.md): the programs that prepare the inputs.
- [Changes](Changes.md): what changed in 2.1.0 and in 2.0.0.
- Data structures and Source files: the code documentation, for developers.

<a id="Citing"></a>
## 5. How to cite

`CITATION.cff` at the repository root lists the authors and the version of the code.

The methods CCEX implements are described in:

- W. Yang and R.-B. Liu, Phys. Rev. B 78, 085315 (2008) and Phys. Rev. B 79, 115320 (2009): the CCE and its ensemble
  form.
- M. Onizhuk et al., PRX Quantum 2, 010311 (2021), and M. Onizhuk and G. Galli, Adv. Theory Simul. 4, 2100254 (2021):
  the generalized CCE.
- M. Onizhuk et al., Phys. Rev. Lett. 132, 250401 (2024): the master-equation CCE.
- H. Seo et al., Nat. Commun. 7, 12935 (2016), and J. Lee, H. Park, H. Seo, npj 2D Mater. Appl. 6, 60 (2022): the
  pseudospin model of nuclear pairs.

[cCCE and gCCE](../methods/cce.md) lists them with their DOIs. A paper on CCEX itself will be added here.
