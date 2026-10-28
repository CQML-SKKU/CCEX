# Overview

CCEX computes the coherence of spin qubits in a bath of nuclear and electron
spins by the cluster-correlation expansion (CCE). It is a C++ and MPI code of the
CQML group at Sungkyunkwan University: <https://github.com/CQML-SKKU/CCEX>.

## What it computes

The input names a qubit (or several), the bath around it and a pulse sequence.
From the positions, gyromagnetic ratios and interaction tensors CCEX builds the
spin Hamiltonian, divides the bath into clusters, solves the dynamics of every
cluster exactly and combines them into the quantity asked for as a function of
time:

- the coherence function `L(t)` of the qubit under the pulse sequence, or with
  gCCE its density matrix (`quantity`);
- by conventional CCE (cCCE) or generalized CCE (gCCE), or their master-equation forms meCCE and megCCE
  (`method`), with partition CCE (pCCE), bundle
  clustering and the disjoint-cluster and inter-bath variants dsj, itb and dsjitb (`cluster.mode`);
- for a Hahn echo, a Ramsey experiment, CPMG (`pulse.count`) or any sequence written
  out pulse by pulse (`pulse.sequence`), with DEER pulses on the bath (`pulse.bath`);
- for one bath configuration or an ensemble of them, and for electron-spin baths
  with sampled bath states (`configuration.s_range`).

Baths can be nuclear spins with point-dipole or DFT tensors (hyperfine and
quadrupole), electron-spin defects such as P1 centres with their own on-site
spins and axes, or both. Beyond the standard clustering by order and distance
CCEX has bundle clustering (`cluster.mode` bundle, `cluster.bundle.max_dim`) and an adaptive
refinement of the clusters (the `cluster.adaptive` block), a CUDA backend for the
cluster propagation (`execution.backend`), and a coordinate frame rotation that lets the
bath file stay in the crystal frame (`frame.rotation`).

<a id="WhatCCEComputes"></a>
## The cluster-correlation expansion in short

The coherence of the qubit in the whole bath is written as a product over clusters $C$ of bath spins (all clusters of
up to `order` spins that the clustering keeps), each with its irreducible contribution $\tilde L_C$:

$$L(t) = \prod_{C} \tilde L_C(t), \qquad \tilde L_C(t) = \frac{L_C(t)}{\prod_{C' \subset C} \tilde L_{C'}(t)},$$

where $L_C$ is the coherence computed exactly with the bath reduced to the spins of $C$ and the product in the
denominator runs over the proper sub-clusters of $C$. Order 1 keeps single spins, order 2 adds pairs, and so on; a
result is trusted when it no longer changes with the order (and with `rbath` and `rdip`, the radii that decide which
spins and which clusters enter).

- **wiDiv and noDiv.** The `_wiDiv` output is this product, the CCE coherence. Written out, the product is a product of the cluster
  coherences $L_C$ raised to positive and negative integer powers; the `_noDiv` output keeps only the factors with a
  positive power (no division). It never diverges but is not a coherence; the analyzer uses it to repair points where
  a division blew up ([Analysis](Analysis.md), [Output files](OutputFiles.md)).
- **Ensemble and sampled states.** Without `configuration.s_range` the bath starts in its unpolarised (ensemble) state.
  With `s_range` [first, last] each state s instead takes a product state of the bath spins drawn from `configuration.seed` (or read
  from a file), one output per state; their mean approaches the ensemble result, and electron-spin baths are usually
  computed this way.
- **Divergence.** A physical coherence has $|L| \le 1$. A CCE product with $|L| > 1$ (or nan) means the expansion did
  not converge for that run: a sub-cluster contribution in a denominator came close to zero, or a strongly coupled
  group of spins was cut apart. Higher order, gCCE, pCCE or bundles, or the analyzer's error correction are the usual
  remedies.

The methods pages give the details: [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
[bundles and the adaptive refinement](../methods/bundle-adaptive.md).

<a id="WhatRunsWhere"></a>
## What runs where

`method` × `cluster.mode` × `execution.backend`. "CPU, GPU" means `cuda` and `auto` may put clusters on the card;
"CPU" means `cuda` is refused and `auto` runs on the CPU; "refused" means the run stops with a message.

| method | `standard` | `pcce` | `bundle` | `dsj`, `itb`, `dsjitb` |
|---|---|---|---|---|
| `cCCE` (default) | CPU, GPU | CPU, GPU | CPU, GPU | CPU |
| `gCCE` | CPU, GPU¹ | CPU | CPU, GPU¹ | refused |
| `meCCE` (cCCE + Lindblad jumps on bath spins) | CPU | CPU | refused | refused |
| `megCCE` (gCCE + jumps on bath spins and the qubit) | CPU | CPU | refused | refused |

¹ with `execution.propagator` `eigen` (the default); `expm` runs on the CPU. Also CPU only: DEER bath pulses
(`pulse.bath`, gCCE and megCCE with `standard` or `pcce`) and an ensemble (no `s_range`) `Manual` sequence. `quantity`
`dm` is for gCCE; meCCE and megCCE compute `coherence` only. The GPU needs the CUDA binary `bin/main_gpu.out`. Details:
[method](InputFile_General.md#tag_method), [cluster.mode](InputFile_Cluster.md#tag_cluster_mode),
[backend](InputFile_Backend.md#tag_backend), [GPU execution](../methods/gpu.md),
[master-equation CCE](../methods/master-equation.md).

Which method and which example fit a given system is the table "What should I use for my case?" on the
[Examples](Examples.md#WhatToUse) page.

## Building

`bash do_compile.sh <preset>` or `cmake -B build && cmake --build build -j 8`; the presets, what is needed and the
options are in [Getting started](GettingStarted.md#Building).

## Running

```
mpirun -n 4 bin/main.out -f ccein.jsonc
```

The input file is JSON (`.json` or `.jsonc`) with `//`, `#` and `!` comments; its keys are
the reference in [Input file](InputFile.md). Every setting of a run is a key of the input
file: the command-line options only name the input (`-f`), the output folder (`-o`, in
place of `output.path`), a longer log (`-v`), the input template (`-e`) and the usage (`-h`);
they are in [Getting started](GettingStarted.md#HowToRun). A run writes the coherence (or density
matrix) per bath configuration into text files ([Output files](OutputFiles.md)), which the analyzer in
`example/analysis/` error-corrects, averages over the ensemble and plots.

## Where to go next

- [Getting started](GettingStarted.md): building, the first calculation, the command-line options, the
  run log, the environment variables and the order in which a run reads its input.
- [Input file](InputFile.md): every key, by group, and the units of every number.
- [Output files](OutputFiles.md): every file a run writes, and its columns.
- [Examples](Examples.md): twelve worked tutorials, each a page of this site with its inputs and expected
  figures, in learning order from a quickstart of seconds to the GPU backend, the automatic convergence and the
  frame rotation.
- [Analysis](Analysis.md) and [Tools](Tools.md): what runs on the output, and what
  makes the inputs (`example/analysis/`, `example/tools/`).
- Methods ([overview](../methods/Methods.md)): [cCCE and gCCE](../methods/cce.md), [pCCE](../methods/pcce.md),
  [bundles and the adaptive refinement](../methods/bundle-adaptive.md),
  [the coordinate frame rotation](../methods/coordinate-frame-rotation.md),
  [DEER bath pulses](../methods/deer.md), [master-equation CCE](../methods/master-equation.md) and
  [GPU execution](../methods/gpu.md).
- [Changes](Changes.md): what changed in 2.0.0.
- Data structures and Source files: the code itself, for developers.

<a id="Citing"></a>
## How to cite / references

The authors and the version of the code are in `CITATION.cff` at the repository root.

References to be added by the CCEX authors.
