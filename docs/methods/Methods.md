# Theory and implementation

Each page explains what one method or part of CCEX computes and how CCEX does it. Read
[cCCE and gCCE](cce.md) first. It describes the expansion, the two output curves and the Hamiltonian, which the
other pages build on. The input keys are in the [input reference](../user/InputFile.md).

1. [Methods](#meth_methods): what each `method` and `cluster.mode` computes.
2. [Other parts of a run](#meth_parts): bath pulses, the frame rotation, the automatic convergence and the GPU.

<a id="meth_methods"></a>
## 1. Methods

| Page | What it describes |
|---|---|
| [Cluster-correlation expansion: cCCE and gCCE](cce.md) | The expansion over clusters, wiDiv and noDiv, the ensemble and sampled bath states, the Hamiltonian with its units and mean fields, knight, hfmedi and delta_A, the defect bath, and dsj, itb and dsjitb. |
| [Partition CCE (pCCE)](pcce.md) | Partitions of fixed size from a k-means search, the dropped rest spins, connections by centres or spins, and the partition cache. |
| [Bundle clustering and adaptive refinement](bundle-adaptive.md) | Bundles from the coupling strengths, their connections (distance, auto, checked), the limits of the expansion, and the merge and regroup refinement with its outputs and exit statuses. |
| [Master-equation CCE](master-equation.md) | meCCE and megCCE, with Lindblad jump operators on the bath and the qubit. |
| [Semiclassical noise](noise.md) | Method noise, the coherence exp(−χ) from the bath's cluster modes or from a model spectrum, and the correlation function and spectrum. |
| [Pseudospin model of nuclear pairs](pseudospin.md) | Method pseudospin, the cut pair Hamiltonian, the channels grouped or pairwise, its accuracy and the channel table. |

<a id="meth_parts"></a>
## 2. Other parts of a run

| Page | What it describes |
|---|---|
| [DEER bath pulses](deer.md) | Pump pulses on bath spins (`pulse.bath`). |
| [Coordinate frame rotation](coordinate-frame-rotation.md) | A bath file written in another frame (`frame.rotation`). |
| [Automatic convergence](auto-convergence.md) | The runner of the `auto_convergence` block, its steps, when two settings agree, and its results and exit codes. |
| [GPU execution](gpu.md) | The CUDA backend, what runs on the GPU, its options and its memory. |
