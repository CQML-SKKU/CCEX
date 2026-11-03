# Methods

How CCEX computes what it computes, one page per method; the keys themselves are in the
[input reference](../user/InputFile.md).

| Page | What it describes |
|---|---|
| [Cluster-correlation expansion: cCCE and gCCE](cce.md) | the expansion, wiDiv and noDiv, ensemble and sampled states, the two methods, the Hamiltonian and its units, the mean fields, knight / hfmedi / delta_A, the defect bath, dsj / itb / dsjitb |
| [Partition CCE (pCCE)](pcce.md) | partitions of fixed size, the k-means search, the dropped rest spins, connections by centres or spins, the partition cache |
| [Bundle clustering and adaptive refinement](bundle-adaptive.md) | bundles from the coupling strengths, the connections (distance, auto, checked), the expansion limits, the merge and regroup refinement, its outputs and exit statuses |
| [Master-equation CCE](master-equation.md) | meCCE and megCCE: Lindblad jump operators on the bath and the qubit |
| [Gaussian phase noise](noise.md) | method noise: the coherence exp(−χ) from the bath's cluster modes or a model spectrum, the correlation function and spectrum |
| [Pseudospin model of nuclear pairs](pseudospin.md) | method pseudospin: the cut pair Hamiltonian, channels grouped or pairwise, its accuracy, the channel table |
| [DEER bath pulses](deer.md) | pump pulses on bath spins (`pulse.bath`) |
| [Coordinate frame rotation](coordinate-frame-rotation.md) | a bath file in another frame (`frame.rotation`) |
| [GPU execution](gpu.md) | the CUDA backend: what runs on the GPU, its options and the memory |
