# Output tags
The `output` block. Its keys at the top level are refused with their place in the block. What each file holds, its
name and its columns are in [Output files](OutputFiles.md).

```
"output": {"path": "results", "file": "ccein", "savemode": {"mode": "coherence", "nodiv": true},
           "cluster_contributions": {"enabled": true, "max_mb": 256}, "cluster_only": false}
```

- [path](#tag_output_path)
- [file](#tag_outfile)
- [savemode](#tag_savemode) (`mode`, `nodiv`)
- [cluster_contributions](#tag_cluster_contributions) (`enabled`, [max_mb](#tag_cluster_contribution_max_mb))
- [cluster_only](#tag_cluster_only)

<a id="tag_output_path"></a>
## output.path
<b>Description:</b><br>
The folder of the outputs. A relative folder is read from the working directory. A missing folder is made. <br>
The folder is `-o DIR` from the command line, else `output.path`, else the working directory. <br>
Refused:

- an absolute [output.file](#tag_outfile) together with `output.path` (or `-o`), since the two would name two places;
- the adaptive refinement, whose steps go to `cluster.adaptive.workdir`;
- the automatic convergence, whose runs go to `auto_convergence.workdir`. The runner refuses `output.path` and
  `output.file`, as it refuses `-o`.

At the start the run writes and removes a file in the output folder (or in the folder given in `output.file`). A folder
that does not exist or cannot be written is refused before anything is calculated. <br>
<b>Data Type:</b> String (a folder)

<b>Default:</b> the working directory

<b>Example:</b> ``` "output" : {"path": "results"} ```

<a id="tag_savemode"></a>
## output.savemode
<b>Description:</b><br>
What is written, and whether `_noDiv` is written beside `_wiDiv`. A string gives the mode, with `nodiv` true. The block
`{"mode": ..., "nodiv": true|false}` gives both; a missing `mode` is `"coherence"` and a missing `nodiv` is `true`. <br>
<b>mode:</b>

- `"coherence"` : The coherence, "file_wiDiv/noDiv". With sampled bath states (`configuration.s_range`) also one file
                 pair per state, "file_state{state}_wiDiv/noDiv". <br>
- `"coherence_tagged"` : Like `"coherence"`, with the same file names. Every time row also carries the zx and zz
                 components of the first bath spin's hyperfine tensor and the first bath file's name, so that results of
                 many bath files can be told apart after collecting them. <br>
- `"cluster_coherence"` : The result and, for each cluster and each state, its coherence in
                 "file_{cluster}_state{state}_wiDiv/noDiv". `{cluster}` is the bath indices of the cluster joined by `_`;
                 `_state{state}` appears only with `configuration.s_range`. <br>
- `"cluster_density_matrix"` : Everything `"cluster_coherence"` writes. With gCCE also, for every cluster, the density
                 matrix of the whole qubit + cluster space at every time step, as `n` × `n` blocks (`n` the full
                 dimension), in "file_{cluster}_full_state{state}_wiDiv/noDiv". The cluster's bath indices are as in
                 `"cluster_coherence"`; the zeroth cluster (the qubits alone) is "file_full_state{state}_...". The
                 `_wiDiv` and `_noDiv` files hold the same matrix. It costs one full-dimension matrix per time step and
                 cluster, in memory and on disk, so it is for small systems
                 (`example/tutorial/06_cCCE_vs_gCCE/gCCE1_Hahn_handwritten/`). cCCE has no qubit + bath density matrix
                 and writes what `"cluster_coherence"` writes. <br>

<b>nodiv:</b> Write `_noDiv` beside `_wiDiv`. `_wiDiv` is the CCE coherence. `_noDiv` is the product of the
positive-power cluster terms only (no division); it never diverges, and it is **not** a coherence. On by default because
these read it:

- the error correction of the group analyzer (`example/analysis/Analyzer`; \|wiDiv − noDiv\| ≥ 2 marks a divergent point);
- the standing checks;
- the paper scripts.

`false` writes `_wiDiv` alone (half the files), for runs whose analysis does not need `_noDiv`. Must be a JSON boolean. <br>

Old values and the old `nodiv` key are refused with a message naming the current form ([old names](InputFile.md#OldNames)).

<b>Data Type:</b> String or object (`mode`: string, `nodiv`: boolean)

<b>Default:</b> "coherence" (`nodiv` true)

<b>Example:</b> ``` "output" : {"savemode": "cluster_coherence"} ```, ``` "output" : {"savemode": {"mode": "coherence", "nodiv": false}} ``` (both in `example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_clusterReport/`)

<a id="tag_outfile"></a>
## output.file
<b>Description:</b><br>
Output file name prefix. <br>
"file_..._wiDiv" is the final result of the coherence function. <br>
"file_..._noDiv" is the result without dividing the final coherence function by those of its sub-clusters. <br>
With [configuration.k_range](InputFile_General.md#tag_configuration_k_range) each configuration k writes
"file_k&lt;k&gt;_...". <br>
A relative name is read from the output folder: `-o DIR`, else [output.path](#tag_output_path), else the working
directory. `-o` names the folder, not this name. <br>
Refused: an absolute name together with `-o` or `output.path`, and any value with the automatic convergence (as
`output.path`). <br>
<b>Data Type:</b> String

<b>Default:</b> the input file's name without `.json` or `.jsonc` (`-f run/ccein.json` writes `ccein_noDiv`, `ccein_wiDiv`, ... in the
output folder).

<b>Example:</b> ``` "output" : {"file": "filepath"} ```

<a id="tag_cluster_only"></a>
## output.cluster_only
<b>Description:</b><br>
`true` writes the cluster plan to `<file>.clusters.json` and stops before the dynamics; nothing else is computed or
written. Use it to see how many clusters an `rbath` and `rdip` give before a long run
(`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_clusterReport/`). <br>
Works with method cCCE or gCCE, and with cCCE `cluster.mode` pcce with explicit connections. Refused:

- a `bundle` input whose connections the run checks (no `rdip`);
- the adaptive probe, factors and except work;
- a `cluster.adaptive` block;
- an `auto_convergence` block.

<b>Data Type:</b> Boolean (a JSON `true` or `false`; anything else is refused)

<b>Default:</b> `false`

<b>Example:</b> ``` "output" : {"cluster_only": true} ```

<a id="tag_cluster_contributions"></a>
## output.cluster_contributions
<b>Description:</b><br>
Record the contribution of every cluster to the coherence over all time points (cCCE and gCCE, `quantity` =
`coherence`). A diagnostic, to find the clusters that dominate the decay or diverge
(`example/tutorial/07_pCCE_and_bundle/cCCE2_Hahn_clusterReport/`). <br>
Files, with columns in [Output files](OutputFiles.md#out_contributions):

- `<output.file>.contributions.state-<s>.tsv` for each bath state: one row per cluster with sites, bundles, production
  power, dimension, the largest |log L̃|, undefined points and worst step;
- `<output.file>.contributions.sites.tsv` once: site, name, spin, gyromagnetic ratio, position.

The raw table per rank is bounded by [max_mb](#tag_cluster_contribution_max_mb). The value is a block
`{"enabled": true|false, "max_mb": n}`; a block without `enabled` is on, and a plain boolean is refused. <br>

<b>Data Type:</b> Object `{"enabled": Boolean, "max_mb": Integer}`

<b>Default:</b> off (absent)

<b>Example:</b> ``` "output" : {"cluster_contributions": {"enabled": true}} ```

<a id="tag_cluster_contribution_max_mb"></a>
## output.cluster_contributions.max_mb
<b>Description:</b><br>
The memory per rank the contribution table of `output.cluster_contributions` may take. The run
stops before the propagation if the table would exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Unit:</b> MiB per rank

<b>Default:</b> `256`

<b>Example:</b> ``` "output" : {"cluster_contributions": {"enabled": true, "max_mb": 512}} ```
