# Output tags

The `output` block says where the outputs go and what is written, and has two keys that inspect the clusters. What
each file holds is in [Output files](OutputFiles.md).

1. [Writing the results](#otag_where): `path`, `file` and `savemode`.
2. [Inspecting the clusters](#otag_diag): `cluster_only` and `cluster_contributions`.

```json
"output": {"path": "results", "file": "ccein", "savemode": {"mode": "coherence", "nodiv": true}}
```

<a id="otag_where"></a>
## 1. Writing the results

<a id="tag_output_path"></a>
### output.path

| Type | Default |
|---|---|
| string (a folder) | the working directory |

The folder of the outputs, made when it is missing. A relative folder is read from the working directory, and `-o DIR`
on the command line takes its place. It is refused with a `cluster.adaptive` or `auto_convergence` block, whose
`workdir` holds the outputs instead.

```json
"output": {"path": "results"}
```

<a id="tag_outfile"></a>
### output.file

| Type | Default |
|---|---|
| string | the input file's name without `.json` or `.jsonc` |

The start of every output name, such as `<file>_wiDiv` and `<file>_noDiv` ([Names](OutputFiles.md#out_names)). An
absolute name cannot go with `-o` or `output.path`. It is refused with an `auto_convergence` block.

```json
"output": {"file": "run1"}
```

<a id="tag_savemode"></a>
### output.savemode

| Type | Default |
|---|---|
| string, or `{"mode": string, "nodiv": boolean}` | `"coherence"`, with `nodiv` true |

What is written. A string gives the `mode`. With `nodiv` false the `_noDiv` file is left out. That file is the
undivided product, which the error correction of the analyzer reads ([Analysis](Analysis.md)).

| mode | Writes |
|---|---|
| `"coherence"` | `<file>_wiDiv` and `_noDiv`, one pair per sampled state. |
| `"coherence_tagged"` | The same, with each line ending in the zx and zz hyperfine of the first bath spin and the name of the bath file ([Tagged rows](OutputFiles.md#out_tagged)). |
| `"cluster_coherence"` | Also one file pair per cluster, with its factor in the product ([One file per cluster](OutputFiles.md#out_clusters)). |
| `"cluster_density_matrix"` | As `"cluster_coherence"`, and with gCCE also the density matrix of every cluster at every time. For small systems only. |

```json
"output": {"savemode": {"mode": "coherence", "nodiv": false}}
```

<a id="otag_diag"></a>
## 2. Inspecting the clusters

These two keys do not change the coherence. They show the clusters of a run, to check its size before a long run
or to find the clusters that dominate the decay or make it diverge.

<a id="tag_cluster_only"></a>
<a id="tag_cluster_contributions"></a>
<a id="tag_cluster_contribution_max_mb"></a>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `cluster_only` | boolean | `false` | With `true` the run writes the list of clusters, `<file>.clusters.json`, and stops before the dynamics ([The cluster plan](OutputFiles.md#out_plan)). It works with cCCE (`standard`, `pcce`, or `bundle` with explicit `connections`) and with gCCE (`standard`, `bundle`). |
| `cluster_contributions` | `{"enabled", "max_mb"}` | off | Writes how much each cluster contributes, to find the clusters that dominate or diverge ([Cluster contributions](OutputFiles.md#out_contributions)). It works with cCCE and gCCE, `cluster.mode` standard or bundle, and `quantity` coherence. |
| `cluster_contributions.max_mb` | positive integer (MiB per rank) | `256` | The memory the contribution table may take. The run stops before the propagation if it needs more. |

`cluster_only` shows how many clusters an `rbath` and `rdip` give before a long run.
[Advanced 08](../../example/tutorial/advanced/08_analysis/README.md) uses `cluster_contributions`.

```json
"output": {"cluster_contributions": {"enabled": true, "max_mb": 512}}
```
