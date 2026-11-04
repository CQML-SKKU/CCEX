# Output tags

The `output` block: where the outputs go, what is written, and two diagnostics. Each key has a small table (its type,
default and unit), what it does, and an example. What each file holds is in [Output files](OutputFiles.md). The keys
of this block are refused at the top level, with their place in the block.

```json
"output": {"path": "results", "file": "ccein", "savemode": {"mode": "coherence", "nodiv": true}}
```

1. [Where and what](#otag_where): `path`, `file`, `savemode`.
2. [Diagnostics](#otag_diag): `cluster_only`, `cluster_contributions`.

<a id="otag_where"></a>
## 1. Where and what

<a id="tag_output_path"></a>
### output.path

| Type | Default |
|---|---|
| string (a folder) | the working directory |

The folder of the outputs, made when it is missing. A relative folder is read from the working directory. `-o DIR` on
the command line takes its place. At the start the run writes and removes a file there; a folder that cannot be
written is refused before anything is calculated.

Refused:

- together with an absolute [output.file](#tag_outfile);
- with a `cluster.adaptive` block, whose steps go to `cluster.adaptive.workdir`;
- with an `auto_convergence` block, whose runs go to `auto_convergence.workdir`.

```json
"output": {"path": "results"}
```

<a id="tag_outfile"></a>
### output.file

| Type | Default |
|---|---|
| string | the input file's name without `.json` or `.jsonc` |

The start of every output name (`<file>_wiDiv`, `<file>_noDiv`, ...; [Names](OutputFiles.md#out_names)). A relative
name is read from the output folder: `-o DIR`, else [output.path](#tag_output_path), else the working directory. With
`-f run/ccein.json` and no `file`, the outputs are `ccein_wiDiv`, `ccein_noDiv`, ... in the output folder.

Refused: an absolute name together with `-o` or `output.path`, and any value with an `auto_convergence` block.

```json
"output": {"file": "run1"}
```

<a id="tag_savemode"></a>
### output.savemode

| Type | Default |
|---|---|
| string, or `{"mode": string, "nodiv": boolean}` | `"coherence"`, with `nodiv` true |

What is written. A string gives `mode`; the object gives both keys, a missing `mode` being `"coherence"` and a missing
`nodiv` being `true`.

| mode | Writes |
|---|---|
| `"coherence"` | the coherence: `<file>_wiDiv` and `_noDiv`, one pair per sampled state with `s_range` |
| `"coherence_tagged"` | the same files, each line ending with the zx and zz hyperfine components of the first bath spin and the first bath file's name, to tell the results of many bath files apart ([Tagged rows](OutputFiles.md#out_tagged)) |
| `"cluster_coherence"` | also one file pair per cluster and state: the factor of that cluster in the product ([One file per cluster](OutputFiles.md#out_clusters)) |
| `"cluster_density_matrix"` | as `"cluster_coherence"`, and with gCCE also the density matrix of the whole qubit + cluster space of every cluster at every time step (`<file>_<cluster>_full_...`); for small systems only. cCCE writes what `"cluster_coherence"` writes |

`nodiv` says whether `_noDiv` is written beside `_wiDiv`. `_noDiv` is the product of the cluster terms with a positive
power only; it never diverges and is not a coherence. The analyzer's error correction reads it
([Analysis](Analysis.md)). `false` writes `_wiDiv` alone.

```json
"output": {"savemode": "cluster_coherence"}
"output": {"savemode": {"mode": "coherence", "nodiv": false}}
```

<a id="otag_diag"></a>
## 2. Diagnostics

<a id="tag_cluster_only"></a>
### output.cluster_only

| Type | Default |
|---|---|
| boolean | `false` |

`true` writes the list of clusters, `<file>.clusters.json`, and stops before the dynamics
([The cluster plan](OutputFiles.md#out_plan)). Use it to see how many clusters an `rbath` and `rdip` give before a
long run. It works with cCCE and gCCE, and with cCCE `cluster.mode` pcce. Refused with a `bundle` input that leaves
`rdip` out (its connections need the dynamics), with a `cluster.adaptive` block and with an `auto_convergence` block.

```json
"output": {"cluster_only": true}
```

<a id="tag_cluster_contributions"></a>
### output.cluster_contributions

| Type | Default |
|---|---|
| `{"enabled": boolean, "max_mb": integer}` | off |

Writes how much each cluster contributes to the coherence over all times, to find the clusters that dominate the decay
or diverge. cCCE and gCCE with `quantity` coherence. A block without `enabled` is on; a plain boolean is refused. The
files are `<file>.contributions.state-<s>.tsv`, one row per cluster, and `<file>.contributions.sites.tsv`
([Cluster contributions](OutputFiles.md#out_contributions)).
[Advanced 08](../../example/tutorial/advanced/08_analysis/README.md) uses it.

```json
"output": {"cluster_contributions": {"enabled": true}}
```

<a id="tag_cluster_contribution_max_mb"></a>
### output.cluster_contributions.max_mb

| Type | Default | Unit |
|---|---|---|
| positive integer | `256` | MiB per rank |

The memory per rank the contribution table may take. The run stops before the propagation if the table would need
more.

```json
"output": {"cluster_contributions": {"enabled": true, "max_mb": 512}}
```
