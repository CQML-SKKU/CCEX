# check_input_forms/

Checks that the new input forms give the physics of the legacy ones, and that each input option does what it says.
`../check_numerics/` holds the legacy form to v1.1.0 and v1.2.2; this suite holds the new form to the legacy form in
one binary. Every check compares two runs; only `evolution_update` and `tensor_format` read stored references.

```
bash test/code_verification/check_input_forms/run.sh          # every stage in order, stops at the first failure
bash test/code_verification/check_input_forms/run.sh axis2    # one stage
# bin/main.out or CCEX_BIN (bash do_compile.sh wsl-gnu); python3 with numpy, scipy, matplotlib (test/requirements.txt)
```

Everything goes to `_work/`, which git ignores. The exit status is 0 only if every stage passed.

## What it checks

Outputs are compared on the raw `_noDiv` and `_wiDiv`, never on the error-corrected file or the ensemble average.

| Stage | Compared with | Tolerance |
|---|---|---|
| `legacy` | writes the legacy input of each system | — |
| `axis2` | legacy form against new form; unit keys in the unit with factor 1 and in another | same bytes with factor 1; a deviation is reported otherwise |
| `errcorr` | the analyzer's `max\|L\|` before and after correction | reported |
| `gyro` | no gyrofile against a gyrofile with the table's value; a species in neither is refused | same bytes |
| `comments`, `nodiv`, `allfull`, `units`, `gpu_options` | `#`/`!` comments, `savemode` forms, `cluster_density_matrix` against `cluster_coherence` (1 and 3 ranks), `{"value", "unit"}` keys, the `gpu_options` block on a CPU run | same bytes |
| `evolution_update` | explicit `exact` against `../check_numerics` references (skipped on CI) | 1e-9 |
| | `recursive` against `exact`, cCCE2 and gCCE2 | 1e-8 |
| | no key against `recursive`; gCCE `expm` against `exact`; `backend` `auto` against `cpu` | same bytes |
| `column_evolution` | no block against `{"enabled": true, "min_dim": 256, "columns": 64}`; off against a `min_dim` no cluster reaches | same bytes |
| | default against off, one cluster of dimension 1296 | 1e-10 |
| `propagator` | `"eigen"` against the default; `"sparse"` (`min_dim` 4096) against `exact` | same bytes |
| | `sparse` with `min_dim` 1 against eigen `exact` | 1e-8, `_wiDiv` at \|L\| ≤ 1 |
| `tensor_format` | v3 hyperfine file (x 1) against the untagged v1 file (x 0.5), five h-BN cases | same bytes |
| | the VASP-tagged v1 file against `head-20261012-aa883d0d4` (skipped on CI) | 1e-12 relative |
| `cli` | `-e` template, its keys against the readers (`example_keys.py`, `example_autoconv.py`), no trailing comments and aligned `:`, the `-o` and `output.*` rules | as stated |
| `configuration` | each `k_range` configuration against its files alone; two bath files against the merged file; bad forms refused | same bytes |
| `autoconv` | `enabled` false against no block; enabled starts the runner once under `mpirun -n 2`; bad values refused | same bytes |
| `pcce_rest` | pCCE with left-over spins against files without them; with a `CCEX_SANITIZE=1` build in `CCEX_BIN`, a write into deleted spins stops the run | same bytes |
| `dsjitb` | cluster.mode dsjitb at order 2 against the standard CCE2 (the same clusters and curve, also on a bath whose first spin has inter-bath pairs) | every value within 1e-12, the same cluster counts |

The systems are reduced fixtures that run in seconds on one core: S1a/b NV + ¹³C CCE2 (ensemble and sampled), S2a
gCCE2, S2b/c NV + P1 10 ppm CCE2 and gCCE2, S3 P1 with ¹⁴N, S4 h-BN V_B⁻ with DFT tensors over `Afile_vertex`, S4b
over `Afile`, S5 two NV with `intmap`. `gen_variants.py` holds the reductions. Refusals are in `../check_refusals/`,
the frame rotation in `../check_rotation/`, the tensor writers in `example/tools/SpinTensorReader/tests/`.

## Files

| File | Contents |
|---|---|
| `run.sh` | runs the stages in order, or the one named |
| `lib.sh` | shared helpers and the run lock (also used by `../check_refusals/`) |
| `legacy.sh`, `gen_variants.py` | write the legacy inputs and the new-form variants beside them |
| other `*.sh` | one stage each, named as in the table |
| `example_keys.py`, `example_autoconv.py` | hold the `-e` template to the input readers and the convergence runner |
| `inputs/` | the NV + ¹³C and NV + P1 fixtures (cCCE and gCCE) |
