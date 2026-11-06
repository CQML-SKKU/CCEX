# check_numerics/

Runs twenty-seven short calculations cut from the tutorials and compares them with outputs that earlier builds
recorded on the same machine. Does this build still give the numbers an earlier one gave? Everything a tutorial
draws at random (bath states, P1 axes, ¹⁴N sublevels) is a file here, so two binaries can be compared.

```
bash test/code_verification/check_numerics/run.sh [NP]   # bin/main.out or CCEX_BIN (bash do_compile.sh wsl-gnu), mpirun, numpy
make cpu-test                                            # the C++ unit checks below (CMake: --target cpu-test)
```

`run.sh` writes to `_work/`, which git ignores.

## What it checks

| Case | Compared with | Tolerance |
|---|---|---|
| The 27 cases of `cases/` (NV + ¹³C, P1, Vp, h-BN; CCE2 to CCE4, CCE2e2n, pCCE2, gCCE, dsj, noise; Hahn, Ramsey, CPMG-8, `sequence`) | the newest `reference/head-*` by name, which decides | 1e-9 on the largest \|difference\| of `_noDiv` and density-matrix files |
| The same, binary built with MKL (the run log says `MKL for large clusters : on`) | the newest `reference/head-*` | 1e-12 relative per element, 1e-12 absolute where the reference is 0 |
| The same | the tagged releases `v1.1.0/`, `v1.2.2/`, `v2.0.1/` | 1e-4 for v1.1.0 and v1.2.2, 1e-9 for v2.0.1, and cases named in `expected.json` are reported as "expected" |
| The same | the older `head-*` folders | 1e-9, printed only |
| `test_fixed_size.cpp`: fixed-size cCCE types (bath dimension 2, 4), `exact`, 2,000 random pairs | the dynamic path | same bits (sign of zero too) on icpx 2023.2 AVX-512, icpx 2026.1 AVX2, g++ 13 alderlake, elsewhere 1e-12 relative or 1e-15 absolute |
| `test_cce_batch.cpp`: batched CPU evolution, bath dimension 2 and 4, on `../../gpu_verification/check_gpu/data/small_cluster/` | the scalar path | same bits, and on other builds the kernel is off ("scalar only") |
| `test_cce_columns.cpp`: column blocks (`column_evolution`), dimension 300 and 600 | the one-time-point path | 1e-12 relative, `exact` phases same bits, `recursive` exact at anchors (between them 1e-11 of `exact` and 1e-14 of the old path) |
| `test_small_eigensolver.cpp`: the small cCCE Jacobi | Eigen | 1e-13 of max\|H\|, same bits alone and in groups |

The `_wiDiv` coherence is printed but not judged, because it divides by the free evolution and amplifies rounding. A case
the deciding reference does not have fails. Against every other reference a case added after it is "new". A case
with no output, or with a reference file the run did not write, fails whatever `expected.json` says. `compare.py`
reads the run logs to know whether the MKL rule applies. `CCEX_CPU_BATCH=0` turns the batch off in a run.

## Updating a reference

Only the newest `head-*` decides. When a change is meant to move a number, record a new one after it with
`bash make_reference.sh current` (this branch, 2 ranks, `CCEX_BLAS=off`) and give the reason in the commit. The name
carries date and time so references of one day sort in order. A tag's reference is made with
`bash make_reference.sh v2.0.1` (a worktree of the tag under /tmp or `CCEX_WORKTREE_DIR`). `expected.json` records
why each case differs from each older reference, and maps renamed density-matrix files (`renames`). Keep
`head-20261012-aa883d0d4`, because `../check_input_forms/tensor_format.sh` needs it. Do not regenerate a reference to
make the check pass.

## Files

| File | Contents |
|---|---|
| `run.sh` | Runs every case and compares it with every folder in `reference/`. |
| `compare.py` | The comparison and the tables. |
| `expected.json` | The known differences per case and reference, and the file renames. |
| `make_cases.py` | Wrote the cases and their input files from the tutorials. |
| `make_reference.sh` | Records a reference from this branch or from a tag. |
| `cases/` | One `ccein.json` per case (`hbn_gcce2_oh` also has `ccein_tags.json`, its input in the tags' form). |
| `inputs/` | The bath, gyro, tensor, state and axis files the cases read. |
| `reference/` | The recorded outputs, one folder per build. |
| `test_fixed_size.cpp`, `test_cce_batch.cpp`, `test_cce_columns.cpp`, `test_small_eigensolver.cpp` | The C++ unit checks of `make cpu-test`. |
