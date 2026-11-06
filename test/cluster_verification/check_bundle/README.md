# check_bundle/

Inputs for bundle clustering and the records of the 2026-09 experiments. This folder holds no runnable check. A
change to `src/clustering/bundle/` is checked by the bundle unit tests in `../check_pcce/` and by the recorded runs of
`../check_adaptive/`, which use the bundle partition.

```
make cpu-test                                                         # CMake: --target cpu-test; the three bundle unit tests
sh test/cluster_verification/check_adaptive/verify.sh regress         # then inprocess
```

## What it checks

| Check | Where |
|---|---|
| Bundle partition, connections and high-order expansion. | `../check_pcce/test_adaptive_partition.cpp`, `test_adaptive_connections.cpp`, `test_adaptive_high_order.cpp` |
| Bundle runs inside the adaptive policies. | `../check_adaptive/verify.sh` `regress` and `inprocess` |
| Bundle with `rdip` omitted and multi-site bundles must finish (`checked_ensemble.json` with `bath_6`). | the `checked_merged` case of `../../code_verification/check_refusals/run.sh` |

The scripts that wrote the records no longer run with the current build and were removed. They are in git history.

## Files

| File | Contents |
|---|---|
| `simple_ensemble.json` | Bundle input with distance connections. |
| `checked_ensemble.json` | Bundle input with `rdip` omitted. |
| `ensemble.json`, `auto_ensemble.json`, `gcce_single.json` | Further bundle inputs. |
| `bath_6`, `bath_hub_6` | Six-spin baths. |
| `baths-20260912.json`, `integration-20260912.json`, `validation-20260912.json` | Records of the 2026-09 experiments. |
