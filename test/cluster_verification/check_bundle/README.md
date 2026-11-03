# check_bundle/: bundle clustering inputs and records

This folder holds no runnable check. The checks of the bundle and adaptive code are `../check_adaptive/verify.sh`
(`unit`, `regress`, `inprocess`) and the unit tests of `../check_pcce/`.

Inputs: `simple_ensemble.json` (distance connections), `checked_ensemble.json` (`rdip` omitted), `ensemble.json`,
`auto_ensemble.json`, `gcce_single.json`, and the six-spin `bath_6` / `bath_hub_6`. `checked_ensemble.json` with
`bath_6` is also the input of the `checked_merged` stage of `../../code_verification/check_refusals/run.sh`.

Records of the 2026-09 experiments: `baths-20260912.json`, `integration-20260912.json` and `validation-20260912.json`.
The scripts that wrote them no longer run with the current build (old option names, removed binaries) and were
removed. They are in git history:

- `check_auto_groups.py` (`bundle.connections=auto` with `rdip` omitted, against the full graph and independent
  subsets) and `check_baths.py` (bundles on the full NV and PL6 baths across MPI ranks and backends):
  `git show 5583a8067:test/cluster_verification/check_bundle/<name>`.
- `check_adaptive.py`, `check_auto_connections.py` and `check_refinement.py`:
  `git show 8d656004d:test/cluster_verification/check_bundle/<name>`.
