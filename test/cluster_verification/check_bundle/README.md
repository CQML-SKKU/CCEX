# check_bundle/ — the bundle clustering against complete calculations

The records of the 2026-09-12..16 experiments. The current checks of the bundle and adaptive code are
`../check_adaptive/verify.sh` (`unit`, `regress`, `inprocess`) and the unit tests of `../check_pcce/`.

Removed on 2026-11-03 (last in commit 5583a8067: `git show 5583a8067:test/cluster_verification/check_bundle/<name>`), the two
scripts behind the records below: `check_auto_groups.py` (the automatic connection choice, `bundle.connections=auto`
with `rdip` omitted, against the full graph and independent subsets) and `check_baths.py` (bundles on the full NV and
PL6 baths across MPI ranks and backends). They set the option names of that time (`adaptive_cluster_max_dim`, ...;
since D24 an unread key near a known one stops the run), ran the binaries `bin/main.cpu.out` / `bin/main.cuda.out` the
build no longer writes, and `check_auto_groups.py` called the removed oracle `check_adaptive/oracle/adaptive_refine.py`:
neither ran with the current build.

Removed on 2026-10-17 (A96; last in commit 8d656004d: `git show 8d656004d:test/cluster_verification/check_bundle/<name>`):
`check_adaptive.py` (adaptive CCE/gCCE against a baseline of full-bath and subset runs),
`check_auto_connections.py` (`bundle.connections=auto` against the full graph) and `check_refinement.py`
(end-to-end refinement against complete small-bath calculations). Nothing ran them; they set keys refused since A47
(`npulse`, ...), and `../check_adaptive/verify.sh` checks the same code now.

Inputs: `simple_ensemble.json` (distance connections), `checked_ensemble.json`
(`rdip` omitted), `ensemble.json`, `auto_ensemble.json`, `gcce_single.json`, and
the six-spin `bath_6` / `bath_hub_6`. `checked_ensemble.json` with `bath_6` is also the input of the
`checked_merged` stage of `../../code_verification/check_refusals/run.sh` (A96).

Kept from the 2026-09 experiments: the records `baths-20260912.json`, `integration-20260912.json` and `validation-20260912.json`.
