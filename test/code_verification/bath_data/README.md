# bath_data/

The frozen input files the checks read. They are data, not checks. Real bath
geometries and DFT tensors cannot be made on the spot and must not change
between runs, so one copy lives here and the checks point at it
(`check_numerics/` keeps its own `inputs/`, cut from the tutorials together
with its cases).

| Folder | Holds | Read by |
|---|---|---|
| `ensemble_data/` | h-BN and MoS₂ ensemble baths, their defect files, gyromagnetic ratios and DFT tensors (5.6 MB). | `check_rotation`, `check_input_forms` (S4), `cluster_verification/check_pcce/check_baths.py` |
| `hexagonal_data/` | One h-BN hexagonal cell with its `Afile`/`Qfile` tensors, the `_vertex` boundary versions and the structure they came from (`CONTCAR`, `POSCAR_test`) (1.5 MB). | `check_rotation` (`run.sh`, `dft_legacy_check.sh`), `check_input_forms` (S4, S4b) |
| `single_data/` | The P1 1 ppm bath with its defect file, gyromagnetic ratios, and two sets of bath-state and ¹⁴N-state files. | `check_input_forms` (S3), `check_refusals` |
| `FullP1SpinBath/cce.json` | The legacy P1 input `check_input_forms` cuts S3 from. | `check_input_forms/legacy.sh` |
| `NuclearSpinBath/cce.in` | A legacy nuclear-spin-bath input (MoS₂). | `cluster_verification/check_pcce/check_baths.py` |

These came from the original CCEX repository's reproduction data
(`CCE_Reprod/Bath_Data`). The numbered folder names were dropped on 2026-09-24
and the unused `3.semi_class_data` with them.
