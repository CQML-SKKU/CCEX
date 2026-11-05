# Advanced 09: Coordinate rotation

Folder: `example/tutorial/advanced/09_coordinate_rotation/`

CCEX builds the Hamiltonian in a frame whose z is the qubit's axis, and reads `bfield` in that frame. Without a
rotation that frame is the bath file's own, so a bath written in a frame where the defect axis is not z has to be
turned, by hand or with the `frame.rotation` block. The block turns every position about the qubit so that
`qubit_axis` (written in the bath file's frame) becomes z, and `"bfield": 500` is then 500 G along `qubit_axis` in the
lab. The block is on when it is written.

Tensors and Defect species written in the input say with the key `"frame"`
whether they are in the bath file's frame (`"file"`, turned with the positions) or already in the turned one
(`"qubit"`). This example has none, since the NV's zero-field splitting is not in the input.
[The coordinate frame rotation](../../../../docs/methods/coordinate-frame-rotation.md) explains what moves and what
does not. Every key not written in an input has its default (`main.out -e`). Times are for 4 ranks of a desktop CPU.

`bath/bath_Diamond_NV_natab_cubic_1` is the bath of Basic 01 turned so that the NV axis points along [111], as a bath
written in the conventional cubic cell would be. `bath/make_cubic_frame.py` made it from `basic/01_first_run/bath/`
and says how.

1. [cCCE2_Hahn_rotation](#rot_case): the same NV three ways.

<a id="rot_case"></a>
## 1. cCCE2_Hahn_rotation: the same NV three ways

- **What it shows:** the Hahn echo of the NV of Basic 01, computed three ways. The first uses the bath as generated,
  with the NV axis along z (`ccein_aligned.jsonc`). The second uses the cubic-frame bath with the rotation on
  (`qubit_axis` [1, 1, 1], `ccein_cubic_rotated.jsonc`). The third uses the cubic-frame bath without it
  (`ccein_cubic_norotation.jsonc`), so the field is 54.7° off the NV axis.
- **Keys that matter:** `frame.rotation.qubit_axis` [1, 1, 1] and `configuration.bathfile` (the cubic-frame file).
- **How to run:** `cd cCCE2_Hahn_rotation && NP=4 bash run.sh`. It takes 1.8 s and 0.05 GB per rank.
- **Expected result:** `expected.png`. The rotated run reproduces the aligned one (largest difference 1e-10, T2
  1.084 ms). Without the rotation T2 is 1.154 ms and the curve differs by up to 0.20.
- **What to change next:** `bfield` as a vector to tilt the field on purpose, without the rotation (with it a
  transverse field is refused), or a bath with DFT tensors (`hamiltonian.common.hyperfine.frame`,
  `quadrupole.frame`).

![cCCE2_Hahn_rotation: expected result](cCCE2_Hahn_rotation/expected.png)

Inputs: [`ccein_aligned.jsonc`](cCCE2_Hahn_rotation/ccein_aligned.jsonc), [`ccein_cubic_norotation.jsonc`](cCCE2_Hahn_rotation/ccein_cubic_norotation.jsonc), [`ccein_cubic_rotated.jsonc`](cCCE2_Hahn_rotation/ccein_cubic_rotated.jsonc)
