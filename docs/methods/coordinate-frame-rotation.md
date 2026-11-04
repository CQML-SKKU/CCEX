# The coordinate frame rotation

This page describes the `frame.rotation` block of the input file and the `frame` key beside each set of directional
numbers. The short form is the [frame.rotation](../user/InputFile_General.md#tag_frame_rotation) entry of the
input-file reference; [tutorial 10](../../example/tutorial/advanced/05_coordinate_rotation/README.md) runs it.

## Two frames

<b>The file frame</b> is the Cartesian frame the bath file is written in: the lab or crystal frame, e.g. the cubic axes
of diamond. Qubit positions (`Qubit.qubit[].xyz`, the position in a `Qubit.qubit[].file`) and the position columns of
the tensor files are always in this frame. <br>
<b>The computational frame</b> is the frame the Hamiltonian is built in. Its z axis is the quantization axis of the
qubit, and `bfield` is read in it. <br>
Without `frame.rotation` the two are the same frame.

With it, CCEX turns the file frame so that `qubit_axis` becomes the computational z. `qubit_axis` is written in the
file frame, e.g. `[1, 1, 1]` for an NV in a cubic diamond bath. The magnetic field is then along the qubit axis:

- `"bfield": 500` (= `[0, 0, 500]`) is 500 G along `qubit_axis` in the lab;
- `[0, 0, -500]` is the same field reversed;
- a field with a transverse part is refused while the rotation is on.

## The rotation

The rotation R is the passive change of basis with rows

<p align="center">e<sub>z</sub> = q̂, &nbsp;&nbsp; e<sub>x</sub> = (q̂ × ẑ)/|q̂ × ẑ|, &nbsp;&nbsp; e<sub>y</sub> = e<sub>z</sub> × e<sub>x</sub></p>

with q̂ the unit `qubit_axis` and ẑ = [0, 0, 1] of the file frame, so R q̂ = [0, 0, 1].

- Every position, qubits and bath spins alike, becomes r<sub>new</sub> = r<sub>0</sub> + R (r<sub>old</sub> −
  r<sub>0</sub>), with r<sub>0</sub> the position of the first qubit.
- A tensor in the file frame becomes R T R<sup>T</sup>.
- For `qubit_axis` = [1, 1, 1]: e<sub>x</sub> = [1, −1, 0]/√2, e<sub>y</sub> = [1, 1, −2]/√6,
  e<sub>z</sub> = [1, 1, 1]/√3.
- `qubit_axis` = [0, 0, 1] is no rotation. [0, 0, −1] leaves e<sub>x</sub> undefined and is refused.

The run log prints R, the origin, every qubit in both frames, the field both ways round
(R<sup>T</sup> [0, 0, B<sub>z</sub>] = B<sub>z</sub> q̂), and how many tensors were turned and kept. <br>
The bath is read, the `rbath` selection made and the tensor files looked up (by position) before the rotation, in the
file frame. All of these depend only on distances to the qubit or on file-frame positions. The rotation then happens in
one step, before the Defect spins are built and the clusters made.

Old inputs: before 2.0.0 the block had five more keys. They are refused, and the message names what replaced them.

- `bath_axis` took the place of ẑ above; it turned e<sub>x</sub> and e<sub>y</sub> about the qubit axis.
- `reference_qubit` was the origin r<sub>0</sub>; another qubit only translates every position.
- `qubit_position_frame` had one value, `"bath"`.
- `hf_tensor_frame` and `qd_tensor_frame` are now the `frame` of the hyperfine and quadrupole blocks.

`example/tools/InputKeys/convert_input_keys.py` rewrites an old input. <br>
An input that used a `bath_axis` other than [0, 0, 1] now gets the x axis of the default, and its result can change.
For the axially symmetric cases of the test suites (cCCE, an NV with E = 0) the coherence is the same to the last digit
for any `bath_axis`. For gCCE with a qubit tensor that is not symmetric about `qubit_axis`, `bath_axis` moved the
coherence by up to about 1e-2 relative; turning the whole system about z by hand, without the rotation, does the same.

## The frame of each set of numbers

Everything with a direction that is not a position says which frame it is written in, with the key `frame`:

| `frame` | meaning | with the rotation on |
|---|---|---|
| `"file"` | written in the file frame, before the rotation | turned with the positions (R T R<sup>T</sup>, R v) |
| `"qubit"` | written in the computational frame (z = `qubit_axis`) | used as written |

| where | numbers it describes |
|---|---|
| [hamiltonian.common.hyperfine.frame](../user/InputFile_General.md#tag_hyperfine_frame) | the tensor components of `hyperfine.file` (needed with `"dipole": "dft"`) |
| [hamiltonian.common.quadrupole.frame](../user/InputFile_General.md#tag_quadrupole_frame) | the tensor components of `quadrupole.file` and `file_without_defect` |
| [Qubit.intmap[].frame](../user/InputFile_Qubit.md#tag_intmap_frame) | that entry's `tensor` (a qubit's ZFS, in both forms of the `Qubit` block, or a qubit pair's coupling) |
| [Defect.species[].frame](../user/InputFile_Defect.md#tag_df_frame) | the species' `axis`, `nuclear_positions`, `hyperfine`, `efg` and `zfs` |

<b>With the rotation on, each of these that is written needs its `frame`.</b> There is no default, because a tensor
turned twice or not at all still looks like a tensor. With the rotation off a written `frame` is not used, and the run
says so once per key. <br>
For old inputs: before 2.0.0 an absent frame meant `"qubit"` for a Defect species, a single qubit's intmap entry and
`Qubit.zfs` (now an intmap self entry). The converter writes that value out for an old input with the rotation on.

What CCEX computes itself is always in the file frame and always turned:

- the point-dipole hyperfine tensors, also with `"contact": true` (where the file gives only the isotropic term) and
  for a spin a DFT file does not cover;
- the point-dipole coupling between qubits;
- the bath-bath dipolar couplings, which come from the turned positions.

`hyperfine.frame` is therefore not used with `"dipole": "point"` (the log notes it). Scalars do not turn: Defect
`detuning`, and `bfield`, which is already in the computational frame.

## Limits

- The computational x axis is `qubit_axis` × [0, 0, 1], which sets the preparation phase of the qubit states relative
  to the system. With gCCE and a non-axial qubit ZFS (E ≠ 0) the coherence depends on it: with the system turned by
  180°, |L| changed by 4.3e-4 for the free induction and by 3.9e-2 for a Hahn echo. cCCE and an axial ZFS do not
  depend on it (gCCE: about 1e-8, from the non-secular hyperfine).
- No block, or `"enabled": false`, leaves the run numerically identical to a run without the rotation.

- An external hyperfine or quadrupole file with more than one qubit is refused. The readers have always refused it;
  the rotation says so before a file is opened.
- cCCE is single-qubit; several qubits are tested with gCCE.

The checks are `test/code_verification/check_rotation` (the geometry, hand-turned inputs against the rotation, the
frames) and `check_refusals`.
