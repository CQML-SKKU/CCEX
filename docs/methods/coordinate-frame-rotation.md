# The coordinate frame rotation

This page describes the `frame.rotation` block of the input file and the `frame` key beside each set of directional
numbers. The short form is the [frame.rotation](../user/InputFile_General.md#tag_frame_rotation) entry of General tags,
[Advanced 05: Coordinate rotation](../../example/tutorial/advanced/05_coordinate_rotation/README.md) runs it, and
`test/code_verification/check_rotation/` checks it against hand-turned inputs. The keys of earlier versions are in
[Upgrading old inputs](../user/InputFile_OldNames.md).

1. [Two frames](#rot_frames): the file frame and the computational frame.
2. [The rotation](#rot_rotation): the matrix, what it turns, and the order of the steps.
3. [The frame of each set of numbers](#rot_frame_key): the key `frame`.
4. [Limits](#rot_limits): what the rotation does not change.

<a id="rot_frames"></a>
## 1. Two frames

- **The file frame** is the Cartesian frame the bath file is written in, the lab or crystal frame, such as the cubic
  axes of diamond. Qubit positions (`Qubit.qubit[].xyz`, the position in a `Qubit.qubit[].file`) and the position
  columns of the tensor files are always in this frame.
- **The computational frame** is the frame the Hamiltonian is built in. Its z axis is the quantization axis of the
  qubit, and `bfield` is read in it.

Without `frame.rotation` the two are the same frame.

With it, CCEX turns the file frame so that `qubit_axis` becomes the computational z. `qubit_axis` is written in the
file frame, for example `[1, 1, 1]` for an NV in a cubic diamond bath. The magnetic field is then along the qubit axis:

- `"bfield": 500` (= `[0, 0, 500]`) is 500 G along `qubit_axis` in the lab.
- `[0, 0, -500]` is the same field reversed.
- A field with a transverse part is refused while the rotation is on.

<a id="rot_rotation"></a>
## 2. The rotation

The rotation $R$ is the passive change of basis with the rows

$$\mathbf e_z = \hat q, \qquad \mathbf e_x = \frac{\hat q \times \hat z}{\lvert \hat q \times \hat z\rvert}, \qquad
\mathbf e_y = \mathbf e_z \times \mathbf e_x,$$

with $\hat q$ the unit `qubit_axis` and $\hat z = [0, 0, 1]$ of the file frame, so $R\hat q = [0, 0, 1]$.

- Every position, qubits and bath spins alike, becomes $\mathbf r_\text{new} = \mathbf r_0 + R(\mathbf r_\text{old} -
  \mathbf r_0)$, with $\mathbf r_0$ the position of the first qubit.
- A tensor in the file frame becomes $R\,T\,R^T$.
- For `qubit_axis` = [1, 1, 1], $\mathbf e_x = [1, -1, 0]/\sqrt2$, $\mathbf e_y = [1, 1, -2]/\sqrt6$ and
  $\mathbf e_z = [1, 1, 1]/\sqrt3$.
- `qubit_axis` = [0, 0, 1] is no rotation. [0, 0, −1] leaves $\mathbf e_x$ undefined and is refused.

The run log prints $R$, the origin, every qubit in both frames, the field both ways round
($R^T [0, 0, B_z] = B_z \hat q$), and how many tensors were turned and kept.

The bath is read, the `rbath` selection made and the tensor files looked up (by position) before the rotation, in the
file frame. All of these depend only on distances to the qubit or on file-frame positions. The rotation then happens in
one step, before the Defect spins are built and the clusters made.

<a id="rot_frame_key"></a>
## 3. The frame of each set of numbers

Everything with a direction that is not a position says which frame it is written in, with the key `frame`:

| `frame` | Meaning | With the rotation on |
|---|---|---|
| `"file"` | Written in the file frame, before the rotation. | Turned with the positions ($R\,T\,R^T$, $R\mathbf v$). |
| `"qubit"` | Written in the computational frame (z = `qubit_axis`). | Used as written. |

| Where | The numbers it describes |
|---|---|
| [hamiltonian.common.hyperfine.frame](../user/InputFile_General.md#tag_hyperfine_frame) | The tensor components of `hyperfine.file`, needed with `"dipole": "dft"`. |
| [hamiltonian.common.quadrupole.frame](../user/InputFile_General.md#tag_quadrupole_frame) | The tensor components of `quadrupole.file` and `file_without_defect`. |
| [Qubit.intmap[].frame](../user/InputFile_Qubit.md#tag_intmap_frame) | That entry's `tensor`, a qubit's ZFS or a qubit pair's coupling. |
| [Defect.species[].frame](../user/InputFile_Defect.md#tag_df_frame) | The species' `axis`, `nuclear_positions`, `hyperfine`, `efg` and `zfs`. |

**With the rotation on, each of these that is written needs its `frame`.** There is no default, because a tensor
turned twice or not at all still looks like a tensor. With the rotation off a written `frame` is not used, and the run
says so once per key.

What CCEX computes itself is always in the file frame and always turned:

- The point-dipole hyperfine tensors, also with `"contact": true` (where the file gives only the isotropic term) and
  for a spin a DFT file does not cover.
- The point-dipole coupling between qubits.
- The bath-bath dipolar couplings, which come from the turned positions.

`hyperfine.frame` is therefore not used with `"dipole": "point"`, and the log notes it. Scalars do not turn, such as
the Defect `detuning` and `bfield`, which is already in the computational frame.

<a id="rot_limits"></a>
## 4. Limits

- The computational x axis is `qubit_axis` × [0, 0, 1], which sets the preparation phase of the qubit states relative
  to the system. With gCCE and a non-axial qubit ZFS (E ≠ 0) the coherence depends on this choice. cCCE and an axial
  ZFS do not.
- No block, or `"enabled": false`, leaves the run numerically identical to a run without the rotation.
- An external hyperfine or quadrupole file with more than one qubit is refused, and the rotation says so before a file
  is opened.
