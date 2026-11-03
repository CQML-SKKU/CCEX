# frame.rotation tests

    bash test/code_verification/check_rotation/run.sh        # needs bin/main.out (bash do_compile.sh wsl-gnu)

    # legacy baseline for the DFT tensor readers, which nothing else covers
    bash test/code_verification/check_rotation/dft_legacy_check.sh <reference-binary> [<test-binary>]

Exit status 0 only if every check passed. Everything the suite writes goes to `_work/`, which is wiped at the start
of each run and is gitignored.

## Where the frame changes

    read options -> CLI flags -> B FIELD ALIGNMENT CHECK
      -> read bath (source frame) -> bath_shift -> rbath selection (inner, outer)
      -> HF tensor lookup (source frame) -> QD tensor lookup (source frame)
      -> ROTATE positions and bath-frame tensors/Defect templates together
      -> Defect construction -> clusterization

`Config_validateBfieldAlignment` runs where `Config_resolveEvolution` does, once every reader has run. It never
modifies the field ("bfield" is already in the computational frame); it only refuses a transverse component while the
rotation is on.

Selection happens before the rotation on purpose: it is a distance to the reference qubit, and the rotation is taken
about that point, so it comes out the same either way. The tensor lookups have to come first, because those files
are keyed by position in the original frame.

## Which tensors move

Per (bath spin, qubit), not per run. `readHftensorfile` records whether each tensor it stored was computed here from
the source geometry or read from the file. `applyCoordinateFrameRotation` transforms the first kind always, and the
second kind only when `hamiltonian.common.hyperfine.frame` says `file`. The two occur together:

| `hyperfine` | what is stored | transformed? |
|---|---|---|
| `point` | point dipole | always |
| `point` + `contact` | point dipole + isotropic FC (the file's tensor is discarded) | always; `hyperfine.frame` is not used (a note says so) |
| `dft`, matched | DFT tensor from the file | only when `hyperfine.frame = "file"` |
| `dft`, fallback | point dipole (spin not covered by the file) | always |

The quadrupole keeps a single flag: `readQdtensorfile` has no computed-here fallback to distinguish. The structures
`strained` and `bubble` share this post-read path. Their two-file reader is exercised with distinct wiDef and woDef
boundaries; `plain` also has the full hand-transformed versus automatic-rotation tensor/coherence check.

## Why it is split in two

`test_rotation_unit.cpp` decides everything that does not need a calculation: the matrix, the fixed point, what the
transform preserves. It links against `obj/*.o` from the main build, so it always tests the same code the binary runs
(with `CCEX_BUILD_DIR=<build>`, the objects and flags of that CMake build instead; an object older than its source
stops the suite). `main.cpp` is compiled at link time by the top-level Makefile rather than into `obj/`, which is why
the test can define `rank` / `verbosity` / `nprocess` itself without a duplicate-symbol clash.

`run.sh` covers what only a real run can show: that a disabled rotation changes nothing, that the bath order and
labels and sampled states survive, that all bath files get the same transform, that the input files are not touched,
and that the refused combinations really exit non-zero with the message they claim.

`check_geometry.py` compares an unrotated log against a rotated one. It rebuilds R from the convention (ez = q,
ex = q x b, ey = ez x ex, as ROWS) instead of holding a table of expected numbers, so a sign flip in the C++ fails here
rather than being copied into the expected values. Positions come out of the BathArray report at 3 decimals; that,
and nothing else, sets its tolerances.

## The test systems

**Nuclear bath.** One qubit at `[10, 20, 30]` (not the origin, so that "r0 stays fixed" is a real check), 6 `13C` in
`bath_a` and 3 more in `bath_b`, all inside `rbath = 10`. One spin sits at `r0 + (1,1,1)`, exactly on the [111] axis,
so after the rotation it must land at `r0 + (0, 0, sqrt(3))`.

**P1 defect bath.** The same 6 positions as `P1` electron spins, with a `Defect` section: `n_nuclear_spins = 0`,
`n_configurations = 12`, a scalar per-axis `detuning`, and `avaax_p1` (`Defect.axis_file`) pinning the assignments so
no RNG is involved. It is checked twice:

- an **identity** rotation must leave the run byte-identical, which isolates the feature itself from any geometry
  change;
- a **[111]** rotation must move the bath (the coherence has to differ) while `paxes[i]` and the `detuning[iax]`
  table they index into stay exactly the same.

Per-spin detuning is `detuning[paxes[i]]`, so those two together pin it down; CCEX does not print the per-spin value
anywhere reachable at `n_nuclear_spins = 0`.

`n_nuclear_spins = 0` is deliberate: `nuclear_positions`, `hyperfine` and `efg` only exist for added subspins. This
P1 entry has no `axis`. With the rotation on it must state `frame: "qubit"`, and its scalar detuning is
frame-independent. Without the frame and with the rotation on it is refused (`defect_legacy_missing_frame`).

**NV defect bath.** Six spin-1 `NV_0` bath spins use three deterministically selected frozen configurations. One
physical `axis` belongs to `NV_0`; the three `n_configurations` entries are configurations, not three orientations.
The suite describes the same [111] NV five ways: shared `{D,E,unit}`, shared tensor, the old three-entry tensor list,
D in GHz, and an already transformed qubit-frame axis. It also expresses the same indexed detuning as a legacy MHz
list, a unit-aware MHz object, and a GHz object. Their numerical coherence outputs must agree. The suite checks that a
shared tensor is broadcast and that a bath-frame Defect is reported as transformed.

**Unit-aware on-site Defect.** One spin-1/2 `DUT` main spin has one on-site spin-1 `14N`, so all four wrapped
quantities reach a real Hamiltonian: `gyros`, `eqs`, `hypf`, and `efg`. A legacy input is compared with canonical unit
objects, then one unit is changed at a time. The suite requires numerical equivalence for MHz/GHz hyperfine,
radkHz/G versus MHz/T gyro, 1e-30 m^2 versus mbarn quadrupole moment, and Hartree/Bohr^2 versus V/angstrom^2 EFG.
Invalid units, value counts, and missing tensor frames must fail.

`bath_a`, `bath_b` and `bath_p1` are toy configurations written for this suite. They are not physical diamond
lattices and are not meant to be.

## bfield and defect_axis_reference

With the rotation on, the computational z axis is the qubit axis and "bfield" is read in that frame: `[0,0,500]` is
500 G along +qubit_axis, `[0,0,-500]` the same magnitude the other way, and both are accepted. A transverse component
is refused. `check_srcfield.py` reads the field back out of the log in both frames and requires
`R^T [0,0,Bz] = Bz * normalize(qubit_axis)`. The expected value is built from qubit_axis alone, so it does not just
re-derive what the code did.

`defect_axis_reference` is removed and refused ("was removed"), on the P1 system as anywhere. The rotation is
`frame.rotation` (`coordinate_frame_rotation` and its older alias `bath_coordinate_rotation` are refused). A written
block is on unless it says `"enabled": false`; the suite checks it against `"enabled": true` (byte-identical) and
against no block (`[disabled]`). The generators write the old flat names and `configuration_block.py`
(`../input_layout.py`) puts them into the blocks.

`frame.rotation` holds `qubit_axis` and `enabled` only. The computational x axis is `qubit_axis x [0,0,1]`.
`bath_axis`, `reference_qubit` and `qubit_position_frame` are refused as removed. The frames of the numbers are the
key `frame` (`"file"` | `"qubit"`) beside them, and the old names and the old value `"bath"` are refused with the new
ones.

## The A/B equivalence test

The strongest check here, and the only one that can catch a frame-dependent quantity nothing transforms.
`make_manual_frame.py` writes the same physical system out by hand in the qubit-aligned frame: bath positions about
r0, tensor file positions and boundary vertices about the origin (they are relative to the qubit), tensor blocks as
R T R^T. That input is run with no rotation option at all. It has to agree with the original input run with the
rotation on.

Numerical agreement is the criterion, not byte identity: the two reach the same numbers by different floating-point
routes, one parsing what the other computes. Byte identity is reported when it happens. Currently it does, for the
point-dipole system and for the DFT-tensor one alike.

## Multi-qubit

The rotation moves every qubit and every bath spin with one R about one r0, the position of the first qubit, which
therefore does not move. That is what makes the source-frame rbath selection valid for several qubits too: no
relative distance changes.

`Qubit.intmap` tensors carry their own provenance, and the rotation reads it rather than guessing from what kind of
tensor it is:

| provenance | source | transformed? |
|---|---|---|
| `INTMAP_DEFAULT_ZERO` | never written | no (zero either way) |
| `INTMAP_AUTO_SOURCE_GEOMETRY` | `setIntmap_dipAuto`, from the source qubit positions | always |
| `INTMAP_EXPLICIT_FILE` | explicit, `frame: "file"` (`Qubit.intmap[]`, the ZFS of the `Qubit.file` form included) | always |
| `INTMAP_EXPLICIT_QUBIT` | explicit, `frame: "qubit"` | never |
| `INTMAP_EXPLICIT_UNSPECIFIED` | explicit, no `frame`: only with the rotation off (with it on the reader refuses it) | - |

Each explicit `Qubit.intmap` entry may also carry `unit: "Hz"`, `"kHz"`, `"MHz"`, or `"GHz"`. An absent unit keeps the
legacy kHz interpretation. The suite feeds the same non-zero self tensor through all five spellings and requires
identical internal intmap tensors. Unknown and non-string units must fail before the simulation begins.

`validateCoordinateFrameRotationInputs` (reader.cpp) runs in main right after `Config_validateBfieldAlignment` (every
option source in, no reader started) and refuses what cannot be made consistent. The validation tests check the
message AND that neither `Read the Hyperfine interaction from DFT inputfile` nor
`Read the Quadrupole interaction from DFT inputfile` appears in the log. That second part would silently regress if
the check drifted later in the run.

External HF/QD tensor files stay single-qubit: `readHftensorfile` and `readQdtensorfile` refuse more than one qubit.
The rotation reports it before a file is opened.

## Turning about the qubit axis, and the qubit state

A rotation about the computational z axis is a symmetry of the Hamiltonian, but CCEX keeps |α⟩, |β⟩ and
psi0 = (|α⟩ + |β⟩)/norm in the computational frame: turned with the system they would pick up the relative phase
exp(−iφ(m<sub>β</sub> − m<sub>α</sub>)), another preparation phase. With a non-axial ZFS (E ≠ 0; the non-secular
hyperfine too, far more weakly) gCCE depends on that phase, so a run whose system is turned about z and whose
m<sub>s</sub> states are kept is a different experiment and differs (4.3e-4 in |L| for the free induction below; it is
also why the x axis of `frame.rotation`, a convention, matters there). cCCE keeps only the secular part and does not.

`check_zrot_gcce.py` turns the system by hand and the state with it where that phase is real (φ = 180° for +1/0, 90°
for +1/−1 and for two qubits (+1,+1)/(0,0)), as `alpha_state`/`beta_state` vectors, and requires the same |L| within
5e-8 (rounding amplified by the GHz phase; 360° gives the same 6e-10 to 8e-9). Cases: free induction, a Hahn echo (a π
pulse about −x is the one about x up to a sign) and two qubits with the pair dipole from the turned positions.
check_analytic compares the turned run with the exact answer.

## What is still refused

Multi-qubit runs with an external HF or QD tensor file (see above); any directional input without its `frame` (a
`Qubit.intmap` entry, the `Qubit.file` ZFS included, a `Defect` species, a DFT hyperfine file, a quadrupole file); and
a tilted magnetic field alongside the rotation. cCCE remains single-qubit (`simulator_cce.cpp`), so multi-qubit
coherence is tested with gCCE.

## What is read as written, and what is not

Four different things get lumped together as "the ZFS"; they are not the same:

| input | frame | transformed? |
|---|---|---|
| `Defect` species with `frame: "file"` | the bath file frame | axis/nuclear_positions and hyperfine/efg/zfs transformed |
| `Defect` species with `frame: "qubit"` | computational frame | no |
| `Defect.detuning` | scalar frequency | never |
| `bfield` | computational-frame input | never |
| the ZFS of the `Qubit.file` form (a `Qubit.intmap` self entry) | whatever its `frame` says | `"file"` yes, `"qubit"` no |
| `Qubit.intmap` self/ZFS tensor | whatever its `frame` says | `"file"` yes, `"qubit"` no |

With the rotation on each of these states its frame; with it off a written `frame` is not used and is warned about
once. A `Defect` section written in the crystal frame says `frame: "file"`. `applyCoordinateFrameRotation` then
transforms its normalized axis and relative `nuclear_positions` with R, and its `hyperfine`, `efg`, and `zfs` with
R T R^T before Defect construction. Scalar `detuning` is unchanged. `example/tools/InputKeys/convert_input_keys.py`
writes `"qubit"` where an older input with the rotation on relied on the hidden default.
