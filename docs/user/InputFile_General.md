# General tags
Tags in this section are about the simulation parameters: the method, the numbers, the files, and the
`hamiltonian` block. The keys before 2026-10-13 (A47) and their new names are in the table of
[InputFile.md](InputFile.md#OldNames).

- [method](#tag_method)       
- [cluster_mode](InputFile_Cluster.md#tag_cluster_mode)
- [quantity](#tag_quantity)     
- [order](#tag_order)        
- [bfield](#tag_bfield)       
- [rbath](#tag_rbath)        
- [rdip](#tag_rdip)         
- [deltat](#tag_deltat)       
- [nstep](#tag_nstep)        
- [rdip_ladder](#tag_rdip_ladder)
- [nstate](#tag_nstate)        
- [seed](#tag_seed)       
- [gyrofile](#tag_gyrofile)      
- [bathfile](#tag_bathfile)      
- [bath_shift](#tag_bath_shift)
- [bath_state_file](#tag_bath_state_file)
- [coordinate_frame_rotation](#tag_coordinate_frame_rotation)
- [defect_axis_reference](#tag_defect_axis_reference)
- [hamiltonian](#tag_hamiltonian)
    - common: [qubit_bath_mean_field](#tag_qubit_bath_mean_field), [bath_bath_mean_field](#tag_bath_bath_mean_field),
      [hyperfine](#tag_hyperfine) (`dipole`, `contact`, [file](#tag_hyperfine_file), [missing_row](#tag_hyperfine_missing_row)),
      [quadrupole](#tag_quadrupole) (`structure`, [file](#tag_quadrupole_file),
      [file_without_defect](#tag_quadrupole_file_without_defect), [cell](#tag_quadrupole_cell))
    - cCCE: [knight](#tag_knight), [hfmedi](#tag_hfmedi), [delta_A](#tag_delta_A) (meCCE reads them too)
    - meCCE, megCCE: [jump_operators](#tag_jump_operators)
    - gCCE, noise: no settings yet (noise is a reserved method)

The qubit file (`Qubit.file`) is on the [Qubit page](InputFile_Qubit.md#tag_Qubit_file); the files of the defect bath
(`Defect.axis_file`, `Defect.state_file`) are on the [Defect page](InputFile_Defect.md#tag_df_axis_file).


<a id="tag_method"></a>
## method
<b>Description:</b><br>
The CCE method: conventional (cCCE, the qubit projected out into its two branches) or generalized (gCCE, the qubit in the
cluster Hamiltonian), and each of them with Lindblad dissipation by the jump operators of
[hamiltonian.\<method\>.jump_operators](#tag_jump_operators) (meCCE, megCCE; [master-equation.md](../methods/master-equation.md)).
How the bath is cut into clusters is [cluster_mode](InputFile_Cluster.md#tag_cluster_mode)
(`standard`, `pcce`, `bundle`, `dsj`, `itb`, `dsjitb`). `pcce` runs with cCCE, gCCE, meCCE and megCCE; `bundle` with cCCE
and gCCE; `dsj`, `itb` and `dsjitb` run with cCCE only: gCCE, meCCE and megCCE with them are refused as not implemented.

You can give the method with the "-m" command-line option (`-m cCCE`, `-m gCCE`, `-m meCCE`, `-m megCCE`); it keeps the
input's `cluster_mode`. <br>
<b>Data Type:</b> String (case is not significant)

<b>Options:</b><br>
| Option  | Description |
|---------|-------------------------------------|
| `cCCE`  | Conventional CCE |
| `gCCE`  | Generalized CCE |
| `meCCE` | ME-CCE: the cCCE with Lindblad jump operators on the bath spins (PRL 132, 250401 (2024)); qubit jumps are refused (use megCCE). CPU, one qubit, ideal π pulses from `pulse.count` |
| `megCCE` | ME-gCCE: the gCCE with Lindblad jump operators on the bath spins and on the qubit (C<sub>q</sub> ⊗ 1, in every cluster and the 0th). CPU; the qubit sequence as given (Manual and any rotation included) with `pulse.per_qubit` false and `pulse.detuning_factor` 1 (others refused); DEER bath pulses (`pulse.bath`) run with it |
| `noise` | reserved: refused as not implemented yet (2026-10-13) |

The values before 2026-10-13 are refused with their new form: `CCE` is `cCCE`; `pCCE`, `dsj`, `itb` and `dsjitb` are
method `cCCE` with that `cluster_mode` (also on the command line: `-m pCCE` is refused); the adding_MECCE branch's
`pgCCE`, `pmeCCE` and `pmegCCE` are method `gCCE`, `meCCE` and `megCCE` with `"cluster_mode": "pcce"`.

<b>Default:</b> `cCCE`

<b>Example:</b> ``` "method" : "gCCE" ```, ``` "method" : "cCCE", "cluster_mode" : "pcce" ```

<b>Related variable:</b> [Config](#Config)::method (inside: CCE, gCCE, pCCE, dsj, itb, dsjitb, meCCE, megCCE, pgCCE, pmeCCE, pmegCCE, from method and cluster_mode)

<a id="tag_quantity"></a>
## quantity
<b>Description:</b><br>
Quantity to calculate the dynamics of the system.

You can specify this value via the "-q" command-line option. <br>
<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `coherence` | Coherence function             |
| `dm`        | Density matrix                 |

Any other value is refused, `noise` included: the parser used to accept it, but there is no noise calculation.

With `dm` (gCCE), the reduced density matrix of the qubits is written as blocks of rows, one
block per time step, and the cluster expansion is applied element by element. An element the
free evolution keeps at zero (a level the initial state does not touch) is left at zero; when
clusters would put a value into it, the run prints a `[dm]` line with the largest such value,
since the elementwise product cannot carry it. Until 2026-09-24 those elements came out nan. <br>
<b>Unit:</b> arbitrary unit (a.u.)

<b>Default:</b> `coherence`

<b>Example:</b> ``` "quantity" : "coherence" ```

<b>Related variable:</b> [Config](#Config)::quantity

<a id="tag_order"></a>
## order
<b>Description:</b><br>
Order of the CCE method: the largest cluster, in bath spins (`cluster_mode` standard), in partitions (pcce) or in bundles
(bundle).

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "order" : 2 ```
<b>Related variable:</b> [Config](#Config)::order

<a id="tag_bfield"></a>
## bfield
<b>Description:</b><br>
Magnetic field.

You can override the z component with the "-B" command-line option, but the key itself must be
in the input file: `-B` is applied after the file is read and does not stand in for a missing
`bfield`. For no field, write `"bfield": 0`. <br>
<b>Data Type:</b> Float or Float List

<b>Unit:</b> gauss (G)

<b>Default:</b> none (required)

<b>Example:</b>

``` "bfield" : 0.0 ```

``` "bfield" : [0.0, 0.0, 0.0] ```

<b>Related variable:</b> [Config](#Config)::bfield

<a id="tag_rbath"></a>
## rbath
<b>Description:</b><br>
Bath radius: the bath spins within it of the qubit are kept. A number is the outer radius; `[inner, outer]` also removes the
spins closer than `inner` (the inner radius was the key `rbathcut` until 2026-10-13). <br>

<b>Data Type:</b> Float, or a list of two floats `[inner, outer]`

<b>Unit:</b> angstrom

<b>Default:</b> - (required); inner 0

<b>Example:</b> ``` "rbath" : 50.0 ```, ``` "rbath" : [2.0, 50.0] ```

<b>Related variable:</b> [Config](#Config)::rbath, [Config](#Config)::rbathcut (the inner radius)

<a id="tag_rdip"></a>
## rdip
<b>Description:</b><br>
Dipolar radius: two bath spins closer than it are connected and can be in one cluster. A number is the outer radius;
`[inner, outer]` also leaves out the pairs closer than `inner` (the inner radius was the key `rdipcut` until 2026-10-13).
With `cluster_mode` bundle it may be left out (the connections are then chosen by a coherence check, see
[cluster_mode](InputFile_Cluster.md#tag_cluster_mode)); an inner radius then needs the outer one. <br>

<b>Data Type:</b> Float, or a list of two floats `[inner, outer]`

<b>Unit:</b> angstrom

<b>Default:</b> - (required, but with `cluster_mode` bundle); inner 0

<b>Example:</b> ``` "rdip" : 10.0 ```, ``` "rdip" : [1.0, 10.0] ```

<b>Related variable:</b> [Config](#Config)::rdip, [Config](#Config)::rdipcut (the inner radius)

<a id="tag_deltat"></a>
## deltat
<b>Description:</b><br>
Time interval.

<b>Data Type:</b> Float

<b>Unit:</b> ms

<b>Default:</b> -

<b>Example:</b> ``` "deltat" : 0.001 ```

<b>Related variable:</b> [Config](#Config)::deltat

<a id="tag_nstep"></a>
## nstep
<b>Description:</b><br>
Number of time steps.

<b>Data Type:</b> Integer

<b>Default:</b> -

<b>Example:</b> ``` "nstep" : 10 ```

<b>Related variable:</b> [Config](#Config)::nstep

<a id="tag_rdip_ladder"></a>
## rdip_ladder
<b>Description:</b><br>
Extra rdip values (rungs) computed from the same run. With rungs r<sub>1</sub> < r<sub>2</sub> < ... below `rdip`,
the run also writes, for every rung, the coherence that a separate run at `rdip` = r<sub>k</sub> gives
(`<outfile>_rdip<r_k>_state<s>_wiDiv`/`_noDiv`, or `<outfile>_rdip<r_k>_wiDiv`/`_noDiv` for `nstate` 0; r<sub>k</sub> is
printed with `%g`). The main output is the run at `rdip` as before. The automatic convergence uses it for its rdip step.

Order 2 only, of `method` cCCE (`cluster_mode` standard or pcce, with both
[pcce.connections](InputFile_Cluster.md#tag_pcce_connections) rules; the partition does not depend on `rdip`, only which
partition pairs are connected) or gCCE (every path: ensemble, `evolution` vector and matrix), with `quantity` coherence,
`savemode` coherence, every pair within `rdip` as a cluster (no `max_clusters_per_order` cut) and without
`cluster_contributions`. Anything else is refused with `Error: rdip_ladder ...`. Each rung is its separate run up to the
rounding of the order of the products (`test/code_verification/check_rdip_ladder`: equal to within 1e-9; on a curve that
diverges far above 1, a difference in the last bit of the double, as between two rank counts of one run).

<b>Data Type:</b> Array of floats, positive, strictly increasing, each between the inner and the outer `rdip` (1 to 64 values)

<b>Unit:</b> angstrom

<b>Default:</b> - (off)

<b>Example:</b> ``` "rdip_ladder" : [200, 300, 450] ```

<b>Related variable:</b> [Config](#Config)::rdip_ladder

<a id="tag_nstate"></a>
## nstate
<b>Description:</b><br>
Number of state templates to be generated.

You can specify this value via the "-N" command-line option. <br>
<b>Data Type:</b> Integer

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `0`    | No state template will be generated (Ensemble approach) |
| `>0`   | State template will be generated (Single-sample approach) <br> The static field of the other bath spins is added (see [bath_bath_mean_field](#tag_bath_bath_mean_field)) |

<b>Default:</b> 0

<b>Example:</b> ``` "nstate" : 10 ```

<b>Related Tags:</b> [seed](#tag_seed), [bath_state_file](#tag_bath_state_file)

<b>Related variable:</b> [Config](#Config)::nstate

<a id="tag_seed"></a>
## seed
<b>Description:</b><br>
Seed for random number generator. <br>
If you want to generate the same state templates, you should set the same seed value.

<b>Data Type:</b> Integer

<b>Default:</b> generated on rank 0 and reported

<b>Example:</b> ``` "seed" : 1241798 ```

<b>Related Tags:</b> [nstate](#tag_nstate)

<b>Related variable:</b> [Config](#Config)::seed

<a id="tag_gyrofile"></a>
## gyrofile
<b>Description:</b><br>
File name including bath spin information, e.g. gyromagnetic ratio, spin number <br>

<pre>
    File format is the following :
    [Name] [SpinNumber] [Gyromagnetic ratio (rad/ms/G)]
</pre>
<pre>
    example)
    13C 0.5 6.728284
    P1  0.5 -17608.597050
</pre>

Optional. A bath species the gyrofile does not name (or every species, without a gyrofile) is taken
from the built-in isotope table: every nuclear isotope (mass number > 0) with a positive spin of the EasySpin isotope
database (magnetic moments from N. Stone, IAEA INDC(NDS)-0658, 2014), named mass number + symbol
(`1H`, `11B`, `13C`, `14N`, `29Si`, ...), with γ = g<sub>n</sub> μ<sub>N</sub>/ħ (CODATA 2018). The table is
`src/input/isotope_table.cpp`, written by `example/tools/IsotopeTable/make_isotope_table.py`. A gyrofile
entry wins over the table (a value that differs from it by more than 0.1 % is noted in the log);
species that are not nuclear isotopes (`P1`, `Vp`, an electron spin) must be in the gyrofile. A name
found in neither stops the run. The log says where each species came from. <br>

<b>Data Type:</b> String (File path)

<b>Default:</b> - (the built-in isotope table)

<b>Example:</b> ``` "gyrofile" : "filepath" ```

<b>Related variable:</b> [Config](#Config)::gyrofile

<a id="tag_bathfile"></a>
## bathfile
<b>Description:</b><br>
File name of the bath spatial configuration. <br>
You can specify this value via the "-I" command-line option. <br>

<pre>
    File format is the following :
    [x] [y] [z] [Name]
</pre>
<pre>
    example)
    1.0 0.0 0.0 0.0   - Total line number 
    1.0 2.0 3.0 13C   - x,y,z in cartesian, angstrom
    4.0 5.0 6.0  P1   - Spin name should be the same to the name in gyrofile
</pre>

<b>Data Type:</b> String List (file paths). In the input file it is always an array, also for one
file; a bare string is refused. The `-I filepath` command-line option takes a single path.

<b>Default:</b> -

<b>Example:</b> 
    ``` "bathfile" : ["filepath"] ```

``` "bathfile" : ["filepath1", "filepath2", ... ] ```

<b>Related variable:</b> [Config](#Config)::bathfile

<a id="tag_bath_shift"></a>
## bath_shift
<b>Description:</b><br>
Shift each bath file's sites by `[x, y, z]` (`bathadjust` until 2026-10-13). <br>
Usually needed when several bath files describe the surroundings of different qubits. <br>

<b>Data Type:</b> List of Double Lists: one `[x, y, z]` shift per bath file, so exactly as many rows
as `bathfile` has entries, each with three numbers. A flat `[1.0, 2.0, 3.0]` is refused.

<b>Unit:</b> angstrom

<b>Default:</b> `[0, 0, 0]` for every bath file

<b>Example:</b> 
    ``` "bath_shift" : [[1.0, 2.0, 3.0]] ```

``` "bath_shift" : [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], ... [7.0, 8.0, 9.0]] ```

<b>Related variable:</b> [Config](#Config)::bathadjust

<a id="tag_bath_state_file"></a>
## bath_state_file
<b>Description:</b><br>
Prefix of the files of the bath states (`statefile` until 2026-10-13): with [nstate](#tag_nstate) &gt; 0 the state k
is read from `<prefix>k`, one sublevel per bath spin, in the order of the bath file. <br>
Without it (or when the file does not exist) the states are drawn from [seed](#tag_seed). <br>
You can specify this value via the "-s" command-line option. <br>

<pre>
    File format is the following :
    [SpinSubLevel]
</pre>
<pre>
    example)
     5.0   - The number of total line
    -0.5   - Spin sublevel when S=1/2
     0.5   -    :
    -0.5
    -0.5
</pre>

<b>Data Type:</b> String (File path prefix)

<b>Default:</b> -

<b>Example:</b> ``` "bath_state_file" : "./state" ``` (reads `./state1`, `./state2`, ...)

<b>Related tags:</b> [nstate](#tag_nstate), [Defect.state_file](InputFile_Defect.md#tag_df_state_file)

<b>Related variable:</b> [Config](#Config)::statefile

<a id="tag_coordinate_frame_rotation"></a>
## coordinate_frame_rotation
<b>Description:</b><br>
Re-express the qubits and the bath in the frame whose z axis is the reference qubit's
symmetry axis. <br>
Use it when the bath file is written in the crystal frame (e.g. diamond [001]) but the
qubit is quantized along another direction (e.g. an NV along [111]). <br>
More than one qubit is supported: every Qubit.xyz and every bath-spin position receives
the same R and the same r<sub>0</sub>. <br>
The keys of the block are "enabled", "bath_axis", "qubit_axis", "reference_qubit",
"qubit_position_frame", "hf_tensor_frame" and "qd_tensor_frame" (the last two: the basis of the components of
[hyperfine.file](#tag_hyperfine_file) and [quadrupole.file](#tag_quadrupole_file)). What is transformed,
what is refused and what each key means are in [The coordinate frame rotation](../methods/coordinate-frame-rotation.md). <br>

<b>Data Type:</b> Object { Boolean, Double List(3), Double List(3), String x4 }

<b>Default:</b> ``` {"enabled": false, "bath_axis": [0,0,1], "qubit_axis": [0,0,1]} ```
An absent section, or "enabled" : false, leaves the run numerically identical to a build
without this option.

<b>Example:</b>
    ```
    "coordinate_frame_rotation" : {
        "enabled"         : true,
        "bath_axis"       : [0.0, 0.0, 1.0],
        "qubit_axis"      : [1.0, 1.0, 1.0],
        "reference_qubit" : "q0",
        "qubit_position_frame" : "bath",
        "hf_tensor_frame" : "bath",
        "qd_tensor_frame" : "bath"
    }
    ```

<b>Related variable:</b> [Config](#Config)::rot_enabled, [Config](#Config)::rot_bath_axis, [Config](#Config)::rot_qubit_axis, [Config](#Config)::rot_reference_qubit

<a id="tag_defect_axis_reference"></a>
## defect_axis_reference
<b>Description:</b><br>
Name the axis a future Jahn-Teller classification should measure defect axes against. <br>
<b>This key changes no number that this version computes.</b> It is read, validated and
written to the run log, and nothing else: it neither performs JT assignment nor changes
the Defect Hamiltonian. Setting it does not move the coherence, the axis assignments or
the detuning table. This is separate from each Defect species' own
[coordinate_frame](InputFile_Defect.md#tag_df_coordinate_frame): bath-frame Defect vectors and tensors
are still converted by coordinate_frame_rotation regardless of this reference key. <br>

The only accepted value is the string "qubit_axis", meaning the central qubit's symmetry
axis. That axis is defined only by
[coordinate_frame_rotation](#tag_coordinate_frame_rotation), so the key requires that
section to be present and enabled; otherwise it is an error rather than a no-op. <br>
Anything else is refused: another string, an array, a boolean. <br>

With the rotation on, the run log reports what was recorded and says plainly that it is
not applied: <br>
<pre>
      Defect axis reference           : qubit_axis
      Reference axis, source frame    : [ normalize(qubit_axis) ]
      Reference axis, calc frame      : [   0, 0, 1 ]
      Status                          : recorded for future JT-axis classification;
                                        not applied to the Defect Hamiltonian
</pre>

<b>Data Type:</b> String

<b>Default:</b> absent -- no reference axis is recorded.

<b>Example:</b>
    ``` "defect_axis_reference" : "qubit_axis" ```

<b>Related variable:</b> [Config](#Config)::defect_axis_reference

<a id="tag_hamiltonian"></a>
## hamiltonian
<b>Description:</b><br>
The terms of the Hamiltonian that are chosen, not given by the files: `common` holds what every method reads, `cCCE`,
`gCCE`, `meCCE` and `megCCE` what one method reads (meCCE reads `cCCE` too: it builds the cCCE's Hamiltonian). A setting
in the block of another method is refused (a `cCCE` setting in a gCCE or megCCE run, a `meCCE` or `megCCE` section in a
run of another method, also when `-m` changes the method); `meCCE` and `megCCE` hold
[jump_operators](#tag_jump_operators) only; `gCCE` has no settings yet, nor has `noise`, the section of the reserved
method (an empty section is allowed, a key in it is refused). Every key has a default, so the block may be left out. <br>

<b>Data Type:</b> Object `{"common": {...}, "cCCE": {...}, "gCCE": {}, "meCCE": {...}, "megCCE": {...}}`

<b>Example:</b>
    ```
    "hamiltonian" : {
        "common" : {
            "qubit_bath_mean_field" : false,
            "bath_bath_mean_field"  : "z",
            "hyperfine"  : {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3",
                            "missing_row": "stop"},
            "quadrupole" : {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}
        },
        "cCCE" : {"knight": false, "hfmedi": false, "delta_A": true}
    }
    ```

<a id="tag_qubit_bath_mean_field"></a>
## qubit_bath_mean_field
<b>Description:</b><br>
`hamiltonian.common`. With [nstate](#tag_nstate) &gt; 0, add the Overhauser field of the bath on the qubit (the static
hyperfine field of the sampled bath states; `Qubit.overhaus` until 2026-10-13). It enters the qubit Hamiltonian, which only
gCCE builds: with cCCE it has no effect (the run says so). Read with both forms of the [Qubit](InputFile_Qubit.md#tag_Qubit)
block. <br>
<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"common": {"qubit_bath_mean_field": true}} ```

<b>Related variable:</b> [QubitArray](#QubitArray)::overhaus

<a id="tag_bath_bath_mean_field"></a>
## bath_bath_mean_field
<b>Description:</b><br>
`hamiltonian.common`. With [nstate](#tag_nstate) &gt; 0 each cluster feels the static field of every bath spin outside it,
from that spin's sampled state m<sub>j</sub>; this key says how (with `nstate` 0 it is not used). <br>
`"z"` keeps only the z part: J<sub>zz</sub> m<sub>j</sub> I<sub>z</sub> (the behaviour before the key). <br>
`"xyz"` uses the whole vector J·(0, 0, m<sub>j</sub>) = J[:, z] m<sub>j</sub>, so the x and y parts enter too (the key
`vector_mean_field: true` until 2026-10-13); a cluster removes its own members' field (all three parts) and keeps their full
coupling tensor instead. The x and y parts matter at first order for a bath spin whose own z field nearly cancels (a close
neighbour's J<sub>zz</sub> m<sub>j</sub> against its Zeeman term, as in clustered P1 baths) or whose axis the qubit tilts;
there the z-only field can make the cluster expansion diverge ((internal notes), 강건성 6–13). The
transverse field is added and removed only between main bath spins: every mean-field term involving a defect sub-spin
(e.g. a P1's <sup>14</sup>N) stays z only, while the couplings inside a cluster keep their full tensors. Against exact
10-spin subsystems it removed the divergence of a broken state but was not more accurate in every normal region, so it is
not the default (same document, 12–13). <br>
`"off"` (new, 2026-10-13) leaves the field of the spins outside the cluster out altogether: a cluster is then the cluster
alone in its sampled state. `test/code_verification/check_analytic` checks that two spins at order 1 then give the product of
each spin alone in its bath, which `"z"` does not. <br>
<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"off"` | no static field of the other bath spins |
| `"z"` | its z part (default) |
| `"xyz"` | the vector J·m |

<b>Default:</b> `"z"`

<b>Example:</b> ``` "hamiltonian" : {"common": {"bath_bath_mean_field": "xyz"}} ```

<b>Related Tags:</b> [nstate](#tag_nstate)

<b>Related variable:</b> [Config](#Config)::vector_mean_field, [Config](#Config)::bath_mean_field_off

<a id="tag_hyperfine"></a>
## hyperfine
<b>Description:</b><br>
`hamiltonian.common`. Where the hyperfine tensors between the qubit and the bath spins come from (`hf_readmode` 0–3 until
2026-10-13). `dipole` is the anisotropic part: `"point"` computes the point-dipole tensor from the positions, `"dft"`
takes the DFT dipolar tensor of [file](#tag_hyperfine_file). `contact` adds the isotropic Fermi-contact term of `file`.
`file` is required when `dipole` is `"dft"` or `contact` is `true`, and [missing_row](#tag_hyperfine_missing_row) is read
then too. <br>

<b>Data Type:</b> Object `{"dipole": "point"|"dft", "contact": true|false, "file": ..., "missing_row": "stop"|"point_dipole"}`

<b>Options:</b><br>
| dipole | contact | Description | hf_readmode before |
|---|---|---|---|
| `"point"` | `false` | point-dipole tensor (the default) | 0 |
| `"point"` | `true` | Fermi-contact term + point-dipole tensor | 1 |
| `"dft"` | `false` | DFT dipolar tensor | 2 |
| `"dft"` | `true` | Fermi-contact term + DFT dipolar tensor | 3 |

<b>Default:</b> `{"dipole": "point", "contact": false}`

<b>Example:</b> ``` "hamiltonian" : {"common": {"hyperfine": {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}}} ```

<b>Related variable:</b> [Config](#Config)::hf_readmode

<a id="tag_hyperfine_file"></a>
## hyperfine.file
<b>Description:</b><br>
The hyperfine tensor file (`hf_tensorfile` until 2026-10-13), read when [hyperfine](#tag_hyperfine) `dipole` is `"dft"` or
`contact` is `true`. It is written by `example/tools/SpinTensorReader` from a DFT run and holds, after a header, one
line per species with its g-factor and its default tensor (`etc____`, used for atoms outside
the range), the range of relative positions the file covers (`MinDif[A]`, `MaxDif[A]`), and one
row per atom: species, position relative to the defect (Å), the isotropic part (MHz) and the
nine tensor components (MHz), each per unit nuclear g-factor (CCEX multiplies by the species' g-factor). <br>
<b>Format v3 (2026-10-13):</b> line 1 is the tag of the DFT code (`VASP` or `QUANTUMESPRESSO`) and line 2 the provenance
line `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`.
The numbers are then the hyperfine of the spin-S defect whichever code made them (SpinTensorReader divides Quantum
ESPRESSO's S = 1/2 coupling by 2S; VASP's OUTCAR is already per S), and CCEX reads them as written. Refused: a first line
that is not the tag of the source, the provenance line elsewhere than line 2 or twice, an unknown format or source, and an
S other than the qubit spin. <br>
A file without the provenance line (written before format v3) is read by the earlier rule, with a warning in the log:
first line `VASP` -> x 1, anything else -> x 0.5/S with S the qubit spin. A SpinTensorReader `--version v1` file (OUTCAR
x 2) is right only untagged and for S = 1; with a `VASP` line it is read 2x (the h-BN tutorial file until 2026-10-13).
Remake old files with SpinTensorReader. (The input keys `DefectTotSpin`, which had to be the qubit spin, and
`CorrTotSpin`, a factor for old-format VASP files, were removed on 2026-10-13: the qubit spin is used, and an old-format
VASP file that needed `CorrTotSpin` has to be remade as format v3.) <br>
A bath spin is matched to a row by its position relative to the qubit. A spin outside the
range gets a point-dipole tensor; a spin inside the range with no row is decided by
[missing_row](#tag_hyperfine_missing_row). The positions are in the bath frame; the frame of the
components is `hf_tensor_frame` of [coordinate_frame_rotation](#tag_coordinate_frame_rotation). One qubit only. <br>
<pre>
VASP                                         the DFT code (format v3)
# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = VASP OUTCAR, S = 1
AtomName   r_1i_X  r_1i_Y  r_1i_Z  iso[MHz]  Axx Axy Axz Ayx Ayy Ayz Azx Azy Azz
10B  g-factor___:  0.600215 [No-unit]       one line per species
10B  etc____:      ---  ---  ---  0.00000    the tensor used outside the range
MinDif[A]:  -11.31  -6.53  0.00              the range of relative positions the file covers
MaxDif[A]:   10.06   6.53  0.00
B   -10.05626  -4.36351  0.00000  0.000000  0.030490 ... -0.022870   one row per atom
</pre>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required when `dipole` is `"dft"` or `contact` is `true`)

<b>Example:</b> ``` "file" : "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3" ``` (inside `hyperfine`)

<b>Related variable:</b> [Config](#Config)::hf_tensorfile

<a id="tag_hyperfine_missing_row"></a>
## hyperfine.missing_row
<b>Description:</b><br>
What to do with a bath spin whose position relative to the qubit lies inside the range the
hyperfine tensor file covers but has no row in the file (`hf_ignore_oor` 0/1 until 2026-10-13). Spins outside the range
always get a point-dipole tensor. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"stop"` | stop the run and print the position that was not found |
| `"point_dipole"` | give the spin a point-dipole tensor and go on |

<b>Default:</b> `"stop"`

<b>Example:</b> ``` "missing_row" : "point_dipole" ``` (inside `hyperfine`)

<b>Related variable:</b> [Config](#Config)::hf_ignore_oor

<a id="tag_quadrupole"></a>
## quadrupole
<b>Description:</b><br>
`hamiltonian.common`. Where the quadrupole (electric field gradient) tensors of the bath spins come from (`qd_readmode`
until 2026-10-13). `"off"` is no quadrupole interaction. An object reads a DFT [file](#tag_quadrupole_file): a bath spin
inside the range the file covers takes the row at its position relative to the defect (a spin inside the range with no row
stops the run); a spin outside the range takes the species default of the file (`etc`) with `structure` `"plain"`, or, with
`"strained"` and `"bubble"`, the tensor of the structure without the defect from
[file_without_defect](#tag_quadrupole_file_without_defect), with the relative position folded by
[cell](#tag_quadrupole_cell) along every axis whose cell parameter is not zero. One qubit only. <br>

<b>Data Type:</b> `"off"`, or Object `{"structure": ..., "file": ..., "file_without_defect": ..., "cell": [a, b, c]}`

<b>Options:</b><br>
| structure | Description | qd_readmode before |
|---|---|---|
| (`"off"`) | no quadrupole interaction (the default) | 0 |
| `"plain"` | one DFT file computed with the defect, `file` | 2 |
| `"strained"` | strained structure: `file` with the defect and `file_without_defect` without it; `cell` required | 3 |
| `"bubble"` | bubble structure: the same two files and `cell` | 4 |

(`qd_readmode` 1, experimental constants, had been removed before.)

<b>Default:</b> `"off"`

<b>Example:</b> ``` "hamiltonian" : {"common": {"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}}} ```

<b>Related variable:</b> [Config](#Config)::qd_readmode

<a id="tag_quadrupole_file"></a>
## quadrupole.file
<b>Description:</b><br>
The quadrupole tensor file computed with the defect (`qd_tensorfile` until 2026-10-13). The same layout as
[hyperfine.file](#tag_hyperfine_file), with the electric field gradient components (Hartree/Bohr²) in place of the
hyperfine ones and one `eQ___` line per species (its quadrupole moment) in place of the g-factor. Written by
`example/tools/SpinTensorReader`. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required with every `structure`)

<b>Example:</b> ``` "file" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ``` (inside `quadrupole`)

<b>Related variable:</b> [Config](#Config)::qd_tensorfile

<a id="tag_quadrupole_file_without_defect"></a>
## quadrupole.file_without_defect
<b>Description:</b><br>
The quadrupole tensor file of the same structure computed without the defect (`qd_tensorfile_woqubit` until 2026-10-13),
in the layout of [quadrupole.file](#tag_quadrupole_file). Read with `structure` `"strained"` or `"bubble"`; it supplies the
tensors of the spins outside the range of the file with the defect. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required with `"strained"` and `"bubble"`)

<b>Example:</b> ``` "file_without_defect" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ``` (inside `quadrupole`)

<b>Related variable:</b> [Config](#Config)::qd_tensorfile_woqubit

<a id="tag_quadrupole_cell"></a>
## quadrupole.cell
<b>Description:</b><br>
Lattice parameters of the strained or bubble structure (`qd_cellpara` until 2026-10-13). Read with `structure`
`"strained"` or `"bubble"`. <br>

<b>Data Type:</b> Double List (3)

<b>Default:</b> - (required with `"strained"` and `"bubble"`)

<b>Unit:</b> angstrom

<b>Example:</b> ``` "cell" : [1.0, 2.0, 3.0] ``` (inside `quadrupole`)

<b>Related variable:</b> [Config](#Config)::qd_cellpara

<a id="tag_hfmedi"></a>
## hfmedi
<b>Description:</b><br>
`hamiltonian.cCCE`. Add the hyperfine-mediated interaction between bath spins, the second-order
coupling through the qubit (N. Zhao et al., Nature Nanotechnology 6, 242 (2011), Supplementary, Eq. 9).
CCEX adds the full second-order correction, for each qubit branch α on its own:
H_med^(α) = Σ_{β≠α} E_αβ E_αβ† / (ω_α − ω_β), E_αβ = Σ_i ⟨α|S|β⟩·A_i·I_i, where α is the eigenstate of the
qubit Hamiltonian closest to the branch state and β runs over every other eigenstate (every intermediate state).
Expanded, it holds both orders of every pair (i, j) and the i = j terms of each bath spin; it is Hermitian.
(Before 2.0.0 it held one pair order and the other qubit level only.)
`example/tutorial/07_hBN_VB_natab/hyperfine_mediated_term/` runs with and without it. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"hfmedi": true}} ```

<b>Related variable:</b> [Config](#Config)::hfmedi

<a id="tag_knight"></a>
## knight
<b>Description:</b><br>
`hamiltonian.cCCE`. Add the Knight-field shift of the bath spins that the qubit induces under a
transverse magnetic field (same reference as [hfmedi](#tag_hfmedi), SI Eq. 14):
H_K^(α) = −γ_n B·(g − 1)·I = −(3|m_s| − 2) (γ_e/Δ) B·Ã·I for each bath spin, with Δ = 3/2 of the zz entry of
`Qubit.zfs`, Ã the hyperfine tensor with its z row set to zero, and |m_s| = ⟨S_z²⟩ of the branch's eigenstate.
It is zero when the field is along z. (Before 2.0.0 the term used the bare field B instead of −γ_n B and |m_s| = 1
for both branches.)
`example/tutorial/02_Diamond_NV_natab/knight_shift/` runs with and without it in a tilted field. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"knight": true}} ```

<b>Related variable:</b> [Config](#Config)::knight

<a id="tag_delta_A"></a>
## delta_A
<b>Description:</b><br>
`hamiltonian.cCCE` (new, 2026-10-13). `false` replaces, in every cluster of more than one bath spin and for each qubit, the
hyperfine tensors of the cluster's spins by their mean, so the Knight fields of the spins no longer differ (the pair term
Z<sub>jk</sub> = A<sub>j,z</sub> − A<sub>k,z</sub> of Zhao 2011, the reference of [hfmedi](#tag_hfmedi), SI eq. 23;
ΔΩ<sub>nm</sub> of Maze 2008, eq. 8). A cluster of one spin keeps its own tensor, so the single-spin (ESEEM) terms
are unchanged and only the flip-flop contribution of pairs and larger clusters goes: a comparison calculation that shows
how much of the decoherence comes from the hyperfine-field difference. `true` (the default) is the ordinary calculation.
It is defined for a bath of nuclear spins: `false` is refused with Defect species and with
[qubit_bath_mean_field](#tag_qubit_bath_mean_field) true (the Overhauser field is computed from the tensors before they are
averaged). `test/code_verification/check_analytic` checks it on two <sup>13</sup>C on the NV axis: with `false` their
Hahn echo is 1, with `true` it is not, and order 1 is the same either way. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"delta_A": false}} ```

<b>Related variable:</b> [Config](#Config)::delta_A

<a id="tag_jump_operators"></a>
## jump_operators
<b>Description:</b><br>
`hamiltonian.meCCE` or `hamiltonian.megCCE` (new, 2026-10-13; a top-level `jump_operators` in the adding_MECCE branch,
refused with its new place). The Lindblad jump operators of the master-equation methods, one entry per operator:
C = √γ S<sub>op</sub> on every spin named `bath_name`, with γ = `rate` × 2π·10³ rad/ms. `operator` is `-` (S<sub>-</sub>),
`+` (S<sub>+</sub>), `z`, `x` or `y` (S<sub>z</sub>, S<sub>x</sub>, S<sub>y</sub>) of that spin, of any S; the branch's
old names `+-` and `-+` are refused with the new ones. `bath_name` is the name of bath spins (as in the bath file or of
the Defect sub-spins, case not significant), or `"qubit"` / the qubit's name for a jump of the qubit, which megCCE alone
takes (one qubit; meCCE refuses it and names megCCE). A `bath_name` that is neither, checked against the whole bath once
it is read, is refused with the names the bath has (a cluster without that species simply gets no jump from it), and so
is a `bath_name` that is both the qubit's name and a bath spin's (megCCE would apply it to the qubit only). The
section must belong to the run's method (after `-m`), and holds `jump_operators` only. Without the key, or with an empty
list, the run is the rate-0 limit: the cCCE (meCCE) or the gCCE (megCCE). The equations, what is supported and the
checks are in [master-equation.md](../methods/master-equation.md). <br>

<b>Data Type:</b> Array of objects `{"bath_name": string, "operator": string, "rate": number}`

<b>Unit:</b> `rate` in MHz (≥ 0)

<b>Default:</b> none (no dissipation)

<b>Example:</b>
    ```
    "method" : "megCCE",
    "hamiltonian" : {"megCCE": {"jump_operators": [
        {"bath_name": "e",     "operator": "-", "rate": 1.5915e-4},
        {"bath_name": "e",     "operator": "+", "rate": 1.5915e-4},
        {"bath_name": "qubit", "operator": "z", "rate": 1.5915e-4}]}}
    ```

<b>Related variable:</b> [Config](#Config)::jump_operators (rates stored in rad/ms), [Config](#Config)::jump_operators_section
