# General tags
Tags in this section are about the simulation parameters: the method, the numbers, the files, the frame, and the
`hamiltonian` block. Old key names and their current places are in the table of
[old names](InputFile.md#OldNames).

- [method](#tag_method)       
- [cluster.mode](InputFile_Cluster.md#tag_cluster_mode)
- [quantity](#tag_quantity)     
- [order](#tag_order)        
- [bfield](#tag_bfield)       
- [rbath](#tag_rbath)        
- [rdip](#tag_rdip)         
- [deltat](#tag_deltat)       
- [nstep](#tag_nstep)        
- [configuration](#tag_configuration)
    - [path](#tag_configuration_path), [k_range](#tag_configuration_k_range) and [{k}](#tag_configuration_k)
    - [gyrofile](#tag_gyrofile), [bathfile](#tag_bathfile), [bath_shift](#tag_bath_shift),
      [bath_state_file](#tag_bath_state_file)
    - [defect_axis_file](#tag_defect_axis_file), [defect_state_file](#tag_defect_state_file)
    - [s_range](#tag_s_range), [seed](#tag_seed), [files that do not open](#tag_configuration_missing)
- [frame](#tag_frame): [rotation](#tag_frame_rotation)
- [hamiltonian](#tag_hamiltonian)
    - common: [qubit_bath_mean_field](#tag_qubit_bath_mean_field), [bath_bath_mean_field](#tag_bath_bath_mean_field),
      [hyperfine](#tag_hyperfine) (`dipole`, `contact`, [file](#tag_hyperfine_file), [missing_row](#tag_hyperfine_missing_row),
      [frame](#tag_hyperfine_frame)),
      [quadrupole](#tag_quadrupole) (`structure`, [file](#tag_quadrupole_file),
      [file_without_defect](#tag_quadrupole_file_without_defect), [cell](#tag_quadrupole_cell), [frame](#tag_quadrupole_frame))
    - cCCE: [knight](#tag_knight), [hfmedi](#tag_hfmedi), [delta_A](#tag_delta_A) (meCCE reads them too)
    - meCCE, megCCE: [jump_operators](#tag_jump_operators)
    - gCCE: no settings yet
    - noise: [spectrum](#tag_noise_spectrum), [pseudo_secular](#tag_noise_pseudo_secular), [compute](#tag_noise_compute),
      [fft](#tag_noise_fft)
    - pseudospin: [channels](#tag_pseudospin_channels)
- [keys read and not used](#tag_not_used), [turning a part off](#tag_turn_off)

Other files and species are on other pages:
- the qubit's position file (`Qubit.qubit[].file`) is on the [Qubit page](InputFile_Qubit.md#tag_Qubit_file);
- the files of the defect bath are configuration keys ([defect_axis_file](#tag_defect_axis_file),
  [defect_state_file](#tag_defect_state_file)), one per bath file like `configuration.bath_state_file`;
- the defect species are on the [Defect page](InputFile_Defect.md).


<a id="tag_method"></a>
## method
<b>Description:</b><br>
The CCE method.

- `cCCE` (conventional): the qubit is projected out into its two branches.
- `gCCE` (generalized): the qubit is in the cluster Hamiltonian.
- `meCCE` and `megCCE`: cCCE and gCCE with Lindblad dissipation by the jump operators of
  [hamiltonian.\<method\>.jump_operators](#tag_jump_operators) ([master-equation.md](../methods/master-equation.md)).

How the bath is cut into clusters is set by [cluster.mode](InputFile_Cluster.md#tag_cluster_mode)
(`standard`, `pcce`, `bundle`, `dsj`, `itb`, `dsjitb`). Some modes run with some methods only:
- `pcce` runs with cCCE, gCCE, meCCE and megCCE;
- `bundle` runs with cCCE and gCCE;
- `dsj`, `itb` and `dsjitb` run with cCCE only. gCCE, meCCE and megCCE with them are refused as not implemented.
<br>
<b>Data Type:</b> String (case is not significant)

<b>Options:</b><br>
| Option  | Description |
|---------|-------------------------------------|
| `cCCE`  | Conventional CCE |
| `gCCE`  | Generalized CCE |
| `meCCE` | ME-CCE: cCCE with Lindblad jump operators on the bath spins (PRL 132, 250401 (2024)) |
| `megCCE` | ME-gCCE: gCCE with Lindblad jump operators on the bath spins and on the qubit (C<sub>q</sub> ⊗ 1, in every cluster and the 0th) |
| `noise` | the Gaussian (second-cumulant) coherence exp(−χ) of the qubit's phase noise, from a model spectrum or from the bath ([spectrum](#tag_noise_spectrum)) |
| `pseudospin` | the pair clusters of cCCE in the pseudospin model of nuclear flip-flops, alone (see below) |

Limits of each method:
- `meCCE`: CPU, one qubit, ideal π pulses from `pulse.count`. Qubit jumps are refused (use megCCE).
- `megCCE`: CPU. The qubit sequence is used as given (Manual and any rotation included), with `pulse.per_qubit` false
  and `pulse.detuning_factor` 1; other values are refused. DEER bath pulses (`pulse.bath`) run with it.
- `noise`: the semi-classical CCE of `semi-cce`. cCCE clusters (`cluster.mode` standard only), one qubit, ideal π
  pulses, quantity `coherence`. With `execution.backend` `cuda` (`bin/main_gpu.out`) the sums over the bath's modes are
  computed on the card and the eigensystems stay on the CPU. `auto` and a model spectrum run on the CPU.

`pseudospin` is for reading a decay rather than for the best curve. It follows Seo et al., Nat. Commun. 7, 12935 (2016),
Supplementary eqs. 20-24 and Fig. 8, and Lee, Park, Seo, npj 2D Mater. Appl. 6, 60 (2022). The treatment of several
channels is CCEX's own.
- The output is Π<sub>pairs</sub> P, the pairs' terms ([the pseudospin model](#tag_pseudospin_model)), without the
  single spins. The single spins' ESEEM is what cCCE order 1 gives.
- Each pair's Hamiltonian in each qubit state is written in the product basis of the two spins' own eigenstates in
  that state (Zeeman, hyperfine, the whole quadrupole). It is cut there to its diagonal and to every nonzero element
  between product states that differ in both spins. These elements are the channels, for any nuclear spin:
  - the flip-flops \|m<sub>1</sub>, m<sub>2</sub>⟩ ↔ \|m<sub>1</sub>−1, m<sub>2</sub>+1⟩;
  - the flip-flips \|m<sub>1</sub>, m<sub>2</sub>⟩ ↔ \|m<sub>1</sub>+1, m<sub>2</sub>+1⟩;
  - with a quadrupole that mixes the \|m⟩, the double-quantum ones such as
    \|m<sub>1</sub>, m<sub>2</sub>⟩ ↔ \|m<sub>1</sub>−2, m<sub>2</sub>+2⟩.
- Every channel's Ω<sub>α</sub>, Ω<sub>β</sub> and level gaps Δ<sub>α</sub>, Δ<sub>β</sub> (kHz) are written to
  `{output.file}_pseudospin_channels`.
- Limits: one qubit, ensemble, order 2 (the default here; another order is refused),
  [channels](#tag_pseudospin_channels) grouped or pairwise, `cluster.mode` standard, the field along the qubit axis,
  ideal π pulses from `pulse.count`, CPU.

Old values are refused with their current form ([old names](InputFile.md#OldNames)).

<b>Default:</b> `cCCE`

<b>Example:</b> ``` "method" : "gCCE" ```, ``` "method" : "cCCE", "cluster" : {"mode": "pcce"} ```

<a id="tag_quantity"></a>
## quantity
<b>Description:</b><br>
What the run computes.
- `coherence`: the decay curve L(t) between |α⟩ and |β⟩. This is the usual choice.
- `dm` (gCCE): the whole qubit density matrix, for populations or several qubit levels
  (`example/tutorial/advanced/01_P1/gCCE2_Ramsey_population/`).

meCCE and megCCE compute `coherence` only.

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `coherence` | Coherence function             |
| `dm`        | Density matrix                 |

Any other value is refused, `noise` included: the noise calculation is [method](#tag_method) `noise`.

With `dm` (gCCE), the reduced density matrix of the qubits is written as blocks of rows, one block per time step. The
cluster expansion is applied element by element. An element that the free evolution keeps at zero (a level the initial
state does not touch) is left at zero. When clusters would put a value into such an element, the run prints a `[dm]`
line with the largest such value, since the elementwise product cannot carry it. <br>
<b>Unit:</b> arbitrary unit (a.u.)

<b>Default:</b> `coherence`

<b>Example:</b> ``` "quantity" : "coherence" ```

<a id="tag_order"></a>
## order
<b>Description:</b><br>
Order of the CCE method: the size of the largest cluster, counted in bath spins (`cluster.mode` standard), in
partitions (pcce) or in bundles (bundle).

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "order" : 2 ```

<a id="tag_bfield"></a>
## bfield
<b>Description:</b><br>
Magnetic field, read in the computational frame. A number here is Bz.

The key must be in the input file; for no field, write `"bfield": 0`. <br>
With [frame.rotation](#tag_frame_rotation) on, the z axis of the computational frame is `qubit_axis`. Then
`[0, 0, B]` (or the number B) is a field of B along `qubit_axis` in the lab (the frame of the bath file), and a
transverse part is refused. <br>
<b>Data Type:</b> Float or Float List, or `{"value": ..., "unit": ...}` with the unit `G`, `mT` or `T` ([Units](InputFile.md#Units))

<b>Unit:</b> gauss (G)

<b>Default:</b> none (required)

<b>Example:</b>

``` "bfield" : 0.0 ```

``` "bfield" : [0.0, 0.0, 0.0] ```

``` "bfield" : {"value": 50, "unit": "mT"} ``` (= 500)

<a id="tag_rbath"></a>
## rbath
<b>Description:</b><br>
Bath radius: the bath spins within this distance of the qubit are kept. A number is the outer radius; `[inner, outer]`
also removes the spins closer than `inner`. <br>

<b>Data Type:</b> Float, or a list of two floats `[inner, outer]`, or `{"value": ..., "unit": ...}` with the unit `A` (angstrom) or `nm` ([Units](InputFile.md#Units))

<b>Unit:</b> angstrom

<b>Default:</b> - (required); inner 0

<b>Example:</b> ``` "rbath" : 50.0 ```, ``` "rbath" : [2.0, 50.0] ```

<a id="tag_rdip"></a>
## rdip
<b>Description:</b><br>
Dipolar radius: two bath spins closer than this are connected and can be in one cluster. A number is the outer
radius; `[inner, outer]` also leaves out the pairs closer than `inner`.

With `cluster.mode` bundle the key may be left out. The connections are then chosen by a coherence check (see
[cluster.mode](InputFile_Cluster.md#tag_cluster_mode)), and an inner radius needs the outer one. <br>

<b>Data Type:</b> Float, or a list of two floats `[inner, outer]`, or `{"value": ..., "unit": ...}` with the unit `A` (angstrom) or `nm` ([Units](InputFile.md#Units))

<b>Unit:</b> angstrom

<b>Default:</b> - (required, but with `cluster.mode` bundle); inner 0

<b>Example:</b> ``` "rdip" : 10.0 ```, ``` "rdip" : [1.0, 10.0] ```, ``` "rdip" : {"value": 1.0, "unit": "nm"} ```

<a id="tag_deltat"></a>
## deltat
<b>Description:</b><br>
Time interval: the step of the time axis. The outputs hold the times 0, `deltat`, ..., (`nstep` − 1) × `deltat`.
Each time is the total time of the whole pulse sequence: a Hahn echo of time t has its π pulse at t/2.

<b>Data Type:</b> Float, or `{"value": ..., "unit": ...}` with the unit `s`, `ms`, `us` or `ns` ([Units](InputFile.md#Units))

<b>Unit:</b> ms

<b>Default:</b> - (required)

<b>Example:</b> ``` "deltat" : 0.001 ```, ``` "deltat" : {"value": 1, "unit": "us"} ``` (= 0.001, bit for bit)

<a id="tag_nstep"></a>
## nstep
<b>Description:</b><br>
Number of time steps: the number of time points of the outputs, from t = 0 (see [deltat](#tag_deltat)).

<b>Data Type:</b> Integer

<b>Default:</b> - (required)

<b>Example:</b> ``` "nstep" : 10 ```

<a id="tag_rdip_ladder"></a>
## rdip_ladder (internal)
An internal key that the [automatic convergence](InputFile_AutoConvergence.md) writes into the inputs of its runs. It
is not listed among the input keys, and a hand-written input does not need it.

It gives extra rdip values (rungs) below `rdip`, computed from the same order-2 run. The rdip step of the automatic
convergence then reads several rdip values from one run (`<output.file>_rdip<r_k>_...`).

CCEX reads it at the top level and checks that the run has:
- order 2;
- cCCE with `cluster.mode` standard or pcce, or gCCE with `cluster.mode` standard;
- `quantity` and `output.savemode` coherence;
- no `cluster.selection.max_per_order` cut and no `output.cluster_contributions`.

`test/code_verification/check_rdip_ladder` holds its rungs to separate runs.

<a id="tag_configuration"></a>
## configuration
<b>Description:</b><br>
The bath configuration in one block: the files of the bath, its states and the random seed. <br>
```
"configuration": {
    "path": "data",                  # optional: the folder of the file names below
    "k_range": [1, 10],              # optional: one configuration per k, run one after another
    "gyrofile": "gyro_13C",
    "bathfile": ["config_{k}"],
    "bath_shift": [[0, 0, 0]],
    "bath_state_file": ["cfg{k}_state{s}"],
    "defect_axis_file": ["cfg{k}_axis"],      # with Defect species in the bath
    "defect_state_file": ["cfg{k}_onsite{s}"],
    "s_range": [1, 5],                        # optional: sampled states; absent: the ensemble
    "seed": 1
}
```

The keys of the block:
- [path](#tag_configuration_path), [k_range](#tag_configuration_k_range) (and [{k}](#tag_configuration_k))
- [gyrofile](#tag_gyrofile), [bathfile](#tag_bathfile), [bath_shift](#tag_bath_shift),
  [bath_state_file](#tag_bath_state_file)
- [defect_axis_file](#tag_defect_axis_file), [defect_state_file](#tag_defect_state_file)
- [s_range](#tag_s_range), [seed](#tag_seed)

Every file the block names and the run reads has to open ([files that do not open](#tag_configuration_missing)).

In their old places these keys are refused, with a message naming their place in the block
([old names](InputFile.md#OldNames)).
`python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites an old input (comments kept).

<b>Default:</b> none: the block names the bath, and [bathfile](#tag_bathfile) has no default (the other keys of the
block have their own).

<b>Related tags:</b> [Defect.species](InputFile_Defect.md#tag_df_species) (the defect model)

<a id="tag_configuration_path"></a>
## configuration.path
<b>Description:</b><br>
The folder of the block's file names. It is put in front of `gyrofile` and of every entry of `bathfile`,
`bath_state_file`, `defect_axis_file` and `defect_state_file`.
- It is relative to the working directory, or absolute.
- A file name that is itself absolute is used as it is.
- It does not apply to files of other blocks (`Qubit.qubit[].file`, the tensor files of `hamiltonian`). Those stay
  relative to the working directory.

<b>Data Type:</b> String (folder)

<b>Default:</b> `"."` (the working directory)

<b>Example:</b> ``` "path" : "../bath/14N_1ppm" ``` with ``` "bathfile" : ["bath_DiaP1_1ppm_1"] ``` reads
`../bath/14N_1ppm/bath_DiaP1_1ppm_1`

<a id="tag_configuration_k_range"></a>
## configuration.k_range
<b>Description:</b><br>
Several bath configurations in one run. Every k from `first` to `last` (both included) is a configuration of its own.
The configurations are calculated one after another in the same CCEX execution, with [{k}](#tag_configuration_k) in
the file names replaced by k.

Each configuration writes its own outputs, `<output.file>_k<k>_...` (e.g. `out_k3_wiDiv`, `out_k3_state1_noDiv`). No
average over configurations is written. The log has a "Configuration k = ..." header before each one. <br>

Seeds:
- Configuration k uses the seed + (k - first). A given [seed](#tag_seed) number is the first configuration's seed; a
  generated seed is generated once.
- A [seed](#tag_seed) list gives each configuration its own seed.
- Each configuration draws its own states. A run of its files alone with its seed S gives it byte for byte. Without
  `seed`, the log prints each seed: "[seed] configuration k=...: seed S".
- Before 2.0.0 a given seed number was every configuration's seed. The converter writes such an input as a list of
  the same number.

Refused, before anything is calculated:
- `k_range` with no `{k}` in any of the names it applies to;
- a file of any k that cannot be read (every bath file and axis file with `{k}`, and the files of the states of
  [s_range](#tag_s_range) of every state file name with `{k}`);
- a `k_range` that is not two nonnegative integers with first &le; last;
- an enabled pCCE partition cache ([pcce.cache](InputFile_Cluster.md) `mode` read or readwrite) whose `file` has no
  `{k}` (one cache file per configuration, e.g. `"partition_{k}.cache"`);
- `k_range` together with a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block (refine one configuration at a
  time).

<b>Data Type:</b> `[first, last]`, two nonnegative integers

<b>Default:</b> - (one configuration; `{k}` is not replaced)

<b>Example:</b> ``` "k_range" : [1, 10] ``` with ``` "bathfile" : ["config_{k}"] ``` runs `config_1` ... `config_10`

<a id="tag_configuration_k"></a>
## {k} in file names
<b>Description:</b><br>
With [k_range](#tag_configuration_k_range), `{k}` is replaced by the configuration's k in the entries of `bathfile`,
`bath_state_file`, `defect_axis_file`, `defect_state_file` and in `pcce.cache.file`.
- An entry without `{k}` is the same for every k. A bath file without `{k}` is merged into every configuration:
  `["nv1_{k}", "common"]` is, for each k, the bath of `nv1_<k>` and `common` together.
- `{{k}}` gives a literal `{k}`.
- Without `k_range` nothing is replaced: a file actually named `config_{k}` is read as it is.

`{s}` in the entries of [bath_state_file](#tag_bath_state_file) and [defect_state_file](#tag_defect_state_file) is the
state s of [s_range](#tag_s_range), by the same rule (`{{s}}` gives a literal `{s}`). A name without `{s}` has the
state appended: `"state_"` reads `state_1`, `state_2`, .... <br>
<b>Example:</b> ``` "bathfile" : ["config_{k}"], "bath_state_file" : ["cfg{k}_state{s}"] ``` with `k_range` [1, 2] and
`s_range` [1, 3] reads `cfg1_state1` ... `cfg2_state3`

<a id="tag_gyrofile"></a>
## configuration.gyrofile
<b>Description:</b><br>
The file of the bath spin species: name, spin number and gyromagnetic ratio. <br>

<pre>
    File format is the following :
    [Name] [SpinNumber] [Gyromagnetic ratio (rad/ms/G)]
</pre>
<pre>
    example)
    13C 0.5 6.728284
    P1  0.5 -17608.597050
</pre>

The key is optional. A bath species that the gyrofile does not name (or every species, without a gyrofile) is taken
from the built-in isotope table.
- The table holds every nuclear isotope (mass number > 0) with a positive spin of the EasySpin isotope database
  (magnetic moments from N. Stone, IAEA INDC(NDS)-0658, 2014), with γ = g<sub>n</sub> μ<sub>N</sub>/ħ (CODATA 2018).
- Its names are mass number + symbol (`1H`, `11B`, `13C`, `14N`, `29Si`, ...).
- The table is `src/input/isotope_table.cpp`, written by `example/tools/IsotopeTable/make_isotope_table.py`.

Which value is used:
- A gyrofile entry wins over the table. A value that differs from the table by more than 0.1 % is noted in the log.
- Species that are not nuclear isotopes (`P1`, `Vp`, an electron spin) must be in the gyrofile.
- A name found in neither stops the run.
- The log says where each species came from. <br>

<b>Data Type:</b> String (File path; in [configuration.path](#tag_configuration_path))

<b>Default:</b> - (the built-in isotope table)

<b>Example:</b> ``` "configuration" : {"gyrofile" : "filepath", ...} ```

<a id="tag_bathfile"></a>
## configuration.bathfile
<b>Description:</b><br>
Files of the bath spatial configuration. Several files are merged into one bath, in their order. The log names the
merge ("Bath = fileA (n spins) + fileB (m spins)").

Two spins of different files close to each other (after [bath_shift](#tag_bath_shift)):
- Within 0.5 &Aring;: the log warns, naming the files and positions, and the calculation goes on. These are two spins
  at one site: species files generated independently on one lattice that share a site, or two configurations merged
  by mistake (give each its own run, [k_range](#tag_configuration_k_range)).
- At the same position (closer than 10<sup>-6</sup> &Aring;): refused before anything is calculated. Their coupling is
  infinite and every result would be nan. <br>

<pre>
    File format is the following :
    [x] [y] [z] [Name]
</pre>
<pre>
    example)
    2.0 0.0 0.0 0.0   - the number of spins (the lines below)
    1.0 2.0 3.0 13C   - x,y,z in cartesian, angstrom
    4.0 5.0 6.0  P1   - Spin name should be the same to the name in gyrofile
</pre>

<b>Data Type:</b> String List (file paths, in [configuration.path](#tag_configuration_path); `{k}` with
[k_range](#tag_configuration_k_range)). In the input file it is always an array, also for one file; a bare string is
refused.

<b>Default:</b> -

<b>Example:</b> 
    ``` "bathfile" : ["filepath"] ```

``` "bathfile" : ["filepath1", "filepath2", ... ] ```

<a id="tag_bath_shift"></a>
## configuration.bath_shift
<b>Description:</b><br>
Shifts each bath file's sites by `[x, y, z]`. It is usually needed when several bath files describe the surroundings
of different qubits. <br>

<b>Data Type:</b> List of Double Lists: one `[x, y, z]` shift per bath file, so exactly as many rows as `bathfile` has
entries, each with three numbers. A flat `[1.0, 2.0, 3.0]` is refused. Or `{"value": [[...], ...], "unit": ...}` with
the unit `A` (angstrom) or `nm` ([Units](InputFile.md#Units)).

<b>Unit:</b> angstrom

<b>Default:</b> `[0, 0, 0]` for every bath file

<b>Example:</b> 
    ``` "bath_shift" : [[1.0, 2.0, 3.0]] ```

``` "bath_shift" : [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], ... [7.0, 8.0, 9.0]] ```

<a id="tag_bath_state_file"></a>
## configuration.bath_state_file
<b>Description:</b><br>
The files of the bath states, one name per bath file, in the order of `bathfile`. They are read with
[s_range](#tag_s_range) only. Without `s_range` (the ensemble) the key is not used (the run says so).
- The state s of the spins of bath file j is read from name j with `{s}` replaced by s (or s appended,
  [{s}](#tag_configuration_k)).
- Each file has one sublevel per line of bath file j, in its order. The files are read one after another in the
  merged bath's spin order.
- Without the key, the states are drawn from [seed](#tag_seed).
- An entry `""` is a bath file without a state file (e.g. a bath file beside the one the states were made for). The
  whole bath is drawn from the seed first, and the files given then set the states of their own bath files. `""` for
  every bath file is the same as leaving the key out.
- A file that does not open stops the run, naming the key and the file
  ([files that do not open](#tag_configuration_missing)).

The same `""` rule holds for [defect_axis_file](#tag_defect_axis_file) and
[defect_state_file](#tag_defect_state_file). <br>

<a id="tag_configuration_files"></a>
<b>File formats</b> (the three keys): the first line is not read (by convention it is the bath file's own first line,
its line count). Then comes one line per line of the bath file after its first, in its order: line i belongs to the
spin of line i of the bath file.

<pre>
    bath_state_file      m of the spin (one of S, S-1, ..., -S)
    defect_axis_file     the axis state, 1 .. n_axis_states (every line has to be a spin of a Defect species with axis
                         states: CCEX refuses a line of another species)
    defect_state_file    the m of each on-site spin of a Defect spin, space-separated, in the order of the species'
                         types (any number for a line of another species)
</pre>
<pre>
    example) bath_state_file of a 13C bath (S = 1/2)
     5.0   - the first line (not read)
    -0.5   - m of the spin of the bath file's second line
     0.5   -    :
    -0.5
    -0.5
</pre>

`example/tools/StateFiles/make_state_files.py INPUT` writes the files of all three keys that an input names, for every
k of `k_range` and every state. It draws them by a fixed rule ([Tools](Tools.md#tag_tool_statefiles)); the automatic
convergence makes its files with the same rule.

<b>Data Type:</b> String List (file names, one per bath file, in [configuration.path](#tag_configuration_path);
`{k}` with [k_range](#tag_configuration_k_range), `{s}` the state). A list of another length, or a bare string, is
refused.

<b>Default:</b> - (the states are drawn from `seed`)

<b>Example:</b> ``` "bath_state_file" : ["./state{s}"] ``` (reads `./state1`, `./state2`, ... for `s_range` [1, 2]; `["./state"]`
reads the same files); ``` "bath_state_file" : ["./stateA_{s}", "./stateB_{s}"] ``` with two bath files

<b>Related tags:</b> [s_range](#tag_s_range), [defect_state_file](#tag_defect_state_file)

<a id="tag_defect_axis_file"></a>
## configuration.defect_axis_file
<b>Description:</b><br>
File paths of the axis state (the principal axis) of each bath spin that is a
[Defect species](InputFile_Defect.md#tag_df_species). One file per bath file, in the order of
[bathfile](#tag_bathfile): file j holds the lines of bath file j, and the files are read one after another in the
merged bath's spin order.
- Without the key (or, for its bath file alone, with `""`), the axis states of the defect spins are drawn from `seed`.
- A file that does not open stops the run.
- Without Defect species the key is not used (the run says so).
- The names are in [configuration.path](#tag_configuration_path); `{k}` is the configuration with
  [k_range](#tag_configuration_k_range). <br>
<pre>
    File format ([file formats](#tag_configuration_files)): a first line (not read), then per bath-file line the axis state
    example)
    5.0 - the first line (not read)
    1.0 - the axis state of the spin of the bath file's second line
    2.0 -  :
    3.0
    3.0
</pre>

<b>Data Type:</b> String List (file paths, one per bath file; `["filepath"]` for one bath file). A list of another
length, or a bare string, is refused.

<b>Default:</b> -

<b>Example:</b> ``` "configuration" : {"bathfile": ["bath_1"], "defect_axis_file": ["axis_1"], ...} ```

<a id="tag_defect_state_file"></a>
## configuration.defect_state_file
<b>Description:</b><br>
File names (prefixes, when the state is appended) of the states of the defects' on-site nuclear spins. They work as
[bath_state_file](#tag_bath_state_file) does for the main bath spins: one name per bath file, and with
[s_range](#tag_s_range) the state s of the spins of bath file j is read from name j with `{s}` replaced by s (or s
appended).
- Without the key (or with `""`) the states are drawn from `seed`.
- A file that does not open stops the run.
- In the ensemble, or without Defect species, the key is not used.
- The names are in `configuration.path`; `{k}` is the configuration with `k_range`. <br>

<pre>
    File format ([file formats](#tag_configuration_files)): a first line (not read), then per bath-file line the m of
    each on-site spin of that Defect spin, space-separated (any number for a line of another species)
    example) one on-site 14N (I = 1) per P1
     5.0   - the first line (not read)
    -1.0   - m of the 14N of the P1 on the bath file's second line
     0.0   -    :
     1.0
    -1.0
</pre>

<b>Data Type:</b> String List (file names, one per bath file; `{k}`, `{s}`). A list of another length, or a bare string,
is refused.

<b>Default:</b> -

<b>Example:</b> ``` "configuration" : {"bathfile": ["bath_1"], "s_range": [1, 4], "defect_state_file": ["./onsite{s}"], ...} ```

<a id="tag_s_range"></a>
## configuration.s_range
<b>Description:</b><br>
The sampled bath states, `[first, last]`. An old `nstate` N is `[1, N]`.
- Without it, the run averages over the unpolarised bath (the ensemble) and writes one output.
- With it, each state s from first to last is a pure bath state: each spin is in one sublevel, read from
  [bath_state_file](#tag_bath_state_file) or drawn from [seed](#tag_seed). Each state writes its own `_state<s>`
  output; their average is left to the analysis.

Use states for an electron-spin bath (P1), where one bath state differs much from another, or for the spread between
single centres. <br>

States drawn from the seed are drawn for 1 .. last and calculated for first .. last. State s is therefore the same
state as in a run of `[1, last]`: `[11, 20]` and `[1, 10]` together are `[1, 20]`, split over two runs. With state
files, only the files of first .. last are read (and have to be there). <br>
<b>Data Type:</b> `[first, last]`, two integers, 1 &le; first &le; last (a list also for one state: `[3, 3]`)

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| absent | the ensemble: the average over the unpolarised bath, one output |
| `[first, last]` | the sampled states first .. last, each its own `_state<s>` output <br> The static field of the other bath spins is added (see [bath_bath_mean_field](#tag_bath_bath_mean_field)) |

<b>Default:</b> - (the ensemble)

<b>Example:</b> ``` "configuration" : {"s_range" : [1, 10], ...} ```

<b>Related Tags:</b> [seed](#tag_seed), [bath_state_file](#tag_bath_state_file), [{s}](#tag_configuration_k)

<a id="tag_seed"></a>
## configuration.seed
<b>Description:</b><br>
Seed of the random number generator. To draw the same states again, give the same seed.
- A number is the seed of the first configuration of [k_range](#tag_configuration_k_range). Configuration k uses
  seed + (k - first); a generated seed follows the same rule. Without k_range it is the run's seed.
- A list gives one seed per configuration: as many as `k_range` has, or one without it.
- Before 2.0.0 a number was every configuration's seed. The converter writes such an input as a list of the same
  number.

<b>Data Type:</b> Integer, or Integer List (one per configuration)

<b>Default:</b> generated on rank 0 and reported in the log (`[seed] ...`) when the run draws random numbers
([s_range](#tag_s_range), Defect species in the bath, a pCCE partition); otherwise the seed is not used

<b>Example:</b> ``` "configuration" : {"seed" : 1241798, ...} ```, ``` "configuration" : {"k_range": [1, 3], "seed" : [7, 7, 9], ...} ```

<b>Related Tags:</b> [s_range](#tag_s_range)

<a id="tag_configuration_missing"></a>
## Files that do not open
<b>Description:</b><br>
A file that the input names and the run reads has to open. Otherwise the run stops before anything is calculated,
naming the key and the file. The keys:
- `configuration.bathfile`, `gyrofile`;
- `bath_state_file` and `defect_state_file` (every state of `s_range`);
- `defect_axis_file` (with Defect species);
- `Qubit.qubit[].file`;
- `hamiltonian.common.hyperfine.file`, `hamiltonian.common.quadrupole.file` and `file_without_defect`;
- `cluster.pcce.cache.file` in mode `read`.

`""` in a list of per-bath-file names is no file: that bath file's states or axes are drawn from the seed. A file the
run does not read (a state file in the ensemble) is not opened; the key is reported as not used.

<a id="tag_frame"></a>
## frame
<b>Description:</b><br>
The frame block: its one key is [rotation](#tag_frame_rotation). <br>

<b>Example:</b> ``` "frame" : {"rotation": {"qubit_axis": [1, 1, 1]}} ```

<a id="tag_frame_rotation"></a>
## frame.rotation
<b>Description:</b><br>
Turns the frame of the bath file (the lab or crystal frame) so that `qubit_axis` becomes the computational z axis.
Use it when the bath file is written in a crystal frame (e.g. the cubic axes of diamond) and the qubit is quantized
along another direction (e.g. an NV along [111]). <br>

A given block is on: `enabled` defaults to `true`, and `"enabled": false` turns a written block off. No block is no
rotation. The rotation, the frames and what is turned are explained in
[The coordinate frame rotation](../methods/coordinate-frame-rotation.md). <br>

What is turned:
- Every position (the bath spins and the qubits) is turned about the first qubit. Every tensor written in the file
  frame is turned with them.
- The [bfield](#tag_bfield) is read in the computational frame, so `[0, 0, B]` is a field along `qubit_axis` in the
  lab.
- The qubit states (|α⟩, |β⟩, psi0) and the pulse axes stay in the computational frame. The computational x is
  `qubit_axis` × [0, 0, 1], and this choice sets the preparation phase relative to the system.

With gCCE and a non-axial qubit ZFS (E ≠ 0) the coherence depends on that phase. With the check's system turned by
180°, |L| changed by 4.3e-4 for the free induction and by 3.9e-2 for a Hahn echo. cCCE and an axial ZFS do not depend
on it (gCCE: about 1e-8, from the non-secular hyperfine). <br>

Each set of directional numbers says which frame it is written in, with a key `frame` beside it:
- `"file"`: written in the bath file's frame, before the rotation; turned with the positions.
- `"qubit"`: written in the computational frame (z = `qubit_axis`); used as written.

With the rotation on, `frame` is required wherever such numbers are written:
[hyperfine.frame](#tag_hyperfine_frame) (with `"dipole": "dft"`), [quadrupole.frame](#tag_quadrupole_frame),
[Qubit.intmap[].frame](InputFile_Qubit.md#tag_intmap_frame) (a qubit's ZFS is a self entry there)
and [Defect.species[].frame](InputFile_Defect.md#tag_df_frame) (a species with `nuclear_positions`, `hyperfine`, `efg` or
`zfs`). With the rotation off, a written `frame` is not used (reported, [keys not used](#tag_not_used)). <br>

Old keys of the block are refused ([old names](InputFile.md#OldNames)):
- `bath_axis` (it turned x and y about `qubit_axis`; x is now `qubit_axis` × [0, 0, 1], the old default);
- `reference_qubit` (the origin, now the first qubit);
- `qubit_position_frame` (one value);
- `hf_tensor_frame`, `qd_tensor_frame` (now [hyperfine.frame](#tag_hyperfine_frame),
  [quadrupole.frame](#tag_quadrupole_frame)). <br>

| Key | Type | Values | Default |
|---|---|---|---|
| `enabled` | Boolean | `true`, `false` | `true` (in a given block) |
| `qubit_axis` | Double List(3) | the qubit's quantization axis in the bath file's frame, any nonzero vector except along [0, 0, -1] | [0, 0, 1] |

<b>Data Type:</b> Object { Boolean, Double List(3) }

<b>Default:</b> no rotation (no block). No block, or "enabled" : false, leaves the run numerically identical to a build
without this option.

<b>Example:</b> an NV along [111] in a cubic diamond bath, with a DFT hyperfine file written in the same frame
    ```
    "frame"       : {"rotation": {"qubit_axis": [1.0, 1.0, 1.0]}},
    "bfield"      : 500,
    "hamiltonian" : {"common": {"hyperfine": {"dipole": "dft", "file": "Afile", "frame": "file"}}}
    ```

<a id="tag_hamiltonian"></a>
## hamiltonian
<b>Description:</b><br>
The terms of the Hamiltonian that are chosen, not given by the files. Every key has a default, so the block may be
left out.

| Section | Read by | Holds |
|---|---|---|
| `common` | every method | the mean fields, `hyperfine`, `quadrupole` |
| `cCCE` | cCCE, and meCCE (it builds the cCCE's Hamiltonian) | `knight`, `hfmedi`, `delta_A` |
| `gCCE` | | no settings yet: an empty section is allowed, a key in it is refused |
| `meCCE`, `megCCE` | that method | [jump_operators](#tag_jump_operators) only |
| `noise` | method `noise` | [spectrum](#tag_noise_spectrum) and the other noise keys |

A setting in the section of another method is refused: a `cCCE` setting in a gCCE or megCCE run, a `meCCE` or
`megCCE` section in a run of another method, a `noise` section in a run of another method. <br>

<b>Data Type:</b> Object `{"common": {...}, "cCCE": {...}, "gCCE": {}, "meCCE": {...}, "megCCE": {...}}`

<b>Default:</b> no block: every key at its default (no mean field on the qubit, the z mean field between bath spins,
point-dipole hyperfine, no quadrupole, no knight, hfmedi or delta_A, no jump operators).

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
        "cCCE" : {"knight": false, "hfmedi": false, "delta_A": false}
    }
    ```

<a id="tag_qubit_bath_mean_field"></a>
## qubit_bath_mean_field
<b>Description:</b><br>
`hamiltonian.common`. With [s_range](#tag_s_range) (sampled states), adds the Overhauser field of the bath on the
qubit: the static hyperfine field of the sampled bath states. In each cluster this is the field of every bath spin
outside the cluster, since the cluster's own spins act on the qubit through their coupling. With `false` the qubit
feels the spins of the cluster only.

It enters the qubit Hamiltonian, which gCCE and megCCE build. With cCCE and meCCE, and in the ensemble (no sampled
state, the field is zero), it is not used (reported, [keys not used](#tag_not_used)). <br>
<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"common": {"qubit_bath_mean_field": true}} ```

<a id="tag_bath_bath_mean_field"></a>
## bath_bath_mean_field
<b>Description:</b><br>
`hamiltonian.common`. With [s_range](#tag_s_range), each cluster feels the static field of every bath spin outside it,
from that spin's sampled state m<sub>j</sub>. This key says how. In the ensemble it is not used, and a written one is
reported.
- `"z"` keeps only the z part: J<sub>zz</sub> m<sub>j</sub> I<sub>z</sub>.
- `"xyz"` uses the whole vector J·(0, 0, m<sub>j</sub>) = J[:, z] m<sub>j</sub>, so the x and y parts enter too. A
  cluster removes its own members' field (all three parts) and keeps their full coupling tensor instead.
- `"off"` leaves out the field of the spins outside the cluster. A cluster is then the cluster alone in its sampled
  state.

The x and y parts matter at first order for a bath spin whose own z field nearly cancels (a close neighbour's
J<sub>zz</sub> m<sub>j</sub> against its Zeeman term, as in clustered P1 baths), or whose axis the qubit tilts. There
the z-only field can make the cluster expansion diverge. Against exact 10-spin subsystems, `"xyz"` removed the
divergence of a broken state but was not more accurate in every normal region, so it is not the default.

The transverse field is added and removed only between main bath spins. Every mean-field term involving a defect
sub-spin (e.g. a P1's <sup>14</sup>N) stays z only, while the couplings inside a cluster keep their full tensors.

- `example/tutorial/basic/06_defect_bath/cCCE2_Hahn_meanField/` runs `"off"`, `"z"` and `"xyz"` on a P1 bath.
- `test/code_verification/check_analytic` checks that with `"off"` two spins at order 1 give the product of each spin
  alone in its bath, which `"z"` does not. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"off"` | no static field of the other bath spins |
| `"z"` | its z part (default) |
| `"xyz"` | the vector J·m |

<b>Default:</b> `"z"`

<b>Example:</b> ``` "hamiltonian" : {"common": {"bath_bath_mean_field": "xyz"}} ```

<b>Related Tags:</b> [s_range](#tag_s_range)

<a id="tag_hyperfine"></a>
## hyperfine
<b>Description:</b><br>
`hamiltonian.common`. Where the hyperfine tensors between the qubit and the bath spins come from.
- `dipole` is the anisotropic part. `"point"` computes the point-dipole tensor from the positions; `"dft"` takes the
  DFT dipolar tensor of [file](#tag_hyperfine_file).
- `contact` adds the isotropic Fermi-contact term of `file`.
- `file` is required when `dipole` is `"dft"` or `contact` is `true`. [missing_row](#tag_hyperfine_missing_row) is
  then read too.
- With [frame.rotation](#tag_frame_rotation) on, [frame](#tag_hyperfine_frame) gives the frame of the file's tensor
  components. <br>

<b>Data Type:</b> Object `{"dipole": "point"|"dft", "contact": true|false, "file": ..., "missing_row": "stop"|"point_dipole", "frame": "file"|"qubit"}`

<b>Options:</b><br>
| dipole | contact | Description |
|---|---|---|
| `"point"` | `false` | point-dipole tensor (the default) |
| `"point"` | `true` | Fermi-contact term + point-dipole tensor |
| `"dft"` | `false` | DFT dipolar tensor |
| `"dft"` | `true` | Fermi-contact term + DFT dipolar tensor |

<b>Default:</b> `{"dipole": "point", "contact": false}`

<b>Example:</b> ``` "hamiltonian" : {"common": {"hyperfine": {"dipole": "dft", "contact": true, "file": "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_v3"}}} ```

<a id="tag_hyperfine_file"></a>
## hyperfine.file
<b>Description:</b><br>
The hyperfine tensor file, read when [hyperfine](#tag_hyperfine) `dipole` is `"dft"` or `contact` is `true`.
`example/tools/SpinTensorReader` writes it from a DFT run. After a header it holds:
- one line per species with its g-factor and its default tensor (`etc____`, used for atoms outside the range);
- the range of relative positions the file covers (`MinDif[A]`, `MaxDif[A]`);
- one row per atom: species, position relative to the defect (Å), the isotropic part (MHz) and the nine tensor
  components (MHz).

The values are per unit nuclear g-factor; CCEX multiplies them by the species' g-factor. <br>

Format v3:
- Line 1 is the tag of the DFT code (`VASP` or `QUANTUMESPRESSO`).
- Line 2 is the provenance line
  `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`.
- The numbers are the hyperfine of the spin-S defect whichever code made them, and CCEX reads them as written.
  SpinTensorReader divides Quantum ESPRESSO's S = 1/2 coupling by 2S; VASP's OUTCAR is already per S.
- Refused: a first line that is not the tag of the source, the provenance line elsewhere than line 2 or twice, an
  unknown format or source, and an S other than the qubit spin. <br>

A file without the provenance line (written before format v3) is read by the earlier rule, with a warning in the log:
first line `VASP` -> × 1, anything else -> × 0.5/S with S the qubit spin. A SpinTensorReader `--version v1` file
(OUTCAR × 2) is right only untagged and for S = 1; with a `VASP` line it is read 2×. Remake old files with
SpinTensorReader. <br>

A bath spin is matched to a row by its position relative to the qubit.
- A spin outside the range gets a point-dipole tensor.
- A spin inside the range with no row is decided by [missing_row](#tag_hyperfine_missing_row).

The positions are in the bath file's frame; the frame of the components is [hyperfine.frame](#tag_hyperfine_frame).
One qubit only. <br>
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

<a id="tag_hyperfine_missing_row"></a>
## hyperfine.missing_row
<b>Description:</b><br>
What to do with a bath spin whose position relative to the qubit lies inside the range the hyperfine tensor file
covers but has no row in the file. Spins outside the range always get a point-dipole tensor. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"stop"` | stop the run and print the position that was not found |
| `"point_dipole"` | give the spin a point-dipole tensor and go on |

<b>Default:</b> `"stop"`

<b>Example:</b> ``` "missing_row" : "point_dipole" ``` (inside `hyperfine`)

<a id="tag_hyperfine_frame"></a>
## hyperfine.frame
<b>Description:</b><br>
The frame the tensor components of [hyperfine.file](#tag_hyperfine_file) are written in, read with
[frame.rotation](#tag_frame_rotation) on. It replaces `frame.rotation.hf_tensor_frame`. The position columns of the
file are always in the bath file's frame.
- Required with the rotation on and `dipole` `"dft"`.
- With `dipole` `"point"` the tensors are point dipoles from the positions and are always turned (with `contact`
  `true` the file gives only the isotropic term, which has no direction). A written frame is then not used (the log
  notes it).
- With the rotation off, a written frame is not used (a warning). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"file"` | the components are in the bath file's frame: turned with the positions, R T R<sup>T</sup> (`"bath"` before 2.0.0) |
| `"qubit"` | the components are already in the computational frame (z = `qubit_axis`): used as written; a spin the file does not cover still gets a turned point-dipole tensor |

<b>Default:</b> none (required with the rotation on and `dipole` `"dft"`)

<b>Example:</b> ``` "frame" : "file" ``` (inside `hyperfine`)

<a id="tag_quadrupole"></a>
## quadrupole
<b>Description:</b><br>
`hamiltonian.common`. Where the quadrupole (electric field gradient) tensors of the bath spins come from. `"off"` is
no quadrupole interaction. An object reads a DFT [file](#tag_quadrupole_file):
- A bath spin inside the range the file covers takes the row at its position relative to the defect. A spin inside
  the range with no row stops the run.
- With `structure` `"plain"`, a spin outside the range takes the species default of the file (`etc`).
- With `"strained"` and `"bubble"`, a spin outside the range takes the tensor of the structure without the defect from
  [file_without_defect](#tag_quadrupole_file_without_defect). Its relative position is folded by
  [cell](#tag_quadrupole_cell) along every axis whose cell parameter is not zero.

One qubit only. With [frame.rotation](#tag_frame_rotation) on, [frame](#tag_quadrupole_frame) gives the frame of the
files' tensor components. <br>

<b>Data Type:</b> `"off"`, or Object `{"structure": ..., "file": ..., "file_without_defect": ..., "cell": [a, b, c], "frame": "file"|"qubit"}`

<b>Options:</b><br>
| structure | Description |
|---|---|
| (`"off"`) | no quadrupole interaction (the default) |
| `"plain"` | one DFT file computed with the defect, `file` |
| `"strained"` | strained structure: `file` with the defect and `file_without_defect` without it; `cell` required |
| `"bubble"` | bubble structure: the same two files and `cell` |

<b>Default:</b> `"off"`

<b>Example:</b> ``` "hamiltonian" : {"common": {"quadrupole": {"structure": "plain", "file": "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE"}}} ```

<a id="tag_quadrupole_file"></a>
## quadrupole.file
<b>Description:</b><br>
The quadrupole tensor file computed with the defect. It has the layout of [hyperfine.file](#tag_hyperfine_file),
with two differences:
- the electric field gradient components (Hartree/Bohr²) in place of the hyperfine ones;
- one `eQ___` line per species (its quadrupole moment) in place of the g-factor.

Written by `example/tools/SpinTensorReader`. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required with every `structure`)

<b>Example:</b> ``` "file" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ``` (inside `quadrupole`)

<a id="tag_quadrupole_file_without_defect"></a>
## quadrupole.file_without_defect
<b>Description:</b><br>
The quadrupole tensor file of the same structure computed without the defect, in the layout of
[quadrupole.file](#tag_quadrupole_file). Read with `structure` `"strained"` or `"bubble"`. It supplies the tensors of
the spins outside the range of the file with the defect. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required with `"strained"` and `"bubble"`)

<b>Example:</b> ``` "file_without_defect" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ``` (inside `quadrupole`)

<a id="tag_quadrupole_cell"></a>
## quadrupole.cell
<b>Description:</b><br>
Lattice parameters of the strained or bubble structure. Read with `structure` `"strained"` or `"bubble"`. <br>

<b>Data Type:</b> Double List (3), or `{"value": [a, b, c], "unit": ...}` with the unit `A` or `nm`
([Units](InputFile.md#Units))

<b>Default:</b> - (required with `"strained"` and `"bubble"`)

<b>Unit:</b> angstrom

<b>Example:</b> ``` "cell" : [1.0, 2.0, 3.0] ``` (inside `quadrupole`)

<a id="tag_quadrupole_frame"></a>
## quadrupole.frame
<b>Description:</b><br>
The frame the tensor components of [quadrupole.file](#tag_quadrupole_file) (and
[file_without_defect](#tag_quadrupole_file_without_defect)) are written in, read with
[frame.rotation](#tag_frame_rotation) on. It replaces `frame.rotation.qd_tensor_frame`.
- Required with the rotation on and a quadrupole structure.
- With the rotation off, a written frame is not used (a warning). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `"file"` | the components are in the bath file's frame: turned with the positions (`"bath"` before 2.0.0) |
| `"qubit"` | the components are already in the computational frame: used as written |

<b>Default:</b> none (required with the rotation on and a quadrupole structure)

<b>Example:</b> ``` "frame" : "file" ``` (inside `quadrupole`)

<a id="tag_hfmedi"></a>
## hfmedi
<b>Description:</b><br>
`hamiltonian.cCCE`. Adds the hyperfine-mediated interaction between bath spins: the second-order coupling through the
qubit (N. Zhao et al., Nature Nanotechnology 6, 242 (2011), Supplementary, Eq. 9).

CCEX adds the full second-order correction, for each qubit branch α on its own:
$H_\text{med}^{(\alpha)} = \sum_{\beta\neq\alpha} E_{\alpha\beta} E_{\alpha\beta}^\dagger / (\omega_\alpha - \omega_\beta)$ with $E_{\alpha\beta} = \sum_i \langle\alpha|\mathbf{S}|\beta\rangle\cdot\mathsf{A}_i\cdot\mathbf{I}_i$.
Here α is the eigenstate of the qubit Hamiltonian closest to the branch state, and β runs over every other eigenstate
(every intermediate state). Expanded, the sum holds both orders of every pair (i, j) and the i = j terms of each bath
spin; it is Hermitian.

`example/tutorial/advanced/06_terms_and_models/cCCE2_Hahn_hfMediated/` runs with and without it. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"hfmedi": true}} ```

<a id="tag_knight"></a>
## knight
<b>Description:</b><br>
`hamiltonian.cCCE`. Adds the Knight-field shift of the bath spins that the qubit induces under a transverse magnetic
field (same reference as [hfmedi](#tag_hfmedi), SI Eq. 14). For each bath spin:

$H_K^{(\alpha)} = -\gamma_n \mathbf{B}\cdot(\mathsf{g} - 1)\cdot\mathbf{I} = -(3|m_s| - 2)\,(\gamma_e/\Delta)\,\mathbf{B}\cdot\tilde{\mathsf{A}}\cdot\mathbf{I}$

- Δ = 3/2 of the zz entry of the qubit's zero-field splitting (the qubit's self term in `Qubit.intmap`);
- Ã is the hyperfine tensor with its z row set to zero;
- $|m_s| = \langle S_z^2\rangle$ of the branch's eigenstate.

The term is zero when the field is along z. A qubit without a zero-field splitting (Δ = 0) is refused, since the term
divides by Δ (`Error: hamiltonian.cCCE.knight divides by the zero-field splitting of the qubit ...`).

`example/tutorial/advanced/06_terms_and_models/cCCE2_Hahn_knightShift/` runs with and without it in a tilted field. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"knight": true}} ```

<a id="tag_delta_A"></a>
## delta_A
<b>Description:</b><br>
`hamiltonian.cCCE`. Adds the second-order hyperfine correction δA of each bath spin (J. R. Maze, J. M. Taylor and
M. D. Lukin, PRB 78, 094303 (2008), Appendix, eqs. A3 and A4). It is of the same family as [knight](#tag_knight)
(Maze's δg) and [hfmedi](#tag_hfmedi). For each qubit branch α and each bath spin n, CCEX adds
δA<sub>n</sub>(m<sub>s</sub>)·I<sub>n</sub> to the branch's bath Hamiltonian, with

δA<sub>n</sub>(m<sub>s</sub>) = −m<sub>s</sub> / (2Δ) · (A<sub>x</sub> × A<sub>y</sub>),  (A<sub>x</sub> × A<sub>y</sub>)<sub>k</sub> = Σ<sub>ij</sub> ε<sub>ijk</sub> A<sub>xi</sub> A<sub>yj</sub>,

- A<sub>x</sub>, A<sub>y</sub>: the x and y rows of the spin's hyperfine tensor (row = qubit component, the qubit axis
  is z);
- Δ = 3/2 of the zz entry of the qubit's zero-field splitting (as for knight);
- m<sub>s</sub>: the branch's m<sub>s</sub> (`Qubit.alpha_state`/`beta_state`).

Requirements:
- a single spin-1 qubit with a zero-field splitting (refused otherwise);
- not together with [hfmedi](#tag_hfmedi) (refused, see below);
- not on `execution.propagator` sparse.

When it matters: the term is of size |A|<sup>2</sup>/Δ, so for nuclei close to the defect, with transverse hyperfine
components of a MHz or more (≈ 2.5 rad/ms for a <sup>13</sup>C 2.5 Å from an NV, a phase of order one in a
millisecond). For the bulk of a bath it is small, as Maze notes, because Δ is large.

[hfmedi](#tag_hfmedi) already contains it. Its sum holds every ordered pair i ≠ j and the i = j (same-nucleus) terms:
the complete second order of the same coupling, with the qubit's own level energies. So the two together are refused.
Use delta_A for the same-nucleus δA correction without the full hyperfine-mediated term. For that <sup>13</sup>C at
5 G, hfmedi alone deviates by 3.7e-4 from the exact echo, delta_A by 4.8e-3.
`test/code_verification/check_analytic` compares delta_A with gCCE order 1 (exact for one spin) for that
<sup>13</sup>C: the Hahn-echo deviation from the exact curve goes from 0.81 to 4.8e-3 at 5 G along z.

Where it comes from: for a spin-1 qubit, H = Δ S<sub>z</sub><sup>2</sup> + S·A·I. The transverse part
S<sub>x</sub> a<sub>x</sub> + S<sub>y</sub> a<sub>y</sub> (a<sub>k</sub> = A<sub>k</sub>·I) changes m<sub>s</sub> by ±1.
To second order the branch m<sub>s</sub> = ±1 gets (1/2Δ) V<sub>±</sub>V<sub>∓</sub>, V<sub>±</sub> = a<sub>x</sub> ∓ i a<sub>y</sub>,
and m<sub>s</sub> = 0 gets −(1/2Δ)(V<sub>+</sub>V<sub>−</sub> + V<sub>−</sub>V<sub>+</sub>). Since
[a<sub>x</sub>, a<sub>y</sub>] = i (A<sub>x</sub> × A<sub>y</sub>)·I,
V<sub>±</sub>V<sub>∓</sub> = a<sub>x</sub><sup>2</sup> + a<sub>y</sub><sup>2</sup> ∓ (A<sub>x</sub> × A<sub>y</sub>)·I.
The part linear in I is ∓(A<sub>x</sub> × A<sub>y</sub>)·I/(2Δ) for m<sub>s</sub> = ±1 and cancels for
m<sub>s</sub> = 0, i.e. −m<sub>s</sub>(A<sub>x</sub> × A<sub>y</sub>)/(2Δ). The a<sub>x</sub><sup>2</sup> + a<sub>y</sub><sup>2</sup>
part is a constant for a spin 1/2 and a quadrupole-like term otherwise; it is not added.
The term has to be odd in m<sub>s</sub>: S and I both change sign under time reversal, which takes m<sub>s</sub> to
−m<sub>s</sub>.

This differs from the coefficient (2 − 3|m<sub>s</sub>|) of Maze, PRB 78, 094303 (2008), Eq. (A4). A direct
derivation, time-reversal symmetry and exact diagonalization give −m<sub>s</sub>. Whether this is a misprint or a
difference of convention is not settled. The two agree for m<sub>s</sub> = +1. For m<sub>s</sub> = ±1 one virtual path
fixes the operator order and sign; for m<sub>s</sub> = 0 the two paths cancel. Exact diagonalization of an NV with one
<sup>13</sup>C (|A| up to 1.5 MHz, B = 5 G along z): the nuclear splitting in each branch is off by 2.0 rad/ms without
the term, by 0.01 with −m<sub>s</sub>, and with (2 − 3|m<sub>s</sub>|) by 1.9 (m<sub>s</sub> = 0) and 4.0
(m<sub>s</sub> = −1).

Limits: the derivation takes the qubit levels as the |m<sub>s</sub>⟩ states, with a field along the qubit axis and
γ<sub>e</sub>B ≪ Δ (secular). The qubit Zeeman energy is left out of the denominators, as in knight, so the next
correction is of relative size γ<sub>e</sub>B<sub>z</sub>/Δ (0.5 % at 5 G, 10 % at 100 G). With a field across the
axis the levels mix. CCEX warns when the qubit eigenstate closest to |α⟩ or |β⟩ has less than 0.99 of its weight on
that |m<sub>s</sub>⟩. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"cCCE": {"delta_A": true}} ```

<a id="tag_pseudospin_model"></a>
<a id="tag_pseudospin_channels"></a>
## channels (the pseudospin model of a pair)
<b>Description:</b><br>
`hamiltonian.pseudospin`. How [method](#tag_method) pseudospin gives a pair's term P, and how its channels are evolved.

The model:
- The α and β product states of one index are taken as one state. Their overlap, the pair's share of the ESEEM, is
  left out.
- With the nuclear state I/d, each of the d product states weighs 1/d, d = (2I<sub>1</sub>+1)(2I<sub>2</sub>+1).
- The channels replace the diagonal values of their states by a coherence that includes the channels' elements.
- A state can sit in two channels (spins 1, no quadrupole: \|1,−1⟩ ↔ \|0,0⟩ ↔ \|−1,1⟩). The channels then form a
  chain.

`channels` says how the channels are evolved:

- `"grouped"` (default): the channels join the product states into blocks (a chain is one block). Every block of k
  states evolves whole, and
  L<sub>pair</sub> = L<sub>diag</sub> + Σ<sub>B</sub> (k<sub>B</sub>/d) (l<sub>B</sub> − l<sub>B,diag</sub>).
  Within the cut Hamiltonian this is exact; what it leaves out is the α-β overlap.
- `"pairwise"`: every channel alone is a two-level system with the fields (Ω<sub>α,β</sub>, 0, Δ<sub>α,β</sub>), and
  L<sub>pair</sub> = L<sub>diag</sub> + (2/d) Σ<sub>c</sub> (l<sub>c</sub> − l<sub>c,diag</sub>). For a Hahn echo with
  Ω<sub>α</sub> = Ω<sub>β</sub>,
  l<sub>c</sub> = 1 − 2 \|n̂<sub>α</sub> × n̂<sub>β</sub>\|² sin²(ω<sub>α</sub>t/4) sin²(ω<sub>β</sub>t/4). It misses what
  couples the channels of a chain. It is kept to compare with that independent-channel picture.

L<sub>diag</sub>, l<sub>B,diag</sub> and l<sub>c,diag</sub> are without the channels' elements. The two agree when no
two channels share a state. Both are CCEX's own generalization of the papers' single flip-flop to any spins and
several channels, not a formula of either paper.

Then P = L<sub>pair</sub> / (L<sub>1</sub> L<sub>2</sub>), with the single spins from their eigenvalues alone. For two
spins 1/2 of one species with the secular hyperfine and the flip-flop alone, both are Seo's eq. 23,
P = 1 − K sin²(√(ΔA² + D²) t/4) sin²(D t/4), K = ΔA²/(ΔA² + D²).

The division by the single spins is 1 under a Hahn echo or CPMG. Without a pulse (`pulse.count` 0) the single spins
oscillate and can pass through zero. A time point where one is zero is reported as `[numerical-invalid]`, as cCCE's
own division by a zero single coherence is.

Accuracy. The comparisons use the inputs of `test/code_verification/check_pseudospin` (two spins on the qubit axis,
Hahn, 200 points). The check itself covers spins 3/2 grouped; the other cases were run once each.
- Without a quadrupole, P times cCCE order 1 was within 1e-10 of cCCE order 2 with `"grouped"` for spins 1, 3/2 and
  3, and 0.43, 0.71 and 0.65 off with `"pairwise"`.
- An axial quadrupole (spins 1) splits the chain's levels; pairwise was then 2.9e-6 off.
- A quadrupole that mixes the \|m⟩ (spins 1, η = 0.5) left both 0.74 off. The likely cause is that the α and β bases
  differ there, which was not tested separately.
- On the h-BN case of check_numerics (DFT hyperfine and quadrupole, Hahn, 200 points), both gave P times cCCE order 1
  within 2.0e-3, 2.2e-2 and 4.2e-2 of cCCE order 2 at 3, 0.3 and 0.03 T. <br>

<b>Data Type:</b> String (`"grouped"` or `"pairwise"`)

<b>Default:</b> `"grouped"`

<b>Example:</b> ``` "method" : "pseudospin", "hamiltonian" : {"pseudospin": {"channels": "pairwise"}} ```

<a id="tag_jump_operators"></a>
## jump_operators
<b>Description:</b><br>
`hamiltonian.meCCE` or `hamiltonian.megCCE`. The Lindblad jump operators of the master-equation methods, one entry per
operator: C = √γ S<sub>op</sub> on every spin named `bath_name`, with γ = `rate` × 2π·10³ rad/ms.
- `operator`: `-` (S<sub>-</sub>), `+` (S<sub>+</sub>), `z`, `x` or `y` (S<sub>z</sub>, S<sub>x</sub>,
  S<sub>y</sub>) of that spin, of any S.
- `bath_name`: the name of bath spins (as in the bath file, or of the Defect sub-spins; case not significant), or
  `"qubit"` / the qubit's name for a jump of the qubit. Only megCCE takes a qubit jump (one qubit); meCCE refuses it
  and names megCCE.
- `rate`: in MHz, ≥ 0.

Refused:
- a `bath_name` that is neither, checked against the whole bath once it is read; the message gives the names the
  bath has. (A cluster without that species simply gets no jump from it.)
- a `bath_name` that is both the qubit's name and a bath spin's (megCCE would apply it to the qubit only);
- the section in a run of another method. The section holds `jump_operators` only.

Without the key, or with an empty list, the run is the rate-0 limit: the cCCE (meCCE) or the gCCE (megCCE). The
equations, what is supported and the checks are in [master-equation.md](../methods/master-equation.md).
`example/tutorial/advanced/06_terms_and_models/meCCE2_Hahn_bathRelaxation/` runs S₋ and S₊ on the ¹³C against cCCE. <br>

<b>Data Type:</b> Array of objects `{"bath_name": string, "operator": string, "rate": number}`; `rate` may also be
`{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

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

<a id="tag_noise_spectrum"></a>
## spectrum
<b>Description:</b><br>
`hamiltonian.noise`, method `noise` only. Where the noise of the qubit transition frequency comes from.

The coherence is $W(T) = e^{-\chi(T)}$ with $\chi = \langle\varphi^2\rangle/2$, where $\varphi = \int f(t)\,\delta\omega(t)\,dt$
is the qubit's phase and $f = \pm 1$ flips at each ideal π pulse of the sequence. This is the Gaussian approximation,
which needs many weakly coupled noise sources. It has no back-action of the qubit on the bath and no non-Gaussian
part, so strongly coupled spins need cCCE or gCCE.

- `"bath"` (the default): the bath of the run. Each cCCE cluster gives
  $\chi_C = \frac{1}{2D}\sum_{ab} |\langle a|\delta\omega|b\rangle|^2\, |Y(\omega_{ab})|^2$ with
  $Y(\omega) = \int f(t)\, e^{i\omega t}\, dt$, over the eigenstates of
  $H_\text{bath} + p \sum_i \mathbf{A}_{z,i}\cdot\mathbf{I}_i$, $p = (|m_\alpha| + |m_\beta|)/2$, at infinite
  temperature, and $\delta\omega = (m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$.
  - The $\chi_C$ add over the CCE terms (the cCCE product of the $e^{-\chi_C}$).
  - With `s_range`, a bath state enters only through the mean field of the other bath spins
    ([bath_bath_mean_field](#tag_bath_bath_mean_field)).
  - For a Hahn echo this is `semi-cce`'s `mode_sum`. $p$ and the $\delta\omega$ with $A_{zz}$ alone are its
    conventions ([pseudo_secular](#tag_noise_pseudo_secular) adds $A_{zx} I_x + A_{zy} I_y$ to $\delta\omega$).
- `{"model": "OU", "sigma": σ, "tau_c": τ_c}`: Ornstein–Uhlenbeck noise of the frequency,
  $\langle\delta\nu(t)\delta\nu(0)\rangle = \sigma^2 e^{-|t|/\tau_c}$ (σ in kHz, τ_c in ms; both take a unit). χ is
  exact for any ideal π sequence.
- `{"model": "power", "amplitude": A, "beta": β, "f_min": f_min, "f_max": f_max}`:
  $S(f) = A\,(f/1\,\text{kHz})^{-\beta}$ (kHz²/kHz, two-sided) on [f_min, f_max] (kHz, both take a unit), and
  $\chi = 4\pi^2 \int S(f)\, |Y(2\pi f)|^2\, df$ by Gauss–Legendre. The work grows as $f_\text{max} T$.

With `"bath"`, the clusters' values ($-\chi_C$, or their correlations with [compute](#tag_noise_compute)) are summed
over the CCE terms; cCCE multiplies its coherences instead. So no $e^{-\chi_C}$ of one cluster can underflow. With
method `noise`, `output.savemode` `cluster_coherence` and `cluster_density_matrix`, `rdip_ladder` and
`output.cluster_contributions` are refused.

A model run reads no bath:
- `configuration`, the Defect block and the clusters are not used;
- `order`, `bfield`, `rbath`, `rdip` and `Qubit` may be left out. The run then takes order 1, bfield 0, rbath 0, rdip 0
  and one qubit at the origin, which it does not use. When they are given they are not used either, and the log warns
  that they are not used.

<b>Data Type:</b> String `"bath"` or Object

<b>Default:</b> `"bath"`

<b>Example:</b> a whole input of a model run
    ```
    {"method": "noise", "hamiltonian": {"noise": {"spectrum": {"model": "OU", "sigma": 10, "tau_c": {"value": 50, "unit": "us"}}}},
     "deltat": 0.002, "nstep": 200, "pulse": {"count": 1}}
    ```

<a id="tag_noise_pseudo_secular"></a>
## pseudo_secular
<b>Description:</b><br>
`hamiltonian.noise`, with [spectrum](#tag_noise_spectrum) `"bath"`. Which terms of the hyperfine coupling between the
qubit and a bath spin, $m_s (A_{zx} I_x + A_{zy} I_y + A_{zz} I_z)$, count as the qubit's frequency noise.
- `false`: the secular term alone, $\delta\omega = (m_\alpha - m_\beta) \sum_i A_{zz,i} I_{z,i}$ (`semi-cce`'s).
- `true`: also the pseudo-secular terms, $\delta\omega = (m_\alpha - m_\beta) \sum_i (A_{zx,i} I_{x,i} + A_{zy,i} I_{y,i} + A_{zz,i} I_{z,i})$,
  the coupling cCCE puts into its two Hamiltonians.

The bath Hamiltonian is the same either way. In `semi-cce`'s 13C example (500 G, Hahn echo) the two differ by up to
0.04 in W, in its P1 example by 1e-3. Which of them is the better approximation has not been settled (both are far
from cCCE there). With a model spectrum the key is not used (the run says so). <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hamiltonian" : {"noise": {"spectrum": "bath", "pseudo_secular": true}} ```

<a id="tag_noise_compute"></a>
## compute
<b>Description:</b><br>
`hamiltonian.noise`, with [spectrum](#tag_noise_spectrum) `"bath"`. What the run computes from the bath's noise
$B = \sum_i A_{zz,i} I_{z,i}$ ($\mathbf{A}_{z,i}\cdot\mathbf{I}_i$ with [pseudo_secular](#tag_noise_pseudo_secular)),
and which files it writes:

| Value | Computes | Writes |
|---|---|---|
| `"coherence"` (default) | $\chi_C$ of each cluster from its eigenmodes, the coherence exp(−χ) at the run's times; no time grid of the noise, exact for ideal π pulses | `_wiDiv` (and `_noDiv`) |
| `"spectrum"` | `semi-cce`'s `psd`: the correlation C(t) on the grid of [fft](#tag_noise_fft), its spectrum S(f) by FFT, and $\chi = 2 (m_\alpha - m_\beta)^2 \int S(f)\, F(2\pi f T)/f^2\, df$ by the trapezoid rule | `_noise_correlation`, `_noise_spectrum`, `_wiDiv` |
| `"correlation"` | PyCCE's noise (`CCENoise`): C(t) at the run's times (`deltat`, `nstep`) | `_noise_correlation` |

$C(t) = \operatorname{Re}\operatorname{Tr}[\rho\,\{B(t), B\}]/2$ in kHz², summed over the CCE terms
(`<outfile>[_state<s>]_noise_correlation`: t in ms, C). S(f) is two-sided, in kHz²/kHz (`_noise_spectrum`: f in kHz,
S).

The state ρ:
- `"coherence"` and `"spectrum"`: ρ is 1/D (infinite temperature). A sampled state (`s_range`) enters through the mean
  field alone, as in `semi-cce`.
- `"correlation"`: with `s_range`, ρ is the cluster's pure state, with the mean field of the other spins, as in PyCCE;
  1/D without `s_range`.

So without `s_range` the three agree. `"correlation"` and `"spectrum"` give the same C(t). `"spectrum"` gives the
coherence of `"coherence"` to the accuracy of its grid (a step below the period of the fastest bath frequency, a
window longer than the times asked for). The agreement is about 1e-9 at best and up to 5e-9 on the default grid: the
rounding of C(t) grows through the FFT and the integral, so the result moves by that much with the number of ranks.
With `s_range`, `"correlation"` is another quantity.

A model spectrum computes the coherence alone (another value is refused). <br>

<b>Data Type:</b> String, `"coherence"`, `"spectrum"` or `"correlation"` (case is not significant)

<b>Default:</b> `"coherence"`

<b>Example:</b> ``` "hamiltonian" : {"noise": {"compute": "correlation"}} ```

<a id="tag_noise_fft"></a>
## fft
<b>Description:</b><br>
`hamiltonian.noise`, with [compute](#tag_noise_compute) `"spectrum"` (not used otherwise; every key has a default).
The time grid of C(t) and the FFT, `semi-cce`'s `psd` keys.

| Key | Meaning | Default |
|---|---|---|
| `time_max` | the end of the grid (ms, takes a unit): C(t) at t = 0, Δt, ..., `time_max` | the run's last time, `deltat` × (`nstep` − 1) of the top level |
| `nstep` | its intervals (`nstep` + 1 points, Δt = `time_max`/`nstep`) | from the bath's modes, see below |
| `window` | `"none"`, `"hann"`, `"tukey"` (the last half tapered) or `"exp"` (exp(−3t/`time_max`)) on C(t) | `"none"` |
| `pad_factor` | the FFT length is `pad_factor` × (2 `nstep` + 1), made odd; the padding holds the mean of the last fifth of C(t) | 1 |
| `subtract_mean` | subtract that mean from C(t) first (and pad with 0) | `false` |
| `clamp_negative` | set S(f) < 0 to 0 before the coherence | `false` |

S(f) = Δt Re FFT of C(t) made even, at f = k/(M Δt), k = 0 .. M/2 for an FFT of length M. The highest frequency is
1/(2Δt) and the resolution 1/(M Δt). Δt has to resolve the fastest bath frequency (the ¹³C Larmor period is 1.9 µs at
500 G), and `time_max` has to cover the times of the coherence. <br>

Without `nstep`:
- Δt gives 8 points per period of the fastest modes left after dropping at most 10⁻³ of the weight (the modes' signed
  weights over the CCE terms, in frequency bands). The log prints the values taken.
- The coherence then stays within 5×10⁻⁷ of `compute` `"coherence"` on 1 ppm P1 and 2×10⁻¹⁰ on natural-abundance ¹³C.
- All the modes stay in memory until the end of the run (24 bytes each; a cluster has dimension² of them). Give
  `nstep` if there are too many. The run stops if `nstep` would pass 10⁸.

The default `time_max` is what the coherence needs; it does not resolve S(f). For a finer S(f), or with a `window`
(it would taper C(t) inside the times of the coherence), give a longer `time_max`. <br>

<b>Data Type:</b> Object

<b>Example:</b> ``` "hamiltonian" : {"noise": {"compute": "spectrum", "fft": {"time_max": 3.0, "nstep": 15000, "window": "hann"}}} ```

<a id="tag_not_used"></a>
## Keys read and not used
<b>Description:</b><br>
A key this run reads but does not use is reported once, with the same line for every such key:
`Warning: input key "<key>" is not used by this run: <why>`. The calculation goes on; the result is the one without
the key. The cases:

| Key | Not used when |
|---|---|
| `Qubit.psi0` | cCCE, meCCE (the qubit is projected out: alpha and beta only) |
| `hamiltonian.common.qubit_bath_mean_field` | cCCE, meCCE, and the ensemble (no sampled state) |
| `hamiltonian.common.bath_bath_mean_field` | the ensemble |
| `configuration.bath_state_file`, `defect_state_file` | the ensemble (no `s_range`) |
| `configuration.defect_axis_file`, `defect_state_file` | no Defect species |
| `pulse.duration` | gCCE with `pulse.detuning_factor` 1 (cCCE, meCCE and megCCE refuse it) |
| `frame` (any) | `frame.rotation` off; a Defect species' `frame` without directional numbers |
| `cluster.selection.max_per_order` | an entry 0-2 other than 0, or every entry with dsj, itb, dsjitb |
| `execution.gpu_options` | `main.out` (built without CUDA; `main_gpu.out` uses them) |
| `cluster.adaptive.policy` | always: `merge` is the only policy |

A key that no part of the run reads is reported as `Warning: input key "<key>" is not read by this run and is ignored`.
A key that looks like a typo of a key that is read stops the run.

<a id="tag_turn_off"></a>
## Turning a part off
<b>Description:</b><br>
How each part that can be switched off is switched off. Leaving the key out gives the default in the last column.

| Key | Off | Without the key |
|---|---|---|
| `hamiltonian.common.quadrupole` | `"off"` | off |
| `hamiltonian.common.bath_bath_mean_field` | `"off"` | `"z"` |
| `cluster.pcce.cache` | `{"mode": "off"}` | off |
| `frame.rotation` | `{"enabled": false}` | off |
| `execution.gcce_refinement` | `{"enabled": false}` | off |
| `execution.column_evolution` | `{"enabled": false}` | on |
| `output.cluster_contributions` | `{"enabled": false}` | off |
| `auto_convergence` | `{"enabled": false}` | off |

A block written with no `"enabled"` is on.
