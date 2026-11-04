# Input file
<a id="InputFile"></a>

The input file is JSON and sets every option of a run. `main.out -e` writes a template with every key
([Getting started](GettingStarted.md#OwnInput)). This page shows how an input is written, the keys most runs need,
and the units. Every key is described on the page of its group. An input of an earlier version can be converted with
[Upgrading old inputs](InputFile_OldNames.md).

1. [How to write the input file](#HowToWrite): the format and a whole input.
2. [The keys most runs need](#CommonKeys): one table, with links.
3. [Units](#Units): the unit of every number.
4. [The pages of the keys](#KeyGroups): what each group page covers.

<a id="HowToWrite"></a>
## 1. How to write the input file

The input is a JSON object of keys and values, named `.json` or `.jsonc`. A `//`, `#` or `!` outside a string starts
a comment to the end of the line. Keys are read without regard to case, and a key that is left out takes its default, where it has one.

A whole input, the one of the first tutorial (run from its folder,
`example/tutorial/basic/01_first_run`):

```json
{
    "order"        : 2,
    "bfield"       : 500,                       // G
    "rbath"        : 50,                        // Å
    "rdip"         : 10,                        // Å
    "deltat"       : 0.001,                     // ms
    "nstep"        : 2000,
    "pulse"        : {"count": 1},              // Hahn echo
    "Qubit"        : {"qubit": [{"file": "bath/bath_Diamond_NV_natab_defect"}]},
    "configuration": {"gyrofile": "bath/gyro_13C", "bathfile": ["bath/bath_Diamond_NV_natab_1"]}
}
```

The run checks the input before it reads the bath:

- a value of the wrong type stops the run and names the key (`"nstep": 10.5`, `"knight": "true"`);
- a key one or two letters away from a real key stops the run (`"s_rang"`);
- a key the run does not read is reported in the log and ignored;
- a key of an old version is refused with its current name ([Upgrading old inputs](InputFile_OldNames.md)).

<a id="CommonKeys"></a>
## 2. The keys most runs need

| Key | What it sets | Example |
|---|---|---|
| [method](InputFile_General.md#tag_method) | the CCE method | `"cCCE"`, `"gCCE"` |
| [cluster.mode](InputFile_Cluster.md#tag_cluster_mode) | how the bath is cut into clusters | `"standard"`, `"pcce"` |
| [order](InputFile_General.md#tag_order) | the largest cluster | `2` |
| [bfield](InputFile_General.md#tag_bfield) | the magnetic field (G) | `500`, `[0, 0, 500]` |
| [rbath](InputFile_General.md#tag_rbath), [rdip](InputFile_General.md#tag_rdip) | the bath radius and the pair radius (Å) | `50`, `10` |
| [deltat](InputFile_General.md#tag_deltat), [nstep](InputFile_General.md#tag_nstep) | the time step (ms) and the number of points | `0.001`, `2000` |
| [pulse.count](InputFile_Pulse.md#tag_pulse_count) | the number of π pulses: 0 Ramsey, 1 Hahn echo, n CPMG | `1` |
| [Qubit.qubit](InputFile_Qubit.md#tag_qubit) | the qubit and its position | `[{"file": "defect"}]` |
| [alpha_state](InputFile_Qubit.md#tag_alpha_state), [beta_state](InputFile_Qubit.md#tag_beta_state) | the two qubit levels of the coherence | `{"ms": [1]}`, `{"ms": [0]}` |
| [configuration.bathfile](InputFile_General.md#tag_bathfile) | the bath files | `["bath_1"]` |
| [configuration.s_range](InputFile_General.md#tag_s_range) | sampled bath states (without it, the ensemble) | `[1, 10]` |
| [configuration.k_range](InputFile_General.md#tag_configuration_k_range) | several bath configurations in one run | `[1, 10]` |
| [configuration.seed](InputFile_General.md#tag_seed) | the random seed of the sampled states and the defect axes | `1` |
| [output.path](InputFile_Output.md#tag_output_path), [output.file](InputFile_Output.md#tag_outfile) | the output folder and name | `"results"`, `"ccein"` |
| [execution.backend](InputFile_Backend.md#tag_backend) | CPU or GPU | `"cpu"`, `"auto"` |

<a id="Units"></a>
## 3. Units

A plain number is in the unit of the table. These keys also take `{"value": ..., "unit": "..."}`, with the value in
the plain form (a number, or a list for `bfield`, `rbath`, `rdip`, `bath_shift`, `xyz` and `energy_shift`):

| Key | Unit of a plain number | Other units |
|---|---|---|
| `deltat`, `auto_convergence.time_window.value` | ms | `s`, `ms`, `us`, `ns` |
| `pulse.duration` | ns | `s`, `ms`, `us`, `ns` |
| `bfield` | G | `G`, `mT`, `T` |
| `rbath`, `rdip`, `configuration.bath_shift`, `Qubit.qubit[].xyz`, `hamiltonian.common.quadrupole.cell` | Å | `A`, `nm` |
| `pulse.bath.energy_shift`, `pulse.bath.bandwidth`, `hamiltonian.<method>.jump_operators[].rate` | MHz | `Hz`, `kHz`, `MHz`, `GHz` |
| `hamiltonian.noise.spectrum` `sigma`, `f_min`, `f_max` | kHz | `Hz`, `kHz`, `MHz`, `GHz` |
| `hamiltonian.noise.spectrum.tau_c`, `hamiltonian.noise.fft.time_max` | ms | `s`, `ms`, `us`, `ns` |
| `Qubit.qubit[].detuning` | kHz | `Hz`, `kHz`, `MHz`, `GHz` |
| `cluster.bundle.merge_cutoff` | rad/ms | `rad/s`, `rad/ms`, `rad/us` |
| `Qubit.qubit[].gyro` | rad kHz/G | the units of Defect `gyros` |

```json
"deltat": {"value": 5, "unit": "us"}          // = 0.005
"bfield": {"value": [0, 0, 50], "unit": "mT"} // = [0, 0, 500]
```

Units are case-sensitive (`mT`, not `MT`).

<a id="FixedUnits"></a>
**Watch the two blocks.** The zero-field splitting and the detuning are in **kHz** in the `Qubit` block
(`Qubit.intmap` without a `unit`, `Qubit.qubit[].detuning`) but in **MHz** in the `Defect` block (`zfs`, `detuning`).
The Defect keys and an intmap entry take their own `unit`
([Qubit tags](InputFile_Qubit.md#tag_intmap), [Defect tags](InputFile_Defect.md)).

The files have fixed units: positions in Å, the gyrofile's γ in rad/ms/G, the hyperfine tensor file in MHz, the
quadrupole tensor file in Hartree/Bohr², and the time of the output files in ms.

<a id="KeyGroups"></a>
## 4. The pages of the keys

<a id="GeneralTag"></a>
- [General tags](InputFile_General.md): the method and the numbers of the calculation, the `configuration` block
  (bath files, states, seed), the frame rotation, and the `hamiltonian` block.
<a id="QubitArrayTags"></a>
- [Qubit tags](InputFile_Qubit.md): the qubits, their interactions and the two states of the coherence.
<a id="ClusterTags"></a>
- [Cluster tags](InputFile_Cluster.md): how the bath is cut into clusters (standard, pcce, bundle) and the adaptive
  refinement.
<a id="BackendTags"></a>
- [Execution tags](InputFile_Backend.md): CPU or GPU, the propagator, and the memory check.
<a id="PulseTags"></a>
- [Pulse tags](InputFile_Pulse.md): the qubit's pulse sequence and DEER pulses on the bath.
<a id="DefectArrayTags"></a>
- [Defect tags](InputFile_Defect.md): electron-spin defects in the bath, such as P1 centres.
<a id="OutputTags"></a>
- [Output tags](InputFile_Output.md): where the outputs go and what they hold.
<a id="AutoConvergenceTags"></a>
- [Automatic convergence tags](InputFile_AutoConvergence.md): the `auto_convergence` block.
