/**
 * @file example_input.cpp
 * @brief The example input (ccex_example.jsonc) that main.out -e writes.
 */
// The input template `main.out -e` writes, the one place it is kept. The active part is the calculation of
// example/tutorial/basic/01_first_run (the NV-13C cCCE2 Hahn echo, ensemble, rbath 50 A, 2 ms) without its files: gyrofile is
// commented (the built-in isotope table) and bathfile is a name the user replaces, so the template does not run as it
// is. Every other input key the readers accept (src/input/option.cpp, the adaptive options table, the automatic-
// convergence runner) is a commented line with its default. Each key has its explanation on the comment lines above it
// (no comment after a key), the ':' of one object's keys share a column ("key" : value), and each section starts with
// its docs/user page and tutorials. test/code_verification checks it: with a real bath it runs, no accepted key is
// missing from it, the comment and ':' layout.
#include "../../include/input/example_input.h"

const char* const kExampleInputName = "ccex_example.jsonc";

const char* const kExampleInput = R"CCEX(// ccex_example.jsonc: an input template written by `main.out -e` (CCEX 2.1.0)
//
// The active part is the Hahn echo of an NV centre in diamond with 13C (example/tutorial/basic/01_first_run). It runs once
// "bathfile" names a real bath file, for example example/tutorial/basic/01_first_run/bath/bath_Diamond_NV_natab_1:
//
//     mpirun -n 4 main.out -f ccex_example.jsonc -o results
//
// Every other key is written below, commented out, with its default. Remove the leading // to use one (and mind the
// commas). Comments start with //, # or !. Keys are read without regard to case. A misspelt key stops the run; a key the
// run does not use is reported. Every key is described in docs/user/InputFile.md and at
// https://cqml-skku.github.io/CCEX/.
//
// Units: each number is in the unit its comment names. Keys marked "(unit)" also take {"value": 5, "unit": "us"}:
// time s, ms, us, ns; field G, mT, T; length A, nm; frequency Hz, kHz, MHz, GHz; angular frequency rad/s, rad/ms,
// rad/us. A plain frequency is kHz for the qubit (detuning, intmap) and MHz for a Defect species.
{
    // ---- method -----------------------------------------------------------------------------------------------------
    // docs: InputFile_General.md#tag_method. Tutorials: basic/05_cCCE_vs_gCCE, advanced/13_noise, advanced/14_pseudospin.

    // cCCE (default): the qubit only sets the field each bath spin feels (pure dephasing).
    // gCCE: the qubit is in every cluster (qubit-bath flip-flops, any spin and field).
    // meCCE, megCCE: cCCE and gCCE with Lindblad jump operators (hamiltonian.meCCE).
    // noise: the Gaussian phase-noise coherence from the bath or a model spectrum (hamiltonian.noise).
    // pseudospin: the bath pairs alone in the pseudospin model of nuclear flip-flops (order 2).
    "method"           : "cCCE",

    // coherence (default): L(t) between alpha and beta. dm (gCCE): the qubit density matrix at every time.
    "quantity"         : "coherence",

    // ---- the qubit --------------------------------------------------------------------------------------------------
    // docs: InputFile_Qubit.md. Tutorials: basic/01_first_run, advanced/08_multi_qubits (several qubits).
    "Qubit"            : {
        // One entry per qubit.
        "qubit"       : [
            {
                // Default "q0", "q1", ... by position in the list.
                "name"     : "NV",
                // Default 1.0.
                "spin"     : 1.0,
                // rad kHz/G (unit: kHz/G, MHz/T, ...). Default -17608.597050, the electron.
//              "gyro"     : -17608.597050,
                // kHz (unit). Adds detuning x S_z to the qubit: its levels shift. Default 0.0.
//              "detuning" : 0.0,
                // Instead of "xyz": a file whose first line is the position (x y z, Angstrom).
//              "file"     : "defect_position",
                // Position in Angstrom (unit), in the frame of the bath file. Required, or "file".
                "xyz"      : [73.73895, 75.72808, 73.97033]
            }
        ],
        // Qubit-qubit tensors. A pair of one name twice is that qubit's zero-field splitting; pairs not written get
        // the point dipole of the positions. Here the NV's D = 2.88 GHz as diag(-D/3, -D/3, 2D/3) in kHz.
//      "intmap"      : [
//          {
//              "between" : ["NV", "NV"],
//              "tensor"  : [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]],
                // Hz, kHz (default for a tensor), MHz, GHz.
//              "unit"    : "kHz",
                // With frame.rotation (then required): file (in the bath file's frame) or qubit (as written).
//              "frame"   : "qubit"
//          }
//      ],
        // L(t) is the coherence between these two states. One of {"ms": [m_1, ...]} (one m_s per qubit),
        // {"index": n} (the n-th eigenstate by descending energy; it changes with the field, the log names its ms, and
        // its largest component is made real and positive) or {"vector": [...]} (each component a number or [re, im]). Defaults {"ms": [1, 1, ...]} and {"ms": [0, 0, ...]}; pulse.per_qubit needs {"ms"}.
        "alpha_state" : {"ms": [1]},
        "beta_state"  : {"ms": [0]}
        // The same as a vector, components ordered ms +1, 0, -1.
//      "beta_state"  : {"vector": [0, 1, 0]},
        // gCCE only: the initial qubit state as written. Default (|alpha> + |beta>) normalized.
//      "psi0"        : [0.7071, 0.7071, 0.0]
    },

    // ---- the bath ---------------------------------------------------------------------------------------------------
    // docs: InputFile_General.md#tag_configuration. Tutorials: basic/02_bath_and_states, basic/06_defect_bath.
    "configuration"    : {
        // Folder of the file names of this block. Default ".".
//      "path"              : "data",
        // Several bath configurations, one after another: each k replaces {k} in the file names and writes
        // <output.file>_k<k>_... ({{k}} is a literal {k}). Default: one configuration.
//      "k_range"           : [1, 10],
        // Isotope, spin and gyromagnetic ratio per line. Default: the built-in isotope table (needed for P1 and other
        // spins that are not isotopes).
//      "gyrofile"          : "gyro.txt",
        // FILL IN: one or more bath files (merged): a count line, then x y z isotope per line (Angstrom).
        "bathfile"          : ["my_bath.txt"]
        // Angstrom (unit): one [x, y, z] shift per bath file. Default [0, 0, 0] each.
//      "bath_shift"        : [[0.0, 0.0, 0.0]],
        // Sampled pure bath states [first, last], each its own output _state<s>. Absent (default): the ensemble,
        // the average over the unpolarised bath.
//      "s_range"           : [1, 20],
        // With s_range: the state files, one name per bath file ("" for one drawn from seed); {s} is the state, {k}
        // the configuration. A count line, then one m_s per bath spin. example/tools/StateFiles writes them.
//      "bath_state_file"   : ["cfg{k}_state{s}"],
        // With Defect species: the axis state (1 .. n_axis_states) of each spin. "" draws them from seed.
//      "defect_axis_file"  : ["cfg{k}_axis"],
        // With s_range and on-site spins: their m_s per bath spin. Default: drawn from seed.
//      "defect_state_file" : ["cfg{k}_onsite{s}"],
        // Random seed of the first configuration (k takes seed + k - first), or a list of one per configuration.
        // Default: generated and printed in the log.
//      "seed"              : 1
    },

    // ---- the expansion and the time axis ----------------------------------------------------------------------------
    // docs: InputFile_General.md#tag_order. Required, except where noted.

    // The largest cluster: bath spins (standard), partitions (pcce) or bundles (bundle). Raise it until L stops changing.
    "order"            : 2,
    // Gauss (unit): Bz, or [Bx, By, Bz]. 0 for no field.
    "bfield"           : 500,
    // Angstrom (unit): the bath spins this close to the qubit are kept, or [inner, outer].
    "rbath"            : 50,
    // Angstrom (unit): two bath spins this close can share a cluster, or [inner, outer] (not needed for bundle).
    "rdip"             : 10,
    // ms (unit): the time step. The outputs are at t = 0, deltat, ..., (nstep - 1) deltat; t is the whole sequence.
    "deltat"           : 0.001,
    // Number of time points.
    "nstep"            : 2000,

    // ---- the pulse sequence -----------------------------------------------------------------------------------------
    // docs: InputFile_Pulse.md. Tutorials: basic/04_pulses, basic/06_defect_bath (DEER).
    "pulse"            : {
        // Number of pi pulses: 0 Ramsey, 1 HahnEcho, 2 or more CPMG. Required.
        "count"           : 1
        // Ramsey, HahnEcho, CPMG, WAHUHA (count 4) or Manual (with sequence). Default: from count.
//      "name"            : "Manual",
        // Manual: count + 1 rows [fraction of t before the pulse, axis X|Y|Z|I, angle in degrees]; the fractions sum
        // to 1 and the last row is the time after the last pulse. Pulses other than pi flips need gCCE.
//      "sequence"        : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
        // false (default): the pulses rotate all qubits together; true: qubit by qubit (keep false with gCCE).
//      "per_qubit"       : false,
        // gCCE: ns (unit), with detuning_factor other than 1. Default 0.0.
//      "duration"        : 0.0,
        // gCCE: 1.0 (default) means ideal instantaneous pulses; another value a finite pulse of that detuning.
//      "detuning_factor" : 1.0,
        // DEER: pulses on bath spins of one species (gCCE and megCCE, CPU; docs/methods/deer.md).
//      "bath"            : {
            // Number of bath pulses. Default 0.
//          "count"        : 1,
            // count + 1 rows, as pulse.sequence, on the same time axis.
//          "sequence"     : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
            // The name of the driven bath spins. Required when count > 0.
//          "defect"       : "P1",
            // MHz (unit): the pump tones. A spin is driven when a tone is within bandwidth of its detuning.
//          "energy_shift" : [0.0],
            // MHz (unit): the half width of that window. Default 10 rad kHz (a provisional value).
//          "bandwidth"    : 0.0015915494309189536,
            // The two levels m of the driven spins. Default [0.5, -0.5].
//          "levels"       : [0.5, -0.5]
//      }
    }

    // ---- the Hamiltonian --------------------------------------------------------------------------------------------
    // docs: InputFile_General.md#tag_hamiltonian. Tutorials: advanced/06_hBN, basic/05_cCCE_vs_gCCE.
    // The last active key is pulse above: a key uncommented below needs a ',' after its '}'.
//  "hamiltonian"      : {
        // Read by every method.
//      "common"     : {
            // gCCE, megCCE with s_range: true adds the static (Overhauser) field of the bath outside each cluster.
            // Default false.
//          "qubit_bath_mean_field" : false,
            // With s_range: the static field of the other bath spins on a cluster: z (default, its z part), xyz (the
            // whole vector) or off.
//          "bath_bath_mean_field"  : "z",
            // The qubit-bath hyperfine tensors. Default {"dipole": "point", "contact": false}.
//          "hyperfine"             : {
                // point (default, from the positions) or dft (from file).
//              "dipole"      : "point",
                // true adds the Fermi contact term of file. Default false.
//              "contact"     : false,
                // The DFT tensor file (example/tools/SpinTensorReader). Required with dft or contact.
//              "file"        : "hyperfine_file",
                // A spin without a row in the file: stop (default) or point_dipole.
//              "missing_row" : "stop",
                // With frame.rotation and dft (then required): file or qubit.
//              "frame"       : "file"
//          },
            // "off" (default) or this object: the quadrupole coupling of spins with I > 1/2 from DFT field gradients.
//          "quadrupole"            : {
                // plain: one file computed with the defect; strained, bubble: spins outside its range take the
                // tensors of the structure without the defect.
//              "structure"           : "plain",
                // The tensor file with the defect. Required.
//              "file"                : "efg_file",
                // The file without the defect. Required with strained and bubble.
//              "file_without_defect" : "efg_file_pristine",
                // Angstrom (unit): the lattice parameters. Required with strained and bubble.
//              "cell"                : [10.0, 10.0, 10.0],
                // With frame.rotation (then required): file or qubit.
//              "frame"               : "file"
//          }
//      },
        // cCCE terms (meCCE reads them too).
//      "cCCE"       : {
            // true adds the Knight shift (a field across the qubit axis). Default false.
//          "knight"  : false,
            // true adds the hyperfine-mediated coupling of bath spins through the qubit. Default false.
//          "hfmedi"  : false,
            // true adds the second-order dA term of each bath spin (strong hyperfine; not with hfmedi). Default false.
//          "delta_A" : false
//      },
        // No settings yet.
//      "gCCE"       : {},
        // method noise.
//      "noise"      : {
            // bath (default): the clusters' eigenmodes. Or a model spectrum without a bath (order, bfield, rbath, rdip
            // and Qubit may then be left out):
            //   {"model": "OU", "sigma": 10.0, "tau_c": 0.1}: sigma kHz (amplitude), tau_c ms (correlation time);
            //   {"model": "power", "amplitude": 1.0, "beta": 1.0, "f_min": 0.1, "f_max": 1000.0}:
            //   S(f) = amplitude (f / 1 kHz)^-beta kHz^2/kHz on f_min .. f_max kHz.
//          "spectrum"       : "bath",
            // false (default): the noise of a bath spin is A_zz I_z; true adds A_zx I_x + A_zy I_y.
//          "pseudo_secular" : false,
            // coherence (default); spectrum: also C(t) and S(f) by FFT (files _noise_correlation, _noise_spectrum);
            // correlation: C(t) alone at the run's times.
//          "compute"        : "coherence",
            // compute spectrum only; every key has a default.
//          "fft"            : {
                // ms (unit): the end of the C(t) grid. Default the run's last time.
//              "time_max"       : 3.0,
                // Intervals of that grid. Default from the bath's fastest modes.
//              "nstep"          : 15000,
                // none (default), hann, tukey or exp, applied to C(t).
//              "window"         : "none",
                // FFT length pad_factor x (2 nstep + 1). Default 1.
//              "pad_factor"     : 1,
                // true subtracts the tail mean of C(t) first. Default false.
//              "subtract_mean"  : false,
                // true sets S(f) < 0 to 0. Default false.
//              "clamp_negative" : false
//          }
//      },
        // method pseudospin.
//      "pseudospin" : {
            // grouped (default): the states joined by flip-flop channels evolve as one block; pairwise: every channel
            // alone as a two-level system.
//          "channels" : "grouped"
//      },
        // Or "megCCE": the Lindblad jump operators (docs/methods/master-equation.md).
//      "meCCE"      : {
            // A jump on every spin named bath_name (megCCE: also "qubit"): operator - + z x y, rate in MHz (unit).
            // - and + relax a spin, z dephases it. Default: none.
//          "jump_operators" : [
//              {"bath_name": "13C", "operator": "z", "rate": 0.001}
//          ]
//      }
//  },

    // ---- the coordinate frame ---------------------------------------------------------------------------------------
    // docs: InputFile_General.md#tag_frame_rotation. Tutorial: advanced/09_coordinate_rotation.
//  "frame"            : {
        // For a bath file in the crystal frame with the qubit along another axis: qubits and bath are turned so that
        // qubit_axis is the new z. On when the block is given; then every directional input states its "frame".
//      "rotation" : {
            // Default true; false turns a written block off.
//          "enabled"    : true,
            // The qubit axis in the bath file's frame. Default [0, 0, 1].
//          "qubit_axis" : [1, 1, 1]
//      }
//  },

    // ---- paramagnetic defects in the bath ---------------------------------------------------------------------------
    // docs: InputFile_Defect.md. Tutorial: basic/06_defect_bath.
//  "Defect"           : {
        // One entry per species. A key in the list form [[axis state, ...], ...] holds a value per axis state.
//      "species" : [
//          {
                // The species name in the bath file.
//              "name"               : "P1",
                // Nuclear spins on the defect site (the 14N of a P1), each in every cluster of its electron. Default 0.
//              "n_onsite_spins"     : 1,
                // The indexed states of the defect (e.g. the 4 Jahn-Teller axes of a P1). Default 0.
//              "n_axis_states"      : 4,
                // With frame.rotation, for directional values (then required): file or qubit.
//              "frame"              : "file",
                // Isotope, spin and gyromagnetic ratio of each on-site spin. Gyros: rad kHz/G, or {values, unit} with
                // radkHz/G, rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T, radMHz/T, rad/us/T or rad/s/T.
//              "types"              : ["14N"],
//              "spins"              : [1.0],
//              "gyros"              : {"values": [3.0777], "unit": "MHz/T"},
                // fm2, or {values, unit} with fm^2 fm2 m^2 m2 cm^2 cm2 barn b mbarn millibarn mb. Default 0.
//              "quadrupole_moments" : {"values": [20.44], "unit": "mbarn"},
                // [axis state, type, xyz]: each on-site spin's position from the electron (Angstrom). Default 0.
//              "nuclear_positions"  : [[1, "14N", [0.0, 0.0, 1.0]]],
                // Per axis state, the electron-nucleus tensor: MHz, or {values, unit} with Hz kHz MHz GHz.
//              "hyperfine"          : {"values": [[1, "14N", [81.3, 0, 0, 0, 81.3, 0, 0, 0, 114.0]]], "unit": "MHz"},
                // Per axis state, the field gradient: au, or {values, unit} with Hartree/Bohr^2 Ha/Bohr^2 Ha/a0^2
                // a.u. au atomic_unit V/m^2 V/A^2.
//              "efg"                : {"values": [[1, "14N", [0, 0, 0, 0, 0, 0, 0, 0, 0]]], "unit": "Hartree/Bohr^2"},
                // The electron's zero-field splitting in MHz (spin > 1/2): {D, axis, unit} (D along the axis, no E;
                // write the tensor for an E other than 0), {tensor, unit} or [[axis state, "e", 9 numbers]]. Default 0.
//              "zfs"                : {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"},
                // MHz per axis state, or {values, unit}: a static shift of the electron (for a P1, the 14N hyperfine as
                // A m_I). Default 0.
//              "detuning"           : [[1, "e", 84.9]]
//          }
//      ]
//  },

    // ---- the clusters -----------------------------------------------------------------------------------------------
    // docs: InputFile_Cluster.md. Tutorials: advanced/11_cluster_shapes, basic/05_cCCE_vs_gCCE (dsj).
//  "cluster"          : {
        // standard (default): clusters of up to order bath spins. pcce: of partitions. bundle: of bundles of strongly
        // coupled spins. dsj, itb, dsjitb (cCCE): disjoint and inter-bath variants.
//      "mode"      : "standard",
        // Which clusters are kept.
//      "selection" : {
            // Order 3 and up: true (default) lists every sub-cluster; false leaves them out of the list (same result).
//          "add_subclusters" : true,
            // order + 1 entries: the most clusters kept per order (the strongest). Default 0, all.
//          "max_per_order"   : [0, 0, 0]
//      },
        // mode pcce: the bath cut into partitions by k-means.
//      "pcce"      : {
            // Spins per partition. Required.
//          "partition_size" : 2,
            // k-means trials (the best is kept) and iterations per trial. Default 30000 each.
//          "trials"         : 30000,
//          "max_iterations" : 30000,
            // k-means++ starting centres. Default true.
//          "kmeans_pp"      : true,
            // The number of trials printed. Default 0.
//          "print_trials"   : 0,
            // deterministic (default, follows seed) or legacy.
//          "partition"      : "deterministic",
            // When two partitions are connected: centers (default) or spins (closer to CCE, costs more).
//          "connections"    : "centers",
            // Reuse a stored partition: mode off (default), read or readwrite. With k_range the file needs {k}.
//          "cache"          : {"mode": "off", "file": "partition_cache.json"}
//      },
        // mode bundle: strongly coupled spins tied into bundles.
//      "bundle"    : {
            // Largest Hilbert dimension of a bundle. Required.
//          "max_dim"           : 64,
            // rad/ms (unit): the smallest coupling for which two groups merge. Default 0.
//          "merge_cutoff"      : 0,
            // distance (pairs within rdip, the default with rdip), auto, or checked (no rdip; the default without).
//          "connections"       : "distance",
            // With auto: the allowed score of the pairs left out. Default none.
//          "connection_budget" : 1000,
            // Limits of the expansion. Defaults 200000 clusters and 20000000 work.
//          "max_clusters"      : 200000,
//          "max_work"          : 20000000
//      },
        // The adaptive refinement (mode pcce or bundle): rerun, measure the error, change the bundles, until atol and
        // rtol hold or max_rounds.
//      "adaptive"  : {
            // Folder of the steps and summary.json. Required.
//          "workdir"                : "./adaptive",
            // Node memory of the probes at once. Default 0, automatic.
//          "probe_memory_gb"        : 0,
            // Absolute and relative tolerance of the coherence. Required.
//          "atol"                   : 1e-3,
//          "rtol"                   : 0,
            // Limits: rounds (default 8), patches per round (64), probes (4096), candidates per round (16).
//          "max_rounds"             : 8,
//          "max_patches"            : 64,
//          "max_probes"             : 4096,
//          "max_candidates"         : 16,
            // Neighbouring bundles per patch. Default 8.
//          "neighbors"              : 8,
            // MB for cached probe curves (default 256), seconds per probe (600), a file of probe results kept across
            // runs (default none), clusters per probe batch (default 0, from cache_mb).
//          "cache_mb"               : 256,
//          "timeout"                : 600,
//          "raw_cache"              : "raw_cache",
//          "probe_batch_clusters"   : 0,
            // The order of the merges: gain-per-cost (default) or worst-error.
//          "merge_priority"         : "gain-per-cost",
            // The work limit of the run (D^3 x states x times). Default 8e12.
//          "max_total_work"         : 8e12,
            // Switches, default false: prefer candidates with |L| <= 1; add candidates around a hub bundle (needs
            // prefer_physical); regroup where |L| > 1 instead of merging (cCCE order 2).
//          "prefer_physical"        : false,
//          "hub_candidates"         : false,
//          "regroup"                : false,
            // With regroup: regions per round (64), candidates per region (256), and how factors and whole curves are
            // obtained: regroup_factors probes (default) or factors, regroup_whole probes (default) or except.
//          "regroup_max_regions"    : 64,
//          "regroup_max_candidates" : 256,
//          "regroup_factors"        : "probes",
//          "regroup_whole"          : "probes"
//      }
//  },

    // ---- where and how it runs --------------------------------------------------------------------------------------
    // docs: InputFile_Backend.md. Tutorial: advanced/04_GPU.
//  "execution"        : {
        // cpu (default); cuda (main_gpu.out: every cluster on the GPU); auto (main_gpu.out: the clusters of
        // gpu_options.min_dim and above on the GPU).
//      "backend"          : "cpu",
        // main_gpu.out only. For cCCE with auto or cuda the keys left out are chosen at the start (a [gpu-auto] line
        // in the log): ranks_per_gpu 4 when a node has more than 4 ranks a GPU, min_dim 256 with it,
        // evolution_batch_size 32, cluster_batch_mb two thirds of the card. The defaults below are the other runs'.
//      "gpu_options"      : {
            // auto: the smallest cluster dimension sent to the GPU (positive). Default 768.
//          "min_dim"               : 768,
            // The packed kernel for small cCCE clusters: 0 (off), 2 or 4 (default in main_gpu.out).
//          "packed_max_dim"        : 4,
            // The warp kernel: 0 (off) to 32, default 27 in main_gpu.out.
//          "warp_max_dim"          : 27,
            // Clusters of one dimension evolved together (default 8) and time points per evolution (default 64).
//          "evolution_batch_size"  : 8,
//          "evolution_batch_steps" : 64,
            // MiB per rank: a GPU batch (default 64) and the reusable GPU memory (nonnegative, default 64; 0 disables).
//          "cluster_batch_mb"      : 64,
//          "cache_mb"              : 64,
            // auto and cCCE: only this many ranks per GPU open it; the others stay on the CPU. Default 0, every rank.
//          "ranks_per_gpu"         : 0
//      },
        // eigen (default), sparse (Krylov, for clusters too large to diagonalize), expm (gCCE), or
        // {"method": "sparse", "min_dim": 4096, "tol": 1e-12}.
//      "propagator"       : "eigen",
        // gCCE: vector (default where it can run, one bath state) or matrix (the density matrix).
//      "evolution"        : "vector",
        // How the phases are made: exact or recursive (default where it exists; the same to about 1e-10).
//      "evolution_update" : "recursive",
        // gCCE: correction sweeps on the eigenvectors. Off by default; a written block turns it on.
//      "gcce_refinement"  : {"max_sweeps": 32},
        // One-state cCCE: large clusters evolve several time points at once. On by default with these values;
        // {"enabled": false} turns it off.
//      "column_evolution" : {"min_dim": 256, "columns": 64},
        // Per-dimension counts and timing in the log. Default false.
//      "profile"          : false
//  },

    // ---- output -----------------------------------------------------------------------------------------------------
    // docs: InputFile_Output.md. Tutorial: advanced/11_cluster_shapes (cCCE2_Hahn_clusterReport).
//  "output"           : {
        // Folder of the outputs, made when missing; -o DIR takes its place. Default: the working directory.
//      "path"                  : "results",
        // Output name. Default: the input file name without .jsonc or .json.
//      "file"                  : "ccex_example",
        // Or a string. coherence (default): the result; coherence_tagged: each row tagged; cluster_coherence: also
        // each cluster's curve; cluster_density_matrix: also each cluster's density matrix (gCCE). nodiv true
        // (default) also writes _noDiv (the product without the division, for error correction).
//      "savemode"              : {"mode": "coherence", "nodiv": true},
        // A table of what each cluster contributes (cCCE, gCCE coherence). Default off; max_mb per rank 256.
//      "cluster_contributions" : {"enabled": true, "max_mb": 256},
        // true: write <file>.clusters.json and stop before the dynamics. Default false.
//      "cluster_only"          : false
//  },

    // ---- automatic convergence --------------------------------------------------------------------------------------
    // docs: InputFile_AutoConvergence.md. Tutorial: advanced/01_auto_convergence.
    // With this block main.out -f finds the time window, rdip, rbath and the number of configurations, then runs the
    // production with the mpirun -n N it was started with. A key written twice below is a choice: one line finds the
    // value, the other fixes it.
//  "auto_convergence" : {
        // Default true; false runs the ordinary calculation.
//      "enabled"     : true,
        // continue (default, reuse the workdir) or fresh.
//      "resume"      : "continue",
        // Folder of the runs, decisions and results. Required.
//      "workdir"     : "auto",
        // Found: T_d where the T2 fit reaches decay_level; the production window is production_factor x T_d.
//      "time_window" : {"decay_level": 0.01, "revival_check_factor": 2.5, "production_factor": 1.5},
        // Or fixed, in ms (unit).
//      "time_window" : {"value": {"value": 1.0, "unit": "ms"}},
        // Found among candidates (Angstrom): two agree within threshold (relative difference of the fitted area).
//      "rdip"        : {"candidates": [4, 6, 8], "threshold": 0.02},
        // Or num_candidates found from the bath (default 7).
//      "rdip"        : {"num_candidates": 7, "threshold": 0.02},
        // Or fixed.
//      "rdip"        : {"value": 8},
        // The same for rbath (num_candidates default 4).
//      "rbath"       : {"candidates": [10, 15], "threshold": 0.02},
//      "rbath"       : {"num_candidates": 4, "threshold": 0.02},
//      "rbath"       : {"value": 15},
        // points peaks (default) or all; floor: fit before |L| falls below it; n_range: the allowed exponent.
//      "t2_fit"      : {"points": "peaks", "floor": 0.05, "n_range": [0.3, 4.0]},
        // How many configurations (and, with s_range, states): found from a test, or fixed.
//      "sampling"    : {
            // Configurations of the test (default 10) and the SE(A)/A target (default 0.05).
//          "test_configurations" : 10,
//          "threshold"           : 0.05,
            // With s_range: states per test configuration (5), their SE target and the most states per configuration.
//          "test_states"         : 5,
//          "states_threshold"    : 0.05,
//          "max_states"          : 30,
            // The fewest and the most production configurations.
//          "min_configurations"  : 10,
//          "max_configurations"  : 50
//      },
        // Or the configuration ids of the test and of the production given.
//      "sampling"    : {"test": [1, 2], "production": [3, 4], "threshold": 0.05},
        // Or fixed: configurations, and with s_range states per configuration.
//      "sampling"    : {"configurations": 10, "states": 5}
//  }
}
)CCEX";
