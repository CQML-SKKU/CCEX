/**
 * @file example_input.cpp
 * @brief The example input (ccex_example.jsonc) that main.out -e writes.
 */
// The input template `main.out -e` writes, the one place it is kept. The active part is the calculation of
// example/tutorial/00_quickstart (the NV-13C cCCE2 Hahn echo, ensemble, rbath 50 A, 2 ms) without its files: gyrofile is
// commented (the built-in isotope table) and bathfile is a name the user replaces, so the template does not run as it
// is. Every other input key the readers accept (src/input/option.cpp, the adaptive options table, the automatic-
// convergence runner) is a commented line with its allowed values and default. Comments are // (vim and VS Code show
// a .jsonc file's // as comments), every trailing comment starts at column 41, and the ':' of one object's keys share a
// column ("key" : value). Each section starts with "Details:" (its docs/user page and the same page of the web site)
// and "Example:" (the tutorials that use its keys). test/code_verification checks it: with a real bath it runs, no accepted key is missing from
// it, the comment column and the ':' column.
#include "../../include/input/example_input.h"

const char* const kExampleInputName = "ccex_example.jsonc";

const char* const kExampleInput = R"CCEX(// ccex_example.jsonc - an input template written by `main.out -e` (CCEX)
//
// This is a template: it does not run until "bathfile" below names a real bath file. A ready bath for the qubit
// written here is example/tutorial/00_quickstart/bath/bath_Diamond_NV_natab_1. With a bath file in place, run
//
//     mpirun -n 4 main.out -f ccex_example.jsonc -o results
//
// The uncommented part is the calculation of example/tutorial/00_quickstart: the Hahn echo of an NV centre (spin 1) in
// diamond with 13C, conventional CCE (cCCE) to order 2, the bath in its unpolarised state (the ensemble: no
// configuration.s_range), 13C within 50 Angstrom, 2 ms. It writes ccex_example_noDiv and ccex_example_wiDiv
// (output.file defaults to the input file name without .jsonc or .json, with -o DIR or output.path they go to that
// folder). Every other key the input accepts is a commented line with its allowed values and its default: remove the
// leading // to use one. A //, # or ! outside a string starts a comment to the end of the line. Keys are read without
// regard to case. A key no part of the run reads is reported, as is a key this run reads and does not use, and a likely
// typo of a key stops the run. The keys are described in docs/user/InputFile.md and the pages it links. Each section
// below starts with "Details:", the page that describes its keys (in docs/user, and the same page on the web site), and
// "Example:", the tutorials under example/tutorial that use them.
//
// Units: a number is in the unit its key names below. The keys marked (unit) also take {"value": ..., "unit": "..."}
// (the value a number or a list, as the plain form): time s, ms, us, ns. Field G, mT, T. Length A (Angstrom), nm.
// Frequency Hz, kHz, MHz, GHz. Angular frequency rad/s, rad/ms, rad/us. Units are written with their case (mHz is not
// MHz). "deltat": {"value": 5, "unit": "us"} is "deltat": 0.005. Mind the two frequency defaults: a plain number is kHz
// for the qubit (Qubit.qubit[].detuning, Qubit.intmap) and MHz for a Defect species (its detuning, zfs, hyperfine).
{

    // ---- method and quantity ----------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_General.md#tag_method
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__General.html#tag_method)
    // Example: example/tutorial/06_cCCE_vs_gCCE (cCCE, gCCE, meCCE), example/tutorial/03_P1 (dm in
    //          gCCE2_Ramsey_population)
    "method"           : "cCCE",        // cCCE (default), gCCE, meCCE, megCCE, noise (Gaussian phase noise) cCCE:
                                        // the qubit only sets the field each bath spin feels (pure dephasing). gCCE:
                                        // the qubit is in every cluster's Hamiltonian (qubit-bath flip-flops, any qubit
                                        // spin and field). meCCE, megCCE: cCCE, gCCE with Lindblad jump operators
    "quantity"         : "coherence",   // coherence (default): the qubit coherence L(t) between alpha and beta. dm
                                        // (gCCE): the qubit density matrix at every time, for populations

    // ---- the qubit --------------------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_Qubit.md#tag_Qubit
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Qubit.html#tag_Qubit)
    // Example: example/tutorial/00_quickstart (a qubit from a defect file), example/tutorial/05_multiNV (several NV),
    //          example/tutorial/04_hBN (beta_state in gCCE2_Hahn_states)
    "Qubit"            : {
        "qubit"       : [               // one entry per qubit
            {
                "name"     : "NV",      // default "q0", "q1", ... by position in the list
                "spin"     : 1.0,       // default 1.0
//              "gyro"     : -17608.597050,
                                        // rad kHz / G, or {"value": ..., "unit": ...} with the units of Defect gyros
                                        // (kHz/G, MHz/T, ...), default -17608.597050 (the electron)
//              "detuning" : 0.0,       // kHz (unit; a Defect detuning is MHz), default 0.0: adds detuning x S_z to
                                        // this qubit's Hamiltonian, which shifts its levels: the phase of L turns
                                        // faster (not |L|), and in gCCE its resonances with bath spins move. A plain
                                        // number other than 0 is noted in the log as read in kHz
//              "file"     : "defect_position",
                                        // in place of "xyz": a file whose first line holds the position x y z
                                        // (Angstrom), e.g. a defect file of example/tools/BathGenerator
                                        // (00_quickstart's bath/bath_Diamond_NV_natab_defect)
                "xyz"      : [73.73895, 75.72808, 73.97033]
                                        // position, Angstrom (unit), in the frame of the bath file (required, or
                                        // "file")
            }
        ],
//      "intmap"      : [               // qubit-qubit tensors, default zero self terms, point dipole between qubits
//          {
//              "between" : ["NV", "NV"],
                                        // two qubit names, the same name twice is a self term (zero-field splitting)
//              "tensor"  : [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]],
                                        // 3x3 (or 9 numbers), here the NV's D = 2.88 GHz as diag(-D/3, -D/3, 2D/3) kHz
//              "unit"    : "kHz",      // Hz, kHz (default: a tensor without unit is kHz, a Defect zfs MHz; the log
                                        // notes it), MHz, GHz
//              "frame"   : "qubit"     // with frame.rotation (required then): file (written in the bath file's frame,
                                        // turned with the positions) or qubit (written with qubit_axis as z, used as
                                        // written)
//          }
//      ],
        "alpha_state" : {"ms": [1]},    // L(t) is the coherence between |alpha> and |beta>: here ms +1 and 0 of the NV.
                                        // One of {"ms": [m_1, m_2, ...]} (one m_s per qubit in the order of "qubit"),
                                        // {"index": n} (the n-th eigenstate of the qubit Hamiltonian by descending
                                        // energy, the order changes with the field, the log names the ms it is; its
                                        // largest component is made real and positive) or {"vector": [...]} (each
                                        // component a number or [re, im]); default {"ms": [1, 1, ...]}. pulse.per_qubit
                                        // needs {"ms"}
        "beta_state"  : {"ms": [0]}     // the same for beta, default {"ms": [0, 0, ...]}
//      "beta_state"  : {"vector": [0, 1, 0]},
                                        // (ms 0 of a spin 1 as a vector, the components ordered ms +1, 0, -1)
//      "psi0"        : [0.7071, 0.7071, 0.0]
                                        // gCCE: the initial qubit state, used as written (give it normalized), each
                                        // component a number or [re, im], default (|alpha> + |beta>) normalized. cCCE
                                        // and meCCE use alpha and beta only (a psi0 there is reported as not used)
    },

    // ---- the bath configuration -------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_General.md#tag_configuration
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__General.html#tag_configuration)
    // Example: example/tutorial/01_ensemble_and_states (k_range, s_range, seed, the built-in isotope table, the state
    //          files of example/tools/StateFiles), example/tutorial/05_multiNV (bath_state_file in
    //          gCCE2_HahnRamsey_NV2)
    "configuration"    : {
//      "path"              : "data",   // folder of the file names of this block (relative to the working directory, or
                                        // absolute, used as it is), default "."
//      "k_range"           : [1, 10],  // several bath configurations in one run, one after another: each k from first
                                        // to last replaces {k} in bathfile, bath_state_file, defect_axis_file,
                                        // defect_state_file and pcce.cache.file and writes <output.file>_k<k>_... (no
                                        // average is written). Default: one configuration ({k} kept, {{k}} a literal
                                        // {k})
//      "gyrofile"          : "gyro.txt",
                                        // isotope, spin, gyromagnetic ratio per line, default: the built-in isotope
                                        // table (needed for spins that are not isotopes, e.g. P1)
        "bathfile"          : ["my_bath.txt"]
                                        // FILL IN: one or more bath files (merged), a count line, then x y z isotope
                                        // per line, in the frame of the qubit xyz above. A ready one for this NV:
                                        // example/tutorial/00_quickstart/bath/bath_Diamond_NV_natab_1. With k_range the
                                        // name holds {k}, e.g. ["config_{k}"] reads config_1 ... config_10 for [1, 10]
//      "bath_shift"        : [[0.0, 0.0, 0.0]],
                                        // Angstrom (unit), one [x, y, z] per bath file that moves its spins (e.g. the
                                        // baths of several qubits), default [0, 0, 0] each
//      "s_range"           : [1, 20],  // the sampled bath states [first, last] (a list also for one, [n, n]): each
                                        // state s a pure bath state, read from bath_state_file or drawn from seed, its
                                        // own output _state<s> (average them), with the static field of the other
                                        // spins. Absent (default): the ensemble, the average over the unpolarised bath.
                                        // A run of [11, 20] calculates the states 11 ... 20 of a run of [1, 20] (the
                                        // states before first are drawn and left out)
//      "bath_state_file"   : ["cfg{k}_state{s}"],
                                        // with s_range: the sampled bath states, one name per bath file ("" for a bath
                                        // file whose states are drawn from seed): {s} is the state, {k} the
                                        // configuration (a name without {s} has the state appended, "state_" reads
                                        // state_1, ...). A count line, then one m_s per bath spin. Every file of
                                        // s_range has to be there (a missing one stops the run).
                                        // example/tools/StateFiles writes them
//      "defect_axis_file"  : ["cfg{k}_axis"],
                                        // with Defect species: the axis state (1 .. Defect n_axis_states) of each
                                        // Defect spin, one file per bath file ("": drawn from seed), it has to be there
//      "defect_state_file" : ["cfg{k}_onsite{s}"],
                                        // with s_range and Defect on-site spins: their m_s, one line per bath spin (the
                                        // m of each on-site spin), {s} as bath_state_file, default: drawn from seed
//      "seed"              : 1         // random seed shared by every rank: the first configuration's, configuration k
                                        // of k_range takes seed + (k - first); or a list of one seed per configuration.
                                        // Default: generated on rank 0 and reported in the log when the run draws
                                        // random numbers. Give it to repeat a run
    },

    // ---- the cluster expansion and the time axis --------------------------------------------------------------------
    // Details: docs/user/InputFile_General.md#tag_order
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__General.html#tag_order)
    // Example: example/tutorial/00_quickstart, example/tutorial/03_P1 (order 4 in cCCE4_Hahn_states)
    "order"            : 2,             // the largest cluster: bath spins (standard), partitions (pcce) or bundles
                                        // (bundle), raised until the curve stops changing (required)
    "bfield"           : 500,           // Gauss (unit): Bz, or the list Bx, By, Bz (required, 0 for no field)
    "rbath"            : 50,            // Angstrom (unit): the bath spins this close to the qubit are kept, or the list
                                        // inner, outer (inner default 0) (required)
    "rdip"             : 10,            // Angstrom (unit): two bath spins this close are connected and can be in one
                                        // cluster, or the list inner, outer (inner default 0) (required, except with
                                        // cluster.mode bundle)
    "deltat"           : 0.001,         // ms (unit): the time step, the outputs are at t = 0, deltat, ..., (nstep - 1)
                                        // deltat, t the whole sequence (required)
    "nstep"            : 2000,          // number of time points (required)

    // ---- the pulse sequence -----------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_Pulse.md#tag_pulse
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Pulse.html#tag_pulse)
    // Example: example/tutorial/02_pulses (count, name, sequence), example/tutorial/03_P1 (bath in gCCE2_Hahn_DEER)
    "pulse"            : {
        "count"           : 1           // number of pi pulses: 0 Ramsey, 1 HahnEcho, > 1 CPMG (required)
//      "name"            : "Manual",   // Ramsey (count 0), HahnEcho (count 1), CPMG (count 2 or more), WAHUHA (count
                                        // 4) or Manual (any count, with sequence); a name that does not fit the count
                                        // stops the run. Default: from count
//      "sequence"        : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
                                        // Manual: count + 1 rows [fraction, axis X|Y|Z|I, angle in degrees], the free
                                        // time before each pulse as a fraction of t (the fractions sum to 1), the last
                                        // row the time after the last pulse. Pulses other than pi flips need gCCE
//      "per_qubit"       : false,      // false (default): a pulse rotates the qubits together, true: qubit by qubit
                                        // between the m_s of alpha_state and beta_state {"ms"} (with gCCE pulses it
                                        // stops, a known defect: keep false)
//      "duration"        : 0.0,        // ns (unit), default 0.0: gCCE only, read with detuning_factor other than 1
                                        // (other methods stop on a value other than 0)
//      "detuning_factor" : 1.0,        // gCCE, default 1.0: ideal instantaneous pulses. Another value: a finite pulse
                                        // of duration, driven at detuning_factor x the qubit frequency (other methods
                                        // stop on a value other than 1)
//      "bath"            : {           // DEER: pulses on bath spins of one species during the qubit's sequence (gCCE,
                                        // megCCE, CPU, docs/methods/deer.md)
//          "count"        : 1,         // bath pulses, default 0 (none)
//          "sequence"     : [["1/2", "X", 180.0], ["1/2", "I", 0.0]],
                                        // count + 1 rows, as pulse.sequence, on the same time axis
//          "defect"       : "P1",      // the name of the bath spins the pump drives (required when count > 0)
//          "energy_shift" : [0.0],     // MHz (unit), the pump tone(s): a spin is driven when a tone is within
                                        // bandwidth of its Defect detuning (required when count > 0)
//          "bandwidth"    : 0.0015915494309189536,
                                        // MHz (unit) > 0: the half width of that window, default 10/(2 pi 1000) MHz =
                                        // 10 rad kHz, a provisional value: set it to the pump's excitation width
//          "levels"       : [0.5, -0.5]
                                        // the two levels m of the driven spins (their spin S is their own) the bath
                                        // pulse rotates between, default [0.5, -0.5]
//      }
    }                                   // the last active key: a key uncommented below needs a ',' after this '}'

    // ---- the Hamiltonian --------------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_General.md#tag_hamiltonian
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__General.html#tag_hamiltonian)
    // Example: example/tutorial/04_hBN (hyperfine, quadrupole, hfmedi, qubit_bath_mean_field),
    //          example/tutorial/06_cCCE_vs_gCCE (knight, meCCE jump_operators), example/tutorial/03_P1
    //          (bath_bath_mean_field in cCCE2_Hahn_meanField)
//  "hamiltonian"      : {
//      "common" : {                    // read by every method
//          "qubit_bath_mean_field" : false,
                                        // gCCE, megCCE with s_range: true adds the static (Overhauser) field of the
                                        // sampled bath outside each cluster to the qubit, default false (the qubit
                                        // feels the spins of its cluster only). Not used by cCCE and meCCE, nor in the
                                        // ensemble (reported)
//          "bath_bath_mean_field"  : "z",
                                        // with s_range: the static field of the bath spins outside a cluster on its
                                        // spins. z (default): its z part, xyz: the whole vector (where z alone
                                        // diverges, e.g. clustered P1), off: left out (each cluster alone in its
                                        // sampled state). Not used in the ensemble (reported)
//          "hyperfine"             : { // default {"dipole": "point", "contact": false}
//              "dipole"      : "point",
                                        // the qubit-bath hyperfine tensors: point (default) from the positions (point
                                        // dipole) or dft (the DFT tensors of file)
//              "contact"     : false,  // true adds the Fermi contact term of file, default false
//              "file"        : "hyperfine_file",
                                        // the DFT tensor file (example/tools/SpinTensorReader), required with dipole
                                        // dft or contact true
//              "missing_row" : "stop", // a spin inside the file's range without a row: stop (default) or point_dipole
                                        // (it gets the point-dipole tensor)
//              "frame"       : "file"  // with frame.rotation and dipole dft (required then): the frame of the file's
                                        // tensor components, file (the bath file's, turned with the positions) or qubit
                                        // (as written)
//          },
//          "quadrupole"            : { // "off" (default) or this object: the quadrupole coupling (electric field
                                        // gradient) of the nuclear spins with I > 1/2 from a DFT file
//              "structure"           : "plain",
                                        // plain: one file computed with the defect. strained, bubble: the spins outside
                                        // its range take the tensors of the structure without the defect (folded by
                                        // cell)
//              "file"                : "efg_file",
                                        // the tensor file with the defect, required with every structure
//              "file_without_defect" : "efg_file_pristine",
                                        // the file without the defect, required with strained and bubble
//              "cell"                : [10.0, 10.0, 10.0],
                                        // Angstrom (unit): the lattice parameters, required with strained and bubble
//              "frame"               : "file"
                                        // with frame.rotation (required then): the frame of the files' tensor
                                        // components, file or qubit, as for hyperfine
//          }
//      },
//      "cCCE"   : {                    // cCCE terms (meCCE reads them too)
//          "knight"  : false,          // true adds the Knight shift: in a field with a part transverse to the qubit
                                        // axis the qubit shifts the bath spins' precession, default false
//          "hfmedi"  : false,          // true adds the hyperfine-mediated coupling of bath spins (second order through
                                        // the qubit), default false
//          "delta_A" : false           // true adds the second-order dA term of each bath spin (Maze 2008; strong
                                        // hyperfine near the defect), not with hfmedi, default false
//      },
//      "gCCE"   : {},                  // no settings yet
//      "noise"  : {                    // method noise: where the phase noise comes from
//          "spectrum"       : "bath",  // "bath" (default: the clusters), or a model,
//                                      // {"model": "OU", "sigma": kHz, "tau_c": ms} or
//                                      // {"model": "power", "amplitude": kHz^2/kHz,
//                                      // "beta", "f_min": kHz, "f_max": kHz}
//          "pseudo_secular" : false,   // the noise: A_zz I_z (secular, default), true adds
//                                      // A_zx I_x + A_zy I_y (pseudo-secular)
//          "compute"        : "coherence",
//                                      // "coherence" (default), "spectrum" (C(t), FFT
//                                      // S(f), the coherence) or "correlation" (C(t))
//          "fft"            : {        // compute "spectrum" only
//            "time_max"       : 3.0,   // ms (unit): C(t) on 0 .. time_max
//            "nstep"          : 15000, // intervals
//            "window"         : "none",// "none", "hann", "tukey" or "exp"
//            "pad_factor"     : 1,     // FFT length pad_factor x (2 nstep + 1)
//            "subtract_mean"  : false, // subtract the tail mean of C(t)
//            "clamp_negative" : false  // S(f) < 0 set to 0
//          }
//      },
//      "meCCE"  : {                    // or "megCCE": master-equation dissipators (docs/methods/master-equation.md)
//          "jump_operators" : [
//              {"bath_name": "13C", "operator": "z", "rate": 0.001}
                                        // a jump sqrt(gamma) S_op on every spin named bath_name (megCCE: also "qubit"),
                                        // operator - + z x y, rate MHz (unit) >= 0, gamma = rate x 2 pi 10^3 rad/ms. -
                                        // and + relax a spin (T1), z dephases it. Default: none (the cCCE or gCCE)
//          ]
//      }
//  },

    // ---- the coordinate frame ---------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_General.md#tag_frame_rotation
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__General.html#tag_frame_rotation)
    // Example: example/tutorial/10_coordinate_rotation
//  "frame"            : {
//      "rotation" : {                  // for a bath file written in the crystal frame with the qubit along another
                                        // axis (an NV along 111): the qubits and the bath are turned about the first
                                        // qubit so that qubit_axis is the new z. On when the block is given, absent =
                                        // no rotation. Every tensor, axis or Defect species with directional numbers
                                        // written in the input then states its "frame": file or qubit
//          "enabled"    : true,        // default true, false turns a written block off
//          "qubit_axis" : [1, 1, 1]    // the qubit's axis in the bath file's frame, the new z (bfield [0, 0, B] is
                                        // then along it in the lab, a transverse part is refused), default [0, 0, 1]
//      }
//  },

    // ---- paramagnetic defects in the bath ---------------------------------------------------------------------------
    // Details: docs/user/InputFile_Defect.md#tag_Defect
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Defect.html#tag_Defect)
    // Example: example/tutorial/03_P1 (P1 by its detunings, and with its 14N on site in cCCE2_Hahn_onsite14N)
//  "Defect"           : {
//      "species" : [                   // one entry per species: a bath spin with its name is that defect, an electron
                                        // spin with its on-site nuclear spins and axis states. A key in the list form
                                        // [[axis state, ...], ...] holds a value per axis state, the object form {value
                                        // or values, unit} one value for every axis state (zfs: D or tensor)
//          {
//              "name"               : "P1",
                                        // the species name in the bath file
//              "n_onsite_spins"     : 1,
                                        // nuclear spins on the defect site (the 14N of a P1): each joins the electron
                                        // in every cluster it is in, default 0 (none)
//              "n_axis_states"      : 4,
                                        // the indexed states of the defect (its axes, e.g. the 4 Jahn-Teller axes of a
                                        // P1, or frozen 14N states), one per spin from defect_axis_file or seed,
                                        // default 0
//              "frame"              : "file",
                                        // with frame.rotation, for a species with nuclear_positions, hyperfine, efg or
                                        // zfs (required then): their frame, file (the bath file's, turned) or qubit (as
                                        // written)
//              "types"              : ["14N"],
                                        // isotope of each on-site nuclear spin
//              "spins"              : [1.0],
                                        // spin of each on-site nuclear spin
//              "gyros"              : {"values": [3.0777], "unit": "MHz/T"},
                                        // a list (radkHz/G) or {values (or value), unit}: radkHz/G, rad/ms/G, kHz/G,
                                        // MHz/G, MHz/T, Hz/T, radMHz/T, rad/us/T, rad/s/T
//              "quadrupole_moments" : {"values": [20.44], "unit": "mbarn"},
                                        // a list (fm2) or {values, unit}: fm^2 fm2 m^2 m2 cm^2 cm2 barn b mbarn
                                        // millibarn mb, default 0
//              "nuclear_positions"  : [[1, "14N", [0.0, 0.0, 1.0]]],
                                        // [axis state, type, xyz]: each on-site spin's position from the electron per
                                        // axis state (Angstrom), default [0, 0, 0]
//              "hyperfine"          : {"values": [[1, "14N", [81.3, 0, 0, 0, 81.3, 0, 0, 0, 114.0]]], "unit": "MHz"},
                                        // per axis state, the electron-nucleus hyperfine tensor: a list (MHz) or
                                        // {values, unit}: Hz kHz MHz GHz
//              "efg"                : {"values": [[1, "14N", [0, 0, 0, 0, 0, 0, 0, 0, 0]]], "unit": "Hartree/Bohr^2"},
                                        // per axis state, the electric field gradient at each on-site spin: a list (au)
                                        // or {values, unit}: Hartree/Bohr^2 Ha/Bohr^2 Ha/a0^2 a.u. au atomic_unit V/m^2
                                        // V/A^2 ...
//              "zfs"                : {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"},
                                        // the electron's zero-field splitting (a spin > 1/2), MHz (a qubit's
                                        // Qubit.intmap self term is kHz): {D, axis, unit} (D along the symmetry axis;
                                        // no E, write the tensor for an E other than 0), {tensor, unit} or the list
                                        // [[axis state, "e", 9 numbers]], default zero tensor
//              "detuning"           : [[1, "e", 84.9]]
                                        // MHz per axis state (a qubit's Qubit.qubit[].detuning is kHz): a static shift
                                        // of the electron (the 14N hyperfine of a P1 as A m_I), or {"values": [...],
                                        // "unit": "MHz"}, default 0
//          }
//      ]
//  },

    // ---- the clusters -----------------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_Cluster.md#tag_cluster_mode
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Cluster.html#tag_cluster_mode)
    // Example: example/tutorial/07_pCCE_and_bundle (pcce, bundle, adaptive), example/tutorial/06_cCCE_vs_gCCE (dsj in
    //          cCCE2_Hahn_dsj), example/tutorial/03_P1 (max_per_order in cCCE4_Hahn_states)
//  "cluster"          : {
//      "mode"      : "standard",       // standard (default): clusters of up to order bath spins connected within rdip.
                                        // pcce: of up to order partitions of partition_size spins (cCCE, gCCE, meCCE,
                                        // megCCE). bundle: of up to order bundles of strongly coupled spins (cCCE,
                                        // gCCE). dsj, itb, dsjitb (cCCE): disjoint, inter-bath
//      "selection" : {                 // which clusters are kept (every mode)
//          "add_subclusters" : true,   // order >= 3: true (default) lists every sub-cluster of the highest-order
                                        // clusters, false leaves them out of the list. The results are the same (the
                                        // missing sub-clusters are added when the powers are computed), false can be a
                                        // little faster
//          "max_per_order"   : [0, 0, 0]
                                        // order + 1 entries, the most clusters kept per order (the strongest), default
                                        // 0 = all. Only entries 3 and up are used (standard, pcce), a nonzero entry 0-2
                                        // is reported as not used
//      },
//      "pcce"      : {                 // mode pcce: the bath cut into partitions by k-means
//          "partition_size" : 2,       // spins per partition (required with mode pcce)
//          "trials"         : 30000,   // k-means trials, the best one kept, default 30000
//          "max_iterations" : 30000,   // k-means iterations per trial, default 30000
//          "kmeans_pp"      : true,    // k-means++ starting centres, default true
//          "print_trials"   : 0,       // the number of rank 0's trials printed, default 0 (none)
//          "partition"      : "deterministic",
                                        // deterministic (default, follows seed) or legacy (the original search)
//          "connections"    : "centers",
                                        // when two partitions are connected (within rdip): centers (default, their
                                        // centroids) or spins (any two of their spins, closer to CCE, costs more)
//          "cache"          : {"mode": "off", "file": "partition_cache.json"}
                                        // reuse a stored partition (cCCE): mode off (default), read (the file has to be
                                        // there) or readwrite, with configuration.k_range the file needs {k} (one cache
                                        // per configuration)
//      },
//      "bundle"    : {                 // mode bundle: strongly coupled spins tied into bundles
//          "max_dim"           : 64,   // largest Hilbert dimension of one bundle (required)
//          "merge_cutoff"      : 0,    // rad/ms (unit, rad/s rad/ms rad/us): the smallest coupling for which two
                                        // groups are merged, default 0
//          "connections"       : "distance",
                                        // distance (pairs within rdip; the default with rdip), auto (weak pairs left
                                        // out up to connection_budget) or checked (no rdip: the run checks the links
                                        // itself; the default without rdip). The log names the one used
//          "connection_budget" : 1000, // with auto: the allowed score of the left-out pairs, default: none
//          "max_clusters"      : 200000,
                                        // the most clusters the expansion may list, default 200000
//          "max_work"          : 20000000
                                        // the work limit of the expansion, default 20000000
//      },
//      "adaptive"  : {                 // the adaptive refinement (mode pcce or bundle): the calculation is repeated,
                                        // the error at the bundle boundaries measured and the bundles changed until
                                        // atol and rtol hold or max_rounds. input, binary, cwd, mpiexec and ranks are
                                        // decided by the run. Keys with _ (max_rounds; - is the command line's)
//          "workdir"                : "./adaptive",
                                        // folder of the steps and summary.json (required)
//          "probe_memory_gb"        : 0,
                                        // node memory of the probes running at once, default 0 (automatic)
//          "atol"                   : 1e-3,
                                        // absolute tolerance of the coherence (required)
//          "rtol"                   : 0,
                                        // relative tolerance of the coherence (required)
//          "max_rounds"             : 8,
                                        // the most rounds, default 8
//          "max_patches"            : 64,
                                        // the most patches (bundle neighbourhoods) checked per round, default 64
//          "max_probes"             : 4096,
                                        // the most probe calculations, default 4096
//          "max_candidates"         : 100,
                                        // the most merge candidates per round, default 16
//          "neighbors"              : 8,
                                        // neighbouring bundles a patch takes, default 8
//          "cache_mb"               : 256,
                                        // MB for the cached probe curves, default 256
//          "timeout"                : 600,
                                        // seconds per probe run, default 600
//          "raw_cache"              : "raw_cache",
                                        // file that keeps probe results across runs, default none
//          "probe_batch_clusters"   : 0,
                                        // clusters per probe batch, default 0 (from cache_mb)
//          "merge_priority"         : "gain-per-cost",
                                        // the order of the merges: gain-per-cost (default) or worst-error
//          "max_total_work"         : 8e12,
                                        // the work limit of the run (D^3 x states x times), default 8000000000000.0
//          "prefer_physical"        : false,
                                        // switch: prefer candidates whose curve keeps |L| <= 1, default false
//          "hub_candidates"         : false,
                                        // switch: add candidates around a hub bundle (needs prefer_physical), default
                                        // false
//          "regroup"                : false,
                                        // switch: regroup the bundles where |L| > 1 instead of merging (cCCE order 2),
                                        // default false
//          "regroup_max_regions"    : 64,
                                        // with regroup: regions checked per round, default 64
//          "regroup_max_candidates" : 256,
                                        // with regroup: candidates per region, default 256
//          "regroup_factors"        : "probes",
                                        // with regroup: probes (default) or factors
//          "regroup_whole"          : "probes"
                                        // with regroup: probes (default) or except
//      }
//  },

    // ---- where and how it is computed -------------------------------------------------------------------------------
    // Details: docs/user/InputFile_Backend.md#tag_backend
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Backend.html#tag_backend)
    // Example: example/tutorial/08_GPU (backend, gpu_options, profile)
//  "execution"        : {
//      "backend"          : "cpu",     // cpu (default), cuda (main_gpu.out: the clusters on the GPU) or auto
                                        // (main_gpu.out: the clusters of gpu_options.min_dim and above on the GPU)
//      "gpu_options"      : {          // main_gpu.out (main.out reports them as not used)
//          "min_dim"               : 768,
                                        // positive, auto: smallest cluster dimension sent to the GPU, default 768
//          "packed_max_dim"        : 4,
                                        // 0 (off), 2 or 4, default 4 in main_gpu.out with cuda/auto and cCCE, else 0
//          "warp_max_dim"          : 27,
                                        // 0 (off) to 32, default 27 in main_gpu.out with cuda/auto and cCCE/gCCE, else
                                        // 0
//          "evolution_batch_size"  : 8,
                                        // clusters of one dimension evolved together on the GPU, default 8
//          "evolution_batch_steps" : 64,
                                        // time points one GPU evolution takes together, default 64
//          "cluster_batch_mb"      : 64,
                                        // MiB per rank for the clusters of one GPU batch, default 64
//          "cache_mb"              : 64
                                        // nonnegative, MiB per rank of GPU memory kept for reuse, 0 disables, default
                                        // 64
//      },
//      "propagator"       : "eigen",   // eigen (default, one diagonalization per cluster), sparse (Krylov, for
                                        // clusters too large to diagonalize: one state, cCCE, gCCE vector), expm
                                        // (gCCE), or {"method": "sparse", "min_dim": 4096, "tol": 1e-12}
//      "evolution"        : "vector",  // gCCE: vector (default where it can run: one bath state) or matrix (the
                                        // density matrix, always works). The two agree to rounding
//      "evolution_update" : "recursive",
                                        // how the phases exp(-iE t) are made: exact (sin, cos at every point) or
                                        // recursive (every 16th point, products between), default recursive where it
                                        // exists, exact elsewhere. Faster, agrees to about 1e-10
//      "gcce_refinement"  : {"max_sweeps": 32},
                                        // gCCE: correction sweeps on the eigenvectors, off by default, writing the
                                        // block turns it on ({"enabled": false} keeps a written block off), max_sweeps
                                        // default 32
//      "column_evolution" : {"min_dim": 256, "columns": 64},
                                        // one-state cCCE: large clusters with columns time points at once (faster, the
                                        // same to 1e-10), on by default with these values, {"enabled": false} turns it
                                        // off
//      "profile"          : false      // per-dimension counts and timing, default false
//  },

    // ---- output -----------------------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_Output.md#tag_output_path
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__Output.html#tag_output_path)
    // Example: example/tutorial/07_pCCE_and_bundle (savemode, cluster_contributions and cluster_only in
    //          cCCE2_Hahn_clusterReport)
//  "output"           : {
//      "path"                  : "results",
                                        // the folder of the outputs, made when missing, -o DIR takes its place,
                                        // default: the working directory (not with auto_convergence: its runs go to
                                        // auto_convergence.workdir)
//      "file"                  : "ccex_example",
                                        // output name, default: the input file name without .jsonc or .json (an
                                        // absolute name is refused with -o DIR or path; not with auto_convergence)
//      "savemode"              : {"mode": "coherence", "nodiv": true},
                                        // or a string (the mode). coherence (default): the result. coherence_tagged:
                                        // the same, each row with the first bath spin's tensor and file name.
                                        // cluster_coherence: also each cluster's curve. cluster_density_matrix: also
                                        // (gCCE) each cluster's density matrix. nodiv true (default) writes _noDiv
                                        // beside _wiDiv: the product without the division, not a coherence (error
                                        // correction)
//      "cluster_contributions" : {"enabled": true, "max_mb": 256},
                                        // a table of what each cluster contributes, <file>.contributions.state-<s>.tsv
                                        // (cCCE, gCCE coherence), default off, max_mb (per rank) default 256
//      "cluster_only"          : false // true: write <file>.clusters.json and stop before the dynamics (cCCE, gCCE,
                                        // and cluster.mode pcce plans with cCCE only. Not with cluster.mode bundle
                                        // without rdip (checked connections) or the adaptive probe, factors and except
                                        // work), default false
//  },

    // ---- automatic convergence --------------------------------------------------------------------------------------
    // Details: docs/user/InputFile_AutoConvergence.md#tag_ac_keys
    //          (https://cqml-skku.github.io/CCEX/md_docs_2user_2InputFile__AutoConvergence.html#tag_ac_keys)
    // Example: example/tutorial/09_auto_convergence
//  "auto_convergence" : {              // main.out -f runs the automatic convergence: it chooses the time window, rdip,
                                        // rbath and the number of configurations, then runs the production. Its runs
                                        // use the mpirun -n N main.out was started with. The bath configurations are
                                        // the configuration block's (with k_range those k, without it the k whose files
                                        // exist); configuration.s_range only says that states are sampled (their number
                                        // is sampling's), configuration.seed is one number (absent: 1). A section is
                                        // found unless its value is given. A key written twice below is a choice: one
                                        // line finds, the other fixes
//      "enabled"     : true,           // default true, false runs the ordinary calculation
//      "resume"      : "continue",     // continue (default, reuse the workdir) or fresh (start over)
//      "workdir"     : "auto",         // folder of the runs, decisions and results (required)
//      "time_window" : {"decay_level": 0.01, "revival_check_factor": 2.5, "production_factor": 1.5},
                                        // found: T_d where the T2 fit reaches decay_level (|L|, default 0.01), the
                                        // production window production_factor x T_d
//      "time_window" : {"value": {"value": 1.0, "unit": "ms"}},
                                        // or fixed: value, ms (unit: s, ms, us, ns)
//      "rdip"        : {"candidates": [4, 6, 8], "threshold": 0.02},
                                        // found among candidates (A, the largest the reference): two agree within
                                        // threshold (the relative difference of the fitted area, default 0.02)
//      "rdip"        : {"num_candidates": 7, "threshold": 0.02},
                                        // or num_candidates found from the bath (default 7)
//      "rdip"        : {"value": 8},   // or fixed (A)
//      "rbath"       : {"candidates": [10, 15], "threshold": 0.02},
                                        // the same for rbath (num_candidates default 4)
//      "rbath"       : {"num_candidates": 4, "threshold": 0.02},
//      "rbath"       : {"value": 15},
//      "t2_fit"      : {"points": "peaks", "floor": 0.05, "n_range": [0.3, 4.0]},
                                        // points peaks (default) or all, floor: fit before |L| falls below it, n_range:
                                        // the allowed stretching exponent
//      "sampling"    : {               // how many configurations (and, with s_range, states per configuration): found
                                        // from a test, or fixed
//          "test_configurations" : 10, // configurations of the test, default 10 (chosen from the bath)
//          "threshold"           : 0.05,
                                        // the SE(A)/A target, default 0.05
//          "test_states"         : 5,  // s_range: states per test configuration, default 5
//          "states_threshold"    : 0.05,
                                        // s_range: the SE target of the states
//          "max_states"          : 30, // s_range: the most states per configuration
//          "min_configurations"  : 10, // the fewest production configurations
//          "max_configurations"  : 50  // the most production configurations
//      },
//      "sampling"    : {"test": [1, 2], "production": [3, 4], "threshold": 0.05},
                                        // or the configuration ids of the test and of the production given (in place of
                                        // test_configurations)
//      "sampling"    : {"configurations": 10, "states": 5}
                                        // or fixed: configurations, and with s_range states per configuration
//  }
}
)CCEX";
