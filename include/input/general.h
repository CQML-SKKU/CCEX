/**
 * @file general.h
 * @brief The Config structure: all run settings, with their setters and getters.
 */
#ifndef __CCEX_SIMULATOR_GENERAL_H_
#define __CCEX_SIMULATOR_GENERAL_H_

#include "common/utilities.h"
#include <string>

/**
 * @enum DefectAxisReference
 * @brief Reference axis for a future Jahn-Teller axis classification
 * @details The input key that set it was removed (it is refused), so it
 *          stays NONE. QUBIT_AXIS means the central qubit's
 *          symmetry axis, which only frame.rotation defines.
*/
typedef enum {
    DEFECT_AXIS_REFERENCE_NONE       = 0, /**< key absent (default) */
    DEFECT_AXIS_REFERENCE_QUBIT_AXIS = 1  /**< "qubit_axis" */
} DefectAxisReference;

struct JumpOperatorArray;   // calculation/superoperator.h: the jump operators of meCCE / megCCE

/**
 * @struct Config
 * @brief All run settings: the input keys (and the -f/-o/-v options) as the simulation uses them.
 * @details Filled by the input reader (src/input/option.cpp) with the defaults of Config_init, then read by
 *          every part of the run.
 */
typedef struct {

    // Method
    /**
     * @brief Clustering algorithm : cce (default) | gcce | dsj | itb | dsjitb | pcce
     * @details The following is the description of each algorithm.
     * 
     *  - cce    : Conventional CCE method *using hash algorithm
     *  - gcce   : Generalized CCE method *using hash algorithm
     *  - pcce   : Partition CCE method *using hash algorithm
     *  - dsj    : Disjoint clustering algorithm
     *             , meaning that it consider the disjointed cluster
     *  - itb    : Inter-bathcluster clustering algorithm
     *             , meaning that it consider the pair interaction between disjointed bath cluster
     *  - dsjitb : Disjoint + inter-bathcluster clustering algorithm
     *             , meaning that it consider the disjointed bundle 
     *             and the interaction between disjointed clusters.
     *             (dsjitb with order = 2 gives the same result as normal CCE2)
     * @todo add pcce-related description HS_pcce
     * @ref Cluster::method
     * 
    */
    char method[MAX_CHARARRAY_LENGTH];
    struct JumpOperatorArray* jump_operators;   /**< hamiltonian.meCCE / megCCE.jump_operators; NULL: none (ME methods then run at rate 0) */
    char jump_operators_section[MAX_CHARARRAY_LENGTH];   /**< "meCCE" or "megCCE": the hamiltonian section it came from ("" if none) */
    char quantity[MAX_CHARARRAY_LENGTH];  /**< Measurement : coherence | dm (noise is refused: not implemented) */
    char propagator[MAX_CHARARRAY_LENGTH];/**< free-evolution propagator : eigen | expm (gCCE) | sparse */
    long sparse_min_dim;                  /**< propagator sparse: clusters of this Hilbert dimension and above take it (default 4096) */
    double sparse_tol;                    /**< propagator sparse: the Krylov expansion tolerance (default 1e-12) */
    char evolution[MAX_CHARARRAY_LENGTH]; /**< gCCE state representation : matrix | vector   */
    bool evolution_isdefault;             /**< true until an "evolution" key is actually read */
    char evolution_update[16];            /**< exact | recursive: the phases of the time grid by sin/cos at every point, or by the previous point's times one step's factor; empty until given or resolved (recursive where supported, else exact) */
    char evolution_update_source[160];    /**< why evolution_update has its value (given, or the default and the path; per path in the GPU build), for the run log */
    bool column_evolution;                /**< one-state cCCE clusters of bath dimension >= column_evolution_min_dim evolve their time points as matrix columns (cceEvolveColumns); on by default */
    int column_evolution_min_dim;         /**< the smallest bath dimension that takes it (default 256) */
    int column_evolution_columns;         /**< time points per column block (default 64) */
    char backend[8];                      /**< cpu (default) | cuda | auto */
    int gpu_min_dim;                       /**< auto uses CUDA at/above this Hilbert dimension */
    int gpu_batch_steps;                   /**< maximum time points per CUDA batch */
    bool gcce_refinement;                 /**< compensated eigensystem refinement; off by default */
    int gcce_refinement_max_sweeps;
    bool vector_mean_field;               /**< the static mean field of the other bath spins as the vector J . m_j, not its z part alone (off by default; hamiltonian.common.bath_bath_mean_field "xyz") */
    bool bath_mean_field_off;             /**< no static mean field of the other bath spins at all (hamiltonian.common.bath_bath_mean_field "off") */
    bool hamiltonian_ccce_given;          /**< hamiltonian.cCCE has a setting: refused in a gCCE run (validateHamiltonianForMethod, main.cpp) */
    bool delta_A;                         /**< hamiltonian.cCCE.delta_A: the second-order term delta_A(m_s).I of each bath spin (Maze 2008 eq. A3, A83) */
    bool pseudospin;                      /**< method pseudospin: the cCCE with its pair clusters in the pseudospin model (pseudospin.cpp) */
    bool pseudospin_block;                /**< hamiltonian.pseudospin.pair "block" (the cut pair exactly); false: "channel" */
    int pcce_detail_limit;                /**< rank-zero trial details, 0 disables */
    int gpu_cache_mb;                      /**< reusable device-buffer cache per rank; 0 disables */
    double probe_memory_gb;                /**< node budget for the probes running at once; 0: 80 % of RAM (and of the GPU) */
    int gpu_batch_clusters;                /**< maximum same-dimension clusters per batch */
    int gpu_cluster_batch_mb;              /**< estimated live buffer budget for bundleing clusters */
    int gpu_packed_max_dim;                /**< gpu_options.packed_max_dim: cCCE clusters of this bath dimension or less (2, 4) to the packed kernel; 0 off; -1 until Config_resolvePackedMaxDim */
    char gpu_packed_source[64];            /**< "given", or the default and why (Config_resolvePackedMaxDim) */
    int gpu_warp_max_dim;                  /**< gpu_options.warp_max_dim: clusters of a compiled total dimension up to this value to the warp kernel; 0 off; -1 until Config_resolveWarpMaxDim */
    char gpu_warp_source[64];              /**< "given", or the default and why (Config_resolveWarpMaxDim) */
    bool profile;                         /**< print per-rank, per-dimension timing totals */
    bool cluster_only;                    /**< write OUTFILE.clusters.json and stop before the dynamics (output.cluster_only) */
    bool cluster_contributions;           /**< record irreducible cluster log importance over all input times */
    int cluster_contribution_max_mb;      /**< diagnostic raw-table budget per MPI rank */
    int nrdip_ladder;                     /**< number of rdip_ladder rungs (0 = off) */
    float rdip_ladder[64];                /**< extra rdip rungs assembled from one order-2 coherence run: CCE, gCCE or pCCE (src/auto_convergence/rdip_ladder.cpp) */
    char pcce_partition[16];              /**< deterministic (default) | legacy */
    char pcce_connections[16];            /**< centers (default) | spins: which distance connects two partitions */
    char pcce_cache_mode[16];             /**< off (default) | read | readwrite */
    char pcce_cache_file[MAX_FILEPATH];    /**< explicit partition cache file, relative to cwd */
    char cluster_mode[16];                 /**< internal clustering mode: standard (cluster.mode standard and pcce) | adaptive (cluster.mode bundle) */
    int bundle_max_dim;            /**< required dimension cap for one bundle */
    int expansion_max_clusters;            /**< subset-closed bundle family count limit */
    int expansion_max_work;       /**< graph/subset enumeration work limit */
    double bundle_merge_cutoff;          /**< coupling-norm bound threshold, rad/ms; default 0 */
    char bundle_connections[16];         /**< distance | legacy auto | internal checked when rdip is omitted */
    double bundle_connection_budget;     /**< auto: maximum tmax * sum of omitted bundle coupling bounds */
    bool adaptive_rdip_is_set;             /**< omitted rdip is unlimited only in adaptive auto */
    bool adaptive_options_present;
    char adaptive_work_file[MAX_FILEPATH]; /**< optional request for the refinement driver */

    // General options
    int   order;        /**< Maximum cluster size ( >= 0 ): bath sites (standard), partitions of sK spins (pCCE), bundles (cluster.mode bundle) */
    float bfield[3];    /**< Magnetic field (Bx,Bx,Bz) (Unit G)                  */
    float rbath;        /**< Bath radius (Unit angstrom)                         */
    float rdip;         /**< Dipole radius (Unit angstrom)                       */
    double deltat;      /**< Time interval (Unit ms). A double: as a float, 0.08 was
                             0.0799999982, and every time on the grid carried that. */
    int   nstep;        /**< Time step                                           */
    float rbathcut;     /**< Cutoff to remove spins within rbathcut              */
    float rdipcut;      /**< Cutoff to remove spin pairs within rdipcut          */
    int   nstate;       /**< configuration.s_range last: the sampled states are state_first .. nstate (0: the ensemble) */
    int   state_first;  /**< configuration.s_range first (A94): states 1 .. state_first-1 are drawn and not calculated */
    int   seed;         /**< Fixed random value (if none, time(null) )           */
    bool  seed_used;    /**< the run draws from seed: nstate > 0, Defect species or pCCE (the log reports the seed only then) */

    // Qubit and Bath file options
    char* gyrofile; /**< Gyro file name */
    int nbathfiles; /**< Number of bath files */
    char** bathfiles; /**< Bath file names */
    double** bathadjust; /**< Bath adjust values */
    // A64: one name per bath file (nbathfiles of each, or NULL when the key is absent), read in the order of the bath files
    char** avaaxfile; /**< configuration.defect_axis_file: available principal axis, per bath file */
    char** statefile; /**< configuration.bath_state_file: state file prefixes, per bath file */
    char** exstatefile ; /**< configuration.defect_state_file: extra spin' state file prefixes, per bath file */
    int* bathfile_first; /**< nbathfiles+1: the spins of bath file j are bathfile_first[j] .. bathfile_first[j+1]-1 */
    int _nflines; /**< number of file lines */
    int* _flines; /**< file line number of each bath spins */

    // Bath coordinate rotation
    /**
     * @brief If true, the qubits and the bath are re-expressed in the qubit-aligned frame
     * @details The transform is PASSIVE : it puts rot_qubit_axis on the computational
     *          +z, so R * normalize(rot_qubit_axis) = [0,0,1].
     *
     *          The origin is the first qubit's position. The computational x axis is
     *          normalize(qubit_axis x [0,0,1]) (identity for qubit_axis = +z).
     *
     *          Each place that holds directional numbers states its "frame" (A90):
     *          "file" = written in the bath file's frame, rotated with the positions;
     *          "qubit" = written in the computational frame, used as written. Required
     *          everywhere while the rotation is on (no default).
     *
     *          Transformed:
     *            - every Qubit.xyz and every bath-spin position
     *            - intmap entries marked INTMAP_AUTO_SOURCE_GEOMETRY (the point-dipole
     *              tensor built from the source qubit geometry)
     *            - Qubit.intmap entries with frame = "file"
     *            - stored HF / QD tensors that are in the file frame -- every
     *              point-dipole tensor CCEX computed itself, plus the tensor-file ones
     *              when hamiltonian.common.hyperfine.frame / quadrupole.frame say "file"
     *            - each Defect species whose frame is "file": its physical axis,
     *              nuclear_positions, and hyperfine/efg/zfs tensors
     *
     *          Kept as written:
     *            - Qubit.intmap entries and Defect species with frame = "qubit"
     *            - scalar Defect detuning values (frame independent)
     *            - bfield
     *
     *          bfield is expressed directly in the computational frame. When
     *          frame.rotation is enabled, Config_validateBfieldAlignment
     *          requires it to lie along the computational z axis, which represents
     *          qubit_axis.
     * @ref buildQubitAlignedRotation, applyCoordinateFrameRotation
    */
    bool rot_enabled;
    double rot_qubit_axis[3]; /**< Qubit symmetry axis, in the bath file frame (default [0,0,1]) */
    /**
     * @brief Frame the HF tensor-file COMPONENTS are written in : "file" | "qubit"
     *        (hamiltonian.common.hyperfine.frame; frame.rotation.hf_tensor_frame until A90)
     * @details "file" means the 3x3 components are in the bath file's Cartesian basis
     *          and are transformed as R*T*R^T; "qubit" means they are already in the
     *          computational basis and are left alone. It says nothing about the POSITION
     *          columns of the tensor file, which are always read in the file frame,
     *          because that is what the lookup matches against.
     *          Empty = not given; required when the rotation is on and hamiltonian.common.hyperfine
     *          has "dipole": "dft". With "dipole": "point" the anisotropic part is a point dipole from the
     *          source geometry, so a given frame is not used (only noted in the log).
    */
    char hf_frame[MAX_CHARARRAY_LENGTH];
    char qd_frame[MAX_CHARARRAY_LENGTH]; /**< Same, for the quadrupole tensors (hamiltonian.common.quadrupole.frame) */

    // Defect axis reference
    /**
     * @brief Which axis a future Jahn-Teller classification should measure against
     * @details Recorded and validated only. Nothing in the Defect Hamiltonian, the JT
     *          assignment, zfs, hypf, efg or rxyzs reads it yet, so setting it changes
     *          no number this version produces -- the run log says as much.
     * @ref DefectAxisReference
    */
    DefectAxisReference defect_axis_reference;

    // DFT Hyperfine tensor
    int hf_readmode; /**< Read mode for hyperfine tensor file : 0 : off , 1~3 : fermi & DFT */ 
    char* hf_tensorfile; /**< Hyperfine tensor file name */
    int hf_ignore_oor; /**< A spin inside the hyperfine file's range but with no entry in it: 0 = error, else its point-dipole tensor is used (no spin is removed) */
    

    // DFT Quadrupole tensor
    int qd_readmode; /**< Read mode for quadrupole tensor file : 0 off; 2, 3, 4 read qd_tensorfile (3 and 4 also need qd_tensorfile_woqubit and qd_cellpara); 1 was removed and is refused */
    char* qd_tensorfile; /**< Quadrupole tensor file name */
    char* qd_tensorfile_woqubit ; /**< Quadrupole tensor file name without qubit */
    double qd_cellpara[3]; 

    // HF medi
    bool hfmedi; /**< If true, we will add Hyperfine mediated term (Available only when you use CCE method), ref. N. Zhao, Nature Nanotechnology, 6, 242-246 (2011) (Supplementary)*/
    bool knight; /**< If true, we will add the knight-field effect induced by transverse magnetic field (Available only when you use CCE method), ref. N. Zhao, Nature Nanotechnology, 6, 242-246 (2011) (Supplementary)*/ 

    // semi-classical options
    // int Interval_filter; // Default : 500      //gsl interval in integration
    // int Step_filter; // Default : 100          //filter function time table step
    // double epsabs_filter; // Default : 1e-10
    // double tolerance_filter; // Default : 1e-10
    // double RoundOff_err_filter; // Default : 1e+6
    // int SaveCorr; // Default : 1    //1 : save the correlation function

    //filter function time table
    // std::vector<double> FilterTable;    //[0 ~ 1], step : Step_filter
    
} Config;

Config* Config_init();
void Config_freeAll(Config* cnf);
void Config_report(Config* cnf, bool qubitBathMeanField);

/* Low level ----------------------------------------------------*/

// get
char*   Config_getMethod(Config* cnf);
char*   Config_getQuantity(Config* cnf);
char*   Config_getPropagator(Config* cnf);
char*   Config_getEvolution(Config* cnf);
int     Config_getOrder(Config* cnf);
float*  Config_getBfield(Config* cnf);
float   Config_getRbath(Config* cnf);
float   Config_getRdip(Config* cnf);
double  Config_getDeltat(Config* cnf);
int     Config_getNstep(Config* cnf);
float   Config_getRbathcut(Config* cnf);
float   Config_getRdipcut(Config* cnf);
int     Config_getNstate(Config* cnf);
int     Config_getSeed(Config* cnf);
char*   Config_getGyrofile(Config* cnf);
int     Config_getNbathfiles(Config* cnf);
char*   Config_getBathfiles_i(Config* cnf,int i);
double* Config_getBathadjust_i(Config* cnf,int i);
char**  Config_getAvaaxfile(Config* cnf);     // nbathfiles names, or NULL
char**  Config_getStatefile(Config* cnf);     // nbathfiles prefixes, or NULL
char**  Config_getExstatefile(Config* cnf);   // nbathfiles prefixes, or NULL
int     Config_get_nflines(Config* cnf);
int*    Config_get_flines(Config* cnf);
int     Config_get_flines_i(Config* cnf, int i);
bool    Config_getRot_enabled(Config* cnf);
double* Config_getRot_qubit_axis(Config* cnf);
char*   Config_getHf_frame(Config* cnf);
char*   Config_getQd_frame(Config* cnf);
DefectAxisReference Config_getDefect_axis_reference(Config* cnf);
char*   Config_getHf_tensorfile(Config* cnf);
int     Config_getHf_ignore_oor(Config* cnf);
int     Config_getHf_readmode(Config* cnf);
char*   Config_getQd_tensorfile(Config* cnf);
char*   Config_getQd_tensorfile_woqubit(Config* cnf);
double* Config_getQd_cellpara(Config* cnf);
int     Config_getQd_readmode(Config* cnf);
bool    Config_getHfmedi(Config* cnf);
bool    Config_getKnight(Config* cnf);

// The internal method names (option.cpp internalMethod): which clusters and which kernel. cluster_mode pcce makes the
// partition clusters of pCCE (cCCE), pgCCE, pmeCCE and pmegCCE (gCCE, meCCE, megCCE with cluster_mode pcce); the kernel
// computing a cluster is "cce" (CCE, pCCE, dsj, itb, dsjitb), "gcce" (gCCE, pgCCE), "mecce" (meCCE, pmeCCE) or "megcce"
// (megCCE, pmegCCE).
bool methodUsesPcceClusters(const char* m);
const char* methodKernel(const char* m);
// The method as the input writes it, for messages (A47): CCE, pCCE, dsj, itb, dsjitb are cCCE (the reader refuses them as
// old method names); pgCCE, pmeCCE, pmegCCE are gCCE, meCCE, megCCE (refused as the adding_MECCE branch's names).
const char* inputMethodName(const char* m);
/** A94 (G3): a key the input wrote that this run reads but does not use, reported once per key on rank 0 (stdout):
 *  "Warning: input key "<key>" is not used by this run: <why>" */
void warnNotUsed(const char* key, const char* why);
/** A94 (G14): the file of state s of a configuration.bath_state_file or defect_state_file name: "{s}" in it is s ("{{s}}" a
 *  literal {s}); a name without {s} has s appended (the form before A94) */
std::string stateFileName(const std::string& name, long s);

// set 
void Config_setMethod(Config* cnf, char* method);
void Config_setQuantity(Config* cnf, char* quantity);
void Config_setPropagator(Config* cnf, char* propagator);
void Config_setEvolution(Config* cnf, char* evolution);
void Config_resolveEvolution(Config* cnf);
void Config_validateEvolutionUpdate(Config* cnf);
void Config_validatePropagator(const Config* cnf, int nqubit);
void Config_resolvePackedMaxDim(Config* cnf);   /**< gpu_options.packed_max_dim without the key */
void Config_resolveWarpMaxDim(Config* cnf);     /**< gpu_options.warp_max_dim without the key */
void Config_setOrder(Config* cnf, int order);
void Config_setBfield(Config* cnf, float* bfield);
void Config_setBfield_z(Config* cnf, float bz);
void Config_setRbath(Config* cnf, float rbath);
void Config_setRdip(Config* cnf, float rdip);
void Config_setDeltat(Config* cnf, double deltat);
void Config_setNstep(Config* cnf, int nstep);
void Config_setRbathcut(Config* cnf, float rbathcut);
void Config_setRdipcut(Config* cnf, float rdipcut);
void Config_setNstate(Config* cnf, int nstate);
void Config_setSeed(Config* cnf, int seed);
void Config_validatePcceCache(const Config* cnf);
void Config_setGyrofile(Config* cnf, char* gyrofile);
void Config_setNbathfiles(Config* cnf, int nbathfiles);
void Config_setBathfiles_i(Config* cnf, char* bathfiles, int i);
void Config_setBathadjust_i(Config* cnf, double* bathadjust, int i);
void Config_setAvaaxfile_i(Config* cnf, char* avaaxfile, int i);
void Config_setStatefile_i(Config* cnf, char* statefile, int i);
void Config_setExstatefile_i(Config* cnf, char* exstatefile, int i);
void Config_set_nflines(Config* cnf, int nflines);
void Config_set_flines_i(Config* cnf, int fline, int i);
void Config_setRot_enabled(Config* cnf, bool rot_enabled);
void Config_setRot_qubit_axis(Config* cnf, double* rot_qubit_axis);
void Config_setHf_frame(Config* cnf, const char* frame);
void Config_setQd_frame(Config* cnf, const char* frame);
void Config_setDefect_axis_reference(Config* cnf, DefectAxisReference reference);

// With the rotation on, the computational z axis IS the qubit axis, so the field has to
// lie along it. Checked once every option source (input file, then CLI) is in.
void Config_validateBfieldAlignment(Config* cnf);

void Config_setHf_tensorfile(Config* cnf, char* hf_tensorfile);
void Config_setHf_ignore_oor(Config* cnf, int hf_ignore_oor);
void Config_setHf_readmode(Config* cnf, int hf_readmode);
void Config_setQd_tensorfile(Config* cnf, char* qd_tensorfile);
void Config_setQd_tensorfile_woqubit(Config* cnf, char* qd_tensorfile_woqubit);
void Config_setQd_readmode(Config* cnf, int qd_readmode);
void Config_setQd_cellpara(Config* cnf, double* qd_cellpara);
void Config_setHfmedi(Config* cnf, bool hfmedi);
void Config_setKnight(Config* cnf, bool knight);

// alloc
void Config_allocBathfiles(Config* cnf);
void Config_reallocBathfiles(Config* cnf, int oldsize, int newsize);
void Config_allocBathadjust(Config* cnf);
void Config_reallocBathadjust(Config* cnf, int oldsize, int newsize);
void Config_allocGyrofile(Config* cnf);
void Config_allocAvaaxfile(Config* cnf);
void Config_allocStatefile(Config* cnf);
void Config_allocExstatefile(Config* cnf);
void Config_alloc_flines(Config* cnf, int size);
void Config_realloc_flines(Config* cnf, int oldsize, int newsize);
void Config_allocHf_tensorfile(Config* cnf);
void Config_allocQd_tensorfile(Config* cnf);
void Config_allocQd_tensorfile_woqubit(Config* cnf);

// free
void Config_freeBathfiles(Config* cnf);
void Config_freeBathadjust(Config* cnf);
void Config_freeGyrofile(Config* cnf);
void Config_freeAvaaxfile(Config* cnf);
void Config_freeStatefile(Config* cnf);
void Config_freeExstatefile(Config* cnf);
void Config_free_flines(Config* cnf);
void Config_freeHf_tensorfile(Config* cnf);
void Config_freeQd_tensorfile(Config* cnf);
void Config_freeQd_tensorfile_woqubit(Config* cnf);

#endif // __CCEX_SIMULATOR_GENERAL_H_
