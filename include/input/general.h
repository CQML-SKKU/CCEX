#ifndef __CCEX_SIMULATOR_GENERAL_H_
#define __CCEX_SIMULATOR_GENERAL_H_

#include "common/utilities.h"

/**
 * @enum DefectAxisReference
 * @brief Reference axis for a future Jahn-Teller axis classification
 * @details Set from the top-level "defect_axis_reference" key. QUBIT_AXIS means the
 *          central qubit's symmetry axis, which only coordinate_frame_rotation defines,
 *          so the key requires that section to be enabled.
*/
typedef enum {
    DEFECT_AXIS_REFERENCE_NONE       = 0, /**< key absent (default) */
    DEFECT_AXIS_REFERENCE_QUBIT_AXIS = 1  /**< "qubit_axis" */
} DefectAxisReference;

struct JumpOperatorArray;   // calculation/superoperator.h: the jump operators of meCCE / megCCE

/**
 * @struct Config
 * @brief This structure contains the parameters for the simulation.
 * @details Each parameter is read from the input file or options.
 *          and those parameters are used in the simulation directly.
 */
typedef struct {

    // Method
    /**
     * @brief Clustering algorithm : cce (default) | gcce | dsj | itb | dsjitb | pcce
     * @details The following is the description of each algorithm.
     * 
     *  - cce    : Conventional CCE method *using hash algorithm
     *  - gcce   : Generalized CCE method *using hash algorithm
     *  - pcce   : Parition CCE method *using hash algorithm
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
    char propagator[MAX_CHARARRAY_LENGTH];/**< free-evolution propagator : eigen | expm (gCCE) | sparse (A2) */
    long sparse_min_dim;                  /**< propagator sparse: clusters of this Hilbert dimension and above take it (default 4096) */
    double sparse_tol;                    /**< propagator sparse: the Krylov expansion tolerance (default 1e-12) */
    char evolution[MAX_CHARARRAY_LENGTH]; /**< gCCE state representation : matrix | vector   */
    bool evolution_isdefault;             /**< true until an "evolution" key is actually read */
    char evolution_update[16];            /**< exact | recursive: the phases of the time grid by sin/cos at every point, or by the previous point's times one step's factor (A16 b); empty until given or resolved (A29: recursive where supported, else exact) */
    char evolution_update_source[160];    /**< why evolution_update has its value (given, or the default and the path; per path in the GPU build, A36), for the run log */
    char backend[8];                      /**< cpu (default) | cuda | auto */
    int gpu_min_dim;                       /**< auto uses CUDA at/above this Hilbert dimension */
    int gpu_batch_steps;                   /**< maximum time points per CUDA batch */
    bool gcce_refinement;                 /**< compensated eigensystem refinement; off by default since 2026-10-02 */
    int gcce_refinement_max_sweeps;
    bool vector_mean_field;               /**< the static mean field of the other bath spins as the vector J . m_j, not its z part alone (2026-10-03, off by default; hamiltonian.common.bath_bath_mean_field "xyz") */
    bool bath_mean_field_off;             /**< no static mean field of the other bath spins at all (hamiltonian.common.bath_bath_mean_field "off", A47) */
    bool hamiltonian_ccce_given;          /**< hamiltonian.cCCE has a setting: refused in a gCCE run (validateHamiltonianForMethod, main.cpp) */
    bool delta_A;                         /**< hamiltonian.cCCE.delta_A: false replaces the hyperfine tensors of a cluster's bath spins by their mean (A47) */
    int pcce_detail_limit;                /**< rank-zero trial details, 0 disables */
    int gpu_cache_mb;                      /**< reusable device-buffer cache per rank; 0 disables */
    double probe_memory_gb;                /**< node budget for the probes running at once; 0: 80 % of RAM (and of the GPU) */
    int gpu_batch_clusters;                /**< maximum same-dimension clusters per batch */
    int gpu_cluster_batch_mb;              /**< estimated live buffer budget for bundleing clusters */
    int gpu_packed_max_dim;                /**< gpu_options.packed_max_dim: cCCE clusters of this bath dimension or less (2, 4) to the packed kernel; 0 off; -1 until Config_resolvePackedMaxDim */
    char gpu_packed_source[64];            /**< "given", or the default and why (Config_resolvePackedMaxDim) */
    int gpu_warp_max_dim;                  /**< gpu_options.warp_max_dim (A40): clusters of a compiled total dimension up to this value to the warp kernel; 0 off; -1 until Config_resolveWarpMaxDim */
    char gpu_warp_source[64];              /**< "given", or the default and why (Config_resolveWarpMaxDim) */
    bool profile;                         /**< print per-rank, per-dimension timing totals */
    bool cluster_contributions;           /**< record irreducible cluster log importance over all input times */
    int cluster_contribution_max_mb;      /**< diagnostic raw-table budget per MPI rank */
    int nrdip_ladder;                     /**< number of rdip_ladder rungs (0 = off) */
    float rdip_ladder[64];                /**< extra rdip rungs assembled from one order-2 coherence run: CCE, gCCE or pCCE (src/auto_convergence/rdip_ladder.cpp) */
    char pcce_partition[16];              /**< deterministic (default) | legacy */
    char pcce_connections[16];            /**< centers (default) | spins: which distance connects two partitions */
    char pcce_cache_mode[16];             /**< off (default) | read | readwrite */
    char pcce_cache_file[MAX_FILEPATH];    /**< explicit partition cache file, relative to cwd */
    char cluster_mode[16];                 /**< standard (default) | adaptive */
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
    int   order;        /**< Maximum cluster size ( >= 0 ): bath sites (standard), partitions of sK spins (pCCE), bundles (cluster_mode adaptive) */
    float bfield[3];    /**< Magnetic field (Bx,Bx,Bz) (Unit G)                  */
    float rbath;        /**< Bath radius (Unit angstrom)                         */
    float rdip;         /**< Dipole radius (Unit angstrom)                       */
    double deltat;      /**< Time interval (Unit ms). A double: as a float, 0.08 was
                             0.0799999982, and every time on the grid carried that. */
    int   nstep;        /**< Time step                                           */
    float rbathcut;     /**< Cutoff to remove spins within rbathcut              */
    float rdipcut;      /**< Cutoff to remove spin pairs within rdipcut          */
    int   nstate;       /**< Number of bath states (single and gcce method only) */
    int   seed;         /**< Fixed random value (if none, time(null) )           */

    // Qubit and Bath file options
    char* qubitfile; /**< Qubit file name */
    char* gyrofile; /**< Gyro file name */
    int nbathfiles; /**< Number of bath files */
    char** bathfiles; /**< Bath file names */
    double** bathadjust; /**< Bath adjust values */
    char* avaaxfile; /**< Available principal axis */
    char* statefile; /**< State file */
    char* exstatefile ; /**< extra spin' state file */
    int _nflines; /**< number of file lines */
    int* _flines; /**< file line number of each bath spins */

    // Bath coordinate rotation
    /**
     * @brief If true, the qubits and the bath are re-expressed in the qubit-aligned frame
     * @details The transform is PASSIVE : it puts rot_qubit_axis on the computational
     *          +z, so R * normalize(rot_qubit_axis) = [0,0,1].
     *
     *          Transformed:
     *            - every Qubit.xyz, not just the reference one
     *            - every bath-spin position
     *            - intmap entries marked INTMAP_AUTO_SOURCE_GEOMETRY (the point-dipole
     *              tensor built from the source qubit geometry)
     *            - explicit Qubit.intmap entries with tensor_frame = "bath", self/ZFS
     *              entries included
     *            - stored HF / QD tensors that are in the bath frame -- every
     *              point-dipole tensor CCEX computed itself, plus the tensor-file ones
     *              when hf_tensor_frame / qd_tensor_frame say "bath"
     *            - each Defect entry whose coordinate_frame is "bath": its physical
     *              axis, rxyzs displacement vectors, and hypf/efg/zfs tensors
     *
     *          Kept as written:
     *            - explicit Qubit.intmap entries with tensor_frame = "qubit"
     *            - the legacy top-level "qzfs", which belongs to the single-qubit
     *              qubitfile path and has no tensor_frame to declare; it is read as a
     *              computational-frame tensor for backward compatibility
     *            - each Defect entry whose coordinate_frame is "qubit"
     *            - scalar Defect detuning values (frame independent)
     *            - bfield
     *
     *          bfield is expressed directly in the computational frame. When
     *          coordinate_frame_rotation is enabled, Config_validateBfieldAlignment
     *          requires it to lie along the computational z axis, which represents
     *          qubit_axis.
     * @ref buildQubitAlignedRotation, applyCoordinateFrameRotation
    */
    bool rot_enabled;
    double rot_bath_axis[3];  /**< Reference axis of the bath file frame; fixes the azimuth only (default [0,0,1]) */
    double rot_qubit_axis[3]; /**< Qubit symmetry axis, in the bath file frame (default [0,0,1]) */
    char rot_reference_qubit[MAX_CHARARRAY_LENGTH]; /**< Qubit.name held fixed by the rotation; empty = Qubit[0] */
    /**
     * @brief Frame every Qubit.xyz is written in : "bath"
     * @details Applies to the whole QubitArray -- per-qubit frames are not a thing. Only
     *          "bath" is supported: the qubit positions use the same Cartesian basis AND
     *          the same origin as the bath files, so one R and one r0 move everything.
     *          A "qubit" value would mean the qubit positions are already computational
     *          while the bath is not, and readBathfiles would then mix frames in its
     *          rbath / rbathcut selection, its mindist and its point-dipole tensors.
     *          Required when nqubit > 1 and the rotation is on; empty means "bath".
    */
    char rot_qubit_position_frame[MAX_CHARARRAY_LENGTH];
    /**
     * @brief Basis the external HF tensor COMPONENTS are written in : "bath" | "qubit"
     * @details "bath" means the 3x3 components are in the original bath Cartesian basis
     *          and must be transformed as R*T*R^T; "qubit" means they are already in the
     *          computational basis and are left alone. It says nothing about the POSITION
     *          columns of the tensor file, which are always read in the original frame,
     *          because that is what the lookup matches against.
     *          Empty = not given; required when the rotation is on and hamiltonian.common.hyperfine
     *          has "dipole": "dft". With "dipole": "point" and "contact": true the anisotropic part is rebuilt as a point dipole from the
     *          source geometry, so a given frame is not used (only noted in the log).
    */
    char rot_hf_tensor_frame[MAX_CHARARRAY_LENGTH];
    char rot_qd_tensor_frame[MAX_CHARARRAY_LENGTH]; /**< Same, for the quadrupole tensors */

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
void Config_report(Config* cnf);

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
char*   Config_getQubitfile(Config* cnf);
char*   Config_getGyrofile(Config* cnf);
int     Config_getNbathfiles(Config* cnf);
char*   Config_getBathfiles_i(Config* cnf,int i);
double* Config_getBathadjust_i(Config* cnf,int i);
char*   Config_getAvaaxfile(Config* cnf);
char*   Config_getStatefile(Config* cnf);
char*   Config_getExstatefile(Config* cnf);
int     Config_get_nflines(Config* cnf);
int*    Config_get_flines(Config* cnf);
int     Config_get_flines_i(Config* cnf, int i);
bool    Config_getRot_enabled(Config* cnf);
double* Config_getRot_bath_axis(Config* cnf);
double* Config_getRot_qubit_axis(Config* cnf);
char*   Config_getRot_reference_qubit(Config* cnf);
char*   Config_getRot_qubit_position_frame(Config* cnf);
char*   Config_getRot_hf_tensor_frame(Config* cnf);
char*   Config_getRot_qd_tensor_frame(Config* cnf);
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

// set 
void Config_setMethod(Config* cnf, char* method);
void Config_setQuantity(Config* cnf, char* quantity);
void Config_setPropagator(Config* cnf, char* propagator);
void Config_setEvolution(Config* cnf, char* evolution);
void Config_resolveEvolution(Config* cnf);
void Config_validateEvolutionUpdate(Config* cnf);
void Config_validatePropagator(const Config* cnf, int nqubit);
void Config_resolvePackedMaxDim(Config* cnf);   /**< gpu_options.packed_max_dim without the key (A26) */
void Config_resolveWarpMaxDim(Config* cnf);     /**< gpu_options.warp_max_dim without the key (A40) */
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
void Config_setQubitfile(Config* cnf, char* qubitfile);
void Config_setGyrofile(Config* cnf, char* gyrofile);
void Config_setNbathfiles(Config* cnf, int nbathfiles);
void Config_setBathfiles_i(Config* cnf, char* bathfiles, int i);
void Config_setBathadjust_i(Config* cnf, double* bathadjust, int i);
void Config_setAvaaxfile(Config* cnf, char* avaaxfile);
void Config_setStatefile(Config* cnf, char* statefile);
void Config_setExstatefile(Config* cnf, char* exstatefile);
void Config_set_nflines(Config* cnf, int nflines);
void Config_set_flines_i(Config* cnf, int fline, int i);
void Config_setRot_enabled(Config* cnf, bool rot_enabled);
void Config_setRot_bath_axis(Config* cnf, double* rot_bath_axis);
void Config_setRot_qubit_axis(Config* cnf, double* rot_qubit_axis);
void Config_setRot_reference_qubit(Config* cnf, char* rot_reference_qubit);
void Config_setRot_qubit_position_frame(Config* cnf, char* frame);
void Config_setRot_hf_tensor_frame(Config* cnf, char* frame);
void Config_setRot_qd_tensor_frame(Config* cnf, char* frame);
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
void Config_allocQubitfile(Config* cnf);
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
void Config_freeQubitfile(Config* cnf);
void Config_freeAvaaxfile(Config* cnf);
void Config_freeStatefile(Config* cnf);
void Config_freeExstatefile(Config* cnf);
void Config_free_flines(Config* cnf);
void Config_freeHf_tensorfile(Config* cnf);
void Config_freeQd_tensorfile(Config* cnf);
void Config_freeQd_tensorfile_woqubit(Config* cnf);

#endif // __CCEX_SIMULATOR_GENERAL_H_
