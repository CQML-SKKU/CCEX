/**
 * @file general.cpp
 * @brief The Config structure: defaults, setters, getters and the checks of combined settings.
 */
#include "../../include/calculation/gpu/accelerator_defaults.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/clustering/bundle/bundle_defaults.h"
#include "../../include/input/general.h"
#include "../../include/calculation/superoperator.h"

#include "../../include/common/memory.h"
#include <string.h>
#include <set>
#include <string>
#include <float.h>  // DBL_EPSILON

// init
Config* Config_init(){
    Config* cnf = (Config*)allocArray1d(1,sizeof(Config));
    cnf->state_first = 1;   // configuration.s_range first (A94)
    cnf->method[0] = '\0';
    cnf->jump_operators = NULL;
    cnf->jump_operators_section[0] = '\0';
    cnf->quantity[0] = '\0';
    cnf->propagator[0] = '\0';
    cnf->sparse_min_dim = 4096;   // A2: propagator = sparse from this cluster dimension on
    cnf->sparse_tol = 1e-12;
    cnf->evolution[0] = '\0';
    cnf->evolution_isdefault = true;
    cnf->evolution_update[0] = '\0';           // not given; Config_validateEvolutionUpdate picks the default (A29)
    cnf->evolution_update_source[0] = '\0';
    cnf->column_evolution = true;              // A62: on by default (job 17492; {"enabled": false} turns it off)
    cnf->column_evolution_min_dim = 256;
    cnf->column_evolution_columns = 64;
    strcpy(cnf->backend, "cpu");
    cnf->gpu_min_dim = 768;                    // A37 (idm4-gpu 699): the card wins from ~670-750 on every path; 256 until 2026-10-13
    cnf->gpu_batch_steps = 64;
    cnf->gpu_cache_mb = ccex::gpuCacheMiB;
    cnf->probe_memory_gb = 0;
    cnf->gcce_refinement=false;   // off by default since 2026-10-02 (A3; idm job 16869: on - off <= 6.7e-7)
    cnf->gcce_refinement_max_sweeps=ccex::refinementSweeps;
    cnf->vector_mean_field=false;
    cnf->bath_mean_field_off=false;
    cnf->hamiltonian_ccce_given=false;
    cnf->delta_A=false;
    cnf->pseudospin=cnf->pseudospin_block=false;
    cnf->pcce_detail_limit=ccex::partitionDetailLimit;
    cnf->gpu_batch_clusters = 8;
    cnf->gpu_cluster_batch_mb = 64;
    cnf->gpu_packed_max_dim = -1;               // not given; Config_resolvePackedMaxDim picks the default (A26)
    cnf->gpu_packed_source[0] = '\0';
    cnf->gpu_warp_max_dim = -1;                 // not given; Config_resolveWarpMaxDim picks the default (A40)
    cnf->gpu_warp_source[0] = '\0';
    cnf->profile = false;
    cnf->cluster_only = false;
    cnf->cluster_contributions = false;
    cnf->cluster_contribution_max_mb = 256;
    cnf->nrdip_ladder = 0;
    strcpy(cnf->pcce_partition, "deterministic");
    strcpy(cnf->pcce_connections, "centers");
    strcpy(cnf->pcce_cache_mode, "off");
    cnf->pcce_cache_file[0] = '\0';
    strcpy(cnf->cluster_mode,"standard");
    cnf->bundle_max_dim=0;
    cnf->expansion_max_clusters=ccex::adaptiveMaxClusters;
    cnf->expansion_max_work=ccex::adaptiveMaxWork;
    cnf->bundle_merge_cutoff=0;
    strcpy(cnf->bundle_connections,"distance");
    cnf->bundle_connection_budget=-1;
    cnf->adaptive_rdip_is_set=false;
    cnf->adaptive_options_present=false;
    cnf->adaptive_work_file[0]='\0';
    cnf->hfmedi = false;
    cnf->knight = false;

    // Bath coordinate rotation : off. Config comes from calloc, so the two axes would
    // otherwise be [0,0,0] -- set them to +z so that "enabled" with no axis given is a
    // well-defined identity instead of a zero-norm error.
    cnf->rot_enabled = false;
    cnf->rot_qubit_axis[0] = 0.0; cnf->rot_qubit_axis[1] = 0.0; cnf->rot_qubit_axis[2] = 1.0;
    cnf->hf_frame[0] = '\0'; // empty : not given
    cnf->qd_frame[0] = '\0';

    cnf->defect_axis_reference = DEFECT_AXIS_REFERENCE_NONE;

    return cnf;
}

// alloc
void Config_allocBathfiles(Config* cnf){
    cnf->bathfiles = allocChar2d(cnf->nbathfiles,MAX_FILEPATH);
}

void Config_reallocBathfiles(Config* cnf, int oldlength, int newlength){
    reallocChar2d(&(cnf->bathfiles),oldlength,newlength,MAX_FILEPATH);
}

void Config_allocBathadjust(Config* cnf){
    cnf->bathadjust = allocDouble2d(cnf->nbathfiles,3);
}

void Config_reallocBathadjust(Config* cnf, int oldrow, int newrow){
    int col = 3;
    reallocDouble2d(&(cnf->bathadjust),oldrow,newrow,col);
}

void Config_allocGyrofile(Config* cnf){
    cnf->gyrofile = allocChar1d(MAX_FILEPATH);
}

void Config_allocHf_tensorfile(Config* cnf){
    cnf->hf_tensorfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocQd_tensorfile(Config* cnf){
    cnf->qd_tensorfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocQd_tensorfile_woqubit(Config* cnf){
    cnf->qd_tensorfile_woqubit = allocChar1d(MAX_FILEPATH);
}

void Config_allocAvaaxfile(Config* cnf){
    cnf->avaaxfile = allocChar2d(cnf->nbathfiles,MAX_FILEPATH);
}

void Config_allocStatefile(Config* cnf){
    cnf->statefile = allocChar2d(cnf->nbathfiles,MAX_FILEPATH);
}

void Config_allocExstatefile(Config* cnf){
    cnf->exstatefile = allocChar2d(cnf->nbathfiles,MAX_FILEPATH);
}

void Config_alloc_flines(Config* cnf, int size){
    cnf->_flines = allocInt1d(size);

}
void Config_realloc_flines(Config* cnf, int oldsize, int newsize){
    reallocInt1d(&(cnf->_flines),newsize);
}

// free
void Config_freeAll(Config* cnf){
    Config_freeBathfiles(cnf);
    Config_freeBathadjust(cnf);
    Config_freeGyrofile(cnf);
    Config_freeAvaaxfile(cnf);
    Config_freeStatefile(cnf);
    Config_freeExstatefile(cnf);
    Config_freeHf_tensorfile(cnf);
    Config_freeQd_tensorfile(cnf);
    Config_freeQd_tensorfile_woqubit(cnf);
    Config_free_flines(cnf);
    freeInt1d(&(cnf->bathfile_first));
    if (cnf->jump_operators) JumpOperatorArray_freeAll(cnf->jump_operators);
    freeArray1d((void**)&(cnf));
}

void Config_freeBathfiles(Config* cnf){
    freeChar2d(&(cnf->bathfiles),cnf->nbathfiles);
}

void Config_freeBathadjust(Config* cnf){
    freeDouble2d(&(cnf->bathadjust),cnf->nbathfiles);
}

void Config_freeGyrofile(Config* cnf){
    freeChar1d(&(cnf->gyrofile));
}

void Config_freeAvaaxfile(Config* cnf){
    freeChar2d(&(cnf->avaaxfile),cnf->nbathfiles);
}

void Config_freeStatefile(Config* cnf){
    freeChar2d(&(cnf->statefile),cnf->nbathfiles);
}

void Config_freeExstatefile(Config* cnf){
    freeChar2d(&(cnf->exstatefile),cnf->nbathfiles);
}

void Config_free_flines(Config* cnf){
    freeInt1d(&(cnf->_flines));
}

void Config_freeHf_tensorfile(Config* cnf){
    freeChar1d(&(cnf->hf_tensorfile));
}

void Config_freeQd_tensorfile(Config* cnf){
    freeChar1d(&(cnf->qd_tensorfile));
}

void Config_freeQd_tensorfile_woqubit(Config* cnf){
    freeChar1d(&(cnf->qd_tensorfile_woqubit));
}

/* Low level ---------------------------------------------------------*/

// get
char* Config_getMethod(Config* cnf){
    return cnf->method;
}

bool Config_getHfmedi(Config* cnf){
    return cnf->hfmedi;
}

bool Config_getKnight(Config* cnf){
    return cnf->knight;
}


char* Config_getQuantity(Config* cnf){
    return cnf->quantity;
}
char* Config_getPropagator(Config* cnf){
    return cnf->propagator;
}
char* Config_getEvolution(Config* cnf){
    return cnf->evolution;
}

int   Config_getOrder(Config* cnf){
    return cnf->order;
}

float* Config_getBfield(Config* cnf){
    return cnf->bfield;
}

float Config_getRbath(Config* cnf){
    return cnf->rbath;
}

float Config_getRdip(Config* cnf){
    return cnf->rdip;
}

double Config_getDeltat(Config* cnf){
    return cnf->deltat;
}

int   Config_getNstep(Config* cnf){
    return cnf->nstep;
}

float Config_getRbathcut(Config* cnf){
    return cnf->rbathcut;
}

float Config_getRdipcut(Config* cnf){
    return cnf->rdipcut;
}

int   Config_getNstate(Config* cnf){
    return cnf->nstate;
}

int   Config_getSeed(Config* cnf){
    return cnf->seed;
}

char*   Config_getGyrofile(Config* cnf){
    return cnf->gyrofile;
}

int     Config_getNbathfiles(Config* cnf){
    return cnf->nbathfiles;
}

char*   Config_getBathfiles_i(Config* cnf,int i){
    return cnf->bathfiles[i];
}

double* Config_getBathadjust_i(Config* cnf,int i){
    return cnf->bathadjust[i];
}

char**  Config_getAvaaxfile(Config* cnf){
    return cnf->avaaxfile;
}

char**  Config_getStatefile(Config* cnf){
    return cnf->statefile;
}

char**  Config_getExstatefile(Config* cnf){
    return cnf->exstatefile;
}

int Config_get_nflines(Config* cnf){
    return cnf->_nflines;
}

int* Config_get_flines(Config* cnf){
    return cnf->_flines;
}

int  Config_get_flines_i(Config* cnf, int i){
    return cnf->_flines[i];
}

bool Config_getRot_enabled(Config* cnf){
    return cnf->rot_enabled;
}

double* Config_getRot_qubit_axis(Config* cnf){
    return cnf->rot_qubit_axis;
}

char* Config_getHf_frame(Config* cnf){
    return cnf->hf_frame;
}

char* Config_getQd_frame(Config* cnf){
    return cnf->qd_frame;
}

DefectAxisReference Config_getDefect_axis_reference(Config* cnf){
    return cnf->defect_axis_reference;
}


char*   Config_getHf_tensorfile(Config* cnf){
    return cnf->hf_tensorfile;
}

int     Config_getHf_ignore_oor(Config* cnf){
    return cnf->hf_ignore_oor;
}

int     Config_getHf_readmode(Config* cnf){
    return cnf->hf_readmode;
}

char*   Config_getQd_tensorfile(Config* cnf){
    return cnf->qd_tensorfile;
}

char*   Config_getQd_tensorfile_woqubit(Config* cnf){
    return cnf->qd_tensorfile_woqubit;
}

double*   Config_getQd_cellpara(Config* cnf){
    return cnf->qd_cellpara;
}

int     Config_getQd_readmode(Config* cnf){
    return cnf->qd_readmode;
}

bool methodUsesPcceClusters(const char* m){
    return !strcasecmp(m,"pcce") || !strcasecmp(m,"pgcce") || !strcasecmp(m,"pmecce") || !strcasecmp(m,"pmegcce");
}

const char* methodKernel(const char* m){
    if (!strcasecmp(m,"gcce")   || !strcasecmp(m,"pgcce"))   return "gcce";
    if (!strcasecmp(m,"mecce")  || !strcasecmp(m,"pmecce"))  return "mecce";
    if (!strcasecmp(m,"megcce") || !strcasecmp(m,"pmegcce")) return "megcce";
    return "cce";
}

const char* inputMethodName(const char* m){
    for (const char* c : {"CCE", "pCCE", "dsj", "itb", "dsjitb"}) if (!strcasecmp(m,c)) return "cCCE";
    if (!strcasecmp(m,"pgcce"))   return "gCCE";
    if (!strcasecmp(m,"pmecce"))  return "meCCE";
    if (!strcasecmp(m,"pmegcce")) return "megCCE";
    return m;
}

void warnNotUsed(const char* key, const char* why){
    static std::set<std::string> warned;
    if (rank != 0 || !warned.insert(key).second) return;
    printf("Warning: input key \"%s\" is not used by this run: %s\n", key, why);
    fflush(stdout);
}

std::string stateFileName(const std::string& name, long s){
    std::string out;
    bool found = false;
    for (size_t i = 0; i < name.size();){
        if (!name.compare(i,5,"{{s}}")){ out += "{s}"; i += 5; }
        else if (!name.compare(i,3,"{s}")){ out += std::to_string(s); i += 3; found = true; }
        else out += name[i++];
    }
    return found ? out : out + std::to_string(s);
}

// set 
void Config_setMethod(Config* cnf, char* method){

    const int opsize = 12;   // 11 names; the empty 12th entry keeps accepting "" as before (opsize was 7 for 6 names)
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"gCCE","CCE","dsj","itb","dsjitb","pCCE","meCCE","megCCE","pgCCE",
                                                  "pmeCCE","pmegCCE"};

    int idx = findIndexCharFix(options,0,opsize-1,method);
    if (idx == -1) {
        ccexRefusalPrint(stderr, "Error: current method options (%s) is not available\n",method);
        ccexRefusalPrint(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            ccexRefusalPrint(stderr, "%s ",options[i]);
        }
        ccexRefusalPrint(stderr, "\n");
        ccexRefusalPrint(stderr, "Please check the input file or change the possible option set\n");
        ccexRefuse();
    }

    strcpy(cnf->method,options[idx]);
}

void Config_setQuantity(Config* cnf, char* quantity){

    // "noise" used to be accepted, but nothing computes it: the run gave the coherence under that name (review 20261013-doxygen D22)
    if (strcasecmp(quantity,"noise") == 0) {
        ccexRefusalPrint(stderr, "Error: quantity \"noise\" is not supported (no noise calculation exists), use coherence or dm\n");
        ccexRefuse();
    }

    const int opsize = 2;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"coherence","dm"};
    int idx = findIndexCharFix(options,0,opsize-1,quantity);
    if (idx == -1) {
        ccexRefusalPrint(stderr, "Error: current quantity options (%s) is not available\n",quantity);
        ccexRefusalPrint(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            ccexRefusalPrint(stderr, "%s ",options[i]);
        }
        ccexRefusalPrint(stderr, "\n");
        ccexRefusalPrint(stderr, "Please check the input file or change the possible option set\n");
        ccexRefuse();
    }

    strcpy(cnf->quantity,quantity);
}

void Config_setPropagator(Config* cnf, char* propagator){

    // Free-evolution propagator.
    //   eigen  : one SelfAdjointEigenSolver per cluster, exp(-i*H*tau) = M*diag(exp(-i*l*tau))*M^dag
    //   expm   : Eigen's general matrix exponential per step, the original CCEX path (gCCE only)
    //   sparse : A2, exp(-i*H*tau) applied to the state by a Krylov expansion, H in compressed rows, for clusters of
    //            dimension sparse_min_dim and above (Config_validatePropagator says where)
    const int opsize = 3;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"eigen","expm","sparse"};
    int idx = findIndexCharFix(options,0,opsize-1,propagator);
    if (idx == -1) {
        ccexRefusalPrint(stderr, "Error: current propagator options (%s) is not available\n",propagator);
        ccexRefusalPrint(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            ccexRefusalPrint(stderr, "%s ",options[i]);
        }
        ccexRefusalPrint(stderr, "\n");
        ccexRefusalPrint(stderr, "Please check the input file or change the possible option set\n");
        ccexRefuse();
    }

    strcpy(cnf->propagator,propagator);
}

void Config_setEvolution(Config* cnf, char* evolution){

    // gCCE state representation.
    //   vector : DEFAULT. rho0 is a pure state whenever nstate > 0 (QubitArray_Rho0 is
    //            always psi0*psi0^dag, and BathArray_Rho0 is too unless it is the
    //            ensemble average), so propagate |psi> instead of rho. Each step becomes
    //            matrix-VECTOR products, O(n^2) where the matrix path is O(n^3), and the
    //            partial trace collapses to reshaping |psi> into a qdim x bdim matrix Psi
    //            and forming Psi*Psi^dag.
    //   matrix : propagate the full density matrix, rho(t) = U*rho0*U^dag, then partial-
    //            trace it. What CCEX did before; works for every configuration.
    // Same physics, different arithmetic: the results agree to rounding but are not
    // guaranteed bit-identical, so a converged series should not switch mid-way.
    // vector needs nstate > 0 (an ensemble rho0 is mixed, so no state vector exists),
    // propagator=eigen (with expm the per-step matrix exponential dominates and there is
    // nothing to win) and method=gCCE. cJSON_readOptionConfig falls the DEFAULT back to
    // matrix wherever those do not hold; an EXPLICIT evolution=vector is never downgraded
    // -- calCoherenceGcce says why it cannot run and stops.
    const int opsize = 2;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"matrix","vector"};
    int idx = findIndexCharFix(options,0,opsize-1,evolution);
    if (idx == -1) {
        ccexRefusalPrint(stderr, "Error: current evolution options (%s) is not available\n",evolution);
        ccexRefusalPrint(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            ccexRefusalPrint(stderr, "%s ",options[i]);
        }
        ccexRefusalPrint(stderr, "\n");
        ccexRefusalPrint(stderr, "Please check the input file or change the possible option set\n");
        ccexRefuse();
    }

    strcpy(cnf->evolution,evolution);
}

/**
 * Pick the evolution path once every option source has been applied.
 *
 * This CANNOT be decided while the input file is being parsed: the readers run one after
 * another, and the method and nstate the decision needs are only settled once all of them
 * have run. Call this
 * after the readers, before Config_report, so the reported value is the one that runs.
 *
 * An explicit "evolution" is never rewritten -- only checked against the one condition
 * that is fatal outside gCCE. calCoherenceGcce enforces the other two.
 */
void Config_resolveEvolution(Config* cnf){

    bool isGCCE = !strcmp(methodKernel(cnf->method),"gcce");   // gCCE, pgCCE

    if (!cnf->evolution_isdefault){
        // Only calCoherenceGcce reads evolution. Any other method would ignore it and run
        // the density-matrix path anyway, so refuse rather than let the input claim
        // something the run does not do.
        if (strcasecmp(cnf->evolution,"matrix")!=0 && !isGCCE){
            ccexRefusalPrint(stderr,"Error : execution.evolution=%s is implemented for method=gCCE only "
                           "(this run has method=%s).\n",cnf->evolution,inputMethodName(cnf->method));
            ccexRefuse();
        }
        return;
    }

    bool vectorFits = isGCCE
                   && (cnf->nstate > 0)
                   && (strcasecmp(cnf->propagator,"expm")!=0);   // eigen, or sparse (A2)
    strcpy(cnf->evolution, vectorFits ? "vector" : "matrix");
}

/**
 * evolution_update = recursive exists on the CPU time evolutions of cCCE (method CCE, pCCE, dsj, itb or
 * dsjitb: the single state, evolveSingleState, and the ensemble nstate = 0, evolveEnsembleRecursiveCce)
 * and of gCCE with propagator=eigen (evolution=vector, propagateStateGcce, and evolution=matrix,
 * calPropagatorGcceRecursive), and in the packed kernel (small_cluster_batch.cpp). The GPU large-cluster path
 * (gpu_backend.cu: every cluster with backend cuda, those of gpu_options.min_dim and above with auto) computes every
 * phase with sin and cos as with exact, and gCCE with propagator=expm has no phases to carry over. Without the key
 * the run takes recursive where it exists and exact elsewhere, silently, per path: in the GPU build with backend auto the clusters on the CPU (below min_dim) and the packed kernel run
 * recursive and the clusters on the GPU exact; with backend cuda the value is recursive when the packed kernel is on
 * (the other clusters exact on the GPU), else exact. In the CPU build backend auto computes every cluster on the CPU.
 * A given recursive is refused where a cluster may reach a path without it (the input would claim something the run
 * does not do): gCCE expm, and a method the GPU takes with backend cuda, or auto in the GPU build.
 * A given exact keeps the bits of the exact phases. The one place this is decided. Call after Config_resolveEvolution and
 * Config_resolvePackedMaxDim (whether the packed kernel is on is part of the decision).
 */
void Config_validateEvolutionUpdate(Config* cnf){

    const bool isGCCE = !strcmp(methodKernel(cnf->method),"gcce");   // gCCE, pgCCE
    const bool eigen = (strcasecmp(cnf->propagator,"eigen")==0);
    char found[96] = "";                    // a CPU path without recursive
    if (isGCCE && !eigen)
        snprintf(found,sizeof(found),"method=gCCE with execution.propagator=%s",cnf->propagator);
    else if (strcasecmp(cnf->propagator,"sparse")==0)   // A2: the Krylov path has no phases to carry over; its clusters
        snprintf(found,sizeof(found),"propagator=sparse");   // below min_dim run exact with it
    else if (!strcmp(methodKernel(cnf->method),"mecce") || !strcmp(methodKernel(cnf->method),"megcce"))   // a matrix exponential per step
        snprintf(found,sizeof(found),"method=%s",inputMethodName(cnf->method));
    // A36: may a cluster reach the GPU large-cluster path? The methods computeBackendInit sends to the GPU, with backend
    // cuda (refused there in the CPU build) or with auto in the GPU build
    const bool cuda = (strcasecmp(cnf->backend,"cuda")==0);
    const bool gpuMethod = computeBackendMethodSupported(cnf);   // gCCE itself: pgCCE runs on the CPU
#ifdef CCEX_USE_CUDA
    const bool gpuLarge = gpuMethod && (cuda || !strcasecmp(cnf->backend,"auto"));
#else
    const bool gpuLarge = gpuMethod && cuda;
#endif
    const bool packed = cnf->gpu_packed_max_dim > 0;
    if (cnf->evolution_update[0] == '\0'){
        const bool allOnGpu = gpuLarge && cuda && !packed;      // nothing on a path with recursive
        char path[128] = "supported path";
        if (found[0] != '\0') snprintf(path,sizeof(path),"%s",found);
        else if (allOnGpu) snprintf(path,sizeof(path),"execution.backend=cuda: every cluster on the GPU, exact there");
        else if (gpuLarge)
            snprintf(path,sizeof(path),"%s: recursive, clusters on the GPU (execution.backend=%s%s): exact",
                     cuda ? "packed kernel" : packed ? "clusters on the CPU and packed kernel" : "clusters on the CPU",
                     cnf->backend, cnf->gpu_warp_max_dim > 0 ? ", warp kernel included" : "");
        strcpy(cnf->evolution_update, found[0] == '\0' && !allOnGpu ? "recursive" : "exact");
        snprintf(cnf->evolution_update_source,sizeof(cnf->evolution_update_source),"default (%s)",path);
        return;
    }
    strcpy(cnf->evolution_update_source,"given");
    if (strcasecmp(cnf->evolution_update,"recursive")!=0) return;
    if (found[0] == '\0' && gpuLarge)
        snprintf(found,sizeof(found),"execution.backend=%s, whose clusters on the GPU are computed exact",cnf->backend);
    if (found[0] == '\0') return;
    ccexRefusalPrint(stderr,"Error : execution.evolution_update=recursive is implemented on the CPU and in the packed kernel for cCCE "
                   "(single state and ensemble) and on the CPU for gCCE with execution.propagator=eigen (evolution=vector or "
                   "matrix), this run has %s. Use execution.evolution_update=exact, or leave the key out (then each path takes "
                   "recursive where it has it).\n",found);
    ccexRefuse();
}

/**
 * propagator: expm exists for gCCE only; sparse for the single-state cCCE (method CCE,
 * pCCE, dsj, itb or dsjitb with nstate > 0) and gCCE with evolution = vector and one qubit, on backend = cpu with
 * evolution_update = exact, without the knight, hfmedi and delta_A terms. Everything else would run another path than the one
 * asked for, so it is refused. Call after Config_resolveEvolution (main.cpp), when the method and nstate are settled.
 */
void Config_validatePropagator(const Config* cnf, int nqubit){

    const bool isGCCE = !strcmp(methodKernel(cnf->method),"gcce");   // gCCE, pgCCE
    if (!isGCCE && strcasecmp(cnf->propagator,"expm")==0){
        ccexRefusalPrint(stderr,"Error : execution.propagator=expm is implemented for method=gCCE only (this run has method=%s), "
                       "cCCE takes eigen or sparse.\n",inputMethodName(cnf->method));
        ccexRefuse();
    }
    if (strcasecmp(cnf->propagator,"sparse")!=0){ return; }
    char found[128] = "";
    if (cnf->nstate == 0)
        snprintf(found,sizeof(found),"the ensemble (no configuration.s_range)");
    else if (isGCCE && strcasecmp(cnf->evolution,"vector")!=0)
        snprintf(found,sizeof(found),"method=gCCE with execution.evolution=%s",cnf->evolution);
    else if (strcasecmp(cnf->backend,"cpu")!=0)
        snprintf(found,sizeof(found),"execution.backend=%s",cnf->backend);
    else if (strcasecmp(cnf->evolution_update,"recursive")==0)   // given; without the key sparse takes exact (A29)
        snprintf(found,sizeof(found),"execution.evolution_update=%s",cnf->evolution_update);
    else if (Config_getKnight((Config*)cnf) || Config_getHfmedi((Config*)cnf) || cnf->delta_A)
        snprintf(found,sizeof(found),"the %s term",Config_getKnight((Config*)cnf) ? "knight" :
                 Config_getHfmedi((Config*)cnf) ? "hfmedi" : "delta_A");
    else if (nqubit != 1)
        snprintf(found,sizeof(found),"nqubit=%d",nqubit);
    else if (!isGCCE && strcasecmp(cnf->quantity,"coherence")!=0)
        snprintf(found,sizeof(found),"quantity=%s",cnf->quantity);
    if (found[0] == '\0'){ return; }
    ccexRefusalPrint(stderr,"Error : execution.propagator=sparse is implemented on execution.backend=cpu with execution.evolution_update=exact for one sampled "
                   "state (configuration.s_range) of the cCCE and of gCCE with evolution=vector, one qubit, without knight, "
                   "hfmedi and delta_A, this run has %s. Use execution.propagator=eigen.\n",found);
    ccexRefuse();
}

/**
 * gpu_options.packed_max_dim without the key: 4 in the GPU build (bin/main_gpu.out) with backend cuda
 * or auto and method CCE - the packed kernel takes the cCCE clusters of bath dimension 2 and 4 - and 0 otherwise: with
 * backend cpu, with another method, and always in the CPU build, where nothing goes to the GPU. A value given in the
 * input is kept as it is (an explicit nonzero one is checked in smallClusterWarmUp). Call once every source is in.
 */
void Config_resolvePackedMaxDim(Config* cnf){
    if (cnf->gpu_packed_max_dim >= 0){ strcpy(cnf->gpu_packed_source,"given"); return; }
#ifdef CCEX_USE_CUDA
    const bool gpuBackend = strcasecmp(cnf->backend,"cpu")!=0;
    const bool on = gpuBackend && strcasecmp(cnf->method,"cce")==0 && strcasecmp(cnf->quantity,"coherence")==0;
    cnf->gpu_packed_max_dim = on ? 4 : 0;
    snprintf(cnf->gpu_packed_source,sizeof(cnf->gpu_packed_source),"default (%s)",
             on ? (strcasecmp(cnf->backend,"cuda") ? "execution.backend=auto, method=CCE" : "execution.backend=cuda, method=CCE")
                : gpuBackend ? "not method=CCE with quantity=coherence" : "backend=cpu");
#else
    cnf->gpu_packed_max_dim = 0;
    strcpy(cnf->gpu_packed_source,"default (CPU build)");
#endif
}

/**
 * gpu_options.warp_max_dim without the key: 27, the largest total dimension the warp kernel has, in
 * the GPU build with backend cuda or auto and method CCE (quantity coherence) or gCCE; 0 otherwise - backend cpu, another
 * method, and always in the CPU build. Where the run is one the kernel does not take (warpClusterRun: gcce_refinement,
 * pulseiter, savemode cluster_density_matrix, ...), the default switches itself off in the run ([gpu-warp] off). A value given in the
 * input is kept (an explicit nonzero one is checked in smallClusterWarmUp). Call once every source is in.
 */
void Config_resolveWarpMaxDim(Config* cnf){
    if (cnf->gpu_warp_max_dim >= 0){ strcpy(cnf->gpu_warp_source,"given"); return; }
#ifdef CCEX_USE_CUDA
    const bool gpuBackend = strcasecmp(cnf->backend,"cpu")!=0;
    const bool cce = strcasecmp(cnf->method,"cce")==0 && strcasecmp(cnf->quantity,"coherence")==0;
    const bool gcce = strcasecmp(cnf->method,"gcce")==0;
    const bool on = gpuBackend && (cce || gcce);
    cnf->gpu_warp_max_dim = on ? 27 : 0;
    snprintf(cnf->gpu_warp_source,sizeof(cnf->gpu_warp_source),"default (%s)",
             on ? (cce ? (strcasecmp(cnf->backend,"cuda") ? "execution.backend=auto, method=CCE" : "execution.backend=cuda, method=CCE")
                       : (strcasecmp(cnf->backend,"cuda") ? "execution.backend=auto, method=gCCE" : "execution.backend=cuda, method=gCCE"))
                : gpuBackend ? "not method=CCE (coherence) or gCCE" : "backend=cpu");
#else
    cnf->gpu_warp_max_dim = 0;
    strcpy(cnf->gpu_warp_source,"default (CPU build)");
#endif
}

void Config_setOrder(Config* cnf, int order){

    if (order < 0) {
        ccexRefusalPrint(stderr, "Error: current order (%d < 0) is not available\n",order);
        ccexRefuse();
    }

    cnf->order = order;
}

void Config_setBfield(Config* cnf, float* bfield){
    copyFloat1d(cnf->bfield,bfield,3);
}

void Config_setBfield_z(Config* cnf, float bz){
    cnf->bfield[2] = bz;
}

void Config_setRbath(Config* cnf, float rbath){

    if (rbath < 0.0) {
        ccexRefusalPrint(stderr, "Error: current rbath (%f < 0) is not available\n",rbath);
        ccexRefuse();
    }
    cnf->rbath = rbath;
}

void Config_setRdip(Config* cnf, float rdip){

    if (rdip < 0.0) {
        ccexRefusalPrint(stderr, "Error: current rdip (%f < 0) is not available\n",rdip);
        ccexRefuse();
    }
    cnf->rdip = rdip;
}

void Config_setDeltat(Config* cnf, double deltat){

    if (deltat < 0.0) {
        ccexRefusalPrint(stderr, "Error: current deltat (%f < 0) is not available\n",deltat);
        ccexRefuse();
    }
    cnf->deltat = deltat;
}

void Config_setNstep(Config* cnf, int nstep){

    // at least one time point: with 0 the result arrays are empty and the output reads result[0] (review 7eed38 R1)
    if (nstep < 1) {
        ccexRefusalPrint(stderr, "Error: current nstep (%d < 1) is not available: at least one time point\n",nstep);
        ccexRefuse();
    }
    cnf->nstep = nstep;
}

void Config_setRbathcut(Config* cnf, float rbathcut){

    if (rbathcut < 0.0) {
        ccexRefusalPrint(stderr, "Error: the inner radius rbath[0] (%f < 0) is not available\n",rbathcut);
        ccexRefuse();
    }
    cnf->rbathcut = rbathcut;
}

void Config_setRdipcut(Config* cnf, float rdipcut){

    if (rdipcut < 0.0) {
        ccexRefusalPrint(stderr, "Error: the inner radius rdip[0] (%f < 0) is not available\n",rdipcut);
        ccexRefuse();
    }
    cnf->rdipcut = rdipcut;
}

void Config_setNstate(Config* cnf, int nstate){

    if (nstate < 0) {
        ccexRefusalPrint(stderr, "Error: current nstate (%d < 0) is not available\n",nstate);
        ccexRefuse();
    }
    cnf->nstate = nstate;
}

void Config_setSeed(Config* cnf, int seed){

    if (seed < 0) {
        ccexRefusalPrint(stderr, "Error: current seed(%d) ( < 0) is not available\n",seed);
        ccexRefuse();
    }
    cnf->seed = seed;
}

void Config_setGyrofile(Config* cnf, char* gyrofile){

    if (cnf->gyrofile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for gyrofile should be allocated\n");
        ccexRefuse();
    }
    copyBounded(cnf->gyrofile,gyrofile,MAX_FILEPATH,"gyrofile");
}

void Config_setNbathfiles(Config* cnf, int nbathfiles){

    if (nbathfiles <= 0) {
        ccexRefusalPrint(stderr, "Error: current nbathfiles(%d) (<= 0) is not available\n",nbathfiles);
        ccexRefuse();
    }
    cnf->nbathfiles = nbathfiles;
}

void Config_setBathfiles_i(Config* cnf, char* bathfiles, int i){

    if (cnf->bathfiles == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for bathfiles should be allocated\n");
        ccexRefuse();
    }
    copyBounded(cnf->bathfiles[i],bathfiles,MAX_FILEPATH,"bathfile");
}

void Config_setBathadjust_i(Config* cnf, double* bathadjust, int i){

    if (cnf->bathadjust == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for bathadjust should be allocated\n");
        ccexRefuse();
    }
    copyDouble1d(cnf->bathadjust[i],bathadjust,3);
}

void Config_setAvaaxfile_i(Config* cnf, char* avaaxfile, int i){

    if (cnf->avaaxfile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for avaaxfile should be allocated\n");
        ccexRefuse();
    }

    copyBounded(cnf->avaaxfile[i],avaaxfile,MAX_FILEPATH,"avaaxfile");
}

void Config_setStatefile_i(Config* cnf, char* statefile, int i){

    if (cnf->statefile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for statefile should be allocated\n");
        ccexRefuse();
    }

    copyBounded(cnf->statefile[i],statefile,MAX_FILEPATH,"statefile");
}

void Config_setExstatefile_i(Config* cnf, char* exstatefile, int i){

    if (cnf->exstatefile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for exstatefile should be allocated\n");
        ccexRefuse();
    }

    copyBounded(cnf->exstatefile[i],exstatefile,MAX_FILEPATH,"exstatefile");
}

void Config_set_nflines(Config* cnf, int nflines){
    cnf->_nflines = nflines;
}

void Config_set_flines_i(Config* cnf, int fline, int i){
    cnf->_flines[i] = fline;
}

// The axes are only copied here. Whether they describe a usable frame -- finite,
// non-zero, not collinear -- is decided in one place, buildQubitAlignedRotation,
// which is also what the rotation itself is built from.
void Config_setRot_enabled(Config* cnf, bool rot_enabled){
    cnf->rot_enabled = rot_enabled;
}

void Config_setRot_qubit_axis(Config* cnf, double* rot_qubit_axis){
    copyDouble1d(cnf->rot_qubit_axis,rot_qubit_axis,3);
}

// "file" or "qubit", checked where it is read (option.cpp readFrameKey)
void Config_setHf_frame(Config* cnf, const char* frame){
    snprintf(cnf->hf_frame,sizeof cnf->hf_frame,"%s",frame);
}

void Config_setQd_frame(Config* cnf, const char* frame){
    snprintf(cnf->qd_frame,sizeof cnf->qd_frame,"%s",frame);
}

void Config_setDefect_axis_reference(Config* cnf, DefectAxisReference reference){
    cnf->defect_axis_reference = reference;
}

/**
 * @brief With the rotation on, require the field to lie along the computational z axis
 * @details Mirrors Config_resolveEvolution in timing: it cannot run during the file
 *          parse: it needs the field and the rotation, which different readers set.
 *
 *          Without the rotation this is a no-op, so an arbitrary field vector keeps
 *          working exactly as it always has. With the rotation on, the computational z
 *          axis IS the qubit axis, and "bfield" is read in that computational frame --
 *          so [0,0,500] means 500 G along +qubit_axis and [0,0,-500] means 500 G the
 *          other way. Both are accepted: the requirement is the LINE, not the sign.
 *          A tilted field is not supported alongside the rotation and is refused here.
 *
 *          The tolerance is tight on purpose. Nothing rotates these numbers -- they are
 *          the input as written -- so anything beyond round-off in the parse is a real
 *          transverse component that the user meant and that this version cannot honour.
*/
void Config_validateBfieldAlignment(Config* cnf){

    if (!cnf->rot_enabled){ return; }

    double bx = cnf->bfield[0], by = cnf->bfield[1], bz = cnf->bfield[2];
    double bperp = sqrt(bx*bx + by*by);
    double bnorm = sqrt(bx*bx + by*by + bz*bz);

    if (bnorm == 0.0){ return; } // no field at all : nothing to align

    double tol = 1.0e-12 * ((bnorm > 1.0) ? bnorm : 1.0);

    if (bperp > tol){
        ccexRefusalPrint(stderr,"Error: frame.rotation requires bfield to be aligned\n");
        ccexRefusalPrint(stderr,"with the computational z axis, which represents qubit_axis.\n");
        ccexRefusalPrint(stderr,"Current bfield: [%g, %g, %g]\n",bx,by,bz);
        ccexRefusalPrint(stderr,"  transverse component |B_perp| = %g  (tolerance %g)\n",bperp,tol);
        ccexRefusalPrint(stderr,"Write the field as [0, 0, Bz], Bz may be negative for the opposite\n");
        ccexRefusalPrint(stderr,"direction along qubit_axis. A tilted field is not supported with the\n");
        ccexRefusalPrint(stderr,"coordinate frame rotation.\n");
        ccexRefuse();
    }
}

void Config_setHf_tensorfile(Config* cnf, char* hf_tensorfile){

    if (cnf->hf_tensorfile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for hf_tensorfile should be allocated\n");
        ccexRefuse();
    }
    copyBounded(cnf->hf_tensorfile,hf_tensorfile,MAX_FILEPATH,"hf_tensorfile");
}

void Config_setHf_ignore_oor(Config* cnf, int hf_ignore_oor){

    if (hf_ignore_oor != 0 && hf_ignore_oor != 1) {
        ccexRefusalPrint(stderr, "Error: possible hf_ignore_oor (current : %d) is 0 or 1\n",hf_ignore_oor);
        ccexRefuse();
    }
    cnf->hf_ignore_oor = hf_ignore_oor;
}

void Config_setHf_readmode(Config* cnf, int hf_readmode){

    // Anything outside 0..3 used to fall through -- the guard read
    // "hf_readmode < 0 && hf_readmode > 3", which is never true -- and readHftensorfile
    // has no else after its `==0` / `==1||2||3` branches, so neither
    // BathArray_setBathHypfs nor the tensor file ever ran. Every bath spin kept an unset
    // hyperfine tensor and the run finished, quietly, with no hyperfine coupling at all.
    if (hf_readmode < 0 || hf_readmode > 3) {
        ccexRefusalPrint(stderr, "Error: current hf_readmode (%d < 0 or > 3) is not available\n",hf_readmode);
        ccexRefuse();
    }

    cnf->hf_readmode = hf_readmode;
}

void Config_setQd_tensorfile(Config* cnf, char* qd_tensorfile){

    if (cnf->qd_tensorfile == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for qd_tensorfile should be allocated\n");
        ccexRefuse();
    }
    copyBounded(cnf->qd_tensorfile,qd_tensorfile,MAX_FILEPATH,"qd_tensorfile");
}

void Config_setQd_tensorfile_woqubit(Config* cnf, char* qd_tensorfile_woqubit){

    if (cnf->qd_tensorfile_woqubit == NULL){
        ccexRefusalPrint(stderr, "Error: the memory for qd_tensorfile_woqubit should be allocated\n");
        ccexRefuse();
    }
    copyBounded(cnf->qd_tensorfile_woqubit,qd_tensorfile_woqubit,MAX_FILEPATH,"qd_tensorfile_woqubit");
}

void Config_setQd_cellpara(Config* cnf, double* qd_cellpara){
    copyDouble1d(cnf->qd_cellpara,qd_cellpara,3);
}

void Config_setQd_readmode(Config* cnf, int qd_readmode){

    // The supported set is not a range: 1 was removed, so 0, 2, 3 and 4 are listed one
    // by one rather than bounded. Anything else used to fall through -- the guard read
    // "qd_readmode < 0 && qd_readmode > 4", which is never true -- and then missed every
    // case of the switch in cJSON_readOptionConfig and every branch of readQdtensorfile,
    // so the run went on quietly with no quadrupole interaction at all.
    if (qd_readmode != 0 &&
        qd_readmode != 2 &&
        qd_readmode != 3 &&
        qd_readmode != 4) {
        ccexRefusalPrint(stderr,
                "Error: current qd_readmode (%d) is not available, "
                "supported values are 0, 2, 3 and 4\n",
                qd_readmode);
        ccexRefuse();
    }
    cnf->qd_readmode = qd_readmode;
}

void Config_setHfmedi(Config* cnf, bool hfmedi){

    cnf->hfmedi = hfmedi;
}

void Config_setKnight(Config* cnf, bool knight){

    cnf->knight = knight;
}


// The rows are named by their input keys (A76). Without -v the report is the run's summary: the rows below that are
// always printed (the checks read evolution_update, propagator, column_evolution and gpu_options, with their "chosen
// by"), those of the cluster mode in use, and the settings that are not their default; -v prints every row.
void Config_report(Config* cnf, bool qubitBathMeanField){

    printTitle("Structure Config");

    const bool pcce = methodUsesPcceClusters(cnf->method);
    const bool bundle = !strcasecmp(cnf->cluster_mode,"adaptive");   // cluster.mode bundle
    printStructElementChar("method",cnf->pseudospin ? (char*)"pseudospin" : Config_getMethod(cnf));   // pseudospin: CCE inside
    if (cnf->pseudospin) printStructElementChar("hamiltonian.pseudospin.pair",(char*)(cnf->pseudospin_block ? "block" : "channel"));
    printStructElementChar("quantity",Config_getQuantity(cnf));
    printStructElementChar("propagator",Config_getPropagator(cnf));
    if (strcasecmp(cnf->propagator,"sparse")==0){
        printStructElementInt("propagator.min_dim",int(cnf->sparse_min_dim));
        printStructElementDouble("propagator.tol",cnf->sparse_tol);
    }
    printStructElementChar("evolution",Config_getEvolution(cnf));
    printStructElementChar("evolution_update",cnf->evolution_update);
    printStructElementChar("  chosen by",cnf->evolution_update_source);
    if (cnf->column_evolution){                // A62: printed only when on
        printStructElementInt("column_evolution.min_dim",cnf->column_evolution_min_dim);
        printStructElementInt("column_evolution.columns",cnf->column_evolution_columns);
    }
    printStructElementChar("backend",cnf->backend);
    printStructElementInt("gpu_options.min_dim",cnf->gpu_min_dim);
    printStructElementInt("gpu_options.evolution_batch_steps",cnf->gpu_batch_steps);
    printStructElementInt("gpu_options.cache_mb",cnf->gpu_cache_mb);
    if (verbosity || bundle)
        printStructElementDouble("cluster.adaptive.probe_memory_gb",cnf->probe_memory_gb);
    if (verbosity || cnf->gcce_refinement){
        printStructElementBool("gcce_refinement.enabled",cnf->gcce_refinement);
        printStructElementInt("gcce_refinement.max_sweeps",cnf->gcce_refinement_max_sweeps);
    }
    if (verbosity || cnf->vector_mean_field || cnf->bath_mean_field_off)
        printStructElementChar("hamiltonian.common.bath_bath_mean_field",
                               cnf->bath_mean_field_off ? (char*)"off" : cnf->vector_mean_field ? (char*)"xyz" : (char*)"z");
    printStructElementInt("gpu_options.evolution_batch_size",cnf->gpu_batch_clusters);
    printStructElementInt("gpu_options.cluster_batch_mb",cnf->gpu_cluster_batch_mb);
    printStructElementInt("gpu_options.packed_max_dim",cnf->gpu_packed_max_dim);
    printStructElementChar("  chosen by",cnf->gpu_packed_source);
    printStructElementInt("gpu_options.warp_max_dim",cnf->gpu_warp_max_dim);
    printStructElementChar("  chosen by",cnf->gpu_warp_source);
    // cluster.mode as the input names it: the method carries pcce, dsj, itb, dsjitb inside (Config_setMethod)
    const char* mode = bundle ? "bundle" : pcce ? "pcce" :
        (!strcasecmp(cnf->method,"dsj") || !strcasecmp(cnf->method,"itb") || !strcasecmp(cnf->method,"dsjitb")) ? cnf->method : "standard";
    printStructElementChar("cluster.mode",(char*)mode);
    if (bundle) printStructElementChar("cluster.bundle.connections",cnf->bundle_connections);   // A94 (G8): the mode chosen
    if (verbosity || pcce){
        printStructElementInt("cluster.pcce.print_trials",cnf->pcce_detail_limit);
        printStructElementChar("cluster.pcce.partition",cnf->pcce_partition);
        printStructElementChar("cluster.pcce.connections",cnf->pcce_connections);
        printStructElementChar("cluster.pcce.cache.mode",cnf->pcce_cache_mode);
        printStructElementChar("cluster.pcce.cache.file",cnf->pcce_cache_file);
    }
    if (bundle) {
        if(cnf->bundle_max_dim)
            printStructElementInt("cluster.bundle.max_dim",cnf->bundle_max_dim);
        else printf("      %-18s : %-21s\n","cluster.bundle.max_dim","auto (from cluster dimension)");
        printStructElementInt("cluster.bundle.max_clusters",cnf->expansion_max_clusters);
        printStructElementInt("cluster.bundle.max_work",cnf->expansion_max_work);
        printStructElementDouble("cluster.bundle.merge_cutoff",cnf->bundle_merge_cutoff);
        printStructElementChar("cluster.bundle.connections",cnf->bundle_connections);
        if(!strcasecmp(cnf->bundle_connections,"auto"))
            printStructElementDouble("cluster.bundle.connection_budget",cnf->bundle_connection_budget);
    }
    printStructElementInt("order",Config_getOrder(cnf));
    printStructElementFloat1d("bfield",Config_getBfield(cnf),3);
    printStructElementFloat("rbath",Config_getRbath(cnf));
    if (verbosity || Config_getRbathcut(cnf)!=0) printStructElementFloat("rbath[0] (inner)",Config_getRbathcut(cnf));
    if(bundle && !strcasecmp(cnf->bundle_connections,"checked"))
        printStructElementChar("rdip","omitted (coherence-checked links)");
    else if(bundle && !strcasecmp(cnf->bundle_connections,"auto") &&
       !cnf->adaptive_rdip_is_set)printStructElementChar("rdip","unlimited");
    else printStructElementFloat("rdip",Config_getRdip(cnf));
    if (verbosity || Config_getRdipcut(cnf)!=0) printStructElementFloat("rdip[0] (inner)",Config_getRdipcut(cnf));
    printStructElementFloat("deltat",Config_getDeltat(cnf));
    printStructElementInt("nstep",Config_getNstep(cnf));
    if (Config_getNstate(cnf) == 0) printStructElementChar("configuration.s_range","(absent: the ensemble)");
    else {
        char range[64];
        snprintf(range,sizeof range,"[%d, %d]",cnf->state_first,Config_getNstate(cnf));
        printStructElementChar("configuration.s_range",range);
    }
    if (cnf->seed_used) printStructElementInt("configuration.seed",Config_getSeed(cnf));
    else printStructElementChar("configuration.seed","(not used)");

    printStructElementChar("configuration.gyrofile",Config_getGyrofile(cnf));
    for (int i = 0; i < Config_getNbathfiles(cnf); i++){
        printStructElementChar("configuration.bathfile",Config_getBathfiles_i(cnf,i));
        const double* shift = Config_getBathadjust_i(cnf,i);
        if (verbosity || shift[0]!=0 || shift[1]!=0 || shift[2]!=0)
            printStructElementDouble1d("configuration.bath_shift",Config_getBathadjust_i(cnf,i),3);
    }

    printStructElementBool("frame.rotation",Config_getRot_enabled(cnf));
    if (Config_getRot_enabled(cnf)){
        printStructElementDouble1d("  qubit_axis",Config_getRot_qubit_axis(cnf),3);
    }

    {
        const char* keys[3] = {"configuration.defect_axis_file", "configuration.bath_state_file", "configuration.defect_state_file"};
        char** names[3] = {Config_getAvaaxfile(cnf), Config_getStatefile(cnf), Config_getExstatefile(cnf)};
        for (int k = 0; k < 3; k++){
            if (!names[k]){ if (verbosity) printStructElementChar((char*)keys[k],(char*)"(null)"); continue; }
            for (int i = 0; i < Config_getNbathfiles(cnf); i++) printStructElementChar((char*)keys[k],names[k][i]);
        }
    }


    // the value the run uses (QubitArray_getOverhaus, set from the input in option.cpp)
    printStructElementBool("hamiltonian.common.qubit_bath_mean_field",qubitBathMeanField);

    // hamiltonian.common.hyperfine {dipole, contact}: hf_readmode 0 point/false, 1 point/true, 2 dft/false, 3 dft/true
    const int hf = Config_getHf_readmode(cnf);
    if (verbosity || hf!=0){
        printStructElementChar("hamiltonian.common.hyperfine.dipole",hf>=2 ? (char*)"dft" : (char*)"point");
        printStructElementBool("hamiltonian.common.hyperfine.contact",hf==1 || hf==3);
    }
    if (hf!=0){
        printStructElementChar("hamiltonian.common.hyperfine.file",Config_getHf_tensorfile(cnf));
        printStructElementChar("hamiltonian.common.hyperfine.missing_row",
                               Config_getHf_ignore_oor(cnf) ? (char*)"point_dipole" : (char*)"stop");
    }
    if (Config_getRot_enabled(cnf) && Config_getHf_frame(cnf)[0])
        printStructElementChar("hamiltonian.common.hyperfine.frame",Config_getHf_frame(cnf));

    // hamiltonian.common.quadrupole: qd_readmode 0 off, 2 plain, 3 strained, 4 bubble
    const int qd = Config_getQd_readmode(cnf);
    if (verbosity || qd!=0)
        printStructElementChar("hamiltonian.common.quadrupole",
                               qd==2 ? (char*)"plain" : qd==3 ? (char*)"strained" : qd==4 ? (char*)"bubble" : (char*)"off");
    if (qd==2){
        printStructElementChar("hamiltonian.common.quadrupole.file",Config_getQd_tensorfile(cnf));
    }else if (qd==3 || qd==4 ){
        printStructElementChar("hamiltonian.common.quadrupole.file",Config_getQd_tensorfile(cnf));
        printStructElementChar("hamiltonian.common.quadrupole.file_without_defect",Config_getQd_tensorfile_woqubit(cnf));
        printStructElementDouble1d("hamiltonian.common.quadrupole.cell",Config_getQd_cellpara(cnf),3);
    }
    if (Config_getRot_enabled(cnf) && Config_getQd_frame(cnf)[0])
        printStructElementChar("hamiltonian.common.quadrupole.frame",Config_getQd_frame(cnf));

    if (verbosity || Config_getHfmedi(cnf)) printStructElementBool("hamiltonian.cCCE.hfmedi",Config_getHfmedi(cnf));
    if (verbosity || Config_getKnight(cnf)) printStructElementBool("hamiltonian.cCCE.knight",Config_getKnight(cnf));
    if (verbosity || cnf->delta_A) printStructElementBool("hamiltonian.cCCE.delta_A",cnf->delta_A);

    if (cnf->jump_operators){
        char section[MAX_CHARARRAY_LENGTH + 16];
        snprintf(section,sizeof(section),"hamiltonian.%s",cnf->jump_operators_section);
        printStructElementChar("jump_operators",section);
        JumpOperatorArray_report(cnf->jump_operators);
    }
}
