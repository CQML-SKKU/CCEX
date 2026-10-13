// ReadOption
#include <fcntl.h>
#include "../include/input/option.h"

// Structures
#include "../include/physical/qubit.h"
#include "../include/physical/bath.h"
#include "../include/physical/defect.h"
#include "../include/clustering/cluster.h"
#include "../include/physical/pulse.h"
#include "../include/physical/pulse_bath.h"
#include "../include/common/output.h"
#include "../include/input/general.h"

// Simulators
#include "../include/calculation/simulator.h"
#include "../include/calculation/simulator_gcce_bath.h"
#include "../include/calculation/gpu/compute_backend.h"
#include "../include/clustering/bundle/bundle_cluster.h"
#include "../include/clustering/bundle/bundle_work.h"
#include "../include/clustering/bundle/bundle_connection_check.h"
#include "../include/clustering/bundle/bundle_contributions.h"
#include "../include/auto_convergence/rdip_ladder.h"

// Utilities
#include "../include/common/utilities.h"
#include "../include/common/memory.h"
#include <getopt.h>
#include <cerrno>
#include <climits>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <pwd.h>
#include <grp.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#include <vector>
#include <cmath>
#include <cstdlib>

// Read files
#include "../include/input/reader.h"

//global variables
#include <memory>
#include "../include/clustering/algo_adaptive/base_run.h"
#include "../include/clustering/algo_adaptive/options_table.h"
#include <string>
extern char** environ;   // the outer MPI launch's variables are dropped from it (autoConvEntry)
int rank,nprocess;
bool verbosity = false;

namespace {
// One complete calculation: build the clusters this request asks for, run them,
// and write what they produced. Everything before it - the options and the five
// input files - is read once and reused, which is the whole point of being able
// to call this more than once in a process (plan stage 8).
//
// The per-process state a second call would otherwise inherit is cleared first:
// the work request this process has read, the bundles a previous plan installed,
// the connection check, and the clusters themselves. The random stream goes back
// to where a calculation starts - see the srand below.
void runOneCalculation(Config* cnf, QubitArray* qa, BathArray* ba, DefectArray* dfa,
                       Pulse* pls, Cluster* cls, Output* op, bool clusterPlanOnly) {
    resetAdaptiveWork();
    clearAdaptiveInstalledBundles();
    clearAdaptiveConnectionCheck();
    Cluster_freeClusinfo(cls);
    // Back to where one calculation per process always started: the seed, then
    // past the numbers the random defect axes drew at set-up. Seeding alone put
    // the state sampling at draw 0 instead, which changed every nstate>0 result
    // for a bath with random axes - the P1 regression cases, by up to 0.88
    // (measured against the build before this function existed, 2026-09-19).
    srand(Config_getSeed(cnf));
    for (long i = 0; i < DefectArray_paxesDraws(); ++i)
        rand();
    double time_start_step = 0., time_end = 0.;
    //=======================================================
    // Clusterize
    //=======================================================

    if (rank==0){
        printf("\n\n");
        printf("   ======================================================================\n");
        printf("       Clusterize \n");
        printf("   ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    const int loadedSites = ba->nspin;
    // A38: nothing below reads every pair row unless one of these is on (data/a38-local-cluster-list-20261012/01-analysis.md)
    cls->localRequest = !clusterPlanOnly && !cnf->cluster_contributions && !adaptiveWork(cnf).active;
    Cluster_clusterize(cls,ba,cnf,qa,dfa);
    cls->localRequest = false;
    checkAdaptiveConnections(qa,ba,dfa,cnf,pls,cls,op);
    prepareContributionClusters(cnf,cls,ba);
    rdipladder::prepare(cnf,cls,ba,op);
    if(clusterPlanOnly && rank==0)
        writeClusterInventory(cls,ba,qa,dfa,cnf,op->outfile,loadedSites);
    if (rank==0){
        Cluster_reportClusinfo(cls);

        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("          Wall time = %.5f s\n", time_end - time_start_step);
        printLine();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //=======================================================
    // Main calculation & save results
    //=======================================================

    if (rank==0){
        printf("\n\n");
        printf("   ======================================================================\n");
        printf("       Main calculation \n");
        printf("   ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    // MPI distribution for the clusters
    int order = Cluster_getOrder(cls);
    int*** clusters = Cluster_getClusinfo(cls);
    // A31: the cost split (Cluster_costFunction); after clusterize, which can cut the bath (pCCE).
    const std::function<double(int,const int*)> clusterCost = Cluster_costFunction(cnf,qa,ba,dfa);
    // A38: when clusterizePairs kept only this rank's pair rows, it decided their split first.
    const std::vector<int> pairRows(cls->pairRows, cls->pairRows + cls->npairRows);
    const std::vector<double> pairLoads(cls->pairLoads, cls->pairLoads + (cls->pairLoads ? nprocess : 0));
    int*** localclusters = MPI_getLocalClusters(order,clusters,clusterCost,cls->pairsLocal ? 2 : 0,&pairRows,&pairLoads);
    // int*** localclusters = clusters; // No MPI
    // reportClusinfo(localclusters,order);
    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate the dynamics
    if(!clusterPlanOnly && adaptiveWork(cnf).mode!="plan")calculate(qa,ba,dfa,cnf,pls,cls,op,localclusters);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_freeLocalClusters(order,&localclusters);
}

// Conventional CCE has no qubit in its matrices: it propagates the bath under the
// two conditional Hamiltonians, and a pulse can only mean "swap them" - an ideal pi
// flip. Any other rotation used to be calculated as one without a word; gCCE puts
// the qubit in the Hamiltonian and applies the rotation that was asked for.
void validatePulsesForMethod(Config* cnf, Pulse* pls) {
    if(!strcasecmp(cnf->method,"gcce")) return;
    for(int i=0;i<Pulse_getNpulse(pls);++i) {
        const char axis=pls->pulse_axes[i];
        if((axis=='X' || axis=='Y') && pls->pulse_angles[i]==180.0) continue;
        if(rank==0) fprintf(stderr,"Error: method %s calculates every pulse as an ideal pi flip (X or Y, 180), "
                                   "but pulse %d of %s is %c %g. Use method gCCE for other rotations.\n",
                            cnf->method,i+1,pls->pulsename,axis,pls->pulse_angles[i]);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    // The CUDA ensemble path takes CPMG intervals the way the CPU one did.
    if(cnf->nstate==0 && !strcasecmp(pls->pulsename,"Manual") && strcasecmp(cnf->backend,"cpu")) {
        if(!strcasecmp(cnf->backend,"cuda")) {
            if(rank==0) fprintf(stderr,"Error: an ensemble (nstate 0) Manual sequence is calculated on the CPU only; "
                                       "set backend to cpu\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        if(rank==0) fprintf(stdout,"[backend] auto: an ensemble Manual sequence uses CPU\n");
        strcpy(cnf->backend,"cpu");
    }
}

// DEER bath pulses act on bath spins, so only a method with the bath in its matrices can apply them (gCCE), and only
// the CPU paths expand a pulse beyond the qubit (gpu_backend.cu applyPulse takes qdim x qdim pulses).
void validateBathPulses(Config* cnf, Pulse* pls) {
    if (!Pulse_hasBathPulses(pls)) return;
    auto refuse = [](const char* msg) { if (rank==0) fprintf(stderr, "Error: %s\n", msg); MPI_Abort(MPI_COMM_WORLD,1); };
    if (strcasecmp(cnf->method,"gcce")) refuse("bath pulses (pulse.bath.count > 0) are implemented for method gCCE");
    if (strcasecmp(cnf->cluster_mode,"standard")) refuse("bath pulses need cluster_mode standard");   // bundle is "adaptive" inside
    if (pls->detuning_factor != 1.0)    refuse("bath pulses need pulse.detuning_factor = 1");
    if (Pulse_getPulseiter(pls))        refuse("bath pulses cannot be combined with pulse.per_qubit = true");
    if (!strcasecmp(cnf->propagator,"sparse"))
        refuse("bath pulses need propagator eigen or expm (propagator sparse applies pulses as U_q x 1)");
    // gpu_warp_max_dim is -1 until Config_resolveWarpMaxDim: > 0 here is a value given in the input
    if (cnf->gpu_warp_max_dim > 0) refuse("bath pulses cannot use gpu_options.warp_max_dim > 0 (the warp kernel applies qubit pulses only)");
    if (!strcasecmp(cnf->backend,"cuda")) refuse("bath pulses (pulse.bath.count > 0) run on the CPU only; set backend to cpu");
    if (strcasecmp(cnf->backend,"cpu")) {
        if (rank==0) fprintf(stdout,"[backend] auto: bath pulses use CPU\n");
        strcpy(cnf->backend,"cpu");
    }
}

// hamiltonian.cCCE holds the terms of the cCCE Hamiltonian (A47): in a gCCE run - from the input or from -m - they
// would do nothing, so they are refused. Checked once -m has been applied.
//
// delta_A false (A47) replaces the hyperfine tensors of a cluster's bath spins by their mean (createBathArray), which is
// defined for a bath of nuclear spins: the spins a defect species brings into a cluster and the qubit's mean field
// (the Overhauser field, computed from the tensors before they are averaged) are refused with it.
void validateHamiltonianForMethod(Config* cnf, QubitArray* qa, DefectArray* dfa) {
    if(!strcasecmp(cnf->method,"gcce") && cnf->hamiltonian_ccce_given) {
        if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE; this run is gCCE\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    if(!cnf->delta_A && (DefectArray_getNdefect(dfa) > 0 || QubitArray_getOverhaus(qa))) {
        if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE.delta_A false is for a bath of nuclear spins: not with Defect species "
                                   "or hamiltonian.common.qubit_bath_mean_field true\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    setClusterHyperfineMean(!cnf->delta_A);
}

// An "adaptive" object in the input file asks for what --adaptive asks for on the
// command line: its keys are the refinement options with '-' written as '_',
//
//     "adaptive": { "policy": "merge", "output": "./adaptive_out",
//                 "atol": 1e-8, "rtol": 2e-6, "max_rounds": 8 }
//
// (written "Adaptive" until A47; keys are read without regard to case). The block
// belongs to cluster_mode pcce or bundle. A step's input, which the refinement
// driver writes, holds "adaptive": {"work_file": ...} alone: that is the step's
// work request (cJSON_readOptionConfig), not a request for a refinement.
//
// and what comes back is the argument list adaptive_main already takes, so the
// policy code and its own checks are the ones the command line goes through.
// What the run itself decides is refused here rather than overridden: the input
// is this file, the binary is this process, the ranks are mpirun's.
std::vector<std::string> adaptiveArgsFromInput(char* fccein, const char* self) {
    std::vector<std::string> args;
    if(access(fccein, R_OK) != 0) return args;          // the ordinary read reports it
    char* text = cJSON_ReadFccein(fccein);
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> root(cJSON_Parse(text),cJSON_Delete);
    freeChar1d(&text);
    cJSON* block = root ? cJSON_GetObjectItem(root.get(),"adaptive") : nullptr;
    const auto fail = [](const std::string& message) {
        if(rank==0) fprintf(stderr,"Error: %s\n",message.c_str());
        MPI_Finalize(); exit(EXIT_FAILURE);
    };
    // The block was called "Refine" until 2026-09-23. The old name is refused, not read:
    // a run that silently did no refinement would look like a finished one.
    if(root && cJSON_GetObjectItem(root.get(),"Refine"))
        fail("the \"Refine\" block is called \"adaptive\" now; rename it (docs/methods/adaptive-input-block.md)");
    if(!block) return args;
    if(!cJSON_IsObject(block)) fail("adaptive must be an object of refinement options");
    if(cJSON_GetObjectItem(block,"work_file")) {
        if(block->child->next || strcasecmp(block->child->string,"work_file"))
            fail("adaptive.work_file is written by the refinement driver into a step's input, beside no other adaptive key");
        return args;
    }
    cJSON* mode = cJSON_GetObjectItem(root.get(),"cluster_mode");
    const char* shape = cJSON_IsString(mode) ? mode->valuestring : "standard";
    if(strcasecmp(shape,"pcce") && strcasecmp(shape,"bundle"))
        fail(std::string("the adaptive block belongs to cluster_mode pcce or bundle (this input has cluster_mode ")+shape+")");
    static const char* fixed[] = {"input","binary","cwd","mpiexec","ranks"};
    const auto value_text = [&](cJSON* v, const std::string& key) -> std::string {
        if(cJSON_IsString(v)) return v->valuestring;
        if(cJSON_IsNumber(v)) {
            char buffer[64];
            const double x = v->valuedouble;
            if(x == std::floor(x) && std::fabs(x) < 1e15) snprintf(buffer,sizeof buffer,"%lld",(long long)x);
            else snprintf(buffer,sizeof buffer,"%.17g",x);
            return buffer;
        }
        fail("adaptive."+key+" must be a number or a string");
        return "";
    };
    for(cJSON* item=block->child; item; item=item->next) {
        const std::string key = item->string;
        std::string name = key;
        for(char& ch : name) if(ch=='_') ch='-';
        const adaptive::options::Option* option = nullptr;
        for(const adaptive::options::Option& one : adaptive::options::kOptions)
            if(name == one.name) option = &one;
        if(!option) fail("adaptive."+key+" is not a refinement option");
        for(const char* f : fixed)
            if(name == f) fail("adaptive."+key+" is decided by the run itself and cannot be set in the input");
        if(cJSON_IsNull(item)) continue;
        if(option->kind == adaptive::options::Kind::Flag) {
            if(!cJSON_IsBool(item)) fail("adaptive."+key+" is a switch: write true or false");
            if(cJSON_IsTrue(item)) args.push_back("--"+name);
        } else if(cJSON_IsArray(item)) {
            fail("adaptive."+key+" takes one value, not a list");
        } else {
            args.push_back("--"+name); args.push_back(value_text(item,key));
        }
    }
    for(const adaptive::options::Option& one : adaptive::options::kOptions) {
        bool is_fixed = false;
        for(const char* f : fixed) if(std::string(one.name) == f) is_fixed = true;
        if(!one.required || is_fixed) continue;
        bool given = false;
        for(const std::string& arg : args) if(arg == std::string("--")+one.name) given = true;
        if(!given) fail(std::string("adaptive.")+one.dest+" is required");
    }
    args.insert(args.begin(), {"--input", fccein, "--binary", self});
    return args;
}

// The settings a refinement step can change, and nothing else: merge writes the
// user's input with outfile and adaptive.work_file set per step and savemode forced
// to coherence. Nothing the bath depends on ever moves, which is why the bath can be
// read once and kept (plan stage 8). (order, bundle_merge_cutoff, npulse and the
// connection pair were changed only by the policies removed on 2026-10-10; a step
// that changes them now is refused like any other.)
//
// Anything else differing is a case this does not handle, and it says so rather
// than calculating with a bath that no longer matches the request.
void applyStepConfig(const std::string& stepInput, const std::string& firstInput,
                     Config* cnf, Output* op) {
    // savemode is forced into every step input whatever the user wrote. It does not
    // touch the bath that was read once: it is an output layout. Without it here a
    // user input that did not already spell out savemode was refused at the first
    // step (measured 2026-09-19; the recorded inputs all spell it).
    static const char* mutable_keys[] = {"outfile", "adaptive", "savemode"};   // adaptive: {"work_file"} in a step
    // Both files are read as CCEX inputs (a line ends at its first ! or #): the
    // first one is the user's, and may carry comments; the step's is what the
    // policy wrote, and carries none.
    char* firstText = cJSON_ReadFccein(const_cast<char*>(firstInput.c_str()));
    char* stepText = cJSON_ReadFccein(const_cast<char*>(stepInput.c_str()));
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> a(cJSON_Parse(firstText),cJSON_Delete);
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> b(cJSON_Parse(stepText),cJSON_Delete);
    freeChar1d(&firstText); freeChar1d(&stepText);
    if(!a || !b) { fprintf(stderr,"Error: refinement step input is not JSON\n"); exit(EXIT_FAILURE); }
    for(cJSON* item=b->child; item; item=item->next) {
        bool allowed=false;
        for(const char* k : mutable_keys) if(!strcmp(item->string,k)) allowed=true;
        if(allowed) continue;
        cJSON* was=cJSON_GetObjectItem(a.get(),item->string);
        char* x=cJSON_PrintUnformatted(item); char* y=was?cJSON_PrintUnformatted(was):nullptr;
        const bool same = y && !strcmp(x,y);
        free(x); if(y) free(y);
        if(!same) {
            fprintf(stderr,"Error: refinement step changed '%s', which the kept bath depends on\n",
                    item->string);
            exit(EXIT_FAILURE);
        }
    }
    cJSON* v;
    cJSON* stepBlock=cJSON_GetObjectItem(b.get(),"adaptive");
    if(stepBlock && (v=cJSON_GetObjectItem(stepBlock,"work_file")) && cJSON_IsString(v)) {
        strncpy(cnf->adaptive_work_file,v->valuestring,MAX_FILEPATH-1);
        cnf->adaptive_work_file[MAX_FILEPATH-1]=0;
    } else cnf->adaptive_work_file[0]=0;
    if((v=cJSON_GetObjectItem(b.get(),"outfile")) && cJSON_IsString(v)) {
        // A top input without "outfile" (the name then comes from -o, or nothing)
        // leaves the buffer unallocated; the step always names one.
        if(!Output_getOutfile(op)) Output_allocOutfile(op);
        Output_setOutfile(op,v->valuestring);
    }
    if((v=cJSON_GetObjectItem(b.get(),"savemode")))   // the mode; nodiv stays the user input's
        cJSON_readSavemode(v,op,false);
}
}   // namespace

// -N and -B as strictly as the input file reads nstate and bfield (review 78d583 R5): atoi/atof took "nope" as 0
// (-N nope ran an ensemble) and "1.5" as 1.
static int commandLineInt(const char* text, const char* option){
    errno = 0; char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || errno == ERANGE || value < INT_MIN || value > INT_MAX){
        fprintf(stderr, "Error: %s must be an integer, not \"%s\"\n", option, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}
static double commandLineNumber(const char* text, const char* option){
    errno = 0; char* end = nullptr;
    const double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || errno == ERANGE || !std::isfinite(value)){
        fprintf(stderr, "Error: %s must be a finite number, not \"%s\"\n", option, text);
        exit(EXIT_FAILURE);
    }
    return value;
}

// main.out --auto_conv -f input.json ((internal notes), "--auto_conv 첫 판 설계 v2"): before MPI starts, hand the
// input to the runner of example/tools/AutoConvergence (python -m autoconv run), which launches the ordinary calculations
// itself. Under mpirun -n N, rank 0 goes on alone and the runner launches with mpirun -n N when the input gives no launcher
// (plan "mpirun으로 부르면 그 rank 수로", 2026-10-10). The tool is CCEX_AUTOCONV_DIR or the first parent of this binary holding
// example/tools/AutoConvergence; the interpreter is CCEX_PYTHON or python3.
// The input of -f FILE or -fFILE (getopt's two spellings).
static const char* inputArgument(int argc, char* argv[]){
    const char* input = nullptr;
    for (int i = 1; i < argc; ++i){
        if (!strcmp(argv[i], "-f") && i + 1 < argc) input = argv[++i];
        else if (!strncmp(argv[i], "-f", 2) && argv[i][2] != '\0') input = argv[i] + 2;
    }
    return input;
}

static int autoConvEntry(int argc, char* argv[], const char* mode = "run"){
    // Under mpirun -n N (plan "mpirun으로 부르면 그 rank 수로", 2026-10-10): every rank reads N and leaves MPI; rank 0 alone
    // goes on, with the outer launch's MPI variables dropped - Open MPI refuses an mpirun started from inside one otherwise
    // ("mpirun does not support recursive calls", idm4-gpu job 672) - and passes N to the runner as CCEX_AUTOCONV_NP.
    // First, so a refusal below is printed once, by rank 0.
    bool launched = false;
    for (const char* v : {"PMI_RANK", "PMI_SIZE", "OMPI_COMM_WORLD_RANK", "MPI_LOCALRANKID", "PMIX_RANK"}) launched = launched || getenv(v);
    if (launched){
        int r = 0, size = 1;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &r);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Finalize();
        if (r != 0) return 0;
        std::vector<std::string> drop;
        for (char** e = environ; *e; ++e){
            const std::string kv = *e, name = kv.substr(0, kv.find('='));
            if (name == "I_MPI_ROOT" || name == "OPAL_PREFIX") continue;      // where the MPI is installed: kept
            for (const char* p : {"PMI_", "PMIX_", "HYDI_", "HYDRA_", "OMPI_", "PRTE_", "I_MPI_", "OPAL_"})
                if (!name.compare(0, strlen(p), p)){ drop.push_back(name); break; }
            if (name == "MPI_LOCALRANKID" || name == "MPI_LOCALNRANKS") drop.push_back(name);
        }
        for (const std::string& name : drop) unsetenv(name.c_str());
        setenv("CCEX_AUTOCONV_NP", std::to_string(size).c_str(), 1);
    }
    // the runner reads only the input: an option such as -P (plan only) or -B (field) would be dropped without a word, so it
    // is refused (review 20261008-autoconv-autostart #1)
    for (int i = 1; i < argc; ++i){
        if (!strcmp(argv[i], "--auto_conv") || !strcmp(argv[i], "--auto_status")) continue;
        if (!strcmp(argv[i], "-f") && i + 1 < argc){ ++i; continue; }
        if (!strncmp(argv[i], "-f", 2) && argv[i][2] != '\0') continue;
        fprintf(stderr, "Error: the automatic convergence takes only -f FILE; \"%s\" is not used by it (set it in the input, or set "
                        "auto_convergence.enabled false for an ordinary run)\n", argv[i]);
        return 2;
    }
    const char* input = inputArgument(argc, argv);
    if (!input){ fprintf(stderr, "Error: --auto_conv needs -f input.json\n"); return 2; }
    // The binary's own path from the system, never from argv[0] (whoever starts it chooses that, and it decides which owner
    // is trusted below; commit security review 2026-10-09): /proc/self/exe on Linux, _NSGetExecutablePath on macOS
    char exe[PATH_MAX];
    const ssize_t n = readlink("/proc/self/exe", exe, sizeof exe - 1);
    if (n > 0) exe[n] = '\0';
    else {
#ifdef __APPLE__
        char raw[PATH_MAX]; uint32_t size = sizeof raw;
        if (_NSGetExecutablePath(raw, &size) != 0 || !realpath(raw, exe))
#endif
        { fprintf(stderr, "Error: --auto_conv cannot find its own path from the system\n"); return 1; }
    }
    // Only the checkout the binary was built in: bin/main.out (repository = two levels up) or build/bin/main.out (three).
    // What it puts first on PYTHONPATH must belong to this user, root or the owner of this binary (a shared installation,
    // user 2026-10-09) and be writable by nobody else - the folders, every .py file and an existing __pycache__ - so a
    // directory or file planted next to a binary in a shared location cannot run code (security review of 5e3b9d439).
    // The folders above it too, from its real path up to /: whoever can write one of them can put another folder in its
    // place between this check and Python's import (commit security review, TOCTOU, 2026-10-09).
    std::string tool = getenv("CCEX_AUTOCONV_DIR") ? getenv("CCEX_AUTOCONV_DIR") : "";
    std::string dir = exe;
    for (int up = 1; tool.empty() && up <= 3 && dir.rfind('/') != std::string::npos && dir.rfind('/') > 0; ++up){
        dir.erase(dir.rfind('/'));
        const std::string candidate = dir + "/example/tools/AutoConvergence";
        if (up >= 2 && access((candidate + "/autoconv/runner.py").c_str(), R_OK) == 0) tool = candidate;
    }
    if (tool.empty()){ fprintf(stderr, "Error: --auto_conv cannot find example/tools/AutoConvergence (set CCEX_AUTOCONV_DIR)\n"); return 1; }
    char real[PATH_MAX];
    if (!realpath(tool.c_str(), real)){ fprintf(stderr, "Error: --auto_conv cannot resolve %s: %s\n", tool.c_str(), strerror(errno)); return 1; }
    tool = real;
    struct stat exe_st;
    const uid_t installer = stat(exe, &exe_st) == 0 ? exe_st.st_uid : getuid();
    // writable by someone else: others without the sticky bit, or a group that is not the owner's own private group (no
    // other members; idm makes every folder group-writable with umask 0002 and such a group)
    auto othersWrite = [](const struct stat& st, bool above){
        if (above && S_ISDIR(st.st_mode) && (st.st_mode & S_ISVTX)) return false;   // /tmp-like folders above the tool only
        if (st.st_mode & S_IWOTH) return true;
        if (!(st.st_mode & S_IWGRP)) return false;
        const passwd* owner = getpwuid(st.st_uid);
        if (!owner || owner->pw_gid != st.st_gid) return true;
        const group* g = getgrgid(st.st_gid);
        if (!g || (g->gr_mem != nullptr && g->gr_mem[0] != nullptr)) return true;      // supplementary members
        bool shared = false;                          // another account with it as its primary group (not in gr_mem; review)
        const uid_t ownerUid = st.st_uid;
        setpwent();
        while (const passwd* pw = getpwent())
            if (pw->pw_gid == st.st_gid && pw->pw_uid != ownerUid){ shared = true; break; }
        endpwent();
        return shared;
    };
    auto trusted = [&](uid_t u){ return u == getuid() || u == 0 || u == installer; };
    std::vector<std::string> paths = {tool, tool + "/autoconv"};
    for (const std::string& folder : {tool + "/autoconv", tool + "/autoconv/__pycache__"}){
        DIR* d = opendir(folder.c_str());
        if (!d) continue;
        if (folder != paths[1]) paths.push_back(folder);
        while (const dirent* e = readdir(d)){
            const std::string name = e->d_name;
            const bool code = (name.size() > 3 && name.compare(name.size() - 3, 3, ".py") == 0) ||
                              (name.size() > 4 && name.compare(name.size() - 4, 4, ".pyc") == 0);
            if (code) paths.push_back(folder + "/" + name);
        }
        closedir(d);
    }
    const size_t own = paths.size();                    // the tool's own folders and files; the folders above follow
    for (std::string up = tool; up.rfind('/') != std::string::npos; ){
        up.erase(up.rfind('/'));
        paths.push_back(up.empty() ? "/" : up);
        if (up.empty()) break;
    }
    for (size_t i = 0; i < paths.size(); ++i){
        const std::string& p = paths[i];
        struct stat st;
        // inside the tool (its real path) nothing may be a symbolic link: the folders above a link's target are not
        // checked (review 6548a9 R2)
        if (i > 0 && i < own && (lstat(p.c_str(), &st) != 0 || S_ISLNK(st.st_mode))){
            fprintf(stderr, "Error: the automatic convergence refuses %s: a symbolic link inside the tool folder (copy the files instead)\n", p.c_str());
            return 1;
        }
        if (stat(p.c_str(), &st) != 0 || !trusted(st.st_uid) || othersWrite(st, i >= own)){
            fprintf(stderr, "Error: the automatic convergence refuses %s: not owned by this user, root or the owner of %s, or writable "
                            "by others (fix it with chmod go-w %s)\n", p.c_str(), exe, p.c_str());
            return 1;
        }
    }
    setenv("PYTHONDONTWRITEBYTECODE", "1", 1);   // no bytecode written into the tool folder (shared or not)
    const char* old = getenv("PYTHONPATH");
    setenv("PYTHONPATH", (old && *old ? tool + ":" + old : tool).c_str(), 1);
    const char* python = getenv("CCEX_PYTHON") ? getenv("CCEX_PYTHON") : "python3";
    printf("[auto] runner %s/autoconv (python %s), binary %s\n", tool.c_str(), python, exe);
    fflush(stdout);
    char* args[] = {const_cast<char*>(python), const_cast<char*>("-m"), const_cast<char*>("autoconv"), const_cast<char*>(mode),
                    const_cast<char*>(input), const_cast<char*>("--binary"), exe, nullptr};
    if (!strcmp(mode, "status")) args[5] = nullptr;          // status reads the workdir only
    execvp(python, args);
    fprintf(stderr, "Error: --auto_conv cannot start %s: %s\n", python, strerror(errno));
    return 1;
}

// main.out -f input.json runs the automatic convergence when the input holds an auto_convergence block whose "enabled"
// is true or absent (user 2026-10-08); "enabled": false runs the ordinary calculation and the block is not read. Decided
// before MPI starts: the runner launches the calculations itself (under mpirun, rank 0 alone; autoConvEntry).
static bool inputAsksAutoConv(int argc, char* argv[]){
    const char* input = inputArgument(argc, argv);
    if (!input || access(input, R_OK) != 0) return false;     // the ordinary path reports a missing file
    char* data = cJSON_ReadFccein(const_cast<char*>(input));
    cJSON* root = cJSON_Parse(data);
    freeChar1d(&data);
    bool yes = false;
    if (root){
        cJSON* block = cJSON_GetObjectItem(root, "auto_convergence");
        if (cJSON_IsObject(block)){
            cJSON* on = cJSON_GetObjectItem(block, "enabled");
            yes = on == nullptr || cJSON_IsTrue(on);
            if (on != nullptr && !cJSON_IsBool(on)){
                fprintf(stderr, "Error: auto_convergence.enabled must be true or false\n");
                cJSON_Delete(root);
                exit(2);
            }
        }
        cJSON_Delete(root);
    }
    return yes;
}

int main(int argc, char* argv[]){

    // main.out --auto_status -f input.json: where an automatic convergence is (running / interrupted / stopped / ended), no
    // calculation (plan "끊김 대비·진행 보기" R4)
    for (int i = 1; i < argc; ++i) if (!strcmp(argv[i], "--auto_status")) return autoConvEntry(argc, argv, "status");
    for (int i = 1; i < argc; ++i) if (!strcmp(argv[i], "--auto_conv")) return autoConvEntry(argc, argv);
    if (inputAsksAutoConv(argc, argv)) return autoConvEntry(argc, argv);

    //Set up MPI. FUNNELED: the packed GPU batch (gpu_options.packed_max_dim) runs threads beside the main one, which alone
    //calls MPI; smallClusterWarmUp refuses that batch, before any thread starts, if the MPI provides less (review
    //20261012-233058-452278 P3-03). Nothing else depends on the level provided.
    int threadLevel = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadLevel);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process id
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess); // Get number of processes

    // wall time check
    double time_start, time_end;
    time_start = MPI_Wtime();

    // Print welcome message
    if (rank==0){
        printBanner();
    }

    //=======================================================
    // Variables declaration
    //=======================================================

    // Input file (cce.in)
    char* fccein = NULL; 

    //Simulation objects
    Config*      cnf = Config_init();
    QubitArray*  qa  = QubitArray_init();
    Cluster*     cls = Cluster_init();
    BathArray*   ba  = BathArray_init();
    DefectArray* dfa = DefectArray_init();
    Pulse*       pls = Pulse_init();
    Output*      op  = Output_init();

    //=======================================================
    // Read options
    //=======================================================
    int c;
    //int nbathfiles_fromfccein = 0;
    int nbathfiles_current = 0;

    //  Declare and initialize the variables 
    // that will be used to store the arguments
    char method[MAX_CHARARRAY_LENGTH] = "";
    char quantity[MAX_CHARARRAY_LENGTH] = "";
    char** bathfiles = NULL;
    int nstate = 0;
    double bfield_z = 0.0;
    char* outfile = NULL;
    bool clusterPlanOnly = false;

    // --adaptive takes everything after it and hands it to the refinement tool's
    // own parser, so getopt below never sees those words (plan stage 8).
    std::vector<std::string> adaptiveArgs;
    int adaptiveStatus = 0;
    std::vector<char*> kept;          // argv points into this; it must outlive getopt
    {
        for (int i = 0; i < argc; ++i) {
            if (adaptiveArgs.empty() && i > 0 && !strcmp(argv[i], "--adaptive")) {
                for (int k = i+1; k < argc; ++k) adaptiveArgs.push_back(argv[k]);
                // The tool's --input is the CCEX input this process reads once.
                for (std::size_t k = 0; k+1 < adaptiveArgs.size(); ++k)
                    if (adaptiveArgs[k] == "--input") {
                        kept.push_back(const_cast<char*>("-f"));
                        kept.push_back(const_cast<char*>(adaptiveArgs[k+1].c_str()));
                    }
                break;
            }
            kept.push_back(argv[i]);
        }
        if (!adaptiveArgs.empty()) {
            argc = int(kept.size()); argv = kept.data(); optind = 1;
        }
    }
    // The same request, written in the input file instead (adaptiveArgsFromInput).
    // argv is left alone: -f is there already, and getopt reads it as it always has.
    for (int i = 1; i+1 < argc; ++i) {
        if (strcmp(argv[i], "-f")) continue;
        char self[4096] = {};
        if (readlink("/proc/self/exe", self, sizeof self-1) <= 0)
            snprintf(self, sizeof self, "%s", argv[0]);
        const std::vector<std::string> fromInput = adaptiveArgsFromInput(argv[i+1], self);
        if (fromInput.empty()) break;
        if (!adaptiveArgs.empty()) {
            if (rank == 0) fprintf(stderr,"Error: the input has an adaptive block and the command line has --adaptive; use one\n");
            MPI_Finalize(); exit(EXIT_FAILURE);
        }
        adaptiveArgs = fromInput;
        break;
    }

    // Two passes over the options: the input file (-f) first, then every other option over it, so that the
    // command line wins whatever the order (review b9eb09 R1: "-N 2 -f in.json" used to take nstate from the file,
    // "-f in.json -N 2" from the command line). The first pass reads -f only and stays quiet about the rest.
    // --help is -h (review 20261013-doxygen D21); getopt reads short options only
    for (int i = 1; i < argc; ++i) if (!strcmp(argv[i], "--help")) argv[i] = const_cast<char*>("-h");
    for (int pass = 0; pass < 2; ++pass) {
    optind = 1; opterr = pass;
    while ((c = getopt(argc, argv, "hPvf:m:q:I:s:N:S:a:B:o:")) != -1) {
        if ((c == 'f') != (pass == 0)) continue;
        switch (c) {
            case 'f':

                // set the parameter as the input file
                if (rank==0){ 
                    // print message
                    char message[500];
                    sprintf(message, "Option filename : '%s' \n", optarg);
                    printMessage(message);
                }
                fccein = optarg;

                // check if the file exists
                if (access(fccein, F_OK) == -1){
                    fprintf(stderr, "Error: Option file does NOT exist\n\n");
                    exit(EXIT_FAILURE);
                }

                // read the input file
                cJSON_readOptionConfig      (cnf, fccein); // general.h
                cJSON_readOptionQubitArray  (qa,  fccein); // qubit.h
                cJSON_readOptionDefectArray (dfa, fccein); // defect.h
                cJSON_readOptionCluster     (cls, fccein); // cluster.h
                cJSON_readOptionPulse       (pls, fccein); // pulse.h
                cJSON_readOptionOutput      (op,  fccein); // output.h
                cJSON_checkInputKeys        (fccein);      // a key nobody read: a typo stops the run (D24)

				if (rank==0){
					printf("\n	>> Read %s file successfully ..\n",optarg);
				}

				//nbathfiles_fromfccein = Config_getNbathfiles(cnf);
                break;

            case 'm':
                cJSON_applyMethodOption(cnf,cls,optarg);   // cCCE or gCCE; the input's cluster_mode stays (A47)
                break;

            case 'q':
                Config_setQuantity(cnf,optarg);
                break;

            case 'I':
                nbathfiles_current++;
                if (nbathfiles_current == 1){
                    Config_freeBathfiles(cnf);
                    Config_freeBathadjust(cnf);
                    Config_setNbathfiles(cnf,nbathfiles_current);
                    Config_allocBathfiles(cnf);
                    Config_allocBathadjust(cnf);
                    Config_setBathfiles_i(cnf,optarg,nbathfiles_current-1);
                }
                else{
                    Config_setNbathfiles(cnf,nbathfiles_current);
                    Config_reallocBathfiles(cnf, nbathfiles_current-1, nbathfiles_current);
                    Config_reallocBathadjust(cnf, nbathfiles_current-1, nbathfiles_current);
                    Config_setBathfiles_i(cnf,optarg,nbathfiles_current-1);
                }
                break;

            case 's':
                Config_freeStatefile(cnf);
                Config_allocStatefile(cnf);
                Config_setStatefile(cnf,optarg);
                break;

            case 'a':
                Config_freeAvaaxfile(cnf);
                Config_allocAvaaxfile(cnf);
                Config_setAvaaxfile(cnf,optarg);
                break;

            case 'S':
                Config_freeExstatefile(cnf);
                Config_allocExstatefile(cnf);
                Config_setExstatefile(cnf,optarg);
                break;

            case 'N':
                Config_setNstate(cnf,commandLineInt(optarg,"-N"));
                break;

            case 'B':
                // Bz in the computational frame, which with the rotation on is the qubit
                // axis -- so this composes with coordinate_frame_rotation rather than
                // conflicting with it. Config_validateBfieldAlignment checks the result.
                Config_setBfield_z(cnf,commandLineNumber(optarg,"-B"));
                break;

            case 'o':
                Output_freeOutfile(op);
                Output_allocOutfile(op);
                Output_setOutfile(op,optarg);
                break;

            case 'v':
                verbosity = true;
                break;

            case 'P':
                clusterPlanOnly = true;
                break;

            case 'h':
                printf("Usage: %s [options]\n", argv[0]);
                printf("Options (short only; given with -f they override the input file, whatever their order):\n");
                printf("  -h, --help\t\t\t\tPrint this message\n");
                printf("  -f FILE\t\t\t\tInput file (cce.in)\n");
                printf("  -I FILE\t\t\t\tInput file name (it can be more than one)\n");
                printf("  -N N\t\t\t\t\tNumber of bath states\n");
                printf("  -m METHOD\t\t\t\tMethod (cCCE | gCCE; meCCE and noise reserved, not implemented yet); the cluster shape is the input's cluster_mode\n");
                printf("  -q QUANTITY\t\t\t\tQuantity (coherence | dm)\n");
                printf("  -B BZ\t\t\t\t\tMagnetic field\n");
                printf("  -v\t\t\t\t\tLong or short print\n");
                printf("  -P\t\t\t\t\tSave <outfile>.clusters.json and stop before dynamics\n");
                printf("  -o FILE\t\t\t\tOutput file name\n");
                printf("  --auto_conv -f FILE\t\t\tChoose the window, rdip, rbath and sampling by itself; the same as -f FILE with an\n\t\t\t\t\tauto_convergence block (example/tools/AutoConvergence/README.md)\n");
                printf("  --auto_status -f FILE\t\t\tWhere that automatic convergence is (running, interrupted, stopped, ended), no calculation\n");
                printf("JSON execution options (in the -f input):\n");
                printf("  backend: cpu (default), cuda, auto; profile: true for per-dimension counts/timing\n");
                printf("  gpu_options: {min_dim, packed_max_dim (0, 2, 4), warp_max_dim (0-32), evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb}\n");
                printf("  propagator: eigen (default), sparse (one state; cCCE, gCCE vector), expm (gCCE)\n  evolution, gcce_refinement {enabled, max_sweeps} (gCCE)\n");
                printf("  evolution_update: exact or recursive (cCCE on the CPU and in the packed kernel, gCCE with propagator=eigen on the CPU); default recursive there, exact elsewhere (GPU clusters)\n");
                printf("  probe_memory_gb (memory of the adaptive probes running at once on a node)\n");
                printf("  seed: shared RNG seed; omit to generate and report a seed on MPI rank 0\n");
                printf("  hamiltonian.common.bath_bath_mean_field: off, z (default) or xyz, the static field of the other bath spins (nstate > 0)\n");
                printf("  cluster_mode: standard (default), pcce, bundle, dsj, itb, dsjitb; the pcce and bundle blocks hold their settings\n");
                printf("  input keys renamed 2026-10-13: example/tools/InputKeys/convert_input_keys.py rewrites an old input\n");
                printf("See docs/user/InputFile_Backend.md, docs/methods/gpu.md and docs/methods/adaptive-cluster-mode.md for defaults and limits.\n");
                exit(EXIT_SUCCESS);
                break;
            default:
                fprintf(stderr, "Error: unknown option\n");
                exit(EXIT_FAILURE);
        }   
    }
    }
    
    if (fccein == NULL){
        fprintf(stderr, "Error: input file is required\n");
        exit(EXIT_FAILURE);
    }

    // Every option source has been applied by here -- the input file (read from case 'f'
    // above) and the flags that can override it, -m and -N in particular. Only now can
    // the defaulted evolution be settled; doing it during the file parse would decide it
    // against an nstate that -N had not yet supplied.
    Pulse_resolveBathSequence(pls);   // after -f and every flag: the merged sequence is what every later check sees
    validateBathPulses(cnf,pls);      // before the packed/warp defaults: they read backend
    Config_resolveEvolution(cnf);
    Config_resolvePackedMaxDim(cnf);
    Config_resolveWarpMaxDim(cnf);             // A40, before the evolution_update check (its log names the warp kernel)
    Config_validatePropagator(cnf,QubitArray_getNqubit(qa));   // A2: before the evolution_update check, so sparse is named
    Config_validateEvolutionUpdate(cnf);       // after the packed kernel's value (A36)
    validatePulsesForMethod(cnf,pls);
    validateHamiltonianForMethod(cnf,qa,dfa);
    validateAdaptiveConfig(cnf,cls);
    Config_validatePcceCache(cnf);
    if (clusterPlanOnly && ((strcasecmp(cnf->method,"cce") &&
         strcasecmp(cnf->method,"gcce") && strcasecmp(cnf->method,"pcce")) ||
         !strcasecmp(cnf->bundle_connections,"checked") || adaptiveWork(cnf).mode=="probe" ||
         adaptiveWork(cnf).mode=="factors" || adaptiveWork(cnf).mode=="except")) {
        if(rank==0)fprintf(stderr,"Error: -P supports cCCE, gCCE and cluster_mode pcce plans with explicit connections; rdip-omitted checked connections and probe, factors and except work require dynamics\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    computeBackendInit(cnf);

    // Same reason: -B lands after the file is parsed, so the field can only be checked
    // against the rotation once both have been applied.
    Config_validateBfieldAlignment(cnf);

    // Needs the Config AND the QubitArray, which two different parsers build, so it can
    // only run once both are done. Still before every reader, so an unsupported
    // combination stops the run before a tensor file is opened.
    validateCoordinateFrameRotationInputs(qa,cnf);

    MPI_Barrier(MPI_COMM_WORLD);

    //if (nbathfiles_current != nbathfiles_fromfccein){
    //    fprintf(stderr, "Error: number of bath files is not consistent\n");
    //    fprintf(stderr, "nbathfiles_fromfccein = %d\n", nbathfiles_fromfccein);
    //    fprintf(stderr, "nbathfiles_current = %d\n", nbathfiles_current);
    //    exit(EXIT_FAILURE);
    //}

    //=======================================================
    // Report options
    //=======================================================
    if (rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("        Report all structures \n");
        printf("    ======================================================================\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Config_report(cnf); 
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        QubitArray_report(qa);
        printf("\n\n");

        if (DefectArray_getNdefect(dfa) > 0){
            printf("    ----------------------------------------------------------------------\n");
            printTitle("DefectArray");
            for (int idf=0; idf<DefectArray_getNdefect(dfa); idf++){
                DefectArray_reportDefect_idf(dfa,idf);
                printf("\n");
            }
        }

        printf("    ----------------------------------------------------------------------\n");
        Cluster_report(cls);
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Pulse_report(pls);
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Output_report(op);
        printf("\n");

        printLineSection();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //=======================================================
    // Set simulators
    //=======================================================

    // Seed once, from the resolved MPI-shared configuration, before any random
    // defect axes, legacy partition trials or sampled bath states are generated.
    srand(Config_getSeed(cnf));

    if (rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("        Read files \n");
        printf("    ======================================================================\n\n");
    }

    // Qubit
    if (rank == 0){
        printSubTitle("Qubit file");
    }
    readQubitfile(qa,cnf); // set qubit xyz position for nqubit=1
    MPI_Barrier(MPI_COMM_WORLD);

    // Bath : Gyromagnetic ratio
    if (rank == 0){
        printSubTitle("Gyromagnetic ratio file");
    }
    readGyrofile(ba,cnf); // set gyro, spin, name (properties)
    MPI_Barrier(MPI_COMM_WORLD);

    // Bath : Configuration files
    if (rank == 0){
        printSubTitle("Bath configuration files");
    }
    readBathfiles(ba,qa,cnf);  // set bath xyz position and properties
    MPI_Barrier(MPI_COMM_WORLD);

    // Hyperfine tensor. hypf_bathframe records, per (bath spin, qubit), whether the
    // tensor stored was computed from the source geometry or read from the file; the
    // rotation below needs it because a single run can hold both.
    HypfProvenance hypf_bathframe = NULL;
    readHftensorfile(ba,qa,cnf,&hypf_bathframe); // set hyperfine tensor only from file
    MPI_Barrier(MPI_COMM_WORLD);

    // Quadrupole tensor
    readQdtensorfile(ba,qa,cnf); // set quadrupole tensor only from file
    MPI_Barrier(MPI_COMM_WORLD);

    // Coordinate frame rotation, and the reason it sits HERE.
    //
    // The tensor files are keyed by POSITION, and those position columns are written in
    // the source bath frame, so every lookup above has to run against un-moved
    // coordinates. Tensor COMPONENTS are a separate question: they are in whichever basis
    // hf_tensor_frame / qd_tensor_frame declares, and a point-dipole tensor CCEX computed
    // itself is in the bath frame no matter what those say.
    //
    // This step normalizes the two into one: every qubit position and every bath-spin
    // position is rotated by the same R about the same r0, and the tensors that are still
    // in the bath frame -- HF, QD, Defect templates, and Qubit.intmap entries whose
    // provenance says so -- are transformed with them. Tensors already in the
    // computational frame are left untouched. The magnetic field is NOT settled here;
    // that happened right after the CLI was parsed, above.
    applyCoordinateFrameRotation(ba,qa,dfa,cnf,hypf_bathframe);
    freeInt2d(&hypf_bathframe,BathArray_getNspin(ba));
    MPI_Barrier(MPI_COMM_WORLD);

    // Defect
    if (DefectArray_getNdefect(dfa) > 0){

        if (rank == 0){
            printSubTitle("Defect Paxes File");
        }

        int nspin = BathArray_getNspin(ba);
        DefectArray_allocPaxes(dfa, nspin);
        setDefectPaxes(dfa,ba,cnf); // Read or random set

        if (rank==0){
            DefectArray_reportPaxes(dfa);
        }

        DefectArray_allocNaddspins(dfa, nspin);
        DefectArray_setNaddspins(dfa,ba); // Naddspins
        if (rank==0){
            DefectArray_reportNaddspins(dfa);
        }

        int nqubit = QubitArray_getNqubit(qa);
        DefectArray_allocSubbath(dfa,ba,nqubit);
        DefectArray_setSubbath(dfa,ba,qa);
        if (rank==0){
            DefectArray_reportSubbath_props(dfa);
            DefectArray_reportSubbath_hypfs(dfa, nqubit);
            DefectArray_reportSubbath_quads(dfa);
            DefectArray_reportSubbath_hypf_subs(dfa);
        }

        // Update main bathspins : mainspidx, zfs(pax), detuning(pax)
        updateMainSpins_fromDefectArray(dfa,ba);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    validateBathPulseTargets(ba,dfa,pls);   // DEER: the driven spins exist and fit pulse.bath (bath and detunings set)

    if (rank==0){
        printf("    ----------------------------------------------------------------------\n");
        printTitle("QubitArray");
        QubitArray_report(qa);
        printf("\n");

        printf("    ----------------------------------------------------------------------\n");
        printTitle("BathArray");
        BathArray_report(ba);
        printf("\n");

    }

    if (rank==0){
        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("          Wall time = %.5f s\n", time_end - time_start);
        printLine();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (adaptiveArgs.empty()) {
        runOneCalculation(cnf,qa,ba,dfa,pls,cls,op,clusterPlanOnly);
    } else {
        // The refinement loop, in this process (plan stage 8). Rank 0 runs the
        // policy - it reads the clock and the filesystem, so only one rank may -
        // and every rank calculates, because the calculation is collective.
        //
        // The request still travels as files: rank 0 writes the step's
        // request.json and input.json exactly as the standalone tool does, and
        // every rank reads that input.json. Keeping the protocol keeps the
        // recordings and their replay, which is the only check this has.
        const std::string firstInput = fccein ? fccein : "";
        // The steps run on this process's own ranks, so --ranks is not the user's
        // to get wrong: left out, it is the MPI size; given, it has to agree. It
        // used to default to 1 whatever mpirun was told, and a run on two ranks
        // died at its first probe with "missing MPI probe response" (2026-09-19).
        bool ranksGiven = false;
        for (std::size_t k = 0; k+1 < adaptiveArgs.size(); ++k)
            if (adaptiveArgs[k] == "--ranks") {
                ranksGiven = true;
                if (atoi(adaptiveArgs[k+1].c_str()) != nprocess) {
                    if (rank == 0)
                        fprintf(stderr,"Error: --ranks %s, but this process has %d MPI rank(s); "
                                       "leave --ranks out and it follows mpirun -n\n",
                                adaptiveArgs[k+1].c_str(), nprocess);
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
            }
        if (!ranksGiven) {
            adaptiveArgs.push_back("--ranks");
            adaptiveArgs.push_back(std::to_string(nprocess));
        }
        // A calculation runs inside this process and cannot be stopped without
        // stopping the driver with it, so a limit the user asked for is one they
        // are told they are not getting (src/clustering/algo_adaptive/base_run.cpp).
        for (const std::string& arg : adaptiveArgs)
            if (arg == "--timeout" && rank == 0)
                fprintf(stderr,"Note: --timeout is not enforced when the steps run in-process; "
                               "use ccex-adaptive with --binary for a per-step limit\n");
        int status = 0;
        if (rank == 0) {
            adaptive::set_local_runner(
                [&](const std::vector<std::string>& command, const std::string&,
                    const std::string& log_path) {
                    std::string stepInput;
                    for (std::size_t i = 0; i+1 < command.size(); ++i)
                        if (command[i] == "-f") stepInput = command[i+1];
                    int length = int(stepInput.size());
                    MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
                    MPI_Bcast(const_cast<char*>(stepInput.data()),length,MPI_CHAR,0,MPI_COMM_WORLD);
                    applyStepConfig(stepInput,firstInput,cnf,op);
                    // The step's output goes to its run.log and then stdout comes
                    // back to wherever it was - a file, a pipe or a terminal. It
                    // used to come back by reopening /dev/tty, which a batch job
                    // does not have: under Slurm everything the policy printed
                    // after the first step was lost, and with a terminal it
                    // ignored a redirection (measured 2026-09-19).
                    fflush(stdout);
                    const int saved = dup(STDOUT_FILENO);
                    const int log = open(log_path.c_str(), O_WRONLY|O_CREAT|O_TRUNC, 0644);
                    if (saved < 0 || log < 0) {
                        fprintf(stderr,"Error: cannot redirect the step log to '%s'\n",
                                log_path.c_str());
                        ccexFatal();   // D38
                    }
                    dup2(log, STDOUT_FILENO);
                    close(log);
                    runOneCalculation(cnf,qa,ba,dfa,pls,cls,op,false);
                    // A separate process per step flushed its measured times on
                    // exit, and the cost table reads them back out of run.log.
                    computeProfileFlush();
                    fflush(stdout);
                    dup2(saved, STDOUT_FILENO);
                    close(saved);
                    return adaptive::NativeResult{0,false};
                });
            status = adaptive_main(adaptiveArgs);
            int done = -1;
            MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
        } else {
            for (;;) {
                int length = 0;
                MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
                if (length < 0) break;
                std::vector<char> path(std::size_t(length)+1,0);
                MPI_Bcast(path.data(),length,MPI_CHAR,0,MPI_COMM_WORLD);
                applyStepConfig(path.data(),firstInput,cnf,op);
                runOneCalculation(cnf,qa,ba,dfa,pls,cls,op,false);
                computeProfileFlush();
            }
        }
        MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
        adaptiveStatus = status;
    }

    //=======================================================
    // Free memory
    //=======================================================

    clearAdaptiveConnectionCheck();
    QubitArray_freeAll(qa);
    Cluster_freeAll(cls);
    Config_freeAll(cnf);
    BathArray_freeAll(ba);
    DefectArray_freeAll(dfa);
    Pulse_freeAll(pls);
    Output_freeAll(op);

    MPI_Barrier(MPI_COMM_WORLD);
    
    /////////////////////////////////////////////////////////

    //wall time check
    if (rank==0){ 
        time_end = MPI_Wtime();
        char message[500];
        sprintf(message, "JOB DONE : Wall time = %.5f s", time_end - time_start);
        printLineSection();
        printMessage(message);
        printLineSection();
        printf("\n");
    }

    computeBackendFinish();
    MPI_Finalize();

    return adaptiveStatus;
}
