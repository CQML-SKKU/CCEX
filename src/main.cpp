/**
 * @file main.cpp
 * @brief The entry point: read the input, build the physical objects, make the clusters, calculate and write.
 */
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
#include "../include/input/example_input.h"

// Simulators
#include "../include/calculation/simulator.h"
#include "../include/calculation/simulator_gcce_bath.h"
#include "../include/calculation/pseudospin.h"
#include "../include/calculation/memory_estimate.h"
#include "../include/calculation/superoperator.h"
#include "../include/calculation/gpu/compute_backend.h"
#include "../include/calculation/gpu/small_cluster_batch.h"
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
#include <cctype>
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
        printf("    ======================================================================\n");
        printf("    Clusterize\n");
        printf("    ======================================================================\n\n");
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
        if (verbosity) Cluster_reportClusinfo(cls);   // A76: -v only (check_master_equation runs with -v to read its cluster lines)

        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("      Wall time = %.5f s\n", time_end - time_start_step);
        printLine();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    logQuietEnd();   // A68: a configuration of k_range after the first (logQuietBegin in runInput)

    //=======================================================
    // Main calculation & save results
    //=======================================================

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
    // A68: one line with the clusters of every order, over all ranks (each rank holds its own share)
    std::vector<long long> perOrder(order+1,0),allOrders(order+1,0);
    for (int n=1; n<=order; n++) perOrder[n] = localclusters[n][0][0]-1;
    MPI_Reduce(perOrder.data(),allOrders.data(),order+1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank==0){
        std::string line;
        long long total = 0;
        for (int n=1; n<=order; n++){
            line += (n>1 ? ", " : "") + std::to_string(n) + ": " + std::to_string(allOrders[n]);
            total += allOrders[n];
        }
        printf("\n      Clusters per order : %s (total %lld, %d bath spins, %d rank%s)\n", line.c_str(), total,
               BathArray_getNspin(ba), nprocess, nprocess>1 ? "s" : "");
    }
    // A67: the memory of the largest cluster against what a rank may use, before anything is calculated (not for a plan
    // alone, nor for the adaptive probes and factors, which have their own budget)
    const std::string work = adaptiveWork(cnf).mode;
    if(!clusterPlanOnly && work!="plan" && work!="probe" && work!="factors")
        checkClusterMemory(qa,ba,dfa,cnf,pls,op,order,localclusters);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("    Main calculation\n");
        printf("    ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    // Calculate the dynamics
    if (cnf->pseudospin){   // a defect spin with added spins makes an order-1 cluster of two spins, taken for a pair
        for (int idf=0; idf<DefectArray_getNdefect(dfa); ++idf)
            if (DefectArray_getDefect_idf_naddspin(dfa,idf) > 0) { ccexRefusalPrint(stderr,"Error: method pseudospin does not take defects with added spins\n"); ccexRefuse(); }
        pseudospinChannelsBegin(op);
    }
    if(!clusterPlanOnly && adaptiveWork(cnf).mode!="plan")calculate(qa,ba,dfa,cnf,pls,cls,op,localclusters);
    MPI_Barrier(MPI_COMM_WORLD);
    if (cnf->pseudospin) pseudospinChannelsEnd();
    MPI_freeLocalClusters(order,&localclusters);
}

// Conventional CCE has no qubit in its matrices: it propagates the bath under the
// two conditional Hamiltonians, and a pulse can only mean "swap them" - an ideal pi
// flip. Any other rotation used to be calculated as one without a word; gCCE puts
// the qubit in the Hamiltonian and applies the rotation that was asked for.
void validatePulsesForMethod(Config* cnf, Pulse* pls) {
    const char* kernel=methodKernel(cnf->method);   // gCCE, megCCE and their pcce forms (pgCCE, pmegCCE)
    // A94 (R2): a finite pulse (pulse.detuning_factor other than 1, of pulse.duration) is the gCCE kernel's alone; the
    // cCCE kernel ignored both keys without a word (meCCE and megCCE refuse detuning_factor, validateMasterEquation)
    if(!strcmp(kernel,"cce") && (pls->detuning_factor!=1.0 || pls->pulse_time!=0.0)) {
        if(rank==0) fprintf(stderr,"Error: method cCCE needs pulse.detuning_factor = 1 and pulse.duration = 0 (it calculates ideal "
                                   "instantaneous pulses; a finite pulse needs method gCCE)\n");
        ccexRefuse();
    }
    if(!strcmp(kernel,"gcce") && pls->detuning_factor==1.0 && pls->pulse_time!=0.0)
        warnNotUsed("pulse.duration","it is read only with pulse.detuning_factor other than 1 (a finite pulse)");
    if(!strcmp(kernel,"gcce") || !strcmp(kernel,"megcce")) return;   // ME-gCCE applies any rotation too
    for(int i=0;i<Pulse_getNpulse(pls);++i) {
        const char axis=pls->pulse_axes[i];
        if((axis=='X' || axis=='Y') && pls->pulse_angles[i]==180.0) continue;
        if(rank==0) fprintf(stderr,"Error: method %s calculates every pulse as an ideal pi flip (X or Y, 180), "
                                   "but pulse %d of %s is %c %g. Use method %s for other rotations.\n",
                            inputMethodName(cnf->method),i+1,pls->pulsename,axis,pls->pulse_angles[i],
                            !strcmp(kernel,"mecce") ? "megCCE" : "gCCE");   // meCCE, pmeCCE: the ME method that rotates
        ccexRefuse();
    }
    // The CUDA ensemble path takes CPMG intervals the way the CPU one did.
    if(cnf->nstate==0 && !strcasecmp(pls->pulsename,"Manual") && strcasecmp(cnf->backend,"cpu")) {
        if(!strcasecmp(cnf->backend,"cuda")) {
            if(rank==0) fprintf(stderr,"Error: an ensemble (no configuration.s_range) Manual sequence is calculated on the CPU only, "
                                       "set execution.backend to cpu\n");
            ccexRefuse();
        }
        if(rank==0) fprintf(stdout,"[backend] auto: an ensemble Manual sequence uses CPU\n");
        strcpy(cnf->backend,"cpu");
    }
}

// "Error: " and the message - a format taking the method's name, if it names one - from rank 0, then the refusal
static void refuseRun(const char* format, const char* method = "") {
    if (rank==0) { fprintf(stderr, "Error: "); fprintf(stderr, format, method); fprintf(stderr, "\n"); }
    ccexRefuse();
}

// DEER bath pulses act on bath spins, so only a method with the bath in its matrices can apply them (the gCCE and megCCE
// kernels, also with cluster_mode pcce: pgCCE, pmegCCE inside), and only the CPU paths expand a pulse beyond the qubit
// (gpu_backend.cu applyPulse takes qdim x qdim pulses). The method is decided by its kernel and named as the input writes
// it (inputMethodName); cluster_mode is "standard" inside for pcce too and "adaptive" for bundle.
void validateBathPulses(Config* cnf, Pulse* pls) {
    if (!Pulse_hasBathPulses(pls)) return;
    const char* kernel = methodKernel(cnf->method);
    if (strcmp(kernel,"gcce") && strcmp(kernel,"megcce"))
        refuseRun("bath pulses (pulse.bath.count > 0) are implemented for methods gCCE and megCCE "
                  "(cluster.mode standard or pcce), this run is method %s", inputMethodName(cnf->method));
    if (strcasecmp(cnf->cluster_mode,"standard"))   // gCCE + bundle (megCCE + bundle is refused by the reader)
        refuseRun("bath pulses need cluster.mode standard or pcce (this run has cluster.mode bundle)");
    if (pls->detuning_factor != 1.0)    refuseRun("bath pulses need pulse.detuning_factor = 1");
    if (Pulse_getPulseiter(pls))        refuseRun("bath pulses cannot be combined with pulse.per_qubit = true");
    if (!strcasecmp(cnf->propagator,"sparse"))   // expm is gCCE only (Config_validatePropagator)
        refuseRun(!strcmp(kernel,"megcce") ? "bath pulses with method megCCE need execution.propagator eigen (execution.propagator sparse applies pulses as U_q x 1)"
                                           : "bath pulses need execution.propagator eigen or expm (execution.propagator sparse applies pulses as U_q x 1)");
    // gpu_warp_max_dim is -1 until Config_resolveWarpMaxDim: > 0 here is a value given in the input
    if (cnf->gpu_warp_max_dim > 0) refuseRun("bath pulses cannot use execution.gpu_options.warp_max_dim > 0 (the warp kernel applies qubit pulses only)");
    if (!strcasecmp(cnf->backend,"cuda")) refuseRun("bath pulses (pulse.bath.count > 0) run on the CPU only, set execution.backend to cpu");
    if (strcasecmp(cnf->backend,"cpu")) {
        if (rank==0) fprintf(stdout,"[backend] auto: bath pulses use CPU\n");
        strcpy(cnf->backend,"cpu");
    }
}

// method pseudospin (pseudospin.cpp): the cCCE clusters, the pairs cut to the pseudospin Hamiltonian; what it leaves out is refused
void validatePseudospin(Config* cnf, Pulse* pls, QubitArray* qa, Cluster* cls) {
    if (!cnf->pseudospin) return;
    if (QubitArray_getNqubit(qa) != 1)         refuseRun("method pseudospin supports one qubit");
    if (strcasecmp(cnf->quantity,"coherence")) refuseRun("method pseudospin computes quantity coherence only");
    if (Config_getNstate(cnf) != 0)            refuseRun("method pseudospin takes the ensemble only (leave configuration.s_range out)");
    if (Cluster_getOrder(cls) > 2)             refuseRun("method pseudospin is a pair model, order 1 or 2");
    if (Pulse_getPulseiter(pls))               refuseRun("method pseudospin needs pulse.per_qubit = false");
    if (strcasecmp(cnf->propagator,"eigen"))   refuseRun("method pseudospin needs execution.propagator eigen");
    if (Config_getKnight(cnf) || Config_getHfmedi(cnf) || cnf->delta_A)
        refuseRun("method pseudospin does not take hamiltonian.cCCE (knight, hfmedi, delta_A)");
    const float* B = Config_getBfield(cnf);   // the cut keeps m_I of the computational z axis: the nuclear Zeeman axis
    if (B[0] != 0 || B[1] != 0)                refuseRun("method pseudospin needs the field along the qubit axis (bfield [0, 0, B])");
    if (!strcasecmp(cnf->backend,"cuda"))      refuseRun("method pseudospin runs on the CPU only, set execution.backend to cpu");
    strcpy(cnf->backend,"cpu");                // auto: the CPU, and none of the GPU batches
    cnf->gpu_packed_max_dim = cnf->gpu_warp_max_dim = 0;
}

bool isMasterEquationMethod(const char* m) {   // meCCE, megCCE and their pcce forms (pmeCCE, pmegCCE)
    const char* kernel = methodKernel(m);
    return !strcmp(kernel,"mecce") || !strcmp(kernel,"megcce");
}

// meCCE / megCCE (adding_MECCE): Lindblad dissipation in Liouville space on the CPU with a matrix exponential per step.
// What their kernels do not take is refused here, before any reader. The backend (computeBackendInit), evolution
// (Config_resolveEvolution), propagator expm (Config_validatePropagator), evolution_update recursive
// (Config_validateEvolutionUpdate) and output.cluster_only are refused by the whitelists that already exist.
void validateMasterEquation(Config* cnf, Pulse* pls, QubitArray* qa, Output* op) {
    if (!isMasterEquationMethod(cnf->method)) return;   // a jump section in another method's run: validateHamiltonianForMethod
    const bool projected = !strcmp(methodKernel(cnf->method),"mecce");   // meCCE, pmeCCE
    const char* method = inputMethodName(cnf->method);                    // pmeCCE, pmegCCE: meCCE, megCCE
    if (projected && !strcasecmp(pls->pulsename,"Manual"))
        refuseRun("method %s takes its pulses from pulse.count (Ramsey, HahnEcho, CPMG), a Manual sequence is not supported", method);
    if (Pulse_getPulseiter(pls))       refuseRun("method %s needs pulse.per_qubit = false", method);
    if (pls->detuning_factor != 1.0)   refuseRun("method %s needs pulse.detuning_factor = 1", method);
    if (pls->pulse_time != 0.0)        refuseRun("method %s needs pulse.duration = 0 (it calculates ideal instantaneous pulses)", method);   // A94 R2
    if (strcasecmp(cnf->propagator,"eigen")) refuseRun("method %s needs execution.propagator eigen (its own matrix exponential)", method);
    if (!strcasecmp(op->savemode,"cluster_density_matrix")) refuseRun("method %s does not write output.savemode cluster_density_matrix", method);
    const int nqubit = QubitArray_getNqubit(qa);
    if (projected && nqubit > 1) refuseRun("method %s supports one qubit", method);
    if (cnf->jump_operators)
        for (int i = 0; i < cnf->jump_operators->njump; ++i) {
            const char* name = cnf->jump_operators->ops[i]->bath_name;
            if (strcasecmp(name,"qubit") && QubitArray_getQubitIdx_fromName(qa,name) < 0) continue;
            // meCCE projects the qubit out; the branch's per-cluster scalar entered once per cluster power (M1)
            if (projected) refuseRun("method %s does not take qubit jump operators (the qubit is projected out), use megCCE", method);
            if (nqubit > 1) refuseRun("a qubit jump operator needs one qubit (method %s)", method);
        }
    // the coherence only: the branch's megCCE density-matrix output is not ported
    if (strcasecmp(cnf->quantity,"coherence")) refuseRun("method %s computes quantity coherence only (dm is not implemented)", method);
    // the ladder's own check (rdipladder::prepare) comes after the clusters
    if (cnf->nrdip_ladder > 0)   // the ladder's internal CCE, gCCE, pCCE (rdip_ladder.cpp) in the input's names
        refuseRun("rdip_ladder needs method cCCE (cluster.mode standard or pcce) or gCCE (cluster.mode standard), this run is method %s", method);
}

// meCCE / megCCE, once the bath is read: a jump operator whose bath_name is neither the qubit nor the name of a bath spin
// of the run would add nothing to any cluster - a typo would run at rate 0 without a word - so it is refused, naming the
// spins there are. The whole bath is checked (the bath file and the Defect sub-spins createBathArray adds, as
// validateBathPulseTargets), not a cluster: most clusters lack a species. The names compare as the kernels do (case
// ignored for bath spins; "qubit" or the qubit's name as validateMasterEquation). A name that is both the qubit's and a
// bath spin's is refused too: the megCCE kernel takes it as a qubit jump and drops the bath spins without a word.
void validateJumpOperatorNames(BathArray* ba, DefectArray* dfa, QubitArray* qa, Config* cnf) {
    if (!isMasterEquationMethod(cnf->method) || !cnf->jump_operators) return;
    std::vector<const char*> names;   // the distinct bath spin names (case ignored), in the order met
    auto add = [&names](const char* n) { for (const char* m : names) if (!strcasecmp(m,n)) return; names.push_back(n); };
    const int nspin = BathArray_getNspin(ba);
    for (int j = 0; j < nspin; ++j) add(ba->bath[j]->name);
    if (DefectArray_getNdefect(dfa) > 0)
        for (int j = 0; j < nspin; ++j)
            for (int k = 0; k < DefectArray_getNaddspins_i(dfa,j); ++k) add(DefectArray_getSubbath_i_isp(dfa,j,k)->name);
    for (int i = 0; i < cnf->jump_operators->njump; ++i) {
        const char* name = cnf->jump_operators->ops[i]->bath_name;
        const bool isQubit = !strcasecmp(name,"qubit") || QubitArray_getQubitIdx_fromName(qa,name) >= 0;
        bool found = false;
        for (const char* m : names) found = found || !strcasecmp(m,name);
        // the qubit and bath spins by one name: the megCCE kernel would apply the jump to the qubit only
        if (isQubit && found) {
            if (rank==0) fprintf(stderr,"Error: hamiltonian.%s.jump_operators[%d].bath_name \"%s\" names both the qubit and "
                                        "bath spins of this run (the comparison ignores case for bath spins), the jump would "
                                        "act on the qubit only. Rename the qubit or the bath species\n",
                                 cnf->jump_operators_section, i, name);
            ccexRefuse();
        }
        if (isQubit || found) continue;
        std::string list;
        for (const char* m : names) { if (!list.empty()) list += ", "; list += m; }
        if (rank==0) fprintf(stderr,"Error: hamiltonian.%s.jump_operators[%d].bath_name \"%s\" names no bath spin of this run "
                                    "(bath spins, compared ignoring case: %s) and is not the qubit (\"qubit\" or its "
                                    "name)\n", cnf->jump_operators_section, i, name, list.empty() ? "none" : list.c_str());
        ccexRefuse();
    }
}

// hamiltonian.cCCE holds the terms of the cCCE Hamiltonian (A47): in a gCCE or megCCE run
// they would do nothing, so they are refused; meCCE reads them as cCCE does. hamiltonian.meCCE / megCCE.jump_operators
// belong to that method.
//
// delta_A true (A83) adds the second-order term delta_A(m_s).I of each bath spin (HamilQubitBathSecularApp), which needs
// a spin-1 qubit with a zero-field splitting Delta = 3/2 D_zz != 0; hfmedi already holds it (its i = j terms are the
// complete second order of the same coupling), so the two together are refused. knight divides by the same Delta, so it is
// refused when Delta = 0 (no zero-field splitting).
void validateHamiltonianForMethod(Config* cnf, QubitArray* qa) {
    const char* kernel=methodKernel(cnf->method);   // the section of a pcce form is its method's (pmeCCE: meCCE)
    if(cnf->jump_operators && strcasecmp(cnf->jump_operators_section,kernel)) {
        if(rank==0) fprintf(stderr,"Error: hamiltonian.%s (jump_operators) belongs to method %s, this run is %s\n",
                            cnf->jump_operators_section,cnf->jump_operators_section,inputMethodName(cnf->method));
        ccexRefuse();
    }
    if((!strcmp(kernel,"gcce") || !strcmp(kernel,"megcce")) && cnf->hamiltonian_ccce_given) {
        if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE, this run is %s\n",
                            inputMethodName(cnf->method));
        ccexRefuse();
    }
    if(Config_getKnight(cnf)) {
        const double Delta = QubitArray_getIntmap_i_j(qa,0,0)(2,2).real() * 3.0 / 2.0;
        if(Delta == 0.0) {
            if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE.knight divides by the zero-field splitting of the qubit (Delta = "
                                       "3/2 D_zz of the qubit's intmap self term), this run has Delta %g kHz\n",
                                       RADKHZ_TO_KHZ(Delta));
            ccexRefuse();
        }
    }
    if(cnf->delta_A) {
        const double Delta = QubitArray_getIntmap_i_j(qa,0,0)(2,2).real() * 3.0 / 2.0;
        if(QubitArray_getNqubit(qa) != 1 || QubitArray_getQubit_i_spin(qa,0) != 1.0f || Delta == 0.0) {
            if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE.delta_A is the second-order term of a spin-1 qubit with a zero-field "
                                       "splitting (Delta = 3/2 D_zz of the qubit's intmap self term), this run has %d qubit(s), "
                                       "spin %g, Delta %g kHz\n", QubitArray_getNqubit(qa),
                                       (double)QubitArray_getQubit_i_spin(qa,0), RADKHZ_TO_KHZ(Delta));
            ccexRefuse();
        }
        if(Config_getHfmedi(cnf)) {
            if(rank==0) fprintf(stderr,"Error: hamiltonian.cCCE.delta_A is refused with hfmedi: hfmedi already includes the "
                                       "delta-A (same-nucleus second-order) term\n");
            ccexRefuse();
        }
    }
}

// A cluster.adaptive object in the input file asks for the refinement (the --adaptive option
// of the command line until A63, a top-level "adaptive" until A65): its keys are the refinement
// options with '-' written as '_',
//
//     "cluster": { "mode": "bundle",
//                  "adaptive": { "policy": "merge", "output": "./adaptive_out",
//                                "atol": 1e-8, "rtol": 2e-6, "max_rounds": 8 } }
//
// (written "Adaptive" until A47; keys are read without regard to case). The block belongs to cluster.mode pcce or
// bundle. A step's input, which the refinement driver writes, holds
// cluster.adaptive: {"work_file": ...} alone: that is the step's
// work request (cJSON_readOptionConfig), not a request for a refinement.
//
// and what comes back is the argument list adaptive_main already takes, so the
// policy code and its own checks are the ones the standalone tool's command line goes through.
// What the run itself decides is refused here rather than overridden: the input
// is this file, the binary is this process, the ranks are mpirun's.
std::vector<std::string> adaptiveArgsFromInput(char* fccein, const char* self) {
    std::vector<std::string> args;
    if(access(fccein, R_OK) != 0) return args;          // the ordinary read reports it
    char* text = cJSON_ReadFccein(fccein);
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> root(cJSON_Parse(text),cJSON_Delete);
    freeChar1d(&text);
    // cluster.adaptive since A65 (a top-level "adaptive" is refused by the reader, with the new place)
    cJSON* clusterBlock = root ? cJSON_GetObjectItem(root.get(),"cluster") : nullptr;
    cJSON* block = cJSON_IsObject(clusterBlock) ? cJSON_GetObjectItem(clusterBlock,"adaptive") : nullptr;
    const auto fail = [](const std::string& message) {
        ccexRefusalPrint(stderr,"Error: %s\n",message.c_str());
        ccexRefuse();
    };
    // The block was called "Refine" until 2026-09-23. The old name is refused, not read:
    // a run that silently did no refinement would look like a finished one.
    if(root && cJSON_GetObjectItem(root.get(),"Refine"))
        fail("the \"Refine\" block is called \"adaptive\" now, rename it (docs/user/InputFile_Cluster.md)");
    if(!block) return args;
    if(!cJSON_IsObject(block)) fail("cluster.adaptive must be an object of refinement options");
    if(cJSON_GetObjectItem(block,"work_file")) {
        // a step's input: the work request, and the probe memory of the user's block (read by the step's Config)
        for(cJSON* item=block->child; item; item=item->next)
            if(strcasecmp(item->string,"work_file") && strcasecmp(item->string,"probe_memory_gb"))
                fail("cluster.adaptive.work_file is written by the refinement driver into a step's input, beside no other "
                     "cluster.adaptive key but probe_memory_gb");
        return args;
    }
    cJSON* mode = cJSON_IsObject(clusterBlock) ? cJSON_GetObjectItem(clusterBlock,"mode") : nullptr;
    const char* shape = cJSON_IsString(mode) ? mode->valuestring : "standard";
    if(strcasecmp(shape,"pcce") && strcasecmp(shape,"bundle"))
        fail(std::string("the cluster.adaptive block belongs to cluster.mode pcce or bundle (this input has cluster.mode ")+shape+")");
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
        fail("cluster.adaptive."+key+" must be a number or a string");
        return "";
    };
    for(cJSON* item=block->child; item; item=item->next) {
        const std::string key = item->string;
        if(!strcasecmp(key.c_str(),"probe_memory_gb")) continue;   // the steps' setting (A65), read by the Config reader
        // A94 (G23): an input key is spelled with _ as every other key (- is the command line's, --max-rounds)
        if(key.find('-') != std::string::npos) {
            std::string spelled = key;
            for(char& ch : spelled) if(ch == '-') ch = '_';
            fail("cluster.adaptive."+key+" is spelled \""+spelled+"\" in the input (- is the command line's spelling)");
        }
        // A94 (G10): the folder of the steps is workdir, as auto_convergence.workdir (output until A94: it read like the
        // top-level output block); the driver's option is still --output
        if(!strcasecmp(key.c_str(),"output"))
            fail("cluster.adaptive.output is now cluster.adaptive.workdir (the folder of the steps and summary.json)\n"
                 "Old inputs: python3 example/tools/InputKeys/convert_input_keys.py FILE... rewrites them with the new keys");
        // A94 (G9): merge is the only policy; the key stays readable for old inputs, and is reported as not needed
        if(!strcasecmp(key.c_str(),"policy"))
            warnNotUsed("cluster.adaptive.policy","merge is the only policy (the default), leave the key out");
        std::string name = !strcasecmp(key.c_str(),"workdir") ? std::string("output") : key;
        for(char& ch : name) ch = ch=='_' ? '-' : char(std::tolower((unsigned char)ch));   // keys are read without regard to case
        const adaptive::options::Option* option = nullptr;
        for(const adaptive::options::Option& one : adaptive::options::kOptions)
            if(name == one.name) option = &one;
        if(!option) fail("cluster.adaptive."+key+" is not a refinement option");
        for(const char* f : fixed)
            if(name == f) fail("cluster.adaptive."+key+" is decided by the run itself and cannot be set in the input");
        if(cJSON_IsNull(item)) continue;
        if(option->kind == adaptive::options::Kind::Flag) {
            if(!cJSON_IsBool(item)) fail("cluster.adaptive."+key+" is a switch: write true or false");
            if(cJSON_IsTrue(item)) args.push_back("--"+name);
        } else if(cJSON_IsArray(item)) {
            fail("cluster.adaptive."+key+" takes one value, not a list");
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
        if(!given) fail(std::string("cluster.adaptive.")+(std::string(one.name) == "output" ? "workdir" : one.dest)+" is required");
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
    // output.file, output.savemode and cluster.adaptive ({"work_file"} in a step) - in blocks since A65
    static const char* mutable_keys[][2] = {{"output", "file"}, {"output", "savemode"}, {"cluster", "adaptive"}};
    // Both files are read as CCEX inputs (a //, # or ! outside a string starts a comment): the
    // first one is the user's, and may carry comments; the step's is what the
    // policy wrote, and carries none.
    char* firstText = cJSON_ReadFccein(const_cast<char*>(firstInput.c_str()));
    char* stepText = cJSON_ReadFccein(const_cast<char*>(stepInput.c_str()));
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> a(cJSON_Parse(firstText),cJSON_Delete);
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> b(cJSON_Parse(stepText),cJSON_Delete);
    freeChar1d(&firstText); freeChar1d(&stepText);
    if(!a || !b) { ccexRefusalPrint(stderr,"Error: refinement step input is not JSON\n"); ccexRefuse(); }
    // the rest is compared without the keys a step may change; an output or cluster block left empty is no block
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> aKept(cJSON_Duplicate(a.get(),1),cJSON_Delete);
    std::unique_ptr<cJSON,decltype(&cJSON_Delete)> bKept(cJSON_Duplicate(b.get(),1),cJSON_Delete);
    for(cJSON* root : {aKept.get(), bKept.get()})
        for(const auto& k : mutable_keys)
            if(cJSON* block=cJSON_GetObjectItem(root,k[0])) {
                cJSON_DeleteItemFromObject(block,k[1]);
                if(cJSON_IsObject(block) && !block->child) cJSON_DeleteItemFromObject(root,k[0]);
            }
    for(cJSON* item=bKept->child; item; item=item->next) {
        cJSON* was=cJSON_GetObjectItem(aKept.get(),item->string);
        char* x=cJSON_PrintUnformatted(item); char* y=was?cJSON_PrintUnformatted(was):nullptr;
        const bool same = y && !strcmp(x,y);
        free(x); if(y) free(y);
        if(!same) {
            ccexRefusalPrint(stderr,"Error: refinement step changed '%s', which the kept bath depends on\n",
                    item->string);
            ccexRefuse();
        }
    }
    cJSON* v;
    cJSON* stepCluster=cJSON_GetObjectItem(b.get(),"cluster");
    cJSON* stepOutput=cJSON_GetObjectItem(b.get(),"output");
    cJSON* stepBlock=stepCluster ? cJSON_GetObjectItem(stepCluster,"adaptive") : nullptr;
    if(stepBlock && (v=cJSON_GetObjectItem(stepBlock,"work_file")) && cJSON_IsString(v)) {
        strncpy(cnf->adaptive_work_file,v->valuestring,MAX_FILEPATH-1);
        cnf->adaptive_work_file[MAX_FILEPATH-1]=0;
    } else cnf->adaptive_work_file[0]=0;
    if(stepOutput && (v=cJSON_GetObjectItem(stepOutput,"file")) && cJSON_IsString(v)) {
        // resolveOutputName has allocated the buffer before any step; the step always names one.
        Output_setOutfile(op,v->valuestring);
    }
    if(stepOutput && (v=cJSON_GetObjectItem(stepOutput,"savemode")))   // the mode; nodiv stays the user input's
        cJSON_readSavemode(v,op,false);
}
}   // namespace

// main.out -f input.json with an auto_convergence block ((internal notes), "--auto_conv 첫 판 설계
// v2"; the --auto_conv option itself was removed in A63): before MPI starts, hand the input to the runner of
// example/tools/AutoConvergence (python -m autoconv run), which launches the ordinary calculations itself. Under
// mpirun -n N, rank 0 goes on alone and the runner launches with mpirun -n N (without mpirun, on one rank)
// (plan "mpirun으로 부르면 그 rank 수로", 2026-10-10). The tool is CCEX_AUTOCONV_DIR or the first parent of this binary holding
// example/tools/AutoConvergence; the interpreter is CCEX_PYTHON or python3.
// The FILE of -f FILE or -fFILE at argv[i] (getopt's two spellings; i moves onto FILE), or nullptr.
static const char* inputAt(int argc, char* argv[], int& i){
    if (!strcmp(argv[i], "-f") && i + 1 < argc) return argv[++i];
    if (!strncmp(argv[i], "-f", 2) && argv[i][2] != '\0') return argv[i] + 2;
    return nullptr;
}
// The input of -f FILE or -fFILE (the last one given).
static const char* inputArgument(int argc, char* argv[]){
    const char* input = nullptr;
    for (int i = 1; i < argc; ++i) if (const char* f = inputAt(argc, argv, i)) input = f;
    return input;
}

static int autoConvEntry(int argc, char* argv[]){
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
    // the runner reads only the input: an option such as -o (output folder) or -v (verbose) would be dropped without a word, so it
    // is refused (review 20261008-autoconv-autostart #1)
    for (int i = 1; i < argc; ++i){
        if (inputAt(argc, argv, i)) continue;
        fprintf(stderr, "Error: the automatic convergence takes only -f FILE, \"%s\" is not used by it (set it in the input, or set "
                        "auto_convergence.enabled false for an ordinary run)\n", argv[i]);
        return 2;
    }
    const char* input = inputArgument(argc, argv);
    if (!input){ fprintf(stderr, "Error: the automatic convergence needs -f input.json\n"); return 2; }
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
        { fprintf(stderr, "Error: the automatic convergence cannot find its own path from the system\n"); return 1; }
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
    if (tool.empty()){ fprintf(stderr, "Error: the automatic convergence cannot find example/tools/AutoConvergence (set CCEX_AUTOCONV_DIR)\n"); return 1; }
    char real[PATH_MAX];
    if (!realpath(tool.c_str(), real)){ fprintf(stderr, "Error: the automatic convergence cannot resolve %s: %s\n", tool.c_str(), strerror(errno)); return 1; }
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
    // python -m puts the working directory first on sys.path, ahead of the folder checked above: an autoconv/ or a
    // module of the standard library's name there would run instead (review 20261014 F3; Python 3.11 and later)
    setenv("PYTHONSAFEPATH", "1", 1);
    const char* old = getenv("PYTHONPATH");
    setenv("PYTHONPATH", (old && *old ? tool + ":" + old : tool).c_str(), 1);
    const char* python = getenv("CCEX_PYTHON") ? getenv("CCEX_PYTHON") : "python3";
    printf("[auto] runner %s/autoconv (python %s), binary %s\n", tool.c_str(), python, exe);
    fflush(stdout);
    char* args[] = {const_cast<char*>(python), const_cast<char*>("-m"), const_cast<char*>("autoconv"), const_cast<char*>("run"),
                    const_cast<char*>(input), const_cast<char*>("--binary"), exe, nullptr};
    execvp(python, args);
    fprintf(stderr, "Error: the automatic convergence cannot start %s: %s\n", python, strerror(errno));
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

// A63 (user 2026-10-15): the command line is -f FILE, -v, -h/--help, -o DIR and -e. Every other setting is an input
// key; an option that was removed is refused with the key that replaces it.
static const char kShortOptions[] = ":hvef:o:";
static const char kRemovedShortOptions[] = "PN:m:q:B:I:s:S:a:";
static const char* removedShortOption(int c){
    switch (c){
        case 'I': return "-I was removed: give the bath files in the input, \"configuration\": {\"bathfile\": [\"FILE\", ...]}";
        case 'N': return "-N was removed: give the bath states in the input, \"configuration\": {\"s_range\": [1, N]}";
        case 'm': return "-m was removed: give the method in the input, \"method\" (with \"cluster.mode\" for the cluster shape)";
        case 'q': return "-q was removed: give the quantity in the input, \"quantity\"";
        case 'B': return "-B was removed: give the field in the input, \"bfield\" (a number is Bz, or [Bx, By, Bz])";
        case 's': return "-s was removed: give the bath state files in the input, \"configuration\": {\"bath_state_file\": [...]}";
        case 'S': return "-S was removed: give the defect state files in the input, \"configuration\": {\"defect_state_file\": [...]}";
        case 'a': return "-a was removed: give the defect axis files in the input, \"configuration\": {\"defect_axis_file\": [...]}";
        case 'P': return "-P was removed: write \"output.cluster_only\": true in the input";
        default:  return nullptr;
    }
}
static const char* removedLongOption(const char* arg){
    if (!strcmp(arg, "--auto_conv"))
        return "--auto_conv was removed: an \"auto_convergence\" block in the input starts the automatic convergence "
               "(main.out -f FILE, example/tools/AutoConvergence/README.md)";
    if (!strcmp(arg, "--auto_status"))
        return "--auto_status was removed: the automatic convergence tool tells where a run is, "
               "PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status FILE (example/tools/AutoConvergence/README.md)";
    if (!strcmp(arg, "--adaptive"))
        return "--adaptive was removed: write the refinement options as a \"cluster.adaptive\" block in the input "
               "(docs/user/InputFile_Cluster.md)";
    return nullptr;
}
// The first removed option on the command line, or nullptr. Found before MPI starts (the automatic convergence starts
// without it); main refuses it once MPI is up, on every rank.
static const char* findRemovedOption(int argc, char* argv[]){
    for (int i = 1; i < argc; ++i){
        if (!strcmp(argv[i], "--")) break;
        if (const char* m = removedLongOption(argv[i])) return m;
    }
    const std::string all = std::string(kShortOptions) + kRemovedShortOptions;
    std::vector<char*> args(argv, argv + argc);       // getopt may permute: a copy keeps argv for the real pass
    args.push_back(nullptr);
    optind = 1; opterr = 0;
    int c;
    const char* found = nullptr;
    while (!found && (c = getopt(argc, args.data(), all.c_str())) != -1)
        found = removedShortOption(c == ':' || c == '?' ? optopt : c);
    optind = 1;
    return found;
}

// mkdir -p; 0 or the errno of the failure (ENOTDIR when the path is a file)
static int makeFolders(const std::string& path){
    for (size_t i = 1; i <= path.size(); ++i)
        if (i == path.size() || path[i] == '/')
            if (mkdir(path.substr(0, i).c_str(), 0777) != 0 && errno != EEXIST) return errno;
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return errno;
    return S_ISDIR(st.st_mode) ? 0 : ENOTDIR;
}
// -o DIR or output.path: rank 0 creates the folder, every rank learns whether it could
static void makeOutputFolder(const char* dir, const char* given = "-o"){
    int error = rank == 0 ? makeFolders(dir) : 0;
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (error){
        ccexRefusalPrint(stderr, "Error: %s %s: cannot create the output folder (%s)\n", given, dir, strerror(error));
        ccexRefuse();
    }
}
// A73: the folder of the outputs takes a file - rank 0 creates and removes one there, every rank learns the outcome - so a
// folder that cannot be written is refused before anything is calculated, not found when the first output is written
// (with configuration.k_range, after the first configuration)
static void checkOutputWritable(const char* outfile){
    std::string folder = outfile;
    const size_t slash = folder.rfind('/');
    folder = slash == std::string::npos ? "." : slash == 0 ? "/" : folder.substr(0, slash);
    int error = 0;
    if (rank == 0){
        std::string probe = folder + (folder.back() == '/' ? "" : "/") + ".ccex-write-check-XXXXXX";
        const int fd = mkstemp(&probe[0]);
        if (fd < 0) error = errno;
        else { close(fd); unlink(probe.c_str()); }
    }
    MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (error){
        ccexRefusalPrint(stderr, "Error: cannot write the outputs %s* in the folder %s (%s), give a folder that exists and can be "
                                 "written (output.file, output.path or -o DIR)\n", outfile, folder.c_str(), strerror(error));
        ccexRefuse();
    }
}
// The output name (A63): the input's output.file, or the input file's name without .json; inside the output folder - -o DIR,
// else the input's output.path (A65), else the working directory -, where an absolute output.file is refused (the two would
// name two places).
static void resolveOutputName(Output* op, const char* input, const char* option){
    const char* given = Output_getOutfile(op);
    std::string name = given ? given : "";
    const char* dir = option ? option : (op->outpath[0] ? op->outpath : nullptr);
    if (dir && !name.empty() && name[0] == '/'){
        if (option)
            ccexRefusalPrint(stderr, "Error: the input's output.file \"%s\" is an absolute path and -o %s names the output folder, "
                                     "give a relative output.file, or leave -o out\n", name.c_str(), dir);
        else
            ccexRefusalPrint(stderr, "Error: the input's output.file \"%s\" is an absolute path and output.path \"%s\" names the "
                                     "output folder, give a relative output.file, or leave output.path out\n", name.c_str(), dir);
        ccexRefuse();
    }
    if (name.empty()){
        name = input;
        if (name.rfind('/') != std::string::npos) name.erase(0, name.rfind('/') + 1);
        for (const char* ext : {".jsonc", ".json"}){
            const size_t n = strlen(ext);
            if (name.size() > n && name.compare(name.size() - n, n, ext) == 0){ name.erase(name.size() - n); break; }
        }
    }
    if (dir){
        makeOutputFolder(dir, option ? "-o" : "output.path");
        std::string folder = dir;
        if (folder.back() != '/') folder += '/';
        name = folder + name;
    }
    if (!Output_getOutfile(op)) Output_allocOutfile(op);
    Output_setOutfile(op, const_cast<char*>(name.c_str()));
}
// -e: the input template, in the working directory or in DIR; an existing file is not overwritten
static void writeExampleInput(const char* dir){
    if (dir) makeOutputFolder(dir);
    std::string folder = dir ? dir : ".";
    if (folder.back() != '/') folder += '/';
    const std::string path = folder + kExampleInputName;
    int status = 0;                   // 1: the file exists, 2: it could not be written
    if (rank == 0){
        if (access(path.c_str(), F_OK) == 0) status = 1;
        else {
            FILE* fp = fopen(path.c_str(), "wx");
            if (!fp) status = 2;
            else {
                const bool written = fputs(kExampleInput, fp) >= 0;
                if (fclose(fp) != 0 || !written) status = 2;
            }
        }
    }
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status == 1){
        ccexRefusalPrint(stderr, "Error: %s exists, -e does not overwrite it (remove it, or give another folder with -o DIR)\n",
                         path.c_str());
        ccexRefuse();
    }
    if (status){
        ccexRefusalPrint(stderr, "Error: -e cannot write %s\n", path.c_str());
        ccexRefuse();
    }
    if (rank == 0)
        printf("Wrote %s, an input template: name your bath file in its \"bathfile\" before running it\n"
               "(a ready bath and input: example/tutorial/00_quickstart)\n", path.c_str());
}

static void printUsage(const char* self){
    const char* name = strrchr(self, '/') ? strrchr(self, '/') + 1 : self;
    printf("Usage: %s -f FILE [-o DIR] [-v]\n", name);
    printf("       %s -e [-o DIR]\n", name);
    printf("  -f FILE   input file (JSON, // or # comments)\n");
    printf("  -o DIR    output folder (takes the place of output.path)\n");
    printf("  -v        long log\n");
    printf("  -e        write %s, every input key with its default, commented, and stop\n", kExampleInputName);
    printf("  -h        this message\n");
    printf("Input keys: docs/user/InputFile.md   Tutorials: example/tutorial/\n");
}

// One run of the input: read and check it, read its files and calculate. Called once, or once for every configuration k
// of configuration.k_range (A64), each a run of its own from the readers to computeBackendFinish, as if the input had
// only the files of k: the defect axes' count of random draws starts again (DefectArray_resetPaxesDraws). Only the CUDA
// context of a GPU run is kept from one configuration to the next (computeBackendFinish(true)).
static int runInput(char* fccein, const char* outputFolder, std::vector<std::string> adaptiveArgs, bool ranged, long k,
                    long rangeFirst, long rangeLast, double time_start){
    double time_end;
    int adaptiveStatus = 0;
    if (ranged && rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("    Configuration k = %ld (configuration.k_range [%ld, %ld])\n", k, rangeFirst, rangeLast);
        printf("    ======================================================================\n\n");
    }
    // A68: after the first configuration, the read-in, the reports and the clusters repeat what the log holds; only
    // the files read, the spin counts, the tagged lines and the warnings are printed (all of it with -v), until
    // runOneCalculation ends this before the calculation, with the cluster counts
    logQuietBegin(ranged && k > rangeFirst && !verbosity);
    DefectArray_resetPaxesDraws();

    //Simulation objects
    Config*      cnf = Config_init();
    QubitArray*  qa  = QubitArray_init();
    Cluster*     cls = Cluster_init();
    BathArray*   ba  = BathArray_init();
    DefectArray* dfa = DefectArray_init();
    Pulse*       pls = Pulse_init();
    Output*      op  = Output_init();

    // read the input file
    cJSON_readOptionConfig      (cnf, fccein); // general.h
    cJSON_readOptionQubitArray  (qa,  fccein); // qubit.h
    cJSON_readOptionDefectArray (dfa, fccein); // defect.h
    cJSON_readOptionCluster     (cls, fccein); // cluster.h
    cJSON_readOptionPulse       (pls, fccein); // pulse.h
    cJSON_readOptionOutput      (op,  fccein); // output.h
    cJSON_checkInputKeys        (fccein);      // a key nobody read: a typo stops the run (D24)

    if (rank==0){
        printf("\n    >> Read %s file successfully\n",fccein);
    }

    // the refinement writes every step into its own folder (cluster.adaptive.output)
    if (!adaptiveArgs.empty() && (outputFolder || op->outpath[0])){
        ccexRefusalPrint(stderr, "Error: %s does not apply to the adaptive refinement, its folder is cluster.adaptive.output\n",
                         outputFolder ? "-o" : "output.path");
        ccexRefuse();
    }
    resolveOutputName(op, fccein, outputFolder);
    if (adaptiveArgs.empty()) checkOutputWritable(Output_getOutfile(op));   // the refinement writes in its own folder
    if (ranged){   // A64: <outfile>_k<k>_..., one set of outputs per configuration
        const std::string name = std::string(Output_getOutfile(op)) + "_k" + std::to_string(k);
        Output_setOutfile(op, const_cast<char*>(name.c_str()));
    }
    const bool clusterPlanOnly = cnf->cluster_only;

    // Every reader has run by here. Only now can the defaulted evolution be settled: it depends on the method and nstate,
    // which different readers set.
    Pulse_resolveBathSequence(pls);   // after every reader: the merged sequence is what every later check sees
    validateBathPulses(cnf,pls);      // before the packed/warp defaults: they read backend
    Config_resolveEvolution(cnf);
    Config_resolvePackedMaxDim(cnf);
    Config_resolveWarpMaxDim(cnf);             // A40, before the evolution_update check (its log names the warp kernel)
    Config_validatePropagator(cnf,QubitArray_getNqubit(qa));   // A2: before the evolution_update check, so sparse is named
    Config_validateEvolutionUpdate(cnf);       // after the packed kernel's value (A36)
    validatePulsesForMethod(cnf,pls);
    validateHamiltonianForMethod(cnf,qa);
    validatePseudospin(cnf,pls,qa,cls);        // before computeBackendInit: it settles the backend
    validateAdaptiveConfig(cnf,cls);
    Config_validatePcceCache(cnf);
    if (clusterPlanOnly && ((strcasecmp(cnf->method,"cce") &&
         strcasecmp(cnf->method,"gcce") && strcasecmp(cnf->method,"pcce")) ||
         !strcasecmp(cnf->bundle_connections,"checked") || adaptiveWork(cnf).mode=="probe" ||
         adaptiveWork(cnf).mode=="factors" || adaptiveWork(cnf).mode=="except")) {
        if(rank==0)fprintf(stderr,"Error: output.cluster_only supports cCCE, gCCE and cluster.mode pcce plans (cCCE only) with explicit connections, rdip-omitted checked connections and probe, factors and except work require dynamics\n");
        ccexRefuse();
    }
    // the refinement runs its steps with dynamics; the -P it replaced was dropped there without a word
    if (clusterPlanOnly && !adaptiveArgs.empty()) {
        ccexRefusalPrint(stderr,"Error: output.cluster_only does not apply to the adaptive refinement (an adaptive block in the input)\n");
        ccexRefuse();
    }
    computeBackendInit(cnf,qa,pls,Output_getSavemode(op));

    // The field can only be checked against the rotation once every reader has run.
    Config_validateBfieldAlignment(cnf);

    // Needs the Config AND the QubitArray, which two different parsers build, so it can
    // only run once both are done. Still before every reader, so an unsupported
    // combination stops the run before a tensor file is opened.
    validateCoordinateFrameRotationInputs(qa,cnf);
    validateMasterEquation(cnf,pls,qa,op);   // meCCE / megCCE: what their kernels do not take, before any reader

    MPI_Barrier(MPI_COMM_WORLD);

    //=======================================================
    // Report options
    //=======================================================
    if (rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("    Report all structures\n");
        printf("    ======================================================================\n\n");

        logAlignBegin();   // A76: the ':' of each table in one column
        printf("    ----------------------------------------------------------------------\n");
        Config_report(cnf,QubitArray_getOverhaus(qa));
        printf("\n\n");

        // A76: the qubits are reported once, after the files (their positions and tensors are final there)
        if (verbosity && DefectArray_getNdefect(dfa) > 0){
            printf("    ----------------------------------------------------------------------\n");
            printTitle("DefectArray");
            for (int idf=0; idf<DefectArray_getNdefect(dfa); idf++){
                DefectArray_reportDefect_idf(dfa,idf);
                printf("\n");
            }
        }

        if (verbosity){
            printf("    ----------------------------------------------------------------------\n");
            Cluster_report(cls);
            printf("\n\n");
        }

        printf("    ----------------------------------------------------------------------\n");
        Pulse_report(pls);
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Output_report(op);
        printf("\n");
        logAlignEnd();

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
        printf("    Read files\n");
        printf("    ======================================================================\n\n");
    }

    // the qubit positions of Qubit.qubit[].file are read with the Qubit block (A94, G12; readQubitfile until A94)

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
    // hamiltonian.common.hyperfine.frame / quadrupole.frame declares, and a point-dipole tensor CCEX computed
    // itself is in the bath frame no matter what those say.
    //
    // This step normalizes the two into one: every qubit position and every bath-spin
    // position is rotated by the same R about the same r0, and the tensors that are still
    // in the bath frame -- HF, QD, Defect templates, and Qubit.intmap entries whose
    // provenance says so -- are transformed with them. Tensors already in the
    // computational frame are left untouched. The magnetic field is NOT settled here;
    // that happened right after the input was read, above.
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

        if (rank==0 && verbosity){
            logAlignBegin();
            DefectArray_reportPaxes(dfa);
            logAlignEnd();
        }

        DefectArray_allocNaddspins(dfa, nspin);
        DefectArray_setNaddspins(dfa,ba); // Naddspins
        if (rank==0 && verbosity){
            logAlignBegin();
            DefectArray_reportNaddspins(dfa);
            logAlignEnd();
        }

        int nqubit = QubitArray_getNqubit(qa);
        DefectArray_allocSubbath(dfa,ba,nqubit);
        DefectArray_setSubbath(dfa,ba,qa);
        if (rank==0 && verbosity){
            logAlignBegin();
            DefectArray_reportSubbath_props(dfa);
            DefectArray_reportSubbath_hypfs(dfa, nqubit);
            DefectArray_reportSubbath_quads(dfa);
            DefectArray_reportSubbath_hypf_subs(dfa);
            logAlignEnd();
        }

        // Update main bathspins : mainspidx, zfs(pax), detuning(pax)
        updateMainSpins_fromDefectArray(dfa,ba);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    validateBathPulseTargets(ba,dfa,pls);   // DEER: the driven spins exist and fit pulse.bath (bath and detunings set)
    validateJumpOperatorNames(ba,dfa,qa,cnf);   // meCCE / megCCE: every jump's bath_name is the qubit or a bath spin
    smallClusterWarmUpAfterBath(cnf,ba,dfa);   // A55: the packed kernel's card, when the bath says the run reaches its threshold

    if (rank==0 && verbosity){
        printf("    ----------------------------------------------------------------------\n");
        logAlignBegin();
        printTitle("QubitArray");
        QubitArray_report(qa);
        printf("\n");

        printf("    ----------------------------------------------------------------------\n");
        printTitle("BathArray");
        BathArray_report(ba);
        printf("\n");
        logAlignEnd();

    }

    if (rank==0){
        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("      Wall time = %.5f s\n", time_end - time_start);
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
        const std::string firstInput = fccein;   // never null here: main refuses a run without -f
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
                        fprintf(stderr,"Error: --ranks %s, but this process has %d MPI rank(s), "
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
                fprintf(stderr,"Note: --timeout is not enforced when the steps run in-process, "
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

    // A64: the CUDA context stays open for the next configuration of configuration.k_range
    computeBackendFinish(ranged && k < rangeLast);
    if (ranged && verbosity){   // -v: the memory of each rank after each configuration (a leak across the loop would grow it)
        long rss = -1, hwm = -1;
        if (FILE* fp = fopen("/proc/self/status", "r")){
            char line[256];
            while (fgets(line, sizeof line, fp)){
                if (!strncmp(line, "VmRSS:", 6)) rss = atol(line + 6);
                if (!strncmp(line, "VmHWM:", 6)) hwm = atol(line + 6);
            }
            fclose(fp);
        }
        printf("[configuration] k=%ld rank=%d VmRSS_kB=%ld VmHWM_kB=%ld\n", k, rank, rss, hwm);
        fflush(stdout);
    }
    return adaptiveStatus;
}

int main(int argc, char* argv[]){

    // --help is -h (review 20261013-doxygen D21); getopt reads short options only
    for (int i = 1; i < argc; ++i) if (!strcmp(argv[i], "--help")) argv[i] = const_cast<char*>("-h");
    // A removed option is refused before anything else reads the command line, the automatic convergence included
    const char* removedOption = findRemovedOption(argc, argv);
    if (!removedOption && inputAsksAutoConv(argc, argv)) return autoConvEntry(argc, argv);

    //Set up MPI. FUNNELED: the packed GPU batch (gpu_options.packed_max_dim) runs threads beside the main one, which alone
    //calls MPI; smallClusterWarmUp refuses that batch, before any thread starts, if the MPI provides less (review
    //20261012-233058-452278 P3-03). Nothing else depends on the level provided.
    int threadLevel = MPI_THREAD_SINGLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadLevel);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process id
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess); // Get number of processes
    if (removedOption){
        ccexRefusalPrint(stderr, "Error: %s\n", removedOption);
        ccexRefuse();
    }

    // wall time check
    double time_start, time_end;
    time_start = MPI_Wtime();

    //=======================================================
    // Variables declaration
    //=======================================================

    // Input file (cce.in)
    char* fccein = NULL;

    //=======================================================
    // Read options
    //=======================================================
    // getopt finds the input file and the folder in every spelling ("-f in.json", "-fin.json", "-vfin.json"; review
    // 20261013-180049 I01); the last -f and the last -o count.
    const char* outputFolder = NULL;
    bool help = false, example = false;
    int c;
    optind = 1; opterr = 0;
    while ((c = getopt(argc, argv, kShortOptions)) != -1) {
        switch (c) {
            case 'f': fccein = optarg; break;
            case 'o': outputFolder = optarg; break;
            case 'v': verbosity = true; break;
            case 'h': help = true; break;
            case 'e': example = true; break;
            case ':':
                ccexRefusalPrint(stderr, "Error: -%c needs a value (-h lists the options)\n", optopt);
                ccexRefuse();
            default:
                if (optopt) ccexRefusalPrint(stderr, "Error: unknown option -%c (-h lists the options)\n", optopt);
                else ccexRefusalPrint(stderr, "Error: unknown option %s (-h lists the options)\n", argv[optind - 1]);
                ccexRefuse();
        }
    }
    if (optind < argc){
        ccexRefusalPrint(stderr, "Error: \"%s\" is not an option (the input file is -f FILE, -h lists the options)\n", argv[optind]);
        ccexRefuse();
    }
    if (help){
        if (rank == 0) printUsage(argv[0]);
        // Open MPI's mpirun reports a process that leaves without MPI_Finalize as failed (exit 1 for -h on CI)
        fflush(stdout);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    // Print welcome message (not with -h)
    if (rank==0){
        printBanner();
    }
    if (example){
        if (fccein){
            ccexRefusalPrint(stderr, "Error: -e writes an example input and calculates nothing, leave -f out\n");
            ccexRefuse();
        }
        writeExampleInput(outputFolder);
        fflush(stdout);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    if (fccein == NULL){
        ccexRefusalPrint(stderr, "Error: input file is required (-f FILE, -e writes an example)\n");
        ccexRefuse();
    }

    // The adaptive refinement, asked for by an adaptive block in the input (adaptiveArgsFromInput; the --adaptive option
    // was removed in A63)
    std::vector<std::string> adaptiveArgs;
    int adaptiveStatus = 0;
    {
        char self[4096] = {};
        if (readlink("/proc/self/exe", self, sizeof self-1) <= 0)
            snprintf(self, sizeof self, "%s", argv[0]);
        adaptiveArgs = adaptiveArgsFromInput(fccein, self);
    }

    if (rank==0){
        char message[500];
        snprintf(message, sizeof message, "Option filename : '%s'\n", fccein);
        printf("      %s\n", message);
    }
    // check if the file exists
    if (access(fccein, F_OK) == -1){
        ccexRefusalPrint(stderr, "Error: Option file does NOT exist\n\n");
        ccexRefuse();
    }

    // A64: with configuration.k_range every configuration k of it is a run of its own, one after another, each checked
    // before anything is calculated (cJSON_readConfigurationRange)
    long rangeFirst = 0, rangeLast = 0;
    const bool ranged = cJSON_readConfigurationRange(fccein, &rangeFirst, &rangeLast);
    if (ranged && !adaptiveArgs.empty()){
        ccexRefusalPrint(stderr, "Error: configuration.k_range does not apply to the adaptive refinement (an adaptive block in the "
                                 "input), refine one configuration at a time\n");
        ccexRefuse();
    }
    for (long k = rangeFirst; k <= rangeLast && adaptiveStatus == 0; ++k){
        cJSON_setConfigurationIndex(ranged, k);
        adaptiveStatus = runInput(fccein, outputFolder, adaptiveArgs, ranged, k, rangeFirst, rangeLast, time_start);
    }

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

    MPI_Finalize();

    return adaptiveStatus;
}
