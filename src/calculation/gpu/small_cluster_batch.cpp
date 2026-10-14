// A26 step 1: the host side of the batched small-cluster cCCE on the GPU (small_cluster_batch.h says how it is used).
#include "../../../include/calculation/gpu/small_cluster_batch.h"
#include "../../../include/calculation/hamiltonian_check.h"
#include "../../../include/calculation/small_hamiltonian.h"
#include "../../../include/physical/pulse_bath.h"
#include "../../../include/calculation/cce_batch_window.h"
#include <algorithm>
#include <atomic>
#include <climits>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
#include "../../../include/calculation/gpu/small_cluster_gpu.h"
#endif

using smallcluster::C;
using smallcluster::Plan;

bool smallClusterPlan(Config* cnf, Pulse* pls, Plan& plan, bool& ensemble, bool& recursive){
    plan = Plan{};
    const int npulse = Pulse_getNpulse(pls);
    if (npulse+1 > smallcluster::maxIntervals) return false;
    ensemble = Config_getNstate(cnf)==0;
    recursive = !strcasecmp(cnf->evolution_update,"recursive");
    // The ensemble CPMG/Hahn/Ramsey branch of calCoherenceCceInto takes tau = t/(2 npulse) at both ends and twice that in
    // between, Ramsey t; a single state and a Manual sequence take the sequence's fractions.
    const bool manual = !strcasecmp(Pulse_getPulsename(pls),"Manual");
    double** sequence = Pulse_getSequence(pls);
    plan.intervals = npulse+1;
    plan.nstep = Config_getNstep(cnf);
    plan.deltat = (double)Config_getDeltat(cnf);
    for (int p=0; p<=npulse; ++p){
        double f;
        if (ensemble && !manual) f = npulse==0 ? 1.0 : (p==0 || p==npulse) ? 1.0/(2.0*npulse) : 1.0/npulse;
        else f = sequence[p][2];
        int slot = -1;
        for (int s=0; s<plan.distinct; ++s) if (plan.fraction[s]==f){ slot = s; break; }
        if (slot < 0){
            if (plan.distinct == 4) return false;
            slot = plan.distinct++;
            plan.fraction[slot] = f;
        }
        plan.same[p] = slot;
    }
    return true;
}

void smallClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result){
    static const char* path = getenv("CCEX_SMALL_CLUSTER_DUMP");
    // cCCE and pCCE (the same cluster calculation; pCCE since A30's check, 2026-10-13)
    if (!path || rank != 0 || (strcasecmp(cnf->method,"cce") && strcasecmp(cnf->method,"pcce"))) return;
    static long left = getenv("CCEX_SMALL_CLUSTER_DUMP_N") ? atol(getenv("CCEX_SMALL_CLUSTER_DUMP_N")) : 400;
    if (left <= 0) return;
    const int d = BathArray_dim(cluster);
    if (d != 2 && d != 4) return;
    Plan plan; bool ensemble, recursive;
    if (!smallClusterPlan(cnf,pls,plan,ensemble,recursive)) return;
    static FILE* fp = fopen(path,"ab");
    if (!fp){ perror(path); ccexFatal(); }
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,cluster,cnf,Ha,Hb);
    const int head[6] = {d, ensemble, recursive, plan.intervals, plan.distinct, plan.nstep};
    fwrite(head,sizeof(int),6,fp);
    fwrite(&plan.deltat,sizeof(double),1,fp);
    fwrite(plan.fraction,sizeof(double),plan.distinct,fp);
    fwrite(plan.same,sizeof(int),plan.intervals,fp);
    fwrite(Ha.data(),sizeof(doublec),size_t(d)*d,fp);
    fwrite(Hb.data(),sizeof(doublec),size_t(d)*d,fp);
    if (!ensemble){ const MatrixXcd psi = BathArray_Psi0(cluster); fwrite(psi.data(),sizeof(doublec),d,fp); }
    for (int i=0; i<plan.nstep; ++i){ const doublec v = result[i](0,0); fwrite(&v,sizeof(doublec),1,fp); }
    fflush(fp);
    if (--left == 0){ fclose(fp); fp = nullptr; }
}

bool warpClusterRun(QubitArray* qa, Config* cnf, Pulse* pls, warpcluster::Plan& plan, const char*& why){
    plan = warpcluster::Plan{};
    const bool gcce = !strcasecmp(cnf->method,"gcce");
    plan.mode = gcce ? warpcluster::modeGcce : warpcluster::modeCce;
    plan.ensemble = Config_getNstate(cnf)==0;
    plan.qdim = 1;
    if (Pulse_hasBathPulses(pls)){ why = "bath pulses (the kernel applies qubit pulses only)"; return false; }
    if (!gcce){
        if (strcasecmp(cnf->method,"cce")){ why = "the method is not CCE or gCCE"; return false; }
        if (QubitArray_getNqubit(qa)!=1 || strcasecmp(cnf->quantity,"coherence")){ why = "cCCE with one qubit and quantity=coherence only"; return false; }
        bool ensemble, recursive;
        if (!smallClusterPlan(cnf,pls,plan.time,ensemble,recursive)){ why = "more than 4 distinct pulse intervals"; return false; }
        return true;
    }
    // gCCE: what calCoherenceGcceInto does with the cluster, the warp kernel does only for these runs
    if (QubitArray_getNqubit(qa)!=1){ why = "gCCE with more than one qubit"; return false; }
    plan.qdim = QubitArray_dim(qa);
    if (plan.qdim > warpcluster::maxQubitDim){ why = "a qubit of dimension above 4"; return false; }
    if (strcasecmp(cnf->quantity,"coherence") && strcasecmp(cnf->quantity,"dm")){ why = "quantity is not coherence or dm"; return false; }
    plan.dm = !strcasecmp(cnf->quantity,"dm");
    if (strcasecmp(cnf->propagator,"eigen")){ why = "propagator is not eigen"; return false; }
    if (cnf->gcce_refinement){ why = "gcce_refinement"; return false; }
    if (plan.ensemble && strcasecmp(cnf->evolution,"matrix")){ why = "an ensemble with evolution=vector"; return false; }
    if (Pulse_getPulseiter(pls) || pls->detuning_factor != 1.0){ why = "pulseiter or a pulse detuning"; return false; }
    const int npulse = Pulse_getNpulse(pls);
    if (npulse+1 > warpcluster::maxIntervals){ why = "too many pulses"; return false; }
    // the fractions as propagateStateGcce reads them (sequence[p][2]; equal ones share their phases)
    double** sequence = Pulse_getSequence(pls);
    smallcluster::Plan& t = plan.time;
    t.intervals = npulse+1; t.nstep = Config_getNstep(cnf); t.deltat = (double)Config_getDeltat(cnf);
    for (int p=0; p<=npulse; ++p){
        const double f = sequence[p][2];
        int slot = -1;
        for (int k=0; k<t.distinct; ++k) if (t.fraction[k]==f){ slot = k; break; }
        if (slot < 0){
            if (t.distinct == 4){ why = "more than 4 distinct pulse intervals"; return false; }
            slot = t.distinct++;
            t.fraction[slot] = f;
        }
        t.same[p] = slot;
    }
    return true;
}

void warpClusterQubit(QubitArray* qa, Pulse* pls, warpcluster::Plan& plan){
    const int npulse = Pulse_getNpulse(pls), q = plan.qdim;
    std::vector<MatrixXcd> pulses(std::max(1,npulse));
    cachedUpulses(pulses.data(),qa,pls);
    for (int p=0; p<npulse; ++p) std::memcpy(plan.pulse[p],pulses[p].data(),sizeof(C)*q*q);
    const MatrixXcd psia = QubitArray_getPsia(qa), psib = QubitArray_getPsib(qa), psi0 = QubitArray_getPsi0(qa);
    std::memcpy(plan.psia,psia.data(),sizeof(C)*q);
    std::memcpy(plan.psib,psib.data(),sizeof(C)*q);
    std::memcpy(plan.qpsi0,psi0.data(),sizeof(C)*q);
}

void warpClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result){
    static const char* path = getenv("CCEX_WARP_CLUSTER_DUMP");
    if (!path || rank != 0) return;
    static const long most = getenv("CCEX_WARP_CLUSTER_DUMP_N") ? atol(getenv("CCEX_WARP_CLUSTER_DUMP_N")) : 400;
    static long taken[warpcluster::maxDim+1] = {};
    warpcluster::Plan plan; const char* why = nullptr;
    if (!warpClusterRun(qa,cnf,pls,plan,why)) return;
    const int d = (plan.mode==warpcluster::modeGcce ? plan.qdim : 1)*BathArray_dim(cluster);
    if (!warpcluster::compiled(plan.mode,plan.qdim,d) || taken[d] >= most) return;
    ++taken[d];
    static FILE* fp = fopen(path,"ab");
    if (!fp){ perror(path); ccexFatal(); }
    std::vector<C> in(warpcluster::inputSize(plan,d));
    if (plan.mode==warpcluster::modeGcce){
        warpClusterQubit(qa,pls,plan);
        MatrixXcd H, Hq;
        gcceClusterHamiltonian(qa,cluster,cnf,H,Hq);
        std::memcpy(in.data(),H.data(),sizeof(C)*d*d);
        if (!plan.ensemble){
            const MatrixXcd psi = kron(QubitArray_getPsi0(qa),BathArray_Psi0(cluster));
            std::memcpy(in.data()+d*d,psi.data(),sizeof(C)*d);
        }
    } else {
        MatrixXcd Ha, Hb;
        cceClusterHamiltonians(qa,cluster,cnf,Ha,Hb);
        std::memcpy(in.data(),Ha.data(),sizeof(C)*d*d);
        std::memcpy(in.data()+d*d,Hb.data(),sizeof(C)*d*d);
        if (!plan.ensemble){ const MatrixXcd psi = BathArray_Psi0(cluster); std::memcpy(in.data()+2*d*d,psi.data(),sizeof(C)*d); }
    }
    const smallcluster::Plan& t = plan.time;
    const int head[8] = {plan.mode, d, plan.qdim, plan.ensemble, plan.dm, t.intervals, t.distinct, t.nstep};
    fwrite(head,sizeof(int),8,fp);
    fwrite(&t.deltat,sizeof(double),1,fp);
    fwrite(t.fraction,sizeof(double),t.distinct,fp);
    fwrite(t.same,sizeof(int),t.intervals,fp);
    if (plan.mode==warpcluster::modeGcce){
        for (int p=0; p<t.intervals-1; ++p) fwrite(plan.pulse[p],sizeof(C),plan.qdim*plan.qdim,fp);
        fwrite(plan.psia,sizeof(C),plan.qdim,fp); fwrite(plan.psib,sizeof(C),plan.qdim,fp); fwrite(plan.qpsi0,sizeof(C),plan.qdim,fp);
    }
    fwrite(in.data(),sizeof(C),in.size(),fp);
    const int per = warpcluster::outputPer(plan);
    for (int i=0; i<t.nstep; ++i){
        if (per==1){ const doublec v = result[i](0,0); fwrite(&v,sizeof(doublec),1,fp); }
        else fwrite(result[i].data(),sizeof(doublec),per,fp);
    }
    fflush(fp);
}

// A40: the total dimensions the warp kernel takes under gpu_options.warp_max_dim: the compiled ones up to the bound, for
// the run's method and qubit dimension (3 before the qubit is read: the only one compiled for gCCE)
static unsigned long long warpMask(const Config* cnf, int qdim = 3){
    const int mode = strcasecmp(cnf->method,"gcce") ? warpcluster::modeCce : warpcluster::modeGcce;
    unsigned long long mask = 0;
    for (int d=2; d<=warpcluster::maxDim && d<=cnf->gpu_warp_max_dim; ++d)
        if (warpcluster::compiled(mode,qdim,d)) mask |= 1ULL << d;
    return mask;
}

// A40: the input checks of a given gpu_options.warp_max_dim that need no qubit (smallClusterWarmUp; the rest in the
// constructor). The default (Config_resolveWarpMaxDim) is nonzero only where these hold.
static void warpInputCheck(const Config* cnf){
    const bool gcce = !strcasecmp(cnf->method,"gcce");
    if (!gcce && strcasecmp(cnf->method,"cce")){
        if (rank==0) fprintf(stderr,"Error: gpu_options.warp_max_dim works with method=CCE or gCCE (this run: %s%s); 0 switches it off\n",
                             inputMethodName(cnf->method),methodUsesPcceClusters(cnf->method) ? ", cluster_mode pcce" : "");
        ccexFatal();
    }
    if (!strcasecmp(cnf->backend,"cpu")){
        if (rank==0) fprintf(stderr,"Error: gpu_options.warp_max_dim works with backend auto or cuda (this run: backend cpu); 0 switches it off\n");
        ccexFatal();
    }
    if (!warpMask(cnf)){
        if (rank==0) fprintf(stderr,"Error: gpu_options.warp_max_dim %d: the warp kernel has no dimension up to it for method=%s (it has %s)\n",
                             cnf->gpu_warp_max_dim,inputMethodName(cnf->method),gcce ? "6, 9, 12, 18, 21, 24 and 27 with a qubit of dimension 3" : "8 and 16");
        ccexFatal();
    }
}

#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
namespace {
// The state of one slot of one rank, in the node's shared memory: the rank packs and submits, the feeder runs it, the rank
// reads it back. std::atomic<int> is lock-free here, so it works between the processes of the node.
enum SlotState { slotIdle = 0, slotSubmitted = 1, slotRunning = 2, slotDone = 3 };
struct alignas(64) SlotHeader { std::atomic<int> state; long count[3]; };   // dimension 2, 4 (packed), warp kernel
struct alignas(64) FeederHeader { std::atomic<int> state; char message[248]; };   // in node rank 0's part: a FeederState
// -2: the card could not be opened and backend is auto - the batch is off and every cluster goes to the CPU (the auto
// contract, docs/methods/gpu.md; review 20261012-233058-452278 P3-01); -1: any other failure, which stops the run
enum FeederState { feederUnavailable = -2, feederFailed = -1, feederStarting = 0, feederReady = 1 };
static_assert(std::atomic<int>::is_always_lock_free, "the slot flags need lock-free atomics");
constexpr size_t headerBytes = 2*sizeof(SlotHeader)+sizeof(FeederHeader);
size_t roundUp(size_t n, size_t to){ return (n+to-1)/to*to; }

// Window: 2048 clusters, fewer when the outputs of a slot would pass 6 MiB - a budget that is kept (until review
// 20261012-233058-452278 P3-04 a floor of 256 raised it again above nstep 1536); 0 when the outputs of one cluster
// alone pass it, and the batch is then off. The shared memory is allocated, page-locked
// and freed once per run, about 0.25 s per GB each on idm4-gpu; windows of 2048, 4096 and 8192 ran P1 CR60 cCCE2 (2 states)
// in about the same time and 16384 and 32768 slower (A26 step 1, jobs 684-688).
int windowFor(int values){
    return int(std::min<long>(2048,(6L<<20)/(16L*std::max(1,values))));
}

// What a slot holds; the same on every rank (it depends on the input only). A40: the warp clusters' inputs go into an area
// of their own (arena values; a window closes early when it is full), their outputs (outPer values per point) after the
// packed ones.
struct Layout {
    int nstep = 0, windowSize = 0, per2 = 0, per4 = 0, outPer = 1;
    long arena = 0;
    bool operator==(const Layout& o) const {
        return nstep==o.nstep && windowSize==o.windowSize && per2==o.per2 && per4==o.per4 && outPer==o.outPer && arena==o.arena;
    }
};
// packed: the packed kernel runs; warpMask: the dimensions of the warp kernel (0: off). gCCE dm assumes the compiled qubit
// dimension 3 (a run with another one has no warp dimension and no warp kernel).
// cap: no window larger than this many clusters (the largest local cluster count of one order, review P3-04).
Layout layoutFor(const Config* cnf, bool packed, unsigned long long warpMask, long cap = 1L<<30){
    Layout l;
    const bool ensemble = cnf->nstate==0, gcce = !strcasecmp(cnf->method,"gcce");
    l.nstep = cnf->nstep;
    l.outPer = warpMask && gcce && !strcasecmp(cnf->quantity,"dm") ? 9 : 1;
    l.windowSize = int(std::min<long>(windowFor(l.nstep*l.outPer),cap));
    if (packed){ const int extra = ensemble ? 0 : 1; l.per2 = 8+2*extra; l.per4 = 32+4*extra; }
    if (warpMask){
        long most = 0;
        for (int d=2; d<=warpcluster::maxDim; ++d)
            if (warpMask>>d & 1ULL) most = std::max(most,(gcce ? 1L : 2L)*d*d + (ensemble ? 0 : d));
        l.arena = std::max(most,std::min(long(l.windowSize)*most,(6L<<20)/16));
    }
    return l;
}

// The node's shared memory: every rank's part is the headers, then two slots of (in2, in4, inW, out, status, dimW).
struct NodeMemory {
    MPI_Comm node = MPI_COMM_NULL; MPI_Win win = MPI_WIN_NULL;
    int nodeRank = 0, nodeSize = 1;
    Layout layout;
    size_t slotBytes = 0, partBytes = 0;
    std::vector<char*> parts;                  // every rank's part (as this process sees it)
    bool allocated() const { return win != MPI_WIN_NULL; }
    bool matches(const Layout& l) const { return allocated() && layout==l; }
    void allocate(const Layout& l){            // collective
        layout = l;
        const size_t w = size_t(l.windowSize);
        slotBytes = roundUp(sizeof(C)*(w*(l.per2+l.per4+size_t(l.nstep)*l.outPer)+size_t(l.arena))+2*w,64);
        partBytes = roundUp(headerBytes+2*slotBytes,4096);
        MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&nodeRank);
        MPI_Comm_size(node,&nodeSize);
        char* base = nullptr;
        MPI_Win_allocate_shared(MPI_Aint(partBytes),1,MPI_INFO_NULL,node,&base,&win);
        parts.resize(nodeSize);
        for (int r=0; r<nodeSize; ++r){
            MPI_Aint size; int unit; char* p;
            MPI_Win_shared_query(win,r,&size,&unit,&p);
            parts[r] = p;
        }
        char* own = parts[nodeRank];
        for (int s=0; s<2; ++s){ SlotHeader* h = new (own+s*sizeof(SlotHeader)) SlotHeader; h->state.store(slotIdle); h->count[0] = h->count[1] = h->count[2] = 0; }
        if (nodeRank==0){ FeederHeader* fh = new (own+2*sizeof(SlotHeader)) FeederHeader; fh->state.store(0); fh->message[0] = 0; }
        MPI_Barrier(node);
        // one contiguous range (MPI_Win_allocate_shared without alloc_shared_noncontig), page-locked in one call
        for (int r=1; r<nodeSize; ++r)
            if (parts[r] != parts[0]+r*partBytes){
                if (rank==0) fprintf(stderr,"Error: the small-cluster GPU batch: the node's shared memory is not contiguous\n");
                ccexFatal();
            }
    }
    void release(){                                          // collective
        if (!allocated()) return;
        MPI_Win_free(&win);
        MPI_Comm_free(&node);
        parts.clear();
    }
    size_t bytes() const { return partBytes*size_t(nodeSize); }
};

// smallClusterWarmUp: the shared memory allocated early, and on node rank 0 the CUDA context opened and that memory
// page-locked in a thread, while the bath is read and the clusters are made
struct Warm {
    std::thread thread;
    GpuSmallDevice* device = nullptr;
    char description[200] = {};
    double openSeconds = 0, lockSeconds = 0, allocateSeconds = 0;
    bool locked = false;
    NodeMemory memory;
    bool deferred = false;                     // A55: the packed start waits for the bath (smallClusterWarmUpAfterBath)
};
Warm warm;
}

struct SmallClusterBatch::Impl {
    QubitArray* qa; BathArray* ba; DefectArray* dfa; Config* cnf;
    Plan plan{};
    bool ensemble = false, recursive = false;
    int maxDim = 0;                            // gpu_options.packed_max_dim (0: the packed kernel is off)
    bool packed = false;
    bool optional = false;                     // backend auto: a card that cannot be opened turns the batch off
    bool off = false;                          // it could not (feederUnavailable seen): every cluster on the CPU
    // A40: the warp kernel
    bool warp = false, gcce = false;
    unsigned long long warpMask = 0;
    int qdim = 1;
    warpcluster::Plan wplan{};
    std::atomic<int> warpReady{0};             // wplan complete (cCCE: from the start; gCCE: after the first beginState)
    NodeMemory mem;
    bool fromWarm = false;                     // mem is the one smallClusterWarmUp made
    double allocateSeconds = 0;
    char* own = nullptr;
    // this rank's windows; kind of a cluster: 0 the CPU, 2 or 4 the packed kernel, 1 the warp kernel
    struct Window { int first = 0, last = 0; bool submitted = false, waited = false, cpu = false; std::vector<int> where; std::vector<char> dim, kind; long count[3] = {0,0,0}; };
    Window window[2];
    int order = 0, ncluster = 0; int** clusters = nullptr;
    MatrixXcd Ha, Hb, Ht, Hq;
    long long gpuClusters = 0, warpClusters = 0, fallback = 0, cpuStarting = 0;
    double prepareSeconds = 0, waitSeconds = 0;
    // the feeder (node rank 0)
    std::thread feeder;
    std::atomic<bool> stop{false};
    double feederOpen = 0, feederLock = 0, feederWaited = 0; long long feederBatches = 0, feederItems = 0; bool locked = false, early = false;

    const Layout& L() const { return mem.layout; }
    SlotHeader* header(char* part, int slot){ return reinterpret_cast<SlotHeader*>(part)+slot; }
    FeederHeader* feederHeader(){ return reinterpret_cast<FeederHeader*>(mem.parts[0]+2*sizeof(SlotHeader)); }
    C* in2(char* part, int slot){ return reinterpret_cast<C*>(part+headerBytes+slot*mem.slotBytes); }
    C* in4(char* part, int slot){ return in2(part,slot)+size_t(L().windowSize)*L().per2; }
    C* inW(char* part, int slot){ return in4(part,slot)+size_t(L().windowSize)*L().per4; }
    C* out(char* part, int slot){ return inW(part,slot)+size_t(L().arena); }
    unsigned char* status(char* part, int slot){ return reinterpret_cast<unsigned char*>(out(part,slot)+size_t(L().windowSize)*L().nstep*L().outPer); }
    unsigned char* dimW(char* part, int slot){ return status(part,slot)+L().windowSize; }

    void feed();
    void prepare(int slot, int first);
    void wait(int slot);
    // whether the feeder is up (its context opened, the memory page-locked); stops the run if it failed, sets off if the
    // card could not be opened under backend auto
    bool ready(){
        FeederHeader* fh = feederHeader();
        const int state = fh->state.load(std::memory_order_acquire);
        if (state == feederFailed){
            fprintf(stderr,"Error: the packed GPU kernel (feeder on node rank 0): %s\n",fh->message);
            ccexFatal();
        }
        if (state == feederUnavailable) off = true;
        return state == feederReady;
    }
    // While the feeder is starting, the next cpuWhileStarting clusters go to the CPU: a window that is never submitted
    void cpuWindow(int slot, int first){
        constexpr int cpuWhileStarting = 256;
        Window& w = window[slot];
        w.first = first; w.last = std::min(first+cpuWhileStarting, ncluster+1);
        w.where.assign(w.last-w.first,-1); w.dim.assign(w.last-w.first,0); w.kind.assign(w.last-w.first,0);
        w.count[0] = w.count[1] = w.count[2] = 0;
        w.submitted = w.waited = true; w.cpu = true;
        cpuStarting += w.last-w.first;
    }
};

// The feeder: runs on node rank 0 beside its main thread. Takes every submitted slot of the node into one batch (up to
// maxClusters clusters), in one of two lanes, and marks the slots done when the lane's batch is back.
void SmallClusterBatch::Impl::feed(){
    FeederHeader* fh = feederHeader();
    GpuSmallDevice* gpu = nullptr;
    try {
        const auto t0 = std::chrono::steady_clock::now();
        if (warm.thread.joinable()) warm.thread.join();      // the context (and the memory's page-lock) of smallClusterWarmUp
        gpu = warm.device; warm.device = nullptr;
        early = gpu != nullptr;
        feederOpen = warm.openSeconds;
        if (fromWarm && warm.locked){ locked = true; warm.locked = false; feederLock = warm.lockSeconds; }
        if (!gpu){
            char description[200] = {};
            gpu = gpuSmallOpen(0,description,sizeof description);
            if (!gpu && optional){
                snprintf(fh->message,sizeof fh->message,"%s",*warm.description ? warm.description : description);
                fh->state.store(feederUnavailable,std::memory_order_release);
                return;
            }
            if (!gpu) throw std::runtime_error(*warm.description ? warm.description : description);
        }
        const auto t1 = std::chrono::steady_clock::now();
        if (!locked){ locked = gpuSmallRegister(mem.parts[0],mem.bytes()); feederLock = std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count(); }
        feederWaited = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        fh->state.store(feederReady,std::memory_order_release);
        const long maxClusters = 1L<<19;
        std::vector<std::pair<int,int>> inLane[2];           // (rank, slot) of each lane's batch
        std::vector<GpuSmallItem> items;
        bool uploaded = !warp;
        int idle = 0;
        while (true){
            bool progress = false;
            for (int lane=0; lane<2; ++lane){
                if (!inLane[lane].empty() && gpuSmallReady(gpu,lane)){
                    for (auto& rs : inLane[lane]) header(mem.parts[rs.first],rs.second)->state.store(slotDone,std::memory_order_release);
                    inLane[lane].clear();
                    progress = true;
                }
            }
            // A40: nothing is submitted before the warp plan is on the card (gCCE: node rank 0's first beginState)
            if (!uploaded && warpReady.load(std::memory_order_acquire)){ gpuSmallWarpPlan(gpu,wplan); uploaded = true; }
            for (int lane=0; lane<2 && uploaded; ++lane){
                if (!inLane[lane].empty()) continue;
                items.clear();
                long total = 0;
                for (int r=0; r<mem.nodeSize; ++r) for (int s=0; s<2; ++s){
                    SlotHeader* h = header(mem.parts[r],s);
                    if (h->state.load(std::memory_order_acquire) != slotSubmitted) continue;
                    const long n = h->count[0]+h->count[1]+h->count[2];
                    if (!items.empty() && total+n > maxClusters) continue;
                    h->state.store(slotRunning,std::memory_order_relaxed);
                    items.push_back(GpuSmallItem{in2(mem.parts[r],s),in4(mem.parts[r],s),h->count[0],h->count[1],
                                                 out(mem.parts[r],s),status(mem.parts[r],s),
                                                 inW(mem.parts[r],s),dimW(mem.parts[r],s),h->count[2]});
                    inLane[lane].push_back({r,s});
                    total += n;
                }
                if (items.empty()) break;
                gpuSmallSubmit(gpu,lane,plan,ensemble,recursive,items.data(),int(items.size()));
                ++feederBatches; feederItems += long(items.size());
                progress = true;
            }
            if (progress){ idle = 0; continue; }
            if (stop.load(std::memory_order_acquire) && inLane[0].empty() && inLane[1].empty()) break;
            // idle: 20 us between looks, up to 200 us after a while (it shares node rank 0's core)
            ++idle;
            std::this_thread::sleep_for(std::chrono::microseconds(std::min(200,20*(1+idle/64))));
        }
        if (locked) gpuSmallUnregister(mem.parts[0]);
        gpuSmallClose(gpu);
    } catch (const std::exception& e){
        snprintf(fh->message,sizeof fh->message,"%s",e.what());
        fh->state.store(feederFailed,std::memory_order_release);
        if (gpu) gpuSmallClose(gpu);
    }
}

void SmallClusterBatch::Impl::prepare(int slot, int first){
    const double start = MPI_Wtime();
    Window& w = window[slot];
    w.first = first; w.last = std::min(first+L().windowSize, ncluster+1);
    w.where.assign(w.last-w.first,-1); w.dim.assign(w.last-w.first,0); w.kind.assign(w.last-w.first,0);
    w.count[0] = w.count[1] = w.count[2] = 0;
    w.cpu = false;
    // the slot was waited for (or never used): the feeder is done with it
    C* dst[2] = {in2(own,slot), in4(own,slot)};
    C* arena = inW(own,slot);
    unsigned char* dims = dimW(own,slot);
    long arenaUsed = 0;
    const int nqubit = QubitArray_getNqubit(qa);
    for (int ic=w.first; ic<w.last; ++ic){
        BathArray* b = createBathArray(clusters[ic]+1,order,ba,dfa,nqubit);
        const int d = BathArray_dim(b);
        const int dw = (gcce ? qdim : 1)*d;   // the warp kernel's total dimension
        if (packed && (d == 2 || d == 4) && d <= maxDim){
            cceClusterHamiltonians(qa,b,cnf,Ha,Hb);
            if (finiteAndHermitian(Ha) && finiteAndHermitian(Hb)){
                const int k = d==4, per = k ? L().per4 : L().per2;
                C* p = dst[k]+w.count[k]*per;
                std::memcpy(p,Ha.data(),sizeof(C)*d*d);
                std::memcpy(p+d*d,Hb.data(),sizeof(C)*d*d);
                if (!ensemble){ const MatrixXcd psi = BathArray_Psi0(b); std::memcpy(p+2*d*d,psi.data(),sizeof(C)*d); }
                w.where[ic-w.first] = int(w.count[k]++);
                w.dim[ic-w.first] = char(d);
                w.kind[ic-w.first] = char(d);
            } else ++fallback;
        } else if (warp && dw <= warpcluster::maxDim && (warpMask>>dw & 1ULL)){
            const long need = warpcluster::inputSize(wplan,dw);
            if (arenaUsed+need > L().arena){       // the area is full: the window ends before this cluster
                BathArray_freeAll(b);
                w.last = ic;
                break;
            }
            C* p = arena+arenaUsed;
            bool fine;
            if (gcce){
                gcceClusterHamiltonian(qa,b,cnf,Ht,Hq);
                double asym = 0, scale = 0;
                // otherwise calCoherenceGcceInto refuses it, with its message
                fine = gcceHamiltonianFinite(Ht) && gcceHamiltonianHermitian(Ht,asym,scale);
                if (fine){
                    std::memcpy(p,Ht.data(),sizeof(C)*dw*dw);
                    if (!ensemble){ const MatrixXcd psi = kron(QubitArray_getPsi0(qa),BathArray_Psi0(b)); std::memcpy(p+dw*dw,psi.data(),sizeof(C)*dw); }
                }
            } else {
                cceClusterHamiltonians(qa,b,cnf,Ha,Hb);
                fine = finiteAndHermitian(Ha) && finiteAndHermitian(Hb);
                if (fine){
                    std::memcpy(p,Ha.data(),sizeof(C)*d*d);
                    std::memcpy(p+d*d,Hb.data(),sizeof(C)*d*d);
                    if (!ensemble){ const MatrixXcd psi = BathArray_Psi0(b); std::memcpy(p+2*d*d,psi.data(),sizeof(C)*d); }
                }
            }
            if (fine){
                dims[w.count[2]] = (unsigned char)dw;
                w.where[ic-w.first] = int(w.count[2]++);
                w.dim[ic-w.first] = char(d);
                w.kind[ic-w.first] = 1;
                arenaUsed += need;
                ++warpClusters;
            } else ++fallback;
        } else ++fallback;
        BathArray_freeAll(b);
    }
    w.where.resize(w.last-w.first); w.dim.resize(w.last-w.first); w.kind.resize(w.last-w.first);
    SlotHeader* h = header(own,slot);
    h->count[0] = w.count[0]; h->count[1] = w.count[1]; h->count[2] = w.count[2];
    w.submitted = true; w.waited = false;
    const long n = w.count[0]+w.count[1]+w.count[2];
    h->state.store(n > 0 ? slotSubmitted : slotDone,std::memory_order_release);
    gpuClusters += n;
    prepareSeconds += MPI_Wtime()-start;
}

void SmallClusterBatch::Impl::wait(int slot){
    const double start = MPI_Wtime();
    SlotHeader* h = header(own,slot);
    FeederHeader* fh = feederHeader();
    long spins = 0;
    while (h->state.load(std::memory_order_acquire) != slotDone){
        if (fh->state.load(std::memory_order_acquire) == feederFailed){
            fprintf(stderr,"Error: the packed GPU kernel (feeder on node rank 0): %s\n",fh->message);
            ccexFatal();
        }
        if (++spins > 64) std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
    window[slot].waited = true;
    waitSeconds += MPI_Wtime()-start;
}
#else
struct SmallClusterBatch::Impl {};
#endif

#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
static void warmStart(const Config* cnf, bool packedLikely, unsigned long long mask);

// A48: one state, packed_max_dim not given and a run the CPU batch takes - the packed kernel only when the clusters of bath
// dimension 2 and 4 reach this many per rank (the constructor; smallClusterWarmUpAfterBath counts them before the clusters)
constexpr long packedMinPerRankState = 40000;

// A48: whether the CPU batch (A30, CceBatchWindow) takes this run's clusters of the packed dimensions (both when
// packed_max_dim is 4). Collective the first time (ccebatchRunTakes -> ccebatchUsable). The warm-up and the constructor
// can get different answers: ccebatchRunTakes asks computeUseCuda, false in the warm-up (before the card is started).
static bool cpuBatchTakesPacked(QubitArray* qa, Config* cnf, Pulse* pls, bool isGCCE, const char** why){
    ccebatch::Plan cpuPlan; bool cpuDims[2];
    return ccebatchRunTakes(qa,cnf,pls,isGCCE,cpuPlan,cpuDims,why) && cpuDims[0] && (cnf->gpu_packed_max_dim < 4 || cpuDims[1]);
}
#endif

void smallClusterWarmUp(const Config* cnf, QubitArray* qa, Pulse* pls, const char* savemode){
    const bool packed = cnf->gpu_packed_max_dim > 0, warp = cnf->gpu_warp_max_dim > 0;
    if (!packed && !warp) return;
    // refused here, before the bath is read: the packed kernel is cCCE, and it is in the GPU binary only
    if (packed && strcasecmp(cnf->method,"cce")){
        if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim works with method=CCE (this run: %s%s)\n",inputMethodName(cnf->method),
                             methodUsesPcceClusters(cnf->method) ? ", cluster_mode pcce" : "");
        ccexFatal();
    }
    if (warp && !strcmp(cnf->gpu_warp_source,"given")) warpInputCheck(cnf);
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    // the batch runs threads beside the main one (this warm-up, the feeder), which MPI allows from MPI_THREAD_FUNNELED on
    // (main.cpp asks for it); every MPI call stays on the main thread (review 20261012-233058-452278 P3-03)
    int provided = MPI_THREAD_SINGLE;
    MPI_Query_thread(&provided);
    if (provided < MPI_THREAD_FUNNELED){
        if (rank==0) fprintf(stderr,"Error: gpu_options.%s runs a thread beside MPI and needs MPI_THREAD_FUNNELED; "
                                    "this MPI provides level %d. Set gpu_options.packed_max_dim to 0%s.\n",
                             packed ? "packed_max_dim" : "warp_max_dim",provided,warp ? " and warp_max_dim to 0" : "");
        ccexFatal();
    }
    // review 20261013-221715 X7: no early start for a kernel the input already rules out - what the constructor will find
    // without looking at the clusters: the packed kernel outside its runs (one qubit, coherence, at most 4 distinct
    // intervals) or off by A48 (not given, the ensemble, and the CPU batch takes the run); the warp kernel where
    // warpClusterRun refuses the run (gcce_refinement, propagator, pulseiter, ...) or savemode cluster_density_matrix (gCCE),
    // or without a dimension for the qubit. A wrong guess here costs only time: the constructor starts the card itself.
    // Whether a run is short (its clusters) is known only after the bath is read; the early start stays for those
    // (deferring every start to the constructor slowed P1 CR0 cCCE2 8.63 -> 8.92 s, job 707).
    bool packedLikely = packed, warpLikely = warp, deferPacked = false;
    int qdim = 3;
    if (qa && pls){
        Config* c = const_cast<Config*>(cnf);
        if (packedLikely){
            Plan plan{}; bool ensemble = false, recursive = false;
            packedLikely = QubitArray_getNqubit(qa)==1 && !strcasecmp(cnf->quantity,"coherence") &&
                           smallClusterPlan(c,pls,plan,ensemble,recursive);
            const char* cpuWhy = nullptr;
            if (packedLikely && strcmp(cnf->gpu_packed_source,"given") && cpuBatchTakesPacked(qa,c,pls,false,&cpuWhy)){
                // A55 (user 2026-10-14): one state - the constructor's threshold (packedMinPerRankState) counts the
                // clusters of bath dimension 2 and 4, which in a standard cCCE run are the bath spins and the pairs
                // closer than rdip: counted right after the bath is read (smallClusterWarmUpAfterBath), and the card
                // starts there only when they reach it (P1 CR0 30 states: +0.2 s for a card it did not use, job 17445)
                if (ensemble) packedLikely = false;
                else if (!strcasecmp(cnf->method,"cce") && strcasecmp(cnf->cluster_mode,"adaptive") && cnf->nrdip_ladder == 0)
                    deferPacked = true;
            }
        }
        if (warpLikely){
            warpcluster::Plan wp{}; const char* why = nullptr;
            const bool gcce = !strcasecmp(cnf->method,"gcce");
            warpLikely = warpClusterRun(qa,c,pls,wp,why) && !(gcce && savemode && !strcasecmp(savemode,"cluster_density_matrix")) &&
                         warpMask(cnf,wp.qdim);
            if (warpLikely) qdim = wp.qdim;
            // A53 (user 2026-10-13): cCCE with the packed kernel off (packed_max_dim 0, its plan, the A48 ensemble rule) -
            // whether a cluster of a warp dimension (bath 8, 16) exists is known only from the bath, and most cCCE baths
            // have none (P1, h-BN ensembles: 0.2-0.35 s of start for nothing, job 17443), so the card starts in the
            // constructor when some rank has one. With the packed kernel on, the early start stays (job 707).
            if (warpLikely && !gcce && !packedLikely){
                warpLikely = false;
                if (rank==0) printf("[gpu-warp] no early start: cCCE with the packed kernel off - the card starts only if some rank "
                                    "has a cluster of a warp dimension (A53)\n");
            }
        }
        if (!packedLikely && !warpLikely) return;
        if (deferPacked){ warm.deferred = true; return; }
    }
    // A40 follow-up (2026-10-13): with the packed kernel on, the early memory has no warp input region - most cCCE runs have
    // no cluster of a warp dimension (P1 CR0 cCCE2: +4.7 % from that region, job 708); the constructor adds it when some
    // rank has one (a new memory, the early CUDA context kept)
    warmStart(cnf,packedLikely,warpLikely && !packedLikely ? warpMask(cnf,qdim) : 0);
#else
    if (rank==0) fprintf(stderr,"Error: gpu_options.%s needs the GPU binary (bin/main_gpu.out, bash do_compile.sh gpu)\n",
                         packed ? "packed_max_dim" : "warp_max_dim");
    ccexFatal();
#endif
}

#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
// the bath dimension of spin ibs: its own, times its defect's added spins, as createBathArray makes the cluster
static int bathSpinDim(BathArray* ba, DefectArray* dfa, int ibs){
    int d = BathSpin_dim(BathArray_getBath_i(ba,ibs));
    const int idf = DefectArray_findDefectIndex(dfa,BathArray_getBath_i_name(ba,ibs));
    if (idf >= 0)
        for (int isp=0; isp<DefectArray_getDefect_idf_naddspin(dfa,idf); isp++) d *= BathSpin_dim(DefectArray_getSubbath_i_isp(dfa,ibs,isp));
    return d;
}

// collective: the early shared memory, and on node rank 0 the CUDA context and that memory's page-lock in a thread
static void warmStart(const Config* cnf, bool packedLikely, unsigned long long mask){
    if (warm.memory.allocated() || layoutFor(cnf,packedLikely,mask).windowSize == 0) return;
    const double start = MPI_Wtime();
    warm.memory.allocate(layoutFor(cnf,packedLikely,mask));
    warm.allocateSeconds = MPI_Wtime()-start;
    if (warm.memory.nodeRank != 0) return;
    // the range by value: SmallClusterBatch takes warm.memory over (and clears it) while this thread may still run
    char* const base = warm.memory.parts[0];
    const size_t bytes = warm.memory.bytes();
    warm.thread = std::thread([base,bytes]{
        const auto t0 = std::chrono::steady_clock::now();
        warm.device = gpuSmallOpen(0,warm.description,sizeof warm.description);
        const auto t1 = std::chrono::steady_clock::now();
        if (warm.device) warm.locked = gpuSmallRegister(base,bytes);
        warm.openSeconds = std::chrono::duration<double>(t1-t0).count();
        warm.lockSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
    });
}
#endif

void smallClusterWarmUpAfterBath(const Config* cnf, BathArray* ba, DefectArray* dfa){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!warm.deferred) return;
    warm.deferred = false;
    // the constructor's count before the clusters: the spins (order 1) and the pairs rdipcut < r < rdip (order 2,
    // clusterizePairs and BathArray_connectivity) of bath dimension 2 and 4, rows shared between the ranks; a rank stops
    // at the threshold. Above about 2e7 distances per rank (40 ms) the count is not made and the card starts as before.
    const long cap = packedMinPerRankState*long(nprocess);
    const int nspin = BathArray_getNspin(ba), maxDim = cnf->gpu_packed_max_dim;
    const float rdip = cnf->rdip, rdipcut = cnf->rdipcut;
    if (0.5*nspin*double(nspin)/nprocess > 2e7){ warmStart(cnf,true,0); return; }
    std::vector<int> dim(nspin);
    for (int i=0; i<nspin; i++) dim[i] = bathSpinDim(ba,dfa,i);
    long count = 0;
    for (int i=rank; i<nspin && count<cap; i+=nprocess){
        if (cnf->order >= 1 && (dim[i] == 2 || dim[i] == 4) && dim[i] <= maxDim) ++count;
        if (cnf->order < 2 || dim[i] != 2 || maxDim < 4) continue;
        for (int j=i+1; j<nspin; j++){
            if (dim[j] != 2) continue;
            const double r = BathArray_dist_i_j(ba,i,j);
            if (r < rdip && r > rdipcut) ++count;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&count,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    if (count >= cap){ warmStart(cnf,true,0); return; }
    if (rank==0) printf("[gpu-packed] no early start: one state, %.0f clusters of bath dimension 2 and 4 per rank in the bath "
                        "read, fewer than %ld - the card starts only if the constructor finds more (A55)\n",
                        double(count)/nprocess,packedMinPerRankState);
#else
    (void)cnf; (void)ba; (void)dfa;
#endif
}


#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
// collective: the early memory released (its page-lock undone), the early CUDA context kept for the feeder
static void warmReleaseMemory(){
    if (warm.thread.joinable()) warm.thread.join();
    if (warm.locked){ gpuSmallUnregister(warm.memory.parts[0]); warm.locked = false; }
    warm.memory.release();
}
#endif

// The packed and warp kernels stay off (the card is not started) when the cluster evaluations (clusters x states) they
// would take, over all ranks, are fewer than this many per rank (A40 follow-up, 2026-10-13; the reasoning in
// (internal notes)): starting the card cost a full node of 32 ranks about 0.5 s
// (0.3-0.85, c4_ens, jobs 707-708), and it saves about 4.4 rank-us per packed evaluation and 21 per warp one (P1 CR0 cCCE2,
// P1 CR60 gCCE2) - break-even 23-66 k per rank; 20 k keeps on the one measured run near break-even (A26 P1 CR0 ensemble,
// about 25 k per rank, 0.95 of the CPU build). Only when neither kernel was given in gpu_options (a given value starts
// the card whatever the count).
constexpr long batchMinClustersPerRank = 20000;

void smallClusterWarmDown(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (warm.thread.joinable()) warm.thread.join();
    if (warm.locked){ gpuSmallUnregister(warm.memory.parts[0]); warm.locked = false; }
    gpuSmallClose(warm.device);
    warm.device = nullptr;
    warm.deferred = false;
    warm.memory.release();                     // collective: every rank made it in smallClusterWarmUp, or none did
#endif
}

SmallClusterBatch::SmallClusterBatch(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE,
                                     int*** localclusters, int order, const char* savemode){
    if (cnf->gpu_packed_max_dim <= 0 && cnf->gpu_warp_max_dim <= 0) return;
    // the packed kernel (A26), as before
    bool packedOn = false;
    Plan plan{}; bool ensemble = false, recursive = false;
    if (cnf->gpu_packed_max_dim > 0){
        if (isGCCE || strcasecmp(cnf->method,"cce") || QubitArray_getNqubit(qa)!=1 || strcasecmp(cnf->quantity,"coherence")){
            // given in the input: refused; the default (Config_resolvePackedMaxDim) is simply off here
            if (strcmp(cnf->gpu_packed_source,"given")){
                if (rank==0) printf("[gpu-packed] off: the packed kernel is for method=CCE with one qubit and quantity=coherence\n");
            } else {
                if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim supports method=CCE with one qubit and quantity=coherence\n");
                ccexFatal();
            }
        } else if (!smallClusterPlan(cnf,pls,plan,ensemble,recursive)){
            if (rank==0) printf("[gpu-packed] more than 4 distinct pulse intervals (or %d pulses): every cluster on the CPU\n",
                                smallcluster::maxIntervals-1);
        } else packedOn = true;
    }
    // A40: the warp kernel (default 27 in the GPU build, Config_resolveWarpMaxDim); what it cannot take stays on its usual path
    bool warpOn = false;
    warpcluster::Plan wp{};
    unsigned long long mask = 0;
    if (cnf->gpu_warp_max_dim > 0){
        const char* why = nullptr;
        if (!warpClusterRun(qa,cnf,pls,wp,why)){
            if (rank==0) printf("[gpu-warp] off: %s - every cluster on its usual path\n",why);
        } else if (isGCCE && savemode && !strcasecmp(savemode,"cluster_density_matrix")){
            if (rank==0) printf("[gpu-warp] off: savemode cluster_density_matrix - every cluster on its usual path\n");
        } else {
            mask = warpMask(cnf,wp.qdim);
            if (!mask){
                if (rank==0) printf("[gpu-warp] off: no warp dimension for a qubit of dimension %d - every cluster on its usual path\n",wp.qdim);
            } else warpOn = true;
        }
    }
    if (!packedOn && !warpOn) return;            // the early memory, if any, is freed at the end (smallClusterWarmDown)
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    // bath dimension of a spin: the product of the spins' dimensions, a defect's added spins included, as createBathArray
    // makes the cluster
    std::vector<int> spinDim(BathArray_getNspin(ba),0);
    auto dimOf = [&](int ibs){
        int& d = spinDim[ibs];
        if (!d){
            d = BathSpin_dim(BathArray_getBath_i(ba,ibs));
            const int idf = DefectArray_findDefectIndex(dfa,BathArray_getBath_i_name(ba,ibs));
            if (idf >= 0)
                for (int isp=0; isp<DefectArray_getDefect_idf_naddspin(dfa,idf); isp++) d *= BathSpin_dim(DefectArray_getSubbath_i_isp(dfa,ibs,isp));
        }
        return d;
    };
    // A48 (user 2026-10-13): with packed_max_dim not given and the CPU batch (A30, CceBatchWindow) taking this run's
    // clusters of the packed dimensions, the packed kernel only where it is faster than A30 (idm4-gpu, 32 ranks, MPS, jobs
    // 17403/17404/17407: per state and rank the card loses about 0.105 s and saves 2.65 us per evaluation - break-even
    // 39.6 k): one state - on when the clusters of bath dimension 2 and 4 (up to packed_max_dim) are at least
    // packedMinPerRankState per rank (each state evaluates them all), else off; the ensemble - off (no gain at 25-178 k).
    // A run A30 does not take keeps the threshold below. Collective (the same on every rank: the input decides).
    if (packedOn && strcmp(cnf->gpu_packed_source,"given")){
        const char* cpuWhy = nullptr;
        const bool cpuTakes = cpuBatchTakesPacked(qa,cnf,pls,isGCCE,&cpuWhy);
        if (!cpuTakes){
            if (rank==0) printf("[gpu-packed] the CPU batch (A30) does not take this run (%s): the packed kernel unless its "
                                "cluster evaluations (clusters x states) are fewer than %ld per rank (A48)\n",
                                cpuWhy ? cpuWhy : "not every packed dimension",batchMinClustersPerRank);
        } else if (ensemble){
            packedOn = false;
            if (rank==0) printf("[gpu-packed] off: the ensemble - the CPU batch (A30, %s) computes the clusters of bath dimension "
                                "2 and 4 (the card saves nothing there; A48)\n",ccebatch::variant());
        } else {
            const long cap = packedMinPerRankState*long(nprocess);
            long count = 0;
            for (int n=1; n<=order && count<cap; n++)
                for (int ic=1; ic<=localclusters[n][0][0]-1 && count<cap; ic++){
                    long d = 1;
                    for (int k=0; k<n && d <= 4; k++) d *= dimOf(localclusters[n][ic][1+k]);
                    if ((d == 2 || d == 4) && d <= cnf->gpu_packed_max_dim) ++count;
                }
            MPI_Allreduce(MPI_IN_PLACE,&count,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
            packedOn = count >= cap;
            if (rank==0){
                if (packedOn) printf("[gpu-packed] on: one state, at least %ld evaluations of bath dimension 2 and 4 per rank per "
                                     "state (A48 threshold; fewer go to the CPU batch, A30)\n",packedMinPerRankState);
                else printf("[gpu-packed] off: one state, %.0f evaluations of bath dimension 2 and 4 per rank per state, fewer than "
                            "%ld - the CPU batch (A30, %s) computes them (A48)\n",double(count)/nprocess,packedMinPerRankState,
                            ccebatch::variant());
            }
        }
        if (!packedOn && !warpOn){ smallClusterWarmDown(); return; }   // collective: the early start released now
    }
    const char* tag = packedOn ? "[gpu-packed]" : "[gpu-warp]";
    // review 20261012-233058-452278 P3-04: no window (and no memory) without clusters, none beyond the 6 MiB budget;
    // GPU 다음 셋 (1): and no feeder (the early start released) when no rank has a cluster of a dimension the batch takes.
    // A40 follow-up (2026-10-13): the warp input region only when some rank has a warp cluster, and the batch off when the
    // clusters it would take are fewer than the threshold (batchMinClustersPerRank). The count is of cluster evaluations,
    // a cluster once per state (nstate states, or 1 for the ensemble), as the batch's own count and the saving per cluster
    // are. A rank counts its clusters from order 1 up and stops once it has the threshold and has seen a warp cluster.
    const bool given = (packedOn && !strcmp(cnf->gpu_packed_source,"given")) || (warpOn && !strcmp(cnf->gpu_warp_source,"given"));
    const long threshold = given ? 0 : batchMinClustersPerRank*long(nprocess), need = std::max(threshold,1L);
    const long states = std::max(Config_getNstate(cnf),1);
    long most = 0, taken = 0, warpTaken = 0;
    const int qubitFactor = isGCCE ? wp.qdim : 1;
    for (int n=1; n<=order; n++) most = std::max<long>(most, localclusters[n][0][0]-1);
    for (int n=1; n<=order && !(taken*states >= need && (warpTaken || !warpOn)); n++)
        for (int ic=1; ic<=localclusters[n][0][0]-1; ic++){
            long d = 1;
            for (int k=0; k<n && d <= warpcluster::maxDim; k++) d *= dimOf(localclusters[n][ic][1+k]);
            const long dw = qubitFactor*d;
            if (packedOn && (d == 2 || d == 4) && d <= cnf->gpu_packed_max_dim) ++taken;
            else if (warpOn && dw <= warpcluster::maxDim && (mask>>dw & 1ULL)){ ++taken; ++warpTaken; }
            if (taken*states >= need && (warpTaken || !warpOn)) break;
        }
    long counts[2] = {most, warpTaken}, total = std::min(taken*states,need);
    MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&total,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    most = counts[0];
    if (!counts[1]) mask = 0;                  // no rank has a warp cluster: no warp kernel, no warp input region
    const bool warpUsed = warpOn && mask;
    const Layout full = layoutFor(cnf,packedOn,mask);
    if (most <= 0 || !total || total < threshold || full.windowSize == 0){
        if (rank==0){
            if (most <= 0) printf("%s off: no clusters to compute\n",tag);
            else if (!total) printf("%s off: no cluster of a dimension it takes (%s): every cluster on its usual path\n",
                                    tag,packedOn && warpOn ? "packed and warp kernel" : packedOn ? "packed kernel" : "warp kernel");
            else if (total < threshold)
                printf("%s off: %ld cluster evaluations (clusters x states) for the %s, fewer than the threshold %ld (%ld per rank "
                       "x %d ranks; starting the card would cost more than it saves): every cluster on its usual path\n",tag,total,
                       packedOn && warpUsed ? "packed and warp kernels" : packedOn ? "packed kernel" : "warp kernel",
                       threshold,batchMinClustersPerRank,nprocess);
            else printf("%s off: the %d values of one cluster pass the 6 MiB of a window: every cluster on the CPU\n",
                        tag,full.nstep*full.outPer);
        }
        smallClusterWarmDown();                // collective: the early memory and context (if any) released now, not at the end
        return;
    }
    if (warpOn && !warpUsed && rank==0) printf("[gpu-warp] off: no cluster of a dimension it takes - no warp input region\n");
    warpOn = warpUsed;
    Impl* m = new Impl;
    m->qa = qa; m->ba = ba; m->dfa = dfa; m->cnf = cnf;
    m->plan = packedOn ? plan : wp.time; m->ensemble = packedOn ? ensemble : bool(wp.ensemble); m->recursive = recursive;
    m->packed = packedOn; m->maxDim = packedOn ? cnf->gpu_packed_max_dim : 0;
    m->optional = !strcasecmp(cnf->backend,"auto");
    m->warp = warpOn; m->warpMask = mask; m->gcce = isGCCE; m->qdim = wp.qdim; m->wplan = wp;
    if (warpOn && !isGCCE) m->warpReady.store(1,std::memory_order_release);   // the cCCE plan needs no qubit state
    // the memory of smallClusterWarmUp when it fits (the same on every rank: it depends on the input only), else a new one
    // whose window is no larger than any rank's largest number of clusters in one order (the early CUDA context kept)
    if (warm.memory.matches(full)){ m->mem = warm.memory; warm.memory = NodeMemory{}; m->fromWarm = true; m->allocateSeconds = warm.allocateSeconds; }
    else {
        if (warm.memory.allocated()) warmReleaseMemory();   // collective, as the test above is the same on every rank
        const double start = MPI_Wtime();
        m->mem.allocate(layoutFor(cnf,packedOn,mask,most));
        m->allocateSeconds = MPI_Wtime()-start;
    }
    m->own = m->mem.parts[m->mem.nodeRank];
    if (m->mem.nodeRank==0) m->feeder = std::thread([m]{ m->feed(); });
    impl = m;
    if (rank==0 && packedOn) printf("[gpu-packed] on: cCCE clusters of bath dimension %s on the GPU (%s, %s, %d intervals, %d distinct), "
                        "window %d, %.1f MiB shared memory per rank; one feeder per node\n", m->maxDim==4 ? "2 and 4" : "2",
                        ensemble ? "ensemble" : "single state", recursive ? "recursive" : "exact", plan.intervals, plan.distinct,
                        m->mem.layout.windowSize, m->mem.partBytes/1048576.0);
    if (rank==0 && warpOn){
        char dims[160] = ""; size_t at = 0;
        for (int d=2; d<=warpcluster::maxDim; ++d) if (mask>>d & 1ULL) at += snprintf(dims+at,sizeof(dims)-at,"%s%d",at ? ", " : "",d);
        printf("[gpu-warp] on (gpu_options.warp_max_dim %d, %s): %s clusters of total dimension %s in the warp kernel (%s%s, exact, "
               "%d intervals, %d distinct), window %d, %.1f MiB of inputs per window, %.1f MiB shared memory per rank\n",
               cnf->gpu_warp_max_dim, cnf->gpu_warp_source, isGCCE ? "gCCE" : "cCCE", dims, wp.ensemble ? "ensemble" : "single state", wp.dm ? ", dm" : "",
               wp.time.intervals, wp.time.distinct, m->mem.layout.windowSize, m->mem.layout.arena*16.0/1048576.0, m->mem.partBytes/1048576.0);
    }
#else
    (void)ba; (void)dfa; (void)localclusters; (void)order;
    if (rank==0) fprintf(stderr,"Error: gpu_options.%s needs the GPU binary (bin/main_gpu.out, bash do_compile.sh gpu)\n",
                         packedOn ? "packed_max_dim" : "warp_max_dim");
    ccexFatal();
#endif
}

SmallClusterBatch::~SmallClusterBatch(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (impl && impl->feeder.joinable()){ impl->stop.store(true); impl->feeder.join(); }
#endif
    delete impl;
}

void SmallClusterBatch::beginState(QubitArray* qa, Pulse* pls){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl || !impl->warp || !impl->gcce) return;
    Impl& m = *impl;
    warpcluster::Plan p = m.wplan;
    warpClusterQubit(qa,pls,p);
    if (!m.warpReady.load(std::memory_order_acquire)){
        m.wplan = p;
        m.warpReady.store(1,std::memory_order_release);
        return;
    }
    if (std::memcmp(p.pulse,m.wplan.pulse,sizeof p.pulse) || std::memcmp(p.psia,m.wplan.psia,sizeof p.psia) ||
        std::memcmp(p.psib,m.wplan.psib,sizeof p.psib) || std::memcmp(p.qpsi0,m.wplan.qpsi0,sizeof p.qpsi0)){
        fprintf(stderr,"Error: the warp kernel (gpu_options.warp_max_dim): the pulses or the qubit states changed between states\n");
        ccexFatal();
    }
#else
    (void)qa; (void)pls;
#endif
}

// the batch's clusters (packed: bath dimension 2 and 4 up to packed_max_dim; warp: the compiled total dimensions up to warp_max_dim - dim
// is the total dimension, qubit x bath for gCCE): the others keep the ordinary CUDA cluster batching (review P1-05)
bool SmallClusterBatch::packs(int dim) const {
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl || impl->off) return false;
    return (impl->packed && (dim == 2 || dim == 4) && dim <= impl->maxDim) ||
           (impl->warp && dim >= 0 && dim <= warpcluster::maxDim && (impl->warpMask>>dim & 1ULL));
#else
    (void)dim;
    return false;
#endif
}

int SmallClusterBatch::cpuEnd(int ic) const {
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl || impl->off || !impl->packed) return INT_MAX;
    for (const Impl::Window& w : impl->window)
        if (w.submitted && w.cpu && ic >= w.first && ic < w.last) return w.last;
    return ic+1;
#else
    (void)ic;
    return INT_MAX;
#endif
}

void SmallClusterBatch::beginOrder(int n, int** list, int count){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl) return;
    for (int s=0; s<2; ++s){
        if (impl->window[s].submitted && !impl->window[s].waited) impl->wait(s);
        impl->window[s] = Impl::Window{};
    }
    impl->order = n; impl->clusters = list; impl->ncluster = count;
#else
    (void)n; (void)list; (void)count;
#endif
}

const std::complex<double>* SmallClusterBatch::get(int ic, int& dim, int& per){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl || impl->off) return nullptr;
    Impl& m = *impl;
    int s = -1;
    for (int k=0; k<2; ++k) if (m.window[k].submitted && ic >= m.window[k].first && ic < m.window[k].last) s = k;
    if (s < 0){
        // both windows are behind ic (or none was used in this order): the older slot
        s = !m.window[0].submitted ? 0 : !m.window[1].submitted ? 1 : m.window[0].first < m.window[1].first ? 0 : 1;
        if (m.window[s].submitted && !m.window[s].waited) m.wait(s);
        if (m.ready()) m.prepare(s,ic);
        else if (m.off) return nullptr;   // the card could not be opened (backend auto): from here on all on the CPU
        else m.cpuWindow(s,ic);           // the card is not up yet (the start of the run): these on the CPU meanwhile
    }
    Impl::Window& w = m.window[s];
    if (!w.waited){
        const int o = 1-s;
        if (w.last <= m.ncluster && !(m.window[o].submitted && m.window[o].first == w.last) && m.ready()){
            if (m.window[o].submitted && !m.window[o].waited) m.wait(o);
            m.prepare(o,w.last);                                    // the next window first
        }
        m.wait(s);
    }
    const int where = w.where[ic-w.first];
    if (where < 0) return nullptr;
    dim = w.dim[ic-w.first];
    const int kind = w.kind[ic-w.first];
    const long packedCount = w.count[0]+w.count[1];
    const long index = kind==2 ? where : kind==4 ? w.count[0]+where : packedCount+where;
    if (m.status(m.own,s)[index] != smallcluster::statusOk){ ++m.fallback; --m.gpuClusters; if (kind==1) --m.warpClusters; return nullptr; }
    const size_t nstep = size_t(m.plan.nstep);
    per = kind==1 ? warpcluster::outputPer(m.wplan) : 1;
    const size_t offset = kind==1 ? size_t(packedCount)*nstep + size_t(where)*nstep*per : size_t(index)*nstep;
    return reinterpret_cast<const std::complex<double>*>(m.out(m.own,s)+offset);
#else
    (void)ic; (void)dim; (void)per;
    return nullptr;
#endif
}

void SmallClusterBatch::report(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl) return;
    Impl& m = *impl;
    const double releaseStart = MPI_Wtime();
    MPI_Barrier(m.mem.node);                   // every rank has its last window back
    double feeder[4] = {0,0,0,0}; long long feederCounts[2] = {0,0}; int flags[3] = {0,0,0};
    char unavailable[sizeof(FeederHeader::message)] = {};
    if (m.mem.nodeRank==0){
        const double start = MPI_Wtime();
        m.stop.store(true,std::memory_order_release);
        m.feeder.join();
        // the feeder's last word, whatever the ranks saw while they computed (they may have finished before it failed):
        // a failure stops the run here too, so that it does not depend on the timing (review 20261012-233058-452278 P3-01)
        const FeederHeader* fh = m.feederHeader();
        const int state = fh->state.load(std::memory_order_acquire);
        if (state == feederFailed){
            fprintf(stderr,"Error: the packed GPU kernel (feeder on node rank 0): %s\n",fh->message);
            ccexFatal();
        }
        if (state == feederUnavailable){ flags[2] = 1; snprintf(unavailable,sizeof unavailable,"%s",fh->message); }
        feeder[0] = m.feederOpen; feeder[1] = m.feederLock; feeder[2] = m.feederWaited; feeder[3] = MPI_Wtime()-start;
        feederCounts[0] = m.feederBatches; feederCounts[1] = m.feederItems;
        flags[0] = m.locked; flags[1] = m.early;
    }
    const double freeStart = MPI_Wtime();
    m.mem.release();
    const double released[3] = {MPI_Wtime()-releaseStart, MPI_Wtime()-freeStart, m.allocateSeconds};
    long long counts[4] = {m.gpuClusters, m.fallback, m.cpuStarting, m.warpClusters}, total[4] = {0,0,0,0}, batches[2] = {0,0};
    double seconds[2] = {m.prepareSeconds, m.waitSeconds}, most[2] = {0,0}, sum[2] = {0,0}, slowest[4] = {0,0,0,0};
    int anyFlags[3] = {0,0,0};
    MPI_Reduce(counts,total,4,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(feederCounts,batches,2,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(seconds,most,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(seconds,sum,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(feeder,slowest,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(flags,anyFlags,3,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    double releasedMost[3] = {0,0,0};
    MPI_Reduce(released,releasedMost,3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (rank==0 && anyFlags[2])
        printf("%s off: the GPU could not be opened (%s); backend=auto, so every cluster was computed on the CPU\n",
               m.packed ? "[gpu-packed]" : "[gpu-warp]", *unavailable ? unavailable : "on another node");
    if (rank==0) printf("%s clusters on the GPU %lld, on the CPU while the GPU started %lld, back to the CPU %lld; "
                        "per rank: prepare %.3f s (max %.3f), "
                        "wait for the GPU %.3f s (max %.3f); feeder: CUDA start %.3f s, page-lock %.3f s (%s)%s, first batch "
                        "after %.3f s, %lld batches of %lld windows, stop %.3f s; shared memory allocated in %.3f s, release %.3f s (of which freeing it %.3f s)\n",
                        m.packed ? "[gpu-packed]" : "[gpu-warp]", total[0], total[2], total[1], sum[0]/nprocess,
                        most[0], sum[1]/nprocess, most[1], slowest[0], slowest[1], anyFlags[0] ? "done" : "refused, pageable copies",
                        anyFlags[1] ? " early" : "", slowest[2], batches[0], batches[1], slowest[3], releasedMost[2], releasedMost[0], releasedMost[1]);
    if (rank==0 && m.warp) printf("[gpu-warp] of these in the warp kernel %lld\n", total[3]);
    delete impl;
    impl = nullptr;
#endif
}
