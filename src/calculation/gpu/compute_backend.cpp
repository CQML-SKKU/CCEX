/**
 * @file compute_backend.cpp
 * @brief The choice of compute backend (cpu, cuda, auto) per cluster and the GPU start-up.
 */
#include "../../../include/calculation/gpu/compute_backend.h"
#include "../../../include/calculation/gpu/small_cluster_batch.h"
#include "../../../include/calculation/small_eigensolver.h"
#include "../../../include/calculation/simulator.h"
#include "../../../include/clustering/cluster.h"
#include <functional>
#include <set>
#include <climits>
#include <cstdarg>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#ifdef CCEX_USE_CUDA
#include "../../../include/calculation/gpu/gpu_backend.h"
#endif

namespace {
bool cudaReady = false;
bool cudaAttempted = false;
bool cudaSupported = false;
int nodeRank = 0, nodeSize = 1;
bool nodeKnown = false;                  // A73: nodeRank/nodeSize found (collective, once per process on every rank)
int cpuFromDimension = INT_MAX;          // backend auto after an out-of-memory: this dimension and above on the CPU
long long outOfMemoryClusters = 0;
bool recovering = false;                 // X12: since an out-of-memory, released device buffers are freed (gpuOutOfMemoryRecovery)
bool cardKept = false;                   // A64: the previous configuration of configuration.k_range left the CUDA context open
struct Timing { long long count=0; double build=0, setup=0, evolution=0; };
std::map<std::tuple<std::string,std::string,int>,Timing> timings;
struct BatchTiming { long long batches=0, clusters=0; int largest=0; };
std::map<int,BatchTiming> clusterBatches;
// A98: gpu_options.ranks_per_gpu (computeGpuSplitSetup)
struct GpuSplit {
    bool on = false;
    std::vector<int> ranks;                 // the GPU ranks, ascending
    std::vector<char> isGpu;                // per rank
    Config* config = nullptr; Pulse* pls = nullptr;
    std::function<double(int,const int*)> dim;
    std::map<int,double> ratio;             // dimension -> GPU time / CPU time; 0: the GPU did not take it
} split;
// A98: on this rank the GPU may take the dimension: a GPU rank, and not measured slower there than on the CPU
bool splitGpuTakes(int dimension){
    if (!split.isGpu[rank]) return false;
    const auto it = split.ratio.find(dimension);
    return it == split.ratio.end() || (it->second > 0 && it->second < 1);
}
void append(std::string& text, const char* format, ...){
    char line[512];
    va_list args; va_start(args,format);
    vsnprintf(line,sizeof line,format,args);
    va_end(args);
    text += line;
}
void fail(const char* message){
    fprintf(stderr,"Error: GPU backend: %s\n",message);
    MPI_Abort(MPI_COMM_WORLD,1);
}
}

bool computeBackendMethodSupported(const Config* config){
    // pCCE changes the cluster partition, then uses the same conditional-bath
    // coherence calculation as conventional CCE.
    return !strcasecmp(config->method,"cce") ||
           !strcasecmp(config->method,"pcce") ||
           (!strcasecmp(config->method,"gcce") &&
            !strcasecmp(config->propagator,"eigen"));
}

void computeBackendInit(const Config* config, QubitArray* qa, Pulse* pls, const char* savemode){
    // A64: a card the previous configuration kept open (computeBackendFinish(true)) is used as it is - the same input, so
    // the same backend, device and cache limit; everything else starts again as in a run of its own
    const bool kept = cardKept && cudaReady;
    cardKept = false;
    if (!kept){ cudaReady=false; cudaAttempted=false; cudaSupported=false; }
    cpuFromDimension=INT_MAX; outOfMemoryClusters=0; recovering=false;
    smallClusterWarmUp(config,qa,pls,savemode);   // A26: the small-cluster batch's CUDA context, early (nothing unless asked for or ruled out, X7)
    if (!strcasecmp(config->backend,"cpu")) return;
    const bool required = !strcasecmp(config->backend,"cuda");
    if (!computeBackendMethodSupported(config)){
        if (required) fail("execution.backend=cuda supports method cCCE (cluster.mode standard or pcce), or method gCCE with execution.propagator=eigen");
        if (rank==0) fprintf(stdout,"[backend] auto: this method/propagator uses CPU\n");
        return;
    }
    cudaSupported=true;
#ifdef CCEX_USE_CUDA
    // A73: the node split is collective, so every rank does it here, the first time - not after `kept`, which differs
    // between ranks under backend auto (only a rank that met an eligible cluster opened the card): a rank that skipped
    // it left the others waiting in the second configuration of configuration.k_range
    if (!nodeKnown){
        MPI_Comm node;
        MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&nodeRank);
        MPI_Comm_size(node,&nodeSize);
        MPI_Comm_free(&node);
        nodeKnown = true;
    }
#endif
    if (kept) return;
#ifdef CCEX_USE_CUDA
    // Only metadata discovery is collective. Lazy initialization below must not
    // use collectives: ranks can encounter eligible clusters at different times,
    // or have no eligible clusters at all.
    if (!required){
        if (rank==0) fprintf(stdout,"[backend] auto: CUDA initialization deferred until an eligible cluster\n");
        return;
    }
    char description[512] = {};
    cudaAttempted=true;
    cudaReady = gpuInitialize(nodeRank,nodeSize,description,sizeof(description));
    if (cudaReady) gpuSetCacheLimit(size_t(config->gpu_cache_mb)*1024*1024);
    // A required GPU run must fail collectively if any rank cannot initialize.
    int ok = cudaReady ? 1 : 0, allOk=0;
    MPI_Allreduce(&ok,&allOk,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if (required && !allOk) fail(description[0] ? description : "CUDA unavailable on another rank");
    fprintf(stdout,"[backend] rank=%d %s\n",rank,description);
#else
    if (required) fail("this executable has no CUDA support, run bin/main_gpu.out (bash do_compile.sh gpu, or a CMake build with CUDA)");
    if (rank==0) fprintf(stdout,"[backend] auto: CPU-only build\n");
#endif
}

std::size_t computeDeviceMemory(const Config* config){
    if (!cudaSupported || !strcasecmp(config->backend,"cpu")) return 0;
#ifdef CCEX_USE_CUDA
    return gpuDeviceMemorySmallest();
#else
    return 0;
#endif
}

bool computeUseCuda(const Config* config, int dimension){
    if (!cudaSupported || !strcasecmp(config->backend,"cpu")) return false;
    if (split.on && !splitGpuTakes(dimension)) return false;   // A98: not a GPU rank, or measured slower on the GPU
#ifdef CCEX_USE_CUDA
    // X12: a new cluster - whatever the last out-of-memory left is destroyed by now, so the cache works again
    if (recovering){ gpuOutOfMemoryRecovery(false); recovering=false; }
#endif
    if (!strcasecmp(config->backend,"auto") && (dimension<config->gpu_min_dim || dimension>=cpuFromDimension)) return false;
#ifdef CCEX_USE_CUDA
    if (!cudaAttempted){
        char description[512]={};
        cudaAttempted=true;
        cudaReady=gpuInitialize(nodeRank,nodeSize,description,sizeof(description));
        if (cudaReady) gpuSetCacheLimit(size_t(config->gpu_cache_mb)*1024*1024);
        fprintf(stdout,"[backend] rank=%d %s\n",rank,description);
    }
#endif
    return cudaReady;
}

bool computeMayUseCuda(const Config* config, int dimension){
#ifdef CCEX_USE_CUDA
    if (!cudaSupported || !strcasecmp(config->backend,"cpu") || (cudaAttempted && !cudaReady)) return false;
    if (split.on && !splitGpuTakes(dimension)) return false;   // A98: not a GPU rank, or measured slower on the GPU
    return strcasecmp(config->backend,"auto") || (dimension>=config->gpu_min_dim && dimension<cpuFromDimension);
#else
    (void)config; (void)dimension;
    return false;                        // a CPU-only build (computeUseCuda: cudaReady stays false)
#endif
}

bool computeCudaOutOfMemory(const Config* config, int dimension, const char* what){
    if (strcasecmp(config->backend,"auto")) return false;
    ++outOfMemoryClusters;
#ifdef CCEX_USE_CUDA
    // X12 (review 20261013b-gpu-perf I1): the cached buffers, and those the caller's partial work and batch objects
    // release after this (until the next cluster, computeUseCuda), are freed instead of kept for reuse
    gpuOutOfMemoryRecovery(true); recovering=true;
#endif
    if (dimension < cpuFromDimension){
        cpuFromDimension = dimension;
        fprintf(stdout,"[gpu-memory-fallback] rank=%d dim=%d: the GPU is out of memory (%s), backend=auto, so this cluster "
                       "and every later one of dimension %d or more on this rank are computed on the CPU\n",rank,dimension,what,dimension);
        fflush(stdout);
    }
    return true;
}

void computeProfile(const Config* config, const char* method, const char* backend,
                    int dimension, double build, double setup, double evolution){
    if (!config->profile) return;
    if (!strcasecmp(method,"cce") && !strcasecmp(config->method,"pcce")) method="pcce";
    Timing& t = timings[std::make_tuple(method,backend,dimension)];
    ++t.count; t.build += build; t.setup += setup; t.evolution += evolution;
}

// Every rank's records are gathered to rank 0 and printed there, in rank order. Each rank
// used to print its own: on a shared stdout the records of two ranks could interleave
// (check_gpu/README.md), and in the in-process adaptive path only rank 0's stdout is
// redirected to the step's run.log, so the other ranks' records never reached it. Every
// rank calls this at the same points (computeBackendFinish, and after each in-process step).
void computeProfileFlush(){
    const bool profiled=!timings.empty();
    std::string text;
    for (const auto& entry : timings){
        const Timing& t=entry.second;
        append(text,"[profile] rank=%d method=%s backend=%s dim=%d clusters=%lld "
                    "build_s=%.9f setup_s=%.9f evolution_s=%.9f\n",rank,
               std::get<0>(entry.first).c_str(),std::get<1>(entry.first).c_str(),
               std::get<2>(entry.first),t.count,t.build,t.setup,t.evolution);
    }
    timings.clear();
    {   // A42: the small cCCE eigensystems by Jacobi and the matrices that went to Eigen because it had not converged
        long long solved=0, fallbacks=0;
        smalleig::takeCounts(solved,fallbacks);
        if (profiled && solved+fallbacks)
            append(text,"[small-eig] rank=%d jacobi=%lld eigen_fallback=%lld\n",rank,solved,fallbacks);
    }
    for(const auto& entry:clusterBatches)
        append(text,"[gpu-batch] rank=%d dim=%d batches=%lld clusters=%lld max_clusters=%d\n",
               rank,entry.first,entry.second.batches,entry.second.clusters,entry.second.largest);
    clusterBatches.clear();
    if(outOfMemoryClusters)
        append(text,"[gpu-memory-fallback] rank=%d clusters=%lld from_dim=%d\n",rank,outOfMemoryClusters,cpuFromDimension);
    outOfMemoryClusters=0;
#ifdef CCEX_USE_CUDA
    if (cudaReady && profiled){
        const auto stats=gpuMemoryStats();
        append(text,"[gpu-memory] rank=%d allocations=%zu reuses=%zu cached_bytes=%zu peak_cached_bytes=%zu pointer_uploads=%zu\n",
               rank,stats.allocations,stats.reuses,stats.cachedBytes,stats.peakCachedBytes,stats.pointerUploads);
    }
#endif
    int size=1;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    int length=int(text.size());
    std::vector<int> lengths(rank==0 ? size : 0), offsets(rank==0 ? size : 0);
    MPI_Gather(&length,1,MPI_INT,lengths.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    std::string all;
    if (rank==0){
        int total=0;
        for (int r=0;r<size;++r){ offsets[r]=total; total+=lengths[r]; }
        all.resize(size_t(total));
    }
    MPI_Gatherv(text.data(),length,MPI_CHAR,rank==0 ? &all[0] : nullptr,lengths.data(),offsets.data(),
                MPI_CHAR,0,MPI_COMM_WORLD);
    if (rank==0 && !all.empty()){
        fflush(stdout);                 // a record starts a line: nothing of this rank is left half-written
        fputs(all.c_str(),stdout);
        fflush(stdout);
    }
}

void computeBackendFinish(bool keepCard){
    computeProfileFlush();
    smallClusterWarmDown();
    if (keepCard && cudaReady){       // A64: the next configuration of configuration.k_range uses this context
#ifdef CCEX_USE_CUDA
        if (recovering) gpuOutOfMemoryRecovery(false);   // as computeUseCuda does at the next cluster
#endif
        cardKept = true;
        return;
    }
    cardKept = false;
#ifdef CCEX_USE_CUDA
    gpuFinalize();
#endif
    cudaReady=false;
    cudaAttempted=false;
    cudaSupported=false;
}

void computeClusterBatch(const Config* config,int dimension,int count){
    if(!config->profile)return;
    auto& timing=clusterBatches[dimension];
    ++timing.batches; timing.clusters+=count;
    timing.largest=std::max(timing.largest,count);
}

void computeGpuSplitSetup(Config* config, QubitArray* qa, BathArray* ba, DefectArray* dfa, Pulse* pls){
    split = GpuSplit{};
    if (config->gpu_ranks_per_gpu <= 0 || !cudaSupported) return;
#ifdef CCEX_USE_CUDA
    // node rank 0 alone asks the driver, so the other ranks never start CUDA
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    int devices = nodeRank==0 ? gpuVisibleDevices() : 0;
    MPI_Bcast(&devices,1,MPI_INT,0,node);
    MPI_Comm_free(&node);
    int mine = nodeRank < (long long)config->gpu_ranks_per_gpu*devices ? 1 : 0;
    std::vector<int> all(nprocess,0);
    MPI_Allgather(&mine,1,MPI_INT,all.data(),1,MPI_INT,MPI_COMM_WORLD);
    split.isGpu.assign(all.begin(),all.end());
    for (int r=0; r<nprocess; r++) if (all[r]) split.ranks.push_back(r);
    if (rank==0) printf("[gpu-split] ranks_per_gpu=%d: %zu GPU rank(s) of %d\n",config->gpu_ranks_per_gpu,split.ranks.size(),nprocess);
    if (split.ranks.empty()) return;   // no GPU visible: every cluster on the CPU, split as without the key
    split.on = true;
    split.config = config; split.pls = pls;
    split.dim = Cluster_dimFunction(config,qa,ba,dfa);
#else
    (void)qa; (void)ba; (void)dfa; (void)pls;
#endif
}

bool computeGpuSplitOn(){ return split.on; }
const std::vector<int>& computeGpuSplitRanks(){ return split.ranks; }

int computeGpuSplitPending(int n, const int* row){
    const double d = split.dim(n,row);
    return d >= split.config->gpu_min_dim && d <= INT_MAX && !split.ratio.count(int(d)) ? int(d) : 0;
}

double computeGpuSplitWeight(int n, const int* row){
    const double d = split.dim(n,row);
    if (d < split.config->gpu_min_dim || d > INT_MAX) return 0;
    const auto it = split.ratio.find(int(d));
    return it == split.ratio.end() || it->second >= 1 ? 0 : it->second;
}

namespace {
// the seconds of a made-up cluster of this dimension on the CPU or the GPU (cceTimeMadeUp; the fastest of `repeats`: the
// first GPU one starts the card); < 0 when the GPU does not take it. The GPU under backend auto from this dimension on,
// so that running out of its memory leaves the cluster to the CPU (computeCudaOutOfMemory), as in the run
double timeDimension(int dimension, bool gpu, int repeats){
    Config* c = split.config;
    char saved[sizeof c->backend];
    strcpy(saved,c->backend);
    const int minDim = c->gpu_min_dim;
    const bool profile = c->profile;
    c->profile = false;                     // not in the [profile] totals of the run
    strcpy(c->backend,gpu ? "auto" : "cpu");
    c->gpu_min_dim = dimension;
    double best = -1;
    if (!gpu || computeUseCuda(c,dimension)){
        for (int i=0; i<repeats; i++){
            const double dt = cceTimeMadeUp(c,split.pls,dimension,gpu);
            if (dt < 0){ best = -1; break; }
            if (best < 0 || dt < best) best = dt;
        }
    }
    strcpy(c->backend,saved);
    c->gpu_min_dim = minDim;
    c->profile = profile;
    return best;
}
}

void computeGpuSplitCalibrate(const std::set<int>& dimensions){
    if (dimensions.empty()) return;
    const std::vector<int> dims(dimensions.begin(),dimensions.end());
    const int nd = int(dims.size());
    // t[2j], t[2j+1]: CPU and GPU seconds summed over the ranks that timed dimension j; count[2j], count[2j+1]: those ranks
    std::vector<double> t(2*nd,0.0), count(2*nd,0.0);
    // the CPU under the load of the run: every rank at once, rank r on dimension r mod nd (fewer ranks than dimensions:
    // on those j with j mod nprocess = r) - the ranks share the memory bandwidth; on a GH200 node a D = 900 cluster took
    // 1.9 s alone and 5.8 s with 288 ranks
    for (int j=0; j<nd; j++)
        if (nprocess >= nd ? rank % nd == j : j % nprocess == rank){ t[2*j] += timeDimension(dims[j],false,1); count[2*j] += 1; }
    // the GPU as loaded as in the run: one dimension at a time, every GPU rank on it at once from a barrier, so that each
    // card has its ranks_per_gpu ranks on it for the whole timing - taking turns without MPS (each waits about k times as
    // long), running together with it. The time is then what such a cluster holds its rank in the run, with or without
    // MPS. (Dimensions dealt round the GPU ranks finished at different times and were timed nearly alone: Hangang, 1 atm
    // case, D = 480 at 14 ms against 13.6 ms alone with 4 ranks per card.) Twice: the first starts the card
    const bool gpuRank = split.isGpu[rank];
    for (int j=0; j<nd; j++){
        double dt = 0;
        for (int rep=0; rep<2; rep++){
            MPI_Barrier(MPI_COMM_WORLD);
            if (gpuRank) dt = timeDimension(dims[j],true,1);
        }
        if (gpuRank && dt > 0){ t[2*j+1] += dt; count[2*j+1] += 1; }
    }
    MPI_Allreduce(MPI_IN_PLACE,count.data(),2*nd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,t.data(),int(t.size()),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    size_t j = 0;
    for (int d : dimensions){
        const double cpu = t[2*j]/count[2*j], gpu = count[2*j+1] > 0 ? t[2*j+1]/count[2*j+1] : 0;
        split.ratio[d] = cpu > 0 && gpu > 0 ? gpu/cpu : 0;
        if (rank==0) printf("[gpu-split] dim=%d cpu=%.4g s gpu=%.4g s gpu/cpu=%.4g%s\n",d,cpu,gpu,split.ratio[d],
                            gpu <= 0 ? " (the GPU did not take it: every rank, on the CPU)" :
                            split.ratio[d] >= 1 ? " (slower on the GPU: every rank, on the CPU)" : "");
        j++;
    }
    fflush(stdout);
}
