#include "../../../include/calculation/gpu/compute_backend.h"
#include "../../../include/calculation/gpu/small_cluster_batch.h"
#include "../../../include/calculation/small_eigensolver.h"
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
int cpuFromDimension = INT_MAX;          // backend auto after an out-of-memory: this dimension and above on the CPU
long long outOfMemoryClusters = 0;
bool recovering = false;                 // X12: since an out-of-memory, released device buffers are freed (gpuOutOfMemoryRecovery)
struct Timing { long long count=0; double build=0, setup=0, evolution=0; };
std::map<std::tuple<std::string,std::string,int>,Timing> timings;
struct BatchTiming { long long batches=0, clusters=0; int largest=0; };
std::map<int,BatchTiming> clusterBatches;
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
    cudaReady=false; cudaAttempted=false; cudaSupported=false;
    cpuFromDimension=INT_MAX; outOfMemoryClusters=0; recovering=false;
    smallClusterWarmUp(config,qa,pls,savemode);   // A26: the small-cluster batch's CUDA context, early (nothing unless asked for or ruled out, X7)
    if (!strcasecmp(config->backend,"cpu")) return;
    const bool required = !strcasecmp(config->backend,"cuda");
    if (!computeBackendMethodSupported(config)){
        if (required) fail("backend=cuda supports method cCCE (cluster_mode standard or pcce), or method gCCE with propagator=eigen");
        if (rank==0) fprintf(stdout,"[backend] auto: this method/propagator uses CPU\n");
        return;
    }
    cudaSupported=true;
#ifdef CCEX_USE_CUDA
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    MPI_Comm_rank(node,&nodeRank);
    MPI_Comm_size(node,&nodeSize);
    MPI_Comm_free(&node);
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
    if (required) fail("this executable has no CUDA support; run bin/main_gpu.out (bash do_compile.sh gpu, or a CMake build with CUDA)");
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
        fprintf(stdout,"[gpu-memory-fallback] rank=%d dim=%d: the GPU is out of memory (%s); backend=auto, so this cluster "
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

void computeBackendFinish(){
    computeProfileFlush();
    smallClusterWarmDown();
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
