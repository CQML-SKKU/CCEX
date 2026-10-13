#include "../../include/clustering/pcce_partition_cache.h"

#include <cerrno>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {
const std::string magic = "CCEX-PCCE-CACHE-v1\n";
// Bump when trial RNG, assignment, centroid arithmetic, scoring, or tie-breaking
// changes in cluster_pcce.cpp/cluster_pcce.h. CPU/CUDA share this algorithm.
const char* algorithm = "pcce-deterministic-20260912-v3";

[[noreturn]] void abortCache(const std::string& message){
    fprintf(stderr,"Error: pCCE partition cache: %s\n",message.c_str());
    MPI_Abort(MPI_COMM_WORLD,1);
    std::abort();
}

void put64(std::string& out, uint64_t value){
    for (int i=0;i<8;++i) out.push_back(char((value>>(8*i))&255));
}
uint64_t get64(const std::string& in, size_t& offset){
    if (offset>in.size() || in.size()-offset<8) throw std::runtime_error("truncated cache");
    uint64_t value=0;
    for (int i=0;i<8;++i) value|=uint64_t(uint8_t(in[offset++]))<<(8*i);
    return value;
}
void putText(std::string& out,const std::string& text){
    put64(out,text.size()); out+=text;
}
void putDouble(std::string& out,double value){
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
                  "partition cache requires IEEE binary64");
    if (!std::isfinite(value)) throw std::runtime_error("nonfinite partition coordinate or spin");
    uint64_t bits; std::memcpy(&bits,&value,8); put64(out,bits);
}
double getDouble(const std::string& in,size_t& offset){
    uint64_t bits=get64(in,offset); double value; std::memcpy(&value,&bits,8);
    if (!std::isfinite(value)) throw std::runtime_error("nonfinite cached center");
    return value;
}

// Accidental-corruption checksum, not authentication. Cache identity is compared
// byte-for-byte; a hash collision cannot substitute different bath/settings.
uint64_t checksum(const std::string& bytes,size_t length){
    uint64_t value=UINT64_C(14695981039346656037);
    for(size_t i=0;i<length;++i){value^=uint8_t(bytes[i]);value*=UINT64_C(1099511628211);}
    return value;
}

std::string prefix(const std::string& identity){
    std::string out=magic; putText(out,identity); return out;
}
size_t fileSize(const std::string& head,const Partition_info* p){
    const uint64_t count=uint64_t(p->ncenter)*3+uint64_t(p->pcce_nspin)+1;
    if(count>(std::numeric_limits<size_t>::max()-head.size())/8)
        throw std::runtime_error("cache size overflow");
    return head.size()+size_t(count)*8;
}

void validateBundles(const Partition_info* p,const Point* centers,const int* assigned){
    std::vector<int> counts(p->ncenter,0);
    for(int i=0;i<p->pcce_nspin;++i){
        if(assigned[i]<0 || assigned[i]>=p->ncenter)
            throw std::runtime_error("cached bundle index out of range");
        ++counts[assigned[i]];
    }
    for(int i=0;i<p->ncenter;++i){
        if(counts[i]!=p->sK) throw std::runtime_error("cached bundle size differs from sK");
        if(!std::isfinite(centers[i].x) || !std::isfinite(centers[i].y) || !std::isfinite(centers[i].z))
            throw std::runtime_error("nonfinite cached center");
    }
}

struct Descriptor {
    int fd;
    explicit Descriptor(int value):fd(value){}
    ~Descriptor(){if(fd>=0)close(fd);}
};

// Only rank zero accesses the shared file. Length is derived from CURRENT input,
// never from untrusted file counts, so truncated/oversized files cannot allocate
// arbitrary memory. A miss in readwrite mode recomputes the original algorithm.
bool readCache(const char* path,const std::string& head,const Partition_info* p,
               Point* centers,int* assigned,std::string& reason){
    Descriptor input(open(path,O_RDONLY|O_NONBLOCK));
    if(input.fd<0){
        if(errno==ENOENT){reason="missing";return false;}
        throw std::runtime_error(std::string("cannot read cache: ")+strerror(errno));
    }
    struct stat info{};
    if(fstat(input.fd,&info)) throw std::runtime_error("cannot stat cache");
    if(!S_ISREG(info.st_mode)) throw std::runtime_error("cache path must be a regular file");
    const size_t expected=fileSize(head,p);
    if(info.st_size<0 || uint64_t(info.st_size)!=expected){reason="length-or-input-mismatch";return false;}
    std::string bytes(expected,'\0'); size_t offset=0;
    while(offset<bytes.size()){
        ssize_t count=read(input.fd,&bytes[offset],bytes.size()-offset);
        if(count<0 && errno==EINTR)continue;
        if(count<0)throw std::runtime_error(std::string("cannot read cache: ")+strerror(errno));
        if(!count){reason="truncated";return false;}
        offset+=size_t(count);
    }
    if(bytes.compare(0,head.size(),head)!=0){reason="input-or-version-mismatch";return false;}
    offset=bytes.size()-8;
    if(get64(bytes,offset)!=checksum(bytes,bytes.size()-8)){reason="checksum-mismatch";return false;}
    try {
        offset=head.size();
        for(int i=0;i<p->ncenter;++i){
            centers[i].x=getDouble(bytes,offset); centers[i].y=getDouble(bytes,offset); centers[i].z=getDouble(bytes,offset);
        }
        for(int i=0;i<p->pcce_nspin;++i){
            uint64_t bundle=get64(bytes,offset);
            if(bundle>=uint64_t(p->ncenter))throw std::runtime_error("cached bundle index out of range");
            assigned[i]=int(bundle);
        }
        validateBundles(p,centers,assigned);
    }catch(const std::exception& error){reason=error.what();return false;}
    reason="matching-input";return true;
}

void writeCache(const char* path,const std::string& head,const Partition_info* p,
                const Point* centers,const int* assigned){
    validateBundles(p,centers,assigned);
    std::string bytes=head; bytes.reserve(fileSize(head,p));
    for(int i=0;i<p->ncenter;++i){putDouble(bytes,centers[i].x);putDouble(bytes,centers[i].y);putDouble(bytes,centers[i].z);}
    for(int i=0;i<p->pcce_nspin;++i)put64(bytes,uint64_t(assigned[i]));
    put64(bytes,checksum(bytes,bytes.size()));
    // Unique temporary file in the same directory: readers see either complete
    // old data or complete new data, including concurrent MPI jobs on IDM/GPU.
    std::string pattern=std::string(path)+".tmp.XXXXXX";
    std::vector<char> temporary(pattern.begin(),pattern.end()); temporary.push_back('\0');
    Descriptor output(mkstemp(temporary.data()));
    if(output.fd<0)throw std::runtime_error(std::string("cannot create cache temporary file: ")+strerror(errno));
    try {
        size_t offset=0;
        while(offset<bytes.size()){
            ssize_t count=write(output.fd,bytes.data()+offset,bytes.size()-offset);
            if(count<0 && errno==EINTR)continue;
            if(count<=0)throw std::runtime_error(std::string("cannot write cache: ")+strerror(errno));
            offset+=size_t(count);
        }
        if(fsync(output.fd))throw std::runtime_error("cannot flush cache");
        const int descriptor=output.fd; output.fd=-1;
        if(close(descriptor))throw std::runtime_error("cannot close cache");
        if(rename(temporary.data(),path))throw std::runtime_error(std::string("cannot publish cache: ")+strerror(errno));
    }catch(...){unlink(temporary.data());throw;}
}
}

void Config_validatePcceCache(const Config* cnf){
    if(!strcasecmp(cnf->pcce_cache_mode,"off"))return;
    if(strcasecmp(cnf->method,"pcce"))abortCache("enabled cache requires method cCCE with cluster_mode pcce");
    if(strcasecmp(cnf->pcce_partition,"deterministic"))
        abortCache("enabled cache requires pcce.partition deterministic; legacy consumes the bath-state RNG");
    if(!cnf->pcce_cache_file[0])abortCache("enabled cache requires pcce.cache.file");
}

std::string pccePartitionCacheIdentity(BathArray* ba,const Partition_info* p,const Cluster* cls,const Config* config){
    if(!strcasecmp(config->pcce_cache_mode,"off"))return {};
    try {
        std::string key; putText(key,algorithm); putText(key,__VERSION__);
        // Conservative arithmetic compatibility. In particular, a native AVX/FMA
        // build must not silently reuse a partition from different host codegen.
        uint64_t arithmetic=0;
#ifdef __AVX2__
        arithmetic|=1;
#endif
#ifdef __AVX512F__
        arithmetic|=2;
#endif
#ifdef __FMA__
        arithmetic|=4;
#endif
#ifdef __FAST_MATH__
        arithmetic|=8;
#endif
        put64(key,arithmetic);put64(key,FLT_EVAL_METHOD);
#ifdef __GLIBCXX__
        put64(key,__GLIBCXX__);
#else
        put64(key,0);
#endif
        put64(key,uint32_t(config->seed));put64(key,p->sK);
        put64(key,cls->max_trial);put64(key,cls->max_iter);put64(key,cls->kmeans_pp);
        put64(key,ba->nspin);put64(key,p->pcce_nspin);put64(key,p->rest_nspin);
        for(int i=0;i<ba->nspin;++i){
            BathSpin* spin=ba->bath[i];
            putText(key,BathSpin_getName(spin));putDouble(key,BathSpin_getSpin(spin));
            const double* xyz=BathSpin_getXyz(spin);
            for(int axis=0;axis<3;++axis)putDouble(key,xyz[axis]);
        }
        return key;
    }catch(const std::exception& error){abortCache(error.what());}
}

static void searchOrRead(BathSpin** bath,Partition_info* p,Cluster* cls,Point* centers,
                         int* assigned,const Config* config,const std::string& identity){
    if(!strcasecmp(config->pcce_cache_mode,"off")){
        simulator_cluster_partition(bath,p,cls,centers,assigned,config);return;
    }
    const double start=MPI_Wtime();
    // All ranks agree on the full identity before sharing rank zero's file.
    if(identity.size()>INT_MAX || p->ncenter>INT_MAX/3)abortCache("MPI cache broadcast exceeds indexing limit");
    int size=rank==0 ? int(identity.size()) : 0;
    MPI_Bcast(&size,1,MPI_INT,0,MPI_COMM_WORLD);
    std::string shared=rank==0 ? identity : std::string(size,'\0');
    MPI_Bcast(&shared[0],size,MPI_CHAR,0,MPI_COMM_WORLD);
    int equal=shared==identity,allEqual=0;
    MPI_Allreduce(&equal,&allEqual,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(!allEqual)abortCache("partition inputs differ between MPI ranks");
    const std::string head=prefix(identity);
    int hit=0; std::string reason;
    if(rank==0){
        try {
            hit=readCache(config->pcce_cache_file,head,p,centers,assigned,reason) ? 1 : 0;
            if(!hit && !strcasecmp(config->pcce_cache_mode,"read"))
                throw std::runtime_error("read mode requires a valid matching cache ("+reason+")");
        }catch(const std::exception& error){hit=-1;fprintf(stderr,"Error: pCCE partition cache: %s\n",error.what());}
    }
    MPI_Bcast(&hit,1,MPI_INT,0,MPI_COMM_WORLD);
    if(hit<0)abortCache("lookup failed on rank 0");
    const double lookupEnd=MPI_Wtime();
    if(hit){
        static_assert(sizeof(Point)==3*sizeof(double),"Point layout must match MPI broadcast");
        MPI_Bcast(centers,p->ncenter*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(assigned,p->pcce_nspin,MPI_INT,0,MPI_COMM_WORLD);
        if(rank==0){
            printf("[partition] mode=deterministic source=cache total_trials=0\n");
            printf("[partition-cache] status=hit lookup_s=%.9f search_s=0 write_s=0 file=%s\n",MPI_Wtime()-start,config->pcce_cache_file);
        }
        return;
    }
    if(rank==0)printf("[partition-cache] status=miss reason=%s file=%s\n",reason.c_str(),config->pcce_cache_file);
    simulator_cluster_partition(bath,p,cls,centers,assigned,config);
    const double searchEnd=MPI_Wtime(); int saved=1;
    if(rank==0){
        try {writeCache(config->pcce_cache_file,head,p,centers,assigned);}
        catch(const std::exception& error){saved=0;fprintf(stderr,"Error: pCCE partition cache: %s\n",error.what());}
    }
    MPI_Bcast(&saved,1,MPI_INT,0,MPI_COMM_WORLD);
    if(!saved)abortCache("write failed on rank 0");
    if(rank==0)printf("[partition-cache] status=stored lookup_s=%.9f search_s=%.9f write_s=%.9f file=%s\n",
                     lookupEnd-start,searchEnd-lookupEnd,MPI_Wtime()-searchEnd,config->pcce_cache_file);
}

static PccePartition lastPartition;

void pccePartitionCached(BathSpin** bath,Partition_info* p,Cluster* cls,Point* centers,
                         int* assigned,const Config* config,const std::string& identity){
    searchOrRead(bath,p,cls,centers,assigned,config,identity);
    lastPartition.centers.assign(centers,centers+p->ncenter);
    lastPartition.assigned.assign(assigned,assigned+p->pcce_nspin);
}

const PccePartition& pcceLastPartition(){ return lastPartition; }
