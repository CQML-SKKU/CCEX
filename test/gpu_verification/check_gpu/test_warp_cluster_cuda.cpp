// A40: the warp kernel (warp_cluster_gpu.cu, gpu_options.warp_max_dim) against clusters the CPU computed.
//   test_warp_cluster_cuda.out data/warp_cluster/*.bin
// Each file holds clusters dumped by the CPU path (CCEX_WARP_CLUSTER_DUMP, small_cluster_batch.h): the plan (time plan;
// gCCE pulses and qubit states), the input (Hamiltonian(s) and, for a single state, the start) and the CPU result (the
// coherence, or the qdim x qdim reduced density matrix for gCCE quantity dm). The card computes each file through the
// batch API the run uses (gpuSmallSubmit: one batch of two items, in lane 0 and again in lane 1); every value must agree
// with the CPU within the GPU tolerance 1e-8 + 2e-6 |CPU|, every cluster must converge, and the two lanes must give the same
// bits. A wrong device compilation - the sm_120 NVVM optimizer gave wrong Jacobi eigenvalues in the A26 and A40 step-0
// prototypes - fails here; the per-dimension lines say which size went wrong.
//   test_warp_cluster_cuda.out --bench N dump.bin...   also times, per file and dimension, one batch of N clusters (the
//                                                      file's clusters of that dimension repeated; copies included)
#include "calculation/gpu/small_cluster_gpu.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <vector>

using smallcluster::C;
using warpcluster::Plan;

struct Record { Plan plan; int dim; std::vector<C> in, cpu; };

// One record of a dump. readEnd only at a clean end of the file, before the first byte of a header; a header or body cut
// short, a read error or a header that makes no sense is readBroken (review 20261013-180049 G2, as test_small_cluster_cuda:
// a cut file used to pass with the records before the cut).
enum ReadResult { readOk, readEnd, readBroken };
static bool readBytes(FILE* f, void* p, size_t bytes){ return fread(p,1,bytes,f) == bytes; }
static ReadResult readRecord(FILE* f, Record& r){
    int head[8];
    const size_t got = fread(head,1,sizeof head,f);
    if (got == 0 && feof(f) && !ferror(f)) return readEnd;
    if (got != sizeof head) return readBroken;
    r = Record{};
    Plan& p = r.plan;
    p.mode = head[0]; r.dim = head[1]; p.qdim = head[2]; p.ensemble = head[3]; p.dm = head[4];
    p.time.intervals = head[5]; p.time.distinct = head[6]; p.time.nstep = head[7];
    if ((p.mode != warpcluster::modeCce && p.mode != warpcluster::modeGcce) || r.dim < 2 || r.dim > warpcluster::maxDim ||
        p.qdim < 1 || p.qdim > warpcluster::maxQubitDim || p.time.intervals < 1 || p.time.intervals > warpcluster::maxIntervals ||
        p.time.distinct < 1 || p.time.distinct > 4 || p.time.nstep < 1) return readBroken;
    if (!readBytes(f,&p.time.deltat,sizeof(double)) || !readBytes(f,p.time.fraction,sizeof(double)*p.time.distinct) ||
        !readBytes(f,p.time.same,sizeof(int)*p.time.intervals)) return readBroken;
    if (p.mode == warpcluster::modeGcce){
        const size_t q = size_t(p.qdim);
        for (int k=0; k<p.time.intervals-1; ++k) if (!readBytes(f,p.pulse[k],sizeof(C)*q*q)) return readBroken;
        if (!readBytes(f,p.psia,sizeof(C)*q) || !readBytes(f,p.psib,sizeof(C)*q) || !readBytes(f,p.qpsi0,sizeof(C)*q)) return readBroken;
    }
    r.in.resize(warpcluster::inputSize(p,r.dim));
    r.cpu.resize(size_t(p.time.nstep)*warpcluster::outputPer(p));
    if (!readBytes(f,r.in.data(),sizeof(C)*r.in.size()) || !readBytes(f,r.cpu.data(),sizeof(C)*r.cpu.size())) return readBroken;
    return readOk;
}
// every record of the file; false if it is broken anywhere
static bool readAll(FILE* f, std::vector<Record>& recs){
    recs.clear();
    Record r;
    while (true){
        const ReadResult res = readRecord(f,r);
        if (res == readEnd) return !ferror(f);
        if (res == readBroken) return false;
        recs.push_back(r);
    }
}
static size_t recordBytes(const Record& r){
    const Plan& p = r.plan;
    size_t b = 8*sizeof(int) + sizeof(double)*(1+p.time.distinct) + sizeof(int)*p.time.intervals;
    if (p.mode == warpcluster::modeGcce) b += sizeof(C)*(size_t(p.qdim)*p.qdim*(p.time.intervals-1) + 3*size_t(p.qdim));
    return b + sizeof(C)*(r.in.size()+r.cpu.size());
}
// The reader on cut copies of a dump (in memory, no GPU): the whole file and the file cut right after its first record
// must read; cut inside the second header, inside its plan, inside its body or one byte short must not.
static bool readerChecks(const char* path){
    FILE* f = fopen(path,"rb");
    if (!f) return false;
    std::vector<char> bytes;
    char buffer[65536]; size_t n;
    while ((n = fread(buffer,1,sizeof buffer,f)) > 0) bytes.insert(bytes.end(),buffer,buffer+n);
    fclose(f);
    std::vector<Record> recs;
    auto readPrefix = [&](size_t size, std::vector<Record>& out){
        FILE* m = fmemopen(bytes.data(),size,"rb");
        if (!m) return false;
        const bool ok = readAll(m,out);
        fclose(m);
        return ok;
    };
    if (!readPrefix(bytes.size(),recs) || recs.empty()){ printf("FAIL warp-cluster reader %s: the whole file does not read\n",path); return false; }
    bool ok = true;
    if (recs.size() >= 2){
        const size_t first = recordBytes(recs[0]);
        struct Cut { size_t size; bool reads; const char* what; } cuts[] = {
            {first, true, "after the first record"}, {first+3, false, "3 bytes into the second header"},
            {first+20, false, "20 bytes into the second header"}, {first+32+4, false, "inside the second record's plan"},
            {first+recordBytes(recs[1])/2, false, "inside the second body"}, {bytes.size()-1, false, "one byte short"}};
        for (const Cut& c : cuts){
            std::vector<Record> got;
            const bool reads = readPrefix(c.size,got);
            const bool right = reads == c.reads && (!reads || got.size() == 1);
            if (!right) printf("FAIL warp-cluster reader %s cut %s (%zu bytes): %s\n",path,c.what,c.size,reads ? "read" : "refused");
            ok &= right;
        }
    } else {
        std::vector<Record> got;
        if (readPrefix(bytes.size()-1,got)){ printf("FAIL warp-cluster reader %s: one byte short was read\n",path); ok = false; }
    }
    std::vector<char> bad(bytes);
    const int three = 3;
    std::memcpy(bad.data(),&three,sizeof three);           // mode 3 in the first header
    FILE* m = fmemopen(bad.data(),bad.size(),"rb");
    std::vector<Record> got;
    const bool reads = m && readAll(m,got);
    if (m) fclose(m);
    if (reads){ printf("FAIL warp-cluster reader %s: a header of mode 3 was read\n",path); ok = false; }
    printf("%s warp-cluster reader %s: whole file, clean cut and broken cuts\n",ok ? "PASS" : "FAIL",path);
    return ok;
}

static bool samePlan(const Plan& a, const Plan& b){
    if (a.mode != b.mode || a.qdim != b.qdim || a.ensemble != b.ensemble || a.dm != b.dm || a.time.intervals != b.time.intervals ||
        a.time.distinct != b.time.distinct || a.time.nstep != b.time.nstep || a.time.deltat != b.time.deltat ||
        std::memcmp(a.time.fraction,b.time.fraction,sizeof(double)*a.time.distinct) ||
        std::memcmp(a.time.same,b.time.same,sizeof(int)*a.time.intervals)) return false;
    if (a.mode != warpcluster::modeGcce) return true;
    const size_t q = size_t(a.qdim);
    for (int k=0; k<a.time.intervals-1; ++k) if (std::memcmp(a.pulse[k],b.pulse[k],sizeof(C)*q*q)) return false;
    return !std::memcmp(a.psia,b.psia,sizeof(C)*q) && !std::memcmp(a.psib,b.psib,sizeof(C)*q) && !std::memcmp(a.qpsi0,b.qpsi0,sizeof(C)*q);
}

struct Score { double worst = 0; long outside = 0, points = 0, notOk = 0, clusters = 0; };
static void compare(const Record& r, const C* got, int status, Score& s){
    ++s.clusters;
    if (status != smallcluster::statusOk) ++s.notOk;
    for (size_t t=0; t<r.cpu.size(); ++t){
        const double diff = std::hypot(got[t].re-r.cpu[t].re, got[t].im-r.cpu[t].im), ref = std::hypot(r.cpu[t].re, r.cpu[t].im);
        if (!(diff <= 1e-8 + 2e-6*ref)) ++s.outside;     // a NaN counts as outside
        if (!(diff <= s.worst)) s.worst = std::isnan(diff) ? INFINITY : diff;
        ++s.points;
    }
}

int main(int argc, char** argv){
    long bench = 0;
    if (argc > 2 && !std::strcmp(argv[1],"--bench")){ bench = atol(argv[2]); argv += 2; argc -= 2; }
    if (argc < 2){ fprintf(stderr,"usage: %s [--bench N] dump.bin...\n",argv[0]); return 2; }
    char description[256] = {};
    GpuSmallDevice* gpu = gpuSmallOpen(0,description,sizeof description);
    if (!gpu){ printf("FAIL warp-cluster: no GPU (%s)\n",description); return 1; }
    printf("[warp-cluster] %s\n",description);
    bool failed = false;
    for (int a=1; a<argc; ++a){
        FILE* f = fopen(argv[a],"rb");
        if (!f){ printf("FAIL warp-cluster %s: cannot open\n",argv[a]); failed = true; continue; }
        std::vector<Record> recs;
        const bool whole = readAll(f,recs);
        fclose(f);
        if (recs.empty() || !whole){ printf("FAIL warp-cluster %s: unreadable or cut short (%zu records before it)\n",argv[a],recs.size()); failed = true; continue; }
        if (!readerChecks(argv[a])){ failed = true; continue; }
        const Plan& plan = recs[0].plan;
        bool fine = true;
        for (const Record& x : recs) fine = fine && samePlan(x.plan,plan) && warpcluster::compiled(plan.mode,plan.qdim,x.dim);
        if (!fine){ printf("FAIL warp-cluster %s: records with different plans or uncompiled dimensions\n",argv[a]); failed = true; continue; }
        const size_t per = size_t(plan.time.nstep)*warpcluster::outputPer(plan);
        std::map<int,Score> scores;
        bool sameBits = true;
        try {
            gpuSmallWarpPlan(gpu,plan);
            // the file as one batch of two items (the first half and the rest), in lane 0 and again in lane 1
            std::vector<C> in; std::vector<unsigned char> dims;
            for (const Record& x : recs){ in.insert(in.end(),x.in.begin(),x.in.end()); dims.push_back((unsigned char)x.dim); }
            const long n = long(recs.size()), half = n/2;
            size_t firstIn = 0;
            for (long c=0; c<half; ++c) firstIn += recs[c].in.size();
            std::vector<C> previous;
            for (int lane=0; lane<2; ++lane){
                std::vector<C> out(size_t(n)*per);
                std::vector<unsigned char> status(n, 0xff);
                GpuSmallItem items[2] = {
                    GpuSmallItem{nullptr,nullptr,0,0,out.data(),status.data(),in.data(),dims.data(),half},
                    GpuSmallItem{nullptr,nullptr,0,0,out.data()+size_t(half)*per,status.data()+half,in.data()+firstIn,dims.data()+half,n-half}};
                gpuSmallSubmit(gpu,lane,plan.time,plan.ensemble,false,items,2);
                gpuSmallFinish(gpu,lane);
                if (lane == 0) for (long c=0; c<n; ++c) compare(recs[c],out.data()+size_t(c)*per,status[c],scores[recs[c].dim]);
                if (lane == 1) sameBits = out.size() == previous.size() && !std::memcmp(out.data(),previous.data(),out.size()*sizeof(C));
                previous.swap(out);
            }
        } catch (const std::exception& e){
            printf("FAIL warp-cluster %s: %s\n",argv[a],e.what()); failed = true; continue;
        }
        bool fileOk = sameBits;
        for (const auto& kv : scores){
            const Score& s = kv.second;
            const bool ok = s.outside == 0 && s.notOk == 0;
            fileOk = fileOk && ok;
            printf("%s warp-cluster %s: %s %s%s, %d intervals (%d distinct), dimension %d: %ld clusters, %ld values, max |GPU-CPU| %.3e, "
                   "outside %ld, not converged %ld\n", ok ? "PASS" : "FAIL", argv[a], plan.mode==warpcluster::modeGcce ? "gCCE" : "cCCE",
                   plan.ensemble ? "ensemble" : "single", plan.dm ? " dm" : "", plan.time.intervals, plan.time.distinct, kv.first,
                   s.clusters, s.points, s.worst, s.outside, s.notOk);
        }
        printf("%s warp-cluster %s: two lanes %s\n", sameBits ? "PASS" : "FAIL", argv[a], sameBits ? "same bits" : "DIFFER");
        failed |= !fileOk;
        for (const auto& kv : scores){
            if (bench <= 0) break;
            std::vector<C> in; std::vector<unsigned char> dims;
            std::vector<const Record*> mine;
            for (const Record& x : recs) if (x.dim == kv.first) mine.push_back(&x);
            for (long c=0; c<bench; ++c){ const Record& x = *mine[c % mine.size()]; in.insert(in.end(),x.in.begin(),x.in.end()); dims.push_back((unsigned char)x.dim); }
            std::vector<C> out(size_t(bench)*per); std::vector<unsigned char> status(bench);
            GpuSmallItem item{nullptr,nullptr,0,0,out.data(),status.data(),in.data(),dims.data(),bench};
            double best = 1e30;
            for (int k=0; k<4; ++k){                  // the first one warms up
                const auto t0 = std::chrono::steady_clock::now();
                gpuSmallSubmit(gpu,0,plan.time,plan.ensemble,false,&item,1);
                gpuSmallFinish(gpu,0);
                const double t = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
                if (k) best = std::min(best,t);
            }
            printf("BENCH warp-cluster %s: dimension %d, %ld clusters, %d points: %.4f s, %.0f clusters/s, %.2f us per cluster\n",
                   argv[a], kv.first, bench, plan.time.nstep, best, bench/best, best/bench*1e6);
        }
    }
    gpuSmallClose(gpu);
    return failed ? 1 : 0;
}
