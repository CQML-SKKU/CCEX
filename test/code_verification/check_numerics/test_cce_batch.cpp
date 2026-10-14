// A30: the batched small-cluster cCCE evolution on the CPU (cce_batch.h) against the scalar path, bit for bit, on clusters
// CCEX dumped (CCEX_SMALL_CLUSTER_DUMP: Halpha, Hbeta, psi and the time plan of each cluster;
// test/gpu_verification/check_gpu/data/small_cluster/*.bin, the A26 fixtures and the production P1 CR60 sigma4 clusters).
// Both sides start from the dumped Hamiltonians and go through ccebatchBothPaths, the function the A51 self-check runs: the
// scalar path is cceSetupSmallAlone (A42; cceSetup when off) + cceEvolve (what calCoherenceCceInto runs after building
// them), the batched one CceBatchWindow's grouped eigensolves (CceEigenGroups) + ccebatchSetup + ccebatchPackLane + the
// kernel. The objects are main.out's, so the check holds for the build it is compiled in - run it in every build (make cpu-test).
// A dump that is empty, cut short or has a header that makes no sense fails (exit 2), as in test_small_eigensolver.
// Each file is run as dumped and with the other evolution_update; a one-state file also as the ensemble Hahn echo.
// Cases the build's kernel does not cover (lanes() 0, the ensemble CPMG/Ramsey) are counted as scalar and pass; a build
// whose kernel fails the A51 self-check (ccebatchUsable) passes as "off on this machine".
// usage: test_cce_batch.out DUMP...
#include "calculation/cce_batch_window.h"
#include "calculation/small_eigensolver.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// A46: with the MKL option (CCEX_BLAS_DIRECT) the scalar path calls the library only from dimension 21 (zgemm) and 63
// (zheevd), and those calls live in their own file (src/common/blas_option.cpp), so the scalar path of the dimensions
// compared here (2, 4) is the same code as without the option and the bit-for-bit comparison stands in that build too.
#ifdef CCEX_BLAS_DIRECT
static_assert(ccexBlas::gemmMin > 4 && ccexBlas::eigMin > 4, "test_cce_batch compares D 2 and 4 bit for bit");
#endif

int rank=0,nprocess=1,seed=19;
bool verbosity=false;

namespace {
struct Rec {
    int d=0, ens=0, rec=0, intervals=0, distinct=0, nstep=0;
    double deltat=0;
    std::vector<double> fraction; std::vector<int> same;
    MatrixXcd Ha, Hb, psi;
};
// false only at a clean end of the file, before the first byte of a header; anything else that is not a whole record stops
bool readRec(FILE* f, Rec& r){
    int h[6];
    const size_t got = fread(h,1,sizeof h,f);
    if (got == 0 && feof(f) && !ferror(f)) return false;
    if (got != sizeof h){ fprintf(stderr,"truncated dump\n"); exit(2); }
    r.d=h[0]; r.ens=h[1]; r.rec=h[2]; r.intervals=h[3]; r.distinct=h[4]; r.nstep=h[5];
    if ((r.d != 2 && r.d != 4) || (r.ens != 0 && r.ens != 1) || (r.rec != 0 && r.rec != 1) || r.intervals < 1 ||
        r.intervals > 65536 || r.distinct < 1 || r.distinct > 4 || r.nstep < 1 || r.nstep > (1<<26)){
        fprintf(stderr,"broken dump: a record header that makes no sense\n"); exit(2);
    }
    r.fraction.resize(r.distinct); r.same.resize(r.intervals);
    bool ok = fread(&r.deltat,sizeof(double),1,f)==1;
    ok = ok && fread(r.fraction.data(),sizeof(double),r.distinct,f)==size_t(r.distinct);
    ok = ok && fread(r.same.data(),sizeof(int),r.intervals,f)==size_t(r.intervals);
    r.Ha.resize(r.d,r.d); r.Hb.resize(r.d,r.d);
    ok = ok && fread(r.Ha.data(),sizeof(doublec),size_t(r.d)*r.d,f)==size_t(r.d)*r.d;
    ok = ok && fread(r.Hb.data(),sizeof(doublec),size_t(r.d)*r.d,f)==size_t(r.d)*r.d;
    if (!r.ens){ r.psi.resize(r.d,1); ok = ok && fread(r.psi.data(),sizeof(doublec),r.d,f)==size_t(r.d); }
    std::vector<doublec> cpu(r.nstep);
    ok = ok && fread(cpu.data(),sizeof(doublec),r.nstep,f)==size_t(r.nstep);
    if (!ok){ fprintf(stderr,"truncated dump\n"); exit(2); }
    return true;
}

long compared=0, scalarOnly=0, failed=0;

// One run of the clusters cl (all of one file) with the given plan: ensemble, evolution_update, intervals (from the
// file, or the Hahn echo)
void runCase(const std::string& name, const std::vector<Rec>& cl, bool ensemble, bool recursive, bool hahn){
    const Rec& r0 = cl[0];
    Config* cnf = Config_init();
    strcpy(cnf->method,"cce"); strcpy(cnf->quantity,"coherence");
    strcpy(cnf->evolution_update, recursive ? "recursive" : "exact");
    cnf->nstep = r0.nstep; cnf->deltat = r0.deltat; cnf->nstate = ensemble ? 0 : 1; cnf->profile = false;
    // the Pulse: sequence[p][2] and sequence_indices as Pulse builds them (assign_sequence_indices)
    const int intervals = hahn ? 2 : r0.intervals;
    std::vector<double> rows(3*intervals,0.0); std::vector<double*> seq(intervals); std::vector<int> idx(intervals);
    for (int p=0; p<intervals; ++p){
        rows[3*p+2] = hahn ? 0.5 : r0.fraction[r0.same[p]];
        seq[p] = &rows[3*p];
        idx[p] = p;
        for (int q=0; q<p; ++q) if (rows[3*q+2] == rows[3*p+2]){ idx[p] = idx[q]; break; }
    }
    Pulse pulse{};
    pulse.npulse = intervals-1; pulse.sequence = seq.data(); pulse.sequence_indices = idx.data();
    strcpy(pulse.pulsename, intervals == 1 ? "Ramsey" : intervals == 2 ? "HahnEcho" : "CPMG");

    ccebatch::Plan plan;
    const bool planned = ccebatchPlan(cnf,&pulse,plan);
    const int nstep = r0.nstep;
    for (int D : {2,4}){
        long differ=0, same=0, scalar=0;
        std::vector<const Rec*> group;
        for (const Rec& r : cl) if (r.d == D) group.push_back(&r);
        if (group.empty()) continue;
        if (!planned || !ccebatch::covers(plan,D)) scalar += long(group.size())*nstep;
        else {
            std::vector<const MatrixXcd*> Ha, Hb, psi;
            for (const Rec* r : group){ Ha.push_back(&r->Ha); Hb.push_back(&r->Hb); psi.push_back(&r->psi); }
            std::vector<doublec> scalarOut(group.size()*size_t(nstep)), batchOut(group.size()*size_t(nstep));
            std::vector<char> batched;
            ccebatchBothPaths(cnf,&pulse,plan,D,int(group.size()),Ha.data(),Hb.data(),psi.data(),scalarOut.data(),
                              batchOut.data(),batched);
            for (size_t c=0; c<group.size(); ++c){
                if (!batched[c]){ scalar += nstep; continue; }
                for (int i=0; i<nstep; ++i){
                    const doublec a = scalarOut[c*nstep+i], b = batchOut[c*nstep+i];
                    if (!memcmp(&a,&b,sizeof a)){ ++same; continue; }
                    if (differ < 3)
                        printf("  %s D %d cluster %zu step %d: scalar (%a, %a), batched (%a, %a)\n",name.c_str(),D,c,i,
                               a.real(),a.imag(),b.real(),b.imag());
                    ++differ;
                }
            }
        }
        printf("%-56s D %d %9ld same %7ld differ %7ld scalar only\n",name.c_str(),D,same,differ,scalar);
        compared += same+differ; failed += differ; scalarOnly += scalar;
    }
}
}

int main(int argc, char** argv){
    if (argc < 2){ fprintf(stderr,"usage: %s DUMP...\n",argv[0]); return 2; }
    printf("A30 batched small-cluster kernel: %s, %d lanes; A42 small eigensolver %s\n",ccebatch::variant(),ccebatch::lanes(),
           smalleig::on() ? "on" : "off (CCEX_SMALL_EIG=0)");
    // A51: a build whose kernel fails the start-up self-check never uses it (CceBatchWindow stays off), so there is nothing
    // to compare: every cluster goes the scalar path
    const char* why = nullptr;
    if (ccebatch::lanes() > 0 && !ccebatchUsable(&why)){
        printf("PASS off on this machine: the A51 self-check turned the batched kernel off (%s); every cluster goes the scalar path\n",why);
        return 0;
    }
    for (int a=1; a<argc; ++a){
        FILE* f = fopen(argv[a],"rb");
        if (!f){ perror(argv[a]); return 2; }
        std::vector<Rec> cl; Rec r;
        while (readRec(f,r)) cl.push_back(r);
        const bool readError = ferror(f);
        fclose(f);
        if (readError || cl.empty()){ fprintf(stderr,"%s: %s\n",argv[a],readError ? "read error" : "no record"); return 2; }
        std::string base = argv[a];
        base = base.substr(base.find_last_of('/')+1);
        const bool ens = cl[0].ens, rec = cl[0].rec;
        runCase(base+(rec ? " recursive" : " exact"), cl, ens, rec, false);
        runCase(base+(rec ? " as exact" : " as recursive"), cl, ens, !rec, false);
        if (!ens){
            runCase(base+" as ensemble Hahn exact", cl, true, false, true);
            runCase(base+" as ensemble Hahn recursive", cl, true, true, true);
        }
    }
    if (failed){ printf("FAIL %ld of %ld batched outputs differ from the scalar path\n",failed,compared); return 1; }
    printf("PASS %ld batched outputs bit-identical with the scalar path (%ld outputs not batched in this build)\n",
           compared,scalarOnly);
    return 0;
}
