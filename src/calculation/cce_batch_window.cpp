// A30: the window of small cCCE clusters calculate() hands to the batched kernel (cce_batch_window.h).
#include "../../include/calculation/cce_batch_window.h"
#include "../../include/calculation/small_hamiltonian.h"
#include "../../include/calculation/sparse_propagator.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/small_eigensolver.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

bool ccebatchPlan(Config* cnf, Pulse* pls, ccebatch::Plan& plan){
    plan = ccebatch::Plan{};
    plan.ensemble = Config_getNstate(cnf) == 0;
    plan.recursive = strcasecmp(cnf->evolution_update,"recursive") == 0;   // as calCoherenceCceInto decides
    plan.nstep = Config_getNstep(cnf);
    plan.deltat = (double)Config_getDeltat(cnf);
    const int npulse = Pulse_getNpulse(pls);
    if (plan.ensemble){
        // the Hahn echo only: the CPMG/Ramsey/Manual ensemble branches stay scalar
        plan.intervals = npulse+1;
        return npulse == 1 && strcasecmp(Pulse_getPulsename(pls),"Manual") != 0;
    }
    if (npulse+1 > ccebatch::maxIntervals) return false;
    plan.intervals = npulse+1;
    double** sequence = Pulse_getSequence(pls);
    int* same = Pulse_getSequenceIndices(pls);
    for (int p=0; p<=npulse; ++p){ plan.fraction[p] = sequence[p][2]; plan.same[p] = same[p]; }
    return true;
}

bool ccebatchPackLane(const CceWorkspace& ws, bool ensemble, int D, ccebatch::Pack& p, int l){
    bool finite = true;
    auto put = [&](double (*re)[ccebatch::maxLanes], double (*im)[ccebatch::maxLanes], int e, doublec z){
        re[e][l] = z.real(); im[e][l] = z.imag();
        finite = finite && std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    for (int i=0; i<D; ++i)
        for (int k=0; k<D; ++k){
            put(p.Mre,p.Mim,i*D+k,ws.M(i,k));
            if (!ensemble) put(p.Are,p.Aim,i*D+k,ws.Madj(i,k));
        }
    for (int k=0; k<D; ++k){
        if (!ensemble){ put(p.r0re,p.r0im,k,ws.rowStart(k)); put(p.c0re,p.c0im,k,ws.colStart(k)); }
        p.Eu[k][l] = ws.valUp(k).real(); p.Ed[k][l] = ws.valDown(k).real();
        finite = finite && std::isfinite(p.Eu[k][l]) && std::isfinite(p.Ed[k][l]);
    }
    return finite;
}

namespace {
// A51: one made-up case of the self-check: n clusters of bath dimension D through the scalar path and through the window's
// (grouped eigensolves, cceSetupSmall, a lane each, the kernel); false with the first differing value in why
bool selfCheckCase(int D, bool ensemble, bool recursive, int intervals, char* why, size_t whySize){
    constexpr int nstep = 40;                         // past two phase anchors (16, 32) of the recursive update
    const int B = ccebatch::lanes(), n = B+1;         // a full group of lanes and a padded one
    Config* cnf = Config_init();
    strcpy(cnf->method,"cce"); strcpy(cnf->quantity,"coherence");
    strcpy(cnf->evolution_update, recursive ? "recursive" : "exact");
    cnf->nstep = nstep; cnf->deltat = 0.01; cnf->nstate = ensemble ? 0 : 1; cnf->profile = false;
    // the Pulse as Pulse builds it: Ramsey, Hahn echo (1/2, 1/2) or CPMG with two pulses (1/4, 1/2, 1/4; slots 0, 1, 0)
    double rows[3][3] = {}; double* seq[3] = {rows[0],rows[1],rows[2]}; int idx[3] = {0,0,0};
    for (int p=0; p<intervals; ++p) rows[p][2] = intervals == 1 ? 1.0 : intervals == 2 ? 0.5 : p == 1 ? 0.5 : 0.25;
    if (intervals == 3) idx[1] = 1;
    Pulse pulse{};
    pulse.npulse = intervals-1; pulse.sequence = seq; pulse.sequence_indices = idx;
    strcpy(pulse.pulsename, intervals == 1 ? "Ramsey" : intervals == 2 ? "HahnEcho" : "CPMG");
    ccebatch::Plan plan;
    bool same = true;
    if (ccebatchPlan(cnf,&pulse,plan) && ccebatch::covers(plan,D)){
        // Hermitian Hamiltonians of the size of P1's (rad/us), Hbeta a little away from Halpha; a basis state. A fixed
        // generator of its own (the run's random numbers are not touched).
        unsigned long long x = 0x9e3779b97f4a7c15ULL*(unsigned long long)(D*64+ensemble*16+recursive*8+intervals);
        auto uniform = [&x](double a){ x ^= x << 13; x ^= x >> 7; x ^= x << 17; return a*(2.0*double(x >> 11)/9007199254740992.0-1.0); };
        std::vector<MatrixXcd> Ha(n), Hb(n), psi(n);
        for (int c=0; c<n; ++c){
            Ha[c].resize(D,D); Hb[c].resize(D,D);
            for (int i=0; i<D; ++i){
                Ha[c](i,i) = doublec(uniform(30.0),0.0); Hb[c](i,i) = Ha[c](i,i)+doublec(uniform(1.0),0.0);
                for (int k=i+1; k<D; ++k){
                    Ha[c](i,k) = doublec(uniform(3.0),uniform(3.0)); Ha[c](k,i) = std::conj(Ha[c](i,k));
                    Hb[c](i,k) = Ha[c](i,k)+doublec(uniform(0.5),uniform(0.5)); Hb[c](k,i) = std::conj(Hb[c](i,k));
                }
            }
            psi[c] = MatrixXcd::Zero(D,1); psi[c](c % D,0) = 1.0;
        }
        std::vector<const MatrixXcd*> pa(n), pb(n), pp(n);
        for (int c=0; c<n; ++c){ pa[c] = &Ha[c]; pb[c] = &Hb[c]; pp[c] = &psi[c]; }
        std::vector<std::complex<double>> scalarOut(size_t(n)*nstep), batchOut(size_t(n)*nstep);
        std::vector<char> batched;
        ccebatchBothPaths(cnf,&pulse,plan,D,n,pa.data(),pb.data(),pp.data(),scalarOut.data(),batchOut.data(),batched);
        for (int c=0; c<n && same; ++c)
            if (!batched[c]){ snprintf(why,whySize,"a test cluster has a value that is not finite"); same = false; }
        for (size_t i=0; same && i<scalarOut.size(); ++i)
            if (std::memcmp(&scalarOut[i],&batchOut[i],sizeof(doublec))){
                snprintf(why,whySize,"bath dimension %d, %s, %s, %s: test cluster %zu, step %zu differs (scalar %a%+ai, batched %a%+ai)",
                         D, ensemble ? "ensemble" : "one state", recursive ? "recursive" : "exact",
                         intervals == 1 ? "Ramsey" : intervals == 2 ? "Hahn echo" : "CPMG 2", i/nstep, i%nstep,
                         scalarOut[i].real(), scalarOut[i].imag(), batchOut[i].real(), batchOut[i].imag());
                same = false;
            }
    }
    Config_freeAll(cnf);
    return same;
}
}

bool ccebatchUsable(const char** why){
    static int state = -1;                            // -1 not checked yet, 0 off, 1 on
    static char reason[400] = "";
    if (state < 0){
        const bool kernel = ccebatch::lanes() > 0;
        bool ok = kernel;
        long long solved = 0, fallbacks = 0;
        smalleig::takeCounts(solved,fallbacks);       // the made-up clusters' eigensolves are not the run's [profile] count
        if (!kernel) snprintf(reason,sizeof reason,"no kernel for this build (a compiler, version or instruction set it was not checked on)");
        for (int D : {2, 4})
            for (int e=0; e<2 && ok; ++e)
                for (int r=0; r<2 && ok; ++r)
                    for (int iv=1; iv<=3 && ok; ++iv)
                        if (!e || iv == 2) ok = selfCheckCase(D,e,r,iv,reason,sizeof reason);   // the ensemble: Hahn only
        long long checkSolved = 0, checkFallbacks = 0;
        smalleig::takeCounts(checkSolved,checkFallbacks);
        smalleig::count(solved,fallbacks);
        int initialized = 0, finalized = 0;
        MPI_Initialized(&initialized);
        if (initialized) MPI_Finalized(&finalized);
        if (kernel && initialized && !finalized){      // every rank the same decision
            int all = ok;
            MPI_Allreduce(MPI_IN_PLACE,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
            if (ok && !all) snprintf(reason,sizeof reason,"the self-check failed on another rank");
            ok = all;
        }
        state = ok;
        if (kernel && !ok && rank == 0)
            printf("[cpu-batch] off: the batched small-cluster kernel (A30, %s) does not reproduce the scalar path in this "
                   "build (%s) - every cluster goes the scalar way (A51 self-check)\n", ccebatch::variant(), reason);
    }
    if (why) *why = reason;
    return state == 1;
}

bool ccebatchRunTakes(QubitArray* qa, Config* cnf, Pulse* pls, bool isGCCE, ccebatch::Plan& plan, bool dimOk[2], const char** why){
    dimOk[0] = dimOk[1] = false;
    if (ccebatch::lanes() == 0){ *why = "no kernel for this build"; return false; }
    if (isGCCE){ *why = "method gCCE"; return false; }
    if (getenv("CCEX_SMALL_CLUSTER_DUMP")){ *why = "CCEX_SMALL_CLUSTER_DUMP is set"; return false; }
    if (getenv("CCEX_CPU_BATCH") && !strcmp(getenv("CCEX_CPU_BATCH"),"0")){ *why = "CCEX_CPU_BATCH=0"; return false; }
    if (QubitArray_getNqubit(qa) != 1 || strcasecmp(Config_getQuantity(cnf),"coherence")){ *why = "not one qubit with quantity coherence"; return false; }
    if (!strcasecmp(cnf->backend,"cuda")){ *why = "backend cuda"; return false; }   // every cluster on the card
    if (!ccebatchPlan(cnf,pls,plan)){
        *why = Config_getNstate(cnf) == 0 ? "an ensemble other than the Hahn echo" : "more pulse intervals than the kernel has";
        return false;
    }
    for (int k=0; k<2; ++k){
        const int d = k == 0 ? 2 : 4;
        dimOk[k] = ccebatch::covers(plan,d) && !sparsePropagatorFor(cnf,d) && !computeUseCuda(cnf,d);
    }
    if (!dimOk[0] && !dimOk[1]){ *why = "its clusters of bath dimension 2 and 4 go another way (propagator sparse, the CUDA cluster path)"; return false; }
    if (!ccebatchUsable(why)){ dimOk[0] = dimOk[1] = false; return false; }
    return true;
}

CceBatchWindow::CceBatchWindow(QubitArray* qa_, BathArray* ba_, DefectArray* dfa_, Config* cnf_, Pulse* pls, bool isGCCE)
    : qa(qa_), ba(ba_), dfa(dfa_), cnf(cnf_){
    B = ccebatch::lanes();
    const char* why = nullptr;
    if (!ccebatchRunTakes(qa,cnf,pls,isGCCE,plan,dimOk,&why)) return;
    nstep = plan.nstep;
    nqubit = QubitArray_getNqubit(qa);
    // the coherences of a window: 256 clusters, fewer when they would pass 4 MiB; none (every cluster the scalar way, the
    // same bits) when one group of lanes alone would (review 20261013-221715 X5: until then at least 4 groups, whatever nstep)
    W = int(std::min<long>(256, (4L<<20)/(16L*std::max(1,nstep))));
    if (W < B){
        if (rank == 0) printf("[cpu-batch] off: the %d time points of one group of %d lanes pass the 4 MiB of a window - every "
                              "cluster goes the scalar way\n",nstep,B);
        return;
    }
    enabled = true;      // the storage is made by the first window that needs it (allocate)
}

// X5: the storage of a window of n clusters, made when a window is first filled - not on a rank without clusters or in a
// run whose orders do not use the batch - and grown (never past W) when a later window is larger
void CceBatchWindow::allocate(int n){
    if (n <= capacity) return;
    capacity = n;
    results.resize(size_t(n)*nstep);
    dims.assign(n,0);
    made.assign(n,nullptr);              // empty: fill dropped the window before
    Ha.resize(n); Hb.resize(n); psiAt.resize(n); buildAt.assign(n,0.0);
    eigW.resize(size_t(2)*n*4); eigV.resize(size_t(2)*n*16); eigOk.assign(size_t(2)*n,false);
    if (!pack){
        pack = static_cast<ccebatch::Pack*>(std::aligned_alloc(64, 2*sizeof(ccebatch::Pack)));
        std::memset(static_cast<void*>(pack), 0, 2*sizeof(ccebatch::Pack));
    }
}

CceBatchWindow::~CceBatchWindow(){
    drop();
    std::free(pack);
}

void CceBatchWindow::drop(){
    for (BathArray*& b : made) if (b){ BathArray_freeAll(b); b = nullptr; }
    first = 1; last = 0;
}

void CceBatchWindow::beginOrder(int n, int** c, int count, bool u){
    if (!enabled) return;
    drop();
    order = n; clusters = c; ncluster = count; use = u;
}

const std::complex<double>* CceBatchWindow::get(int ic, int& dim, BathArray** out, int end){
    *out = nullptr;
    if (!enabled || !use) return nullptr;
    if (ic < first || ic > last) fill(ic,end);
    const int at = ic-first;
    if (dims[at]){ dim = dims[at]; return &results[size_t(at)*nstep]; }
    *out = made[at]; made[at] = nullptr;
    return nullptr;
}

void CceBatchWindow::fill(int ic, int end){
    drop();
    first = ic; last = std::min({ncluster, ic+W-1, std::max(ic,end-1)});
    allocate(last-first+1);
    filled[0] = filled[1] = 0;
    const bool profile = cnf->profile;
    // the window's clusters of dimension 2 and 4 the batch takes: bath array, Hamiltonians (and the state)
    std::vector<int> mine;
    for (int j=first; j<=last; ++j){
        const int at = j-first;
        dims[at] = 0;
        BathArray* b = createBathArray(clusters[j]+1, order, ba, dfa, nqubit);
        const int d = BathArray_dim(b);
        const int k = d == 2 ? 0 : d == 4 ? 1 : -1;
        if (k < 0 || !dimOk[k]){ made[at] = b; continue; }
        const double t0 = profile ? MPI_Wtime() : 0.0;
        cceClusterHamiltonians(qa,b,cnf,Ha[at],Hb[at]);
        if (profile) buildAt[at] = MPI_Wtime()-t0;
        if (!plan.ensemble) psiAt[at] = BathArray_Psi0(b);
        made[at] = b;                                  // freed below when its setup is finite
        dims[at] = d;
        mine.push_back(at);
    }
    // A42: the eigensystems, smalleig::maxLanes matrices at a time per dimension (the two of a cluster, then the next
    // cluster's), with the bits each matrix gets alone in cceSetup
    const bool jacobi = smalleig::on();
    if (jacobi){
        for (int d : {2, 4}){
            const double t0 = profile ? MPI_Wtime() : 0.0;
            CceEigenGroups groups(d,eigW.data(),eigV.data(),eigOk.data());
            long long matrices = 0;
            for (int at : mine){
                if (dims[at] != d) continue;
                groups.add(2*at,Ha[at]);
                groups.add(2*at+1,Hb[at]);
                matrices += 2;
            }
            groups.flush();
            smalleig::count(groups.solved,groups.fallbacks);
            if (profile && matrices) eigEach[d == 2 ? 0 : 1] = (MPI_Wtime()-t0)/double(matrices/2);
        }
    }
    // the rest of each cluster's setup in order, a lane, the lanes evolved when a group is full
    for (int at : mine){
        const int d = dims[at], k = d == 2 ? 0 : 1;
        const double t1 = profile ? MPI_Wtime() : 0.0;
        ws.Halpha.swap(Ha[at]); ws.Hbeta.swap(Hb[at]);
        const MatrixXcd* psi = plan.ensemble ? nullptr : &psiAt[at];
        ccebatchSetup(ws,psi,jacobi,at,eigW.data(),eigV.data(),eigOk.data());
        ws.Halpha.swap(Ha[at]); ws.Hbeta.swap(Hb[at]);
        const int lane = filled[k];
        if (!ccebatchPackLane(ws,plan.ensemble,d,pack[k],lane)){ dims[at] = 0; continue; }   // the scalar path decides
        BathArray_freeAll(made[at]); made[at] = nullptr;
        slot[k][lane] = at;
        if (profile){ build[k][lane] = buildAt[at]; setup[k][lane] = MPI_Wtime()-t1 + (jacobi ? eigEach[k] : 0.0); }
        if (++filled[k] == B) run(k);
    }
    for (int k=0; k<2; ++k) if (filled[k]) run(k);
}

void CceEigenGroups::add(int m, const MatrixXcd& H){
    h[n] = H.data(); wl[n] = &w[size_t(m)*4]; Vl[n] = &V[size_t(m)*16]; slot[n] = m;
    if (++n == smalleig::maxLanes) flush();
}

void CceEigenGroups::flush(){
    if (!n) return;
    bool converged[smalleig::maxLanes];
    smalleig::solve(D,n,h,wl,Vl,converged);
    for (int g=0; g<n; ++g){
        ok[slot[g]] = converged[g];
        if (converged[g]) ++solved; else ++fallbacks;
    }
    n = 0;
}

void ccebatchSetup(CceWorkspace& ws, const MatrixXcd* psi, bool jacobi, int c,
                   const double* w, const std::complex<double>* V, const char* ok){
    if (!jacobi){ cceSetup(ws,psi); return; }
    const double* wc[2] = {nullptr,nullptr};
    const doublec* Vc[2] = {nullptr,nullptr};
    for (int s=0; s<2; ++s) if (ok[2*c+s]){ wc[s] = &w[size_t(2*c+s)*4]; Vc[s] = &V[size_t(2*c+s)*16]; }
    cceSetupSmall(ws,psi,wc,Vc);
}

void ccebatchBothPaths(Config* cnf, Pulse* pls, const ccebatch::Plan& plan, int D, int n,
                       const MatrixXcd* const* Ha, const MatrixXcd* const* Hb, const MatrixXcd* const* psi,
                       std::complex<double>* scalarOut, std::complex<double>* batchOut, std::vector<char>& batched){
    const int nstep = plan.nstep, B = ccebatch::lanes();
    CceWorkspace ws;
    std::vector<MatrixXcd> res(nstep);
    for (int c=0; c<n; ++c){
        const MatrixXcd* p = plan.ensemble ? nullptr : psi[c];
        ws.Halpha = *Ha[c]; ws.Hbeta = *Hb[c];
        if (!cceSetupSmallAlone(ws,p)) cceSetup(ws,p);
        cceEvolve(cnf,pls,res.data(),ws);
        for (int i=0; i<nstep; ++i) scalarOut[size_t(c)*nstep+i] = res[i](0,0);
    }
    const bool jacobi = smalleig::on();
    std::vector<double> eigW(size_t(n)*8); std::vector<doublec> eigV(size_t(n)*32); std::vector<char> eigOk(size_t(n)*2,0);
    if (jacobi){
        CceEigenGroups groups(D,eigW.data(),eigV.data(),eigOk.data());
        for (int c=0; c<n; ++c){ groups.add(2*c,*Ha[c]); groups.add(2*c+1,*Hb[c]); }
        groups.flush();
    }
    ccebatch::Pack pack{};
    std::complex<double>* out[ccebatch::maxLanes] = {};
    int count = 0;
    batched.assign(n,0);
    for (int c=0; c<n; ++c){
        ws.Halpha = *Ha[c]; ws.Hbeta = *Hb[c];
        ccebatchSetup(ws,plan.ensemble ? nullptr : psi[c],jacobi,c,eigW.data(),eigV.data(),eigOk.data());
        if (!ccebatchPackLane(ws,plan.ensemble,D,pack,count)) continue;
        batched[c] = 1;
        out[count++] = &batchOut[size_t(c)*nstep];
        if (count == B){ ccebatchEvolvePadded(plan,D,pack,count,out); count = 0; }
    }
    if (count) ccebatchEvolvePadded(plan,D,pack,count,out);
}

void ccebatchEvolvePadded(const ccebatch::Plan& plan, int D, ccebatch::Pack& p, int count, std::complex<double>* const* out){
    double (*rows)[ccebatch::maxLanes] = reinterpret_cast<double (*)[ccebatch::maxLanes]>(&p);
    const size_t nrows = sizeof(ccebatch::Pack)/sizeof(rows[0]);
    for (int l=count; l<ccebatch::lanes(); ++l) for (size_t r=0; r<nrows; ++r) rows[r][l] = rows[r][count-1];
    std::complex<double>* lanes[ccebatch::maxLanes] = {};
    for (int l=0; l<count; ++l) lanes[l] = out[l];
    ccebatch::evolve(plan,D,p,lanes);
}

// Evolves the filled lanes of dimension 2 (k 0) or 4 (k 1)
void CceBatchWindow::run(int k){
    const int count = filled[k], d = k == 0 ? 2 : 4;
    std::complex<double>* out[ccebatch::maxLanes] = {};
    for (int l=0; l<count; ++l) out[l] = &results[size_t(slot[k][l])*nstep];
    const double t0 = cnf->profile ? MPI_Wtime() : 0.0;
    ccebatchEvolvePadded(plan,d,pack[k],count,out);
    if (cnf->profile){
        const double each = (MPI_Wtime()-t0)/count;
        for (int l=0; l<count; ++l) computeProfile(cnf,"cce","cpu-batch",d,build[k][l],setup[k][l],each);
    }
    filled[k] = 0;
}
