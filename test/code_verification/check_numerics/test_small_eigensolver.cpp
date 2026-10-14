// A42: the fixed-sweep Jacobi of the small cCCE eigensystems (small_eigensolver.h) against Eigen and against itself.
//   - accuracy: eigenvalues within 1e-13 of max|H| of Eigen's, residual |H v - w v| and |V^dagger V - 1| below 1e-13
//     (relative to max|H|), for random Hermitian matrices, degenerate ones (U diag U^dagger with repeated values, a
//     multiple of 1), diagonal, real symmetric, nearly diagonal, the zero matrix, and the Hamiltonians of the dumped
//     clusters given on the command line (check_gpu/data/small_cluster/*.bin, as test_cce_batch reads them);
//   - lanes: every matrix solved alone (count 1), in full groups of 8 and in groups of 3 (padded) - the same bits;
//   - not converged after the sweeps (the caller then uses Eigen): allowed only for the dense random sets (4 x 4 needs
//     more than 4 sweeps for some of them), never for the structured sets and the dumped clusters; a NaN always is, and
//     so are an infinite part and a matrix outside [1e-130, 1e150] (review 20261013-221715 X1);
//   - the dump reader refuses an empty, cut or broken dump (X2), on cut copies of a dump made in memory.
// usage: test_small_eigensolver.out [DUMP...]
#include "calculation/small_eigensolver.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

typedef std::complex<double> cd;
typedef Eigen::MatrixXcd MatrixXcd;

int rank=0,nprocess=1,seed=19;
bool verbosity=false;

namespace {
struct Set { std::string name; int D; std::vector<MatrixXcd> H; bool mayFallBack = false; };

MatrixXcd randomUnitary(int D, std::mt19937_64& g){
    std::normal_distribution<double> n;
    MatrixXcd A(D,D);
    for (int i=0;i<D*D;++i) A.data()[i] = cd(n(g),n(g));
    return Eigen::HouseholderQR<MatrixXcd>(A).householderQ();
}

std::vector<Set> synthetic(){
    std::mt19937_64 g(20261013);
    std::normal_distribution<double> n;
    std::vector<Set> sets;
    for (int D : {2,4}){
        Set r{"random", D, {}}, deg{"degenerate", D, {}}, dia{"diagonal", D, {}}, rs{"real symmetric", D, {}},
            near{"nearly diagonal", D, {}}, z{"zero and multiples of 1", D, {}};
        for (int k=0;k<500;++k){
            MatrixXcd A(D,D);
            for (int i=0;i<D*D;++i) A.data()[i] = cd(n(g),n(g));
            const double scale = std::pow(10.0, double(k%13)-3.0);              // 1e-3 .. 1e9, as Hamiltonians in rad/us
            r.H.push_back(scale*(A+A.adjoint())/2.0);
            Eigen::VectorXd d(D);
            for (int i=0;i<D;++i) d(i) = double((i*k)%3);                       // values repeated
            const MatrixXcd U = randomUnitary(D,g);
            deg.H.push_back(scale*U*d.cast<cd>().asDiagonal()*U.adjoint());
            Eigen::VectorXd e(D); for (int i=0;i<D;++i) e(i) = scale*n(g);
            dia.H.push_back(e.cast<cd>().asDiagonal());
            Eigen::MatrixXd R(D,D); for (int i=0;i<D*D;++i) R.data()[i] = n(g);
            rs.H.push_back(scale*((R+R.transpose())/2.0).cast<cd>());
            MatrixXcd N = e.cast<cd>().asDiagonal();
            for (int i=0;i<D;++i) for (int j=0;j<i;++j){ const cd x(1e-12*scale*n(g),1e-12*scale*n(g)); N(i,j) = x; N(j,i) = std::conj(x); }
            near.H.push_back(N);
        }
        z.H.push_back(MatrixXcd::Zero(D,D));
        z.H.push_back(MatrixXcd::Identity(D,D)*cd(3.5,0));
        z.H.push_back(MatrixXcd::Identity(D,D)*cd(-1e8,0));
        r.mayFallBack = deg.mayFallBack = rs.mayFallBack = true;
        for (Set* s : {&r,&deg,&dia,&rs,&near,&z}) sets.push_back(*s);
    }
    return sets;
}

// The Hamiltonians of a dump of CCEX_SMALL_CLUSTER_DUMP (record layout: test_cce_batch.cpp). readEnd only at a clean end of
// the file, before the first byte of a header; a header, its metadata or a body cut short, a read error or a header that
// makes no sense is readBroken (review 20261013-221715 X2, as test_warp_cluster_cuda: a cut dump used to pass with the
// records before the cut, a 1-byte one with none).
enum ReadResult { readOk, readEnd, readBroken };
bool readBytes(FILE* f, void* p, size_t bytes){ return fread(p,1,bytes,f) == bytes; }
ReadResult readRecord(FILE* f, MatrixXcd& a, MatrixXcd& b){
    int h[6];
    const size_t got = fread(h,1,sizeof h,f);
    if (got == 0 && feof(f) && !ferror(f)) return readEnd;
    if (got != sizeof h) return readBroken;
    const int d = h[0], ens = h[1], rec = h[2], intervals = h[3], distinct = h[4], nstep = h[5];
    if ((d != 2 && d != 4) || (ens != 0 && ens != 1) || (rec != 0 && rec != 1) || intervals < 1 || intervals > 65536 ||
        distinct < 1 || distinct > 4 || nstep < 1 || nstep > (1<<26)) return readBroken;
    std::vector<char> skip(sizeof(double)*(1+distinct) + sizeof(int)*intervals);
    if (!readBytes(f,skip.data(),skip.size())) return readBroken;
    a.resize(d,d); b.resize(d,d);
    std::vector<cd> rest((ens ? 0 : d) + nstep);
    if (!readBytes(f,a.data(),sizeof(cd)*d*d) || !readBytes(f,b.data(),sizeof(cd)*d*d) ||
        !readBytes(f,rest.data(),sizeof(cd)*rest.size())) return readBroken;
    return readOk;
}
// every record of f into s2/s4 (by dimension); false if it is broken anywhere or has no record
bool readAll(FILE* f, Set& s2, Set& s4){
    MatrixXcd a, b;
    long records = 0;
    while (true){
        const ReadResult res = readRecord(f,a,b);
        if (res == readEnd) return !ferror(f) && records > 0;
        if (res == readBroken) return false;
        Set& s = a.rows() == 2 ? s2 : s4;
        s.H.push_back(a); s.H.push_back(b);
        ++records;
    }
}
void readDump(const char* path, std::vector<Set>& sets){
    FILE* f = fopen(path,"rb");
    if (!f){ perror(path); exit(2); }
    std::string base = path; base = base.substr(base.find_last_of('/')+1);
    Set s2{base,2,{}}, s4{base,4,{}};
    const bool ok = readAll(f,s2,s4);
    fclose(f);
    if (!ok){ fprintf(stderr,"empty, truncated or broken dump %s\n",path); exit(2); }
    for (Set* s : {&s2,&s4}) if (!s->H.empty()) sets.push_back(*s);
}

// The reader on cut copies of a two-record dump made here (in memory): the whole and the first record alone read; empty,
// one byte, cut inside the second header, inside its metadata, inside its body and one byte short do not, nor a header
// of dimension 3.
bool readerChecks(){
    std::vector<char> bytes;
    auto put = [&](const void* p, size_t n){ bytes.insert(bytes.end(),(const char*)p,(const char*)p+n); };
    size_t first = 0;
    for (int r=0; r<2; ++r){
        const int h[6] = {2,0,0,2,1,3};                            // d, ens, rec, intervals, distinct, nstep
        const double meta[2] = {0.01,0.5}; const int same[2] = {0,0};
        const cd H[4] = {cd(1,0),cd(0.5,0.25),cd(0.5,-0.25),cd(-1,0)}, rest[2+3] = {cd(1,0),cd(0,0),cd(1,0),cd(0.9,0),cd(0.8,0)};
        put(h,sizeof h); put(meta,sizeof meta); put(same,sizeof same); put(H,sizeof H); put(H,sizeof H); put(rest,sizeof rest);
        if (r == 0) first = bytes.size();
    }
    auto reads = [&](const std::vector<char>& data, size_t size, long& records){
        Set s2{"",2,{}}, s4{"",4,{}};
        if (size == 0){                                            // fmemopen of 0 bytes is not portable: an empty file
            FILE* e = tmpfile();
            const bool ok = e && readAll(e,s2,s4);
            if (e) fclose(e);
            records = 0; return ok;
        }
        FILE* m = fmemopen(const_cast<char*>(data.data()),size,"rb");
        const bool ok = m && readAll(m,s2,s4);
        if (m) fclose(m);
        records = long(s2.H.size()+s4.H.size())/2;
        return ok;
    };
    struct Cut { size_t size; bool reads; const char* what; } cuts[] = {
        {bytes.size(), true, "whole"}, {first, true, "after the first record"}, {0, false, "empty"}, {1, false, "one byte"},
        {first+10, false, "inside the second header"}, {first+24+4, false, "inside the second record's metadata"},
        {first+24+16+8+20, false, "inside the second body"}, {bytes.size()-1, false, "one byte short"}};
    bool ok = true;
    for (const Cut& c : cuts){
        long records = 0;
        const bool r = reads(bytes,c.size,records);
        const bool right = r == c.reads && (!r || records == (c.size == first ? 1 : 2));
        if (!right) printf("FAIL dump reader, cut %s (%zu bytes): %s\n",c.what,c.size,r ? "read" : "refused");
        ok &= right;
    }
    std::vector<char> bad(bytes);
    const int three = 3;
    std::memcpy(bad.data(),&three,sizeof three);
    long records = 0;
    if (reads(bad,bad.size(),records)){ printf("FAIL dump reader: a header of dimension 3 was read\n"); ok = false; }
    printf("%s dump reader: whole, first record, empty, 1 byte, cut header/metadata/body, 1 byte short, dimension 3\n",
           ok ? "ok" : "FAIL");
    return ok;
}

long failures = 0;

void check(const Set& s){
    const int D = s.D, m = int(s.H.size());
    std::vector<double> w1(size_t(m)*D), w8(size_t(m)*D), w3(size_t(m)*D);
    std::vector<cd> V1(size_t(m)*D*D), V8(size_t(m)*D*D), V3(size_t(m)*D*D);
    std::vector<char> ok1(m), ok8(m), ok3(m);
    auto run = [&](int group, std::vector<double>& w, std::vector<cd>& V, std::vector<char>& ok){
        for (int i=0;i<m;i+=group){
            const int c = std::min(group,m-i);
            const cd* h[smalleig::maxLanes]; double* ww[smalleig::maxLanes]; cd* vv[smalleig::maxLanes]; bool k[smalleig::maxLanes];
            for (int l=0;l<c;++l){ h[l] = s.H[i+l].data(); ww[l] = &w[size_t(i+l)*D]; vv[l] = &V[size_t(i+l)*D*D]; }
            smalleig::solve(D,c,h,ww,vv,k);
            for (int l=0;l<c;++l) ok[i+l] = k[l];
        }
    };
    run(1,w1,V1,ok1); run(smalleig::maxLanes,w8,V8,ok8); run(3,w3,V3,ok3);
    long laneDiffer = 0, notOk = 0;
    double dw = 0, res = 0, orth = 0;
    for (int i=0;i<m;++i){
        laneDiffer += ok1[i] != ok8[i] || ok1[i] != ok3[i]
            || memcmp(&w1[size_t(i)*D],&w8[size_t(i)*D],sizeof(double)*D) || memcmp(&w1[size_t(i)*D],&w3[size_t(i)*D],sizeof(double)*D)
            || memcmp(&V1[size_t(i)*D*D],&V8[size_t(i)*D*D],sizeof(cd)*D*D) || memcmp(&V1[size_t(i)*D*D],&V3[size_t(i)*D*D],sizeof(cd)*D*D);
        if (!ok1[i]){ ++notOk; continue; }
        const MatrixXcd& H = s.H[i];
        const double hmax = std::max(H.cwiseAbs().maxCoeff(), 1e-300);
        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
        const Eigen::Map<const Eigen::VectorXd> w(&w1[size_t(i)*D],D);
        const Eigen::Map<const MatrixXcd> V(&V1[size_t(i)*D*D],D,D);
        dw = std::max(dw, (w-es.eigenvalues()).cwiseAbs().maxCoeff()/hmax);
        res = std::max(res, (H*V - V*w.cast<cd>().asDiagonal()).cwiseAbs().maxCoeff()/hmax);
        orth = std::max(orth, (V.adjoint()*V - MatrixXcd::Identity(D,D)).cwiseAbs().maxCoeff());
        for (int k=1;k<D;++k) if (w(k) < w(k-1)) ++laneDiffer;     // ascending
    }
    const bool pass = laneDiffer == 0 && (notOk == 0 || s.mayFallBack) && dw <= 1e-13 && res <= 1e-13 && orth <= 1e-13;
    printf("%-58s D %d %6d matrices: lanes 1/8/3 differ %ld, not converged %ld, max dw %.1e, residual %.1e, "
           "orthonormality %.1e  %s\n",s.name.c_str(),D,m,laneDiffer,notOk,dw,res,orth,pass ? "ok" : "FAIL");
    if (!pass) ++failures;
}
}

int main(int argc, char** argv){
    std::vector<Set> sets = synthetic();
    for (int a=1;a<argc;++a) readDump(argv[a],sets);
    for (const Set& s : sets) check(s);
    // not finite: reported, so that the caller solves it with Eigen
    for (int D : {2,4}){
        MatrixXcd H = MatrixXcd::Identity(D,D); H(1,0) = cd(std::nan(""),0); H(0,1) = H(1,0);
        double w[4]; cd V[16]; double* wp = w; cd* Vp = V; bool ok = true; const cd* h = H.data();
        smalleig::solve(D,1,&h,&wp,&Vp,&ok);
        printf("NaN matrix D %d: %s\n",D,ok ? "reported converged - FAIL" : "not converged (Eigen takes it) ok");
        if (ok) ++failures;
    }
    // the input range (review 20261013-221715 X1): H = s (e_0 e_1^T + e_1 e_0^T). Inside [1e-130, 1e150] the Jacobi
    // answers +-s; outside it, or with an infinite part, it reports not converged (Eigen, which answers +-s, takes it)
    // - it used to report (0, 0) for 1e-151 and (-inf, inf) for 1e155 as converged
    for (int D : {2,4}){
        for (double v : {1e-151, 1e-160, 1e153, 1e155, 1e200, 1e-120, 1.0, 1e140}){
            MatrixXcd H = MatrixXcd::Zero(D,D); H(0,1) = H(1,0) = v;
            double w[4]; cd V[16]; double* wp = w; cd* Vp = V; bool ok = false; const cd* h = H.data();
            smalleig::solve(D,1,&h,&wp,&Vp,&ok);
            const bool inside = v >= 1e-130 && v <= 1e150;
            Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
            const Eigen::Map<const Eigen::VectorXd> got(w,D);
            const Eigen::VectorXd& want = es.eigenvalues();      // -s, 0.., s
            const bool eigenRight = es.info() == Eigen::Success && std::fabs(want(0)+v) <= 1e-13*v && std::fabs(want(D-1)-v) <= 1e-13*v;
            const bool right = ok == inside && eigenRight && (!ok || (got-want).cwiseAbs().maxCoeff() <= 1e-13*v);
            printf("off-diagonal %.0e D %d: %s, %s\n",v,D,ok ? "converged" : "not converged (Eigen takes it)",right ? "ok" : "FAIL");
            if (!right) ++failures;
        }
        for (int where : {0,1,2}) for (double inf : {INFINITY, -INFINITY}){
            MatrixXcd H = MatrixXcd::Identity(D,D);
            if (where == 0) H(0,0) = inf; else if (where == 1){ H(0,1) = inf; H(1,0) = inf; } else { H(1,0) = cd(0,inf); H(0,1) = cd(0,-inf); }
            double w[4]; cd V[16]; double* wp = w; cd* Vp = V; bool ok = true; const cd* h = H.data();
            smalleig::solve(D,1,&h,&wp,&Vp,&ok);
            printf("%s in %s D %d: %s\n",inf > 0 ? "+inf" : "-inf",where == 0 ? "the diagonal" : where == 1 ? "the off-diagonal" : "an imaginary part",D,
                   ok ? "reported converged - FAIL" : "not converged ok");
            if (ok) ++failures;
        }
    }
    if (!readerChecks()) ++failures;
    if (failures){ printf("FAIL %ld sets\n",failures); return 1; }
    printf("PASS small eigensolver (A42): %zu sets\n",sets.size());
    return 0;
}
