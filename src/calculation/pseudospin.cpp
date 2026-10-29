/**
 * @file pseudospin.cpp
 * @brief method pseudospin: the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops
 * (Seo et al., Nat. Commun. 7, 12935 (2016), SI eqs. 20-24; Lee, Park, Seo, npj 2D Mater. Appl. 6, 60 (2022), SI note 4).
 *
 * A pair's cCCE Hamiltonians H_alpha, H_beta are cut to the pseudospin Hamiltonian: the diagonal (Zeeman, m_s A_zz, the
 * quadrupole's diagonal, J_zz) and the flip-flop elements |m1,m2> <-> |m1-1,m2+1>; the pseudo-secular hyperfine, the
 * quadrupole's off-diagonal and the non-secular dipolar terms are dropped. With L_cut the pair under the cut Hamiltonians:
 *  - pair "block": L_cut, the pair exactly;
 *  - pair "channel": every channel c alone, a two-level system l_c (l_c,diag without its flip-flop), and
 *    L_cut ~ L_diag + (2/d) sum_c (l_c - l_c,diag), d = (2I_1+1)(2I_2+1), L_diag the pair without flip-flops; exact for two
 *    spins 1/2, for larger spins it misses what couples the channels of a chain (a state shared by two channels).
 * The pair's term is P = L_cut / (L_1,cut L_2,cut); the cluster's value is P L_1 L_2, L_k the full single-spin coherence
 * (what the order-1 clusters give), so the run's coherence is prod_k L_k prod_pairs P.
 */
#include "../../include/calculation/pseudospin.h"
#include "../../include/calculation/small_hamiltonian.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {
FILE* channelFile = nullptr;
std::string channelName;

// the ensemble coherence of the bath under Ha, Hb with the run's pulses, at every time point
std::vector<doublec> evolve(const MatrixXcd& Ha, const MatrixXcd& Hb, Config* cnf, Pulse* pls){
    CceWorkspace ws;
    ws.Halpha = Ha;
    ws.Hbeta = Hb;
    cceSetup(ws,nullptr);
    std::vector<MatrixXcd> r(Config_getNstep(cnf));
    cceEvolve(cnf,pls,r.data(),ws);
    std::vector<doublec> v(r.size());
    for (size_t i=0; i<r.size(); ++i) v[i] = r[i](0,0);
    return v;
}

MatrixXcd diagonalOf(const MatrixXcd& H){ return MatrixXcd(H.diagonal().asDiagonal()); }

// spin k of the pair alone: its full coherence and its coherence under the diagonal (cut) Hamiltonians
void single(QubitArray* qa, BathArray* ba, int k, Config* cnf, Pulse* pls, std::vector<doublec>& full, std::vector<doublec>& cut){
    BathArray* one = BathArray_init();
    BathArray_setNspin(one,1);
    BathArray_allocBath(one,QubitArray_getNqubit(qa));
    BathArray_setBath_i(one,BathArray_getBath_i(ba,k),0,QubitArray_getNqubit(qa));
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,one,cnf,Ha,Hb);
    BathArray_freeAll(one);
    full = evolve(Ha,Hb,cnf,pls);
    cut = evolve(diagonalOf(Ha),diagonalOf(Hb),cnf,pls);
}
}

void calCoherencePseudospinInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result){
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,ba,cnf,Ha,Hb);
    // the product basis: spin 0 the slow index, row a d2 + b, m = I - a (Pauli_matrix_Z, addExpandedHamiltonian)
    const int d1 = int(std::lround(2*BathArray_getBath_i_spin(ba,0)))+1, d2 = int(std::lround(2*BathArray_getBath_i_spin(ba,1)))+1;
    const int d = d1*d2;
    std::vector<std::pair<int,int>> channels;   // (|m1,m2>, |m1-1,m2+1>)
    for (int a=0; a+1<d1; ++a)
        for (int b=1; b<d2; ++b) channels.push_back({a*d2+b, (a+1)*d2+b-1});
    MatrixXcd Ca = diagonalOf(Ha), Cb = diagonalOf(Hb);
    for (auto [r,s] : channels){
        Ca(r,s) = Ha(r,s); Ca(s,r) = Ha(s,r);
        Cb(r,s) = Hb(r,s); Cb(s,r) = Hb(s,r);
    }

    std::vector<doublec> pair;
    if (cnf->pseudospin_block) pair = evolve(Ca,Cb,cnf,pls);
    else {
        pair = evolve(diagonalOf(Ca),diagonalOf(Cb),cnf,pls);
        for (auto [r,s] : channels){
            const int idx[2] = {r,s};
            MatrixXcd ca(2,2), cb(2,2);
            for (int i=0; i<2; ++i) for (int j=0; j<2; ++j){ ca(i,j) = Ca(idx[i],idx[j]); cb(i,j) = Cb(idx[i],idx[j]); }
            const std::vector<doublec> l = evolve(ca,cb,cnf,pls), l0 = evolve(diagonalOf(ca),diagonalOf(cb),cnf,pls);
            for (size_t i=0; i<pair.size(); ++i) pair[i] += 2.0/d*(l[i]-l0[i]);
        }
    }
    std::vector<doublec> full1, cut1, full2, cut2;
    single(qa,ba,0,cnf,pls,full1,cut1);
    single(qa,ba,1,cnf,pls,full2,cut2);
    for (size_t i=0; i<pair.size(); ++i)
        result[i] = MatrixXcd::Constant(1,1,pair[i]/(cut1[i]*cut2[i])*full1[i]*full2[i]);

    if (channelFile){   // kHz: the Hamiltonians are in rad/ms
        const double* x1 = BathArray_getBath_i_xyz(ba,0);
        const double* x2 = BathArray_getBath_i_xyz(ba,1);
        const double I1 = 0.5*(d1-1), I2 = 0.5*(d2-1);
        for (auto [r,s] : channels)
            fprintf(channelFile,"%s %.4f %.4f %.4f %s %.4f %.4f %.4f %g %g %.16e %.16e %.16e\n",
                    BathArray_getBath_i_name(ba,0),x1[0],x1[1],x1[2],BathArray_getBath_i_name(ba,1),x2[0],x2[1],x2[2],
                    I1-r/d2,I2-r%d2,2*std::abs(Ca(r,s))/(2*M_PI),
                    (Ca(r,r)-Ca(s,s)).real()/(2*M_PI),(Cb(r,r)-Cb(s,s)).real()/(2*M_PI));
    }
}

void pseudospinChannelsBegin(Output* op){
    channelName = std::string(Output_getOutfile(op))+"_pseudospin_channels";
    channelFile = std::fopen((channelName+".rank"+std::to_string(rank)).c_str(),"w");
    if (!channelFile){
        fprintf(stderr,"Error: cannot write %s.rank%d\n",channelName.c_str(),rank);
        ccexFatal();
    }
}

void pseudospinChannelsEnd(){
    std::fclose(channelFile);
    channelFile = nullptr;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != 0) return;
    FILE* out = std::fopen(channelName.c_str(),"w");
    if (!out){
        fprintf(stderr,"Error: cannot write %s\n",channelName.c_str());
        ccexFatal();
    }
    fprintf(out,"# name1 x1 y1 z1 name2 x2 y2 z2 (A) m1 m2: the channel |m1,m2> <-> |m1-1,m2+1>; Omega (kHz, the flip-flop "
                "splitting 2|<m1,m2|H|m1-1,m2+1>|), Delta_alpha, Delta_beta (kHz, E(m1,m2) - E(m1-1,m2+1) in the qubit's alpha "
                "and beta states)\n");
    for (int r=0; r<nprocess; ++r){
        const std::string part = channelName+".rank"+std::to_string(r);
        if (FILE* in = std::fopen(part.c_str(),"r")){
            char buf[4096];
            size_t n;
            while ((n = std::fread(buf,1,sizeof buf,in)) > 0) std::fwrite(buf,1,n,out);
            std::fclose(in);
            std::remove(part.c_str());
        }
    }
    std::fclose(out);
}
