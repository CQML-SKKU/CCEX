/**
 * @file pseudospin.cpp
 * @brief method pseudospin: the pair clusters of the cCCE in the pseudospin model of nuclear flip-flops
 * (Seo et al., Nat. Commun. 7, 12935 (2016), SI eqs. 20-24; Lee, Park, Seo, npj 2D Mater. Appl. 6, 60 (2022), SI note 4).
 *
 * A pair's cCCE Hamiltonian in each qubit state, H_alpha and H_beta, is written in the product basis of the two spins'
 * eigenstates in that state (each spin's own H_alpha or H_beta: Zeeman, its hyperfine with the qubit, the whole
 * quadrupole), and cut there to the pseudospin Hamiltonian: the diagonal and every element between two product states that
 * differ in both spins, whatever the Hamiltonian gives there (the flip-flops |m1,m2> <-> |m1-1,m2+1>, the flip-flips
 * |m1,m2> <-> |m1+1,m2+1> of the non-secular dipolar coupling, and with a quadrupole that mixes the |m> the double-quantum
 * ones such as |m1,m2> <-> |m1-2,m2+2>). Dropped are the elements of the spins' coupling that change one spin only in that
 * basis; the single spins themselves are exact (their own Hamiltonians are diagonal there).
 *  - pair "block": the cut Hamiltonians, back in the |m> basis, evolved exactly: L_cut.
 *    The cluster's value is L_cut (the single spins under the same cut are exact, P = L_cut / (L_1 L_2)).
 *  - pair "channel", the pseudospin model: the alpha and beta product states of one index are taken as one state (their
 *    overlap, the pair's share of the ESEEM, left out), so that every channel c (two product states that differ in both
 *    spins, its element not zero in alpha or beta) is a two-level system with its Omega and Delta_alpha, Delta_beta:
 *    l_c, and l_c,diag without its element. Then L_ps = L_diag + (2/d) sum_c (l_c - l_c,diag), d = (2I_1+1)(2I_2+1),
 *    L_diag from the diagonals, P = L_ps / (L_1,ps L_2,ps) with L_k,ps spin k from its eigenvalues alone, and the
 *    cluster's value P L_1 L_2, L_k the exact single-spin coherence. Exact for two spins 1/2 with one basis in alpha and
 *    beta; it misses what couples the channels sharing a state (a chain).
 * The run's coherence is prod_k L_k prod_pairs P.
 */
#include "../../include/calculation/pseudospin.h"
#include "../../include/calculation/small_hamiltonian.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {
bool channelOn = false;      // between pseudospinChannelsBegin and pseudospinChannelsEnd
std::string channelText;     // this rank's lines: sent to rank 0 at the end, which alone writes the file
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

// the eigenstates of one spin's Hamiltonian H as columns, ordered by <I_z> from +I down, so that a diagonal H gives the |m>
// basis of the rows (m = I - row) as it is; a degenerate level is split by I_z, each column's largest entry real and positive
MatrixXcd nuclearBasis(const MatrixXcd& H){
    const int d = int(H.rows());
    Eigen::VectorXd Iz(d);
    for (int a=0; a<d; ++a) Iz(a) = 0.5*(d-1) - a;
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(H);
    MatrixXcd V = es.eigenvectors();
    const Eigen::VectorXd& E = es.eigenvalues();
    const double tol = 1e-9*std::max(1.0,E.cwiseAbs().maxCoeff());
    for (int i=0, j; i<d; i=j){
        for (j=i+1; j<d && E(j)-E(i) < tol; ++j) {}
        if (j-i > 1){
            const MatrixXcd B = V.middleCols(i,j-i);
            Eigen::SelfAdjointEigenSolver<MatrixXcd> split(B.adjoint()*Iz.asDiagonal()*B);
            V.middleCols(i,j-i) = B*split.eigenvectors();
        }
    }
    std::vector<double> mz(d);
    std::vector<int> order(d);
    for (int c=0; c<d; ++c){ mz[c] = (V.col(c).adjoint()*Iz.asDiagonal()*V.col(c)).real()(0,0); order[c] = c; }
    std::stable_sort(order.begin(),order.end(),[&](int a, int b){ return mz[a] > mz[b]; });
    MatrixXcd W(d,d);
    for (int c=0; c<d; ++c){
        W.col(c) = V.col(order[c]);
        int k; W.col(c).cwiseAbs().maxCoeff(&k);
        W.col(c) *= std::conj(W(k,c))/std::abs(W(k,c));
    }
    return W;
}

// the |m> of a basis column's largest weight, for the channel file
double label(const MatrixXcd& W, int c){ int k; W.col(c).cwiseAbs().maxCoeff(&k); return 0.5*(W.rows()-1) - k; }

// spin k of the pair alone: its eigenbases in the qubit's alpha and beta states; with pseudo, its full coherence and its
// coherence from its eigenvalues alone, the alpha and beta eigenstates of one index taken as one state (the channel's)
void bases(QubitArray* qa, BathArray* ba, int k, Config* cnf, MatrixXcd& Wa, MatrixXcd& Wb, Pulse* pls = nullptr,
           std::vector<doublec>* full = nullptr, std::vector<doublec>* pseudo = nullptr){
    BathArray* one = BathArray_init();
    BathArray_setNspin(one,1);
    BathArray_allocBath(one,QubitArray_getNqubit(qa));
    BathArray_setBath_i(one,BathArray_getBath_i(ba,k),0,QubitArray_getNqubit(qa));
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,one,cnf,Ha,Hb);
    BathArray_freeAll(one);
    Wa = nuclearBasis(Ha);
    Wb = nuclearBasis(Hb);
    if (pls){
        *full = evolve(Ha,Hb,cnf,pls);
        *pseudo = evolve(diagonalOf(Wa.adjoint()*Ha*Wa),diagonalOf(Wb.adjoint()*Hb*Wb),cnf,pls);
    }
}

// W1 (x) W2: spin 0 the slow index, row a d2 + b (Pauli_matrix_Z, addExpandedHamiltonian)
MatrixXcd product(const MatrixXcd& W1, const MatrixXcd& W2){
    const int d1 = int(W1.rows()), d2 = int(W2.rows());
    MatrixXcd V(d1*d2,d1*d2);
    for (int a=0; a<d1; ++a) for (int b=0; b<d1; ++b) V.block(a*d2,b*d2,d2,d2) = W1(a,b)*W2;
    return V;
}
}

void calCoherencePseudospinInto(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, MatrixXcd* result){
    MatrixXcd W1a, W1b, W2a, W2b;
    std::vector<doublec> full1, pseudo1, full2, pseudo2;
    bases(qa,ba,0,cnf,W1a,W1b,pls,&full1,&pseudo1);
    bases(qa,ba,1,cnf,W2a,W2b,pls,&full2,&pseudo2);
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,ba,cnf,Ha,Hb);
    const int d2 = int(W2a.rows()), d = int(Ha.rows());
    const MatrixXcd Va = product(W1a,W2a), Vb = product(W1b,W2b);
    const MatrixXcd Ta = Va.adjoint()*Ha*Va, Tb = Vb.adjoint()*Hb*Vb;
    const MatrixXcd Da = diagonalOf(Ta), Db = diagonalOf(Tb);
    // the channels: the elements between product states that differ in both spins and are not zero (to rounding)
    const double zero = 1e-12*std::max({1.0,Ta.cwiseAbs().maxCoeff(),Tb.cwiseAbs().maxCoeff()});
    std::vector<std::pair<int,int>> channels;
    for (int r=0; r<d; ++r) for (int s=r+1; s<d; ++s)
        if (r/d2 != s/d2 && r%d2 != s%d2 && (std::abs(Ta(r,s)) > zero || std::abs(Tb(r,s)) > zero)) channels.push_back({r,s});
    auto with = [](MatrixXcd C, const MatrixXcd& T, int r, int s){ C(r,s) = T(r,s); C(s,r) = T(s,r); return C; };
    auto lab = [](const MatrixXcd& V, const MatrixXcd& C){ return MatrixXcd(V*C*V.adjoint()); };   // back in the |m> basis

    if (cnf->pseudospin_block){
        MatrixXcd Ca = Da, Cb = Db;
        for (auto [r,s] : channels){ Ca = with(Ca,Ta,r,s); Cb = with(Cb,Tb,r,s); }
        const std::vector<doublec> pair = evolve(lab(Va,Ca),lab(Vb,Cb),cnf,pls);
        for (size_t i=0; i<pair.size(); ++i) result[i] = MatrixXcd::Constant(1,1,pair[i]);
    } else {
        // the pseudospin model: the alpha and beta product states of one index taken as one state, so every channel is a
        // two-level system of its Omega and Delta_alpha, Delta_beta; the pair from the diagonals and the channels' changes
        std::vector<doublec> pair = evolve(Da,Db,cnf,pls);
        for (auto [r,s] : channels){
            const int idx[2] = {r,s};
            MatrixXcd ca(2,2), cb(2,2);
            for (int i=0; i<2; ++i) for (int j=0; j<2; ++j){ ca(i,j) = Ta(idx[i],idx[j]); cb(i,j) = Tb(idx[i],idx[j]); }
            const std::vector<doublec> l = evolve(ca,cb,cnf,pls), l0 = evolve(diagonalOf(ca),diagonalOf(cb),cnf,pls);
            for (size_t i=0; i<pair.size(); ++i) pair[i] += 2.0/d*(l[i]-l0[i]);
        }
        for (size_t i=0; i<pair.size(); ++i)
            result[i] = MatrixXcd::Constant(1,1,pair[i]/(pseudo1[i]*pseudo2[i])*full1[i]*full2[i]);
    }

    if (channelOn){   // kHz: the Hamiltonians are in rad/ms
        const double* x1 = BathArray_getBath_i_xyz(ba,0);
        const double* x2 = BathArray_getBath_i_xyz(ba,1);
        char line[512];
        for (auto [r,s] : channels){
            snprintf(line,sizeof line,"%s %.4f %.4f %.4f %s %.4f %.4f %.4f %g %g %g %g %.16e %.16e %.16e %.16e\n",
                     BathArray_getBath_i_name(ba,0),x1[0],x1[1],x1[2],BathArray_getBath_i_name(ba,1),x2[0],x2[1],x2[2],
                     label(W1a,r/d2),label(W2a,r%d2),label(W1a,s/d2),label(W2a,s%d2),
                     2*std::abs(Ta(r,s))/(2*M_PI),2*std::abs(Tb(r,s))/(2*M_PI),
                     (Ta(r,r)-Ta(s,s)).real()/(2*M_PI),(Tb(r,r)-Tb(s,s)).real()/(2*M_PI));
            channelText += line;
        }
    }
}

void pseudospinChannelsBegin(Output* op){
    channelName = std::string(Output_getOutfile(op))+"_pseudospin_channels";
    channelText.clear();
    channelOn = true;
}

// Rank 0 writes the header, its own lines and those of the other ranks in rank order, sent by MPI (no rank file, so the
// ranks need no shared file system, and only rank 0 opens a file); in pieces of at most 1 GiB (an int count)
void pseudospinChannelsEnd(){
    channelOn = false;
    FILE* out = nullptr;
    if (rank == 0){
        out = std::fopen(channelName.c_str(),"w");
        if (!out){
            fprintf(stderr,"Error: cannot write %s\n",channelName.c_str());
            ccexFatal();
        }
        fprintf(out,"# name1 x1 y1 z1 name2 x2 y2 z2 (A) m1 m2 m1' m2': the channel |m1,m2> <-> |m1',m2'> between two product "
                    "states of the spins' eigenstates in each qubit state, labelled by the |m> of the largest weight in alpha; Omega_alpha, "
                    "Omega_beta (kHz, the splitting 2|<m1,m2|H|m1',m2'>|), Delta_alpha, Delta_beta (kHz, E(m1,m2) - E(m1',m2')) in "
                    "the qubit's alpha and beta states\n");
        std::fwrite(channelText.data(),1,channelText.size(),out);
    }
    constexpr long long piece = 1LL << 30;
    for (int r=1; r<nprocess; ++r){
        if (rank == r){
            long long size = (long long)channelText.size();
            MPI_Send(&size,1,MPI_LONG_LONG,0,0,MPI_COMM_WORLD);
            for (long long at=0; at<size; at+=piece)
                MPI_Send(channelText.data()+at,int(std::min(piece,size-at)),MPI_CHAR,0,1,MPI_COMM_WORLD);
        } else if (rank == 0){
            long long size = 0;
            MPI_Recv(&size,1,MPI_LONG_LONG,r,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            std::string text(size,'\0');
            for (long long at=0; at<size; at+=piece)
                MPI_Recv(&text[at],int(std::min(piece,size-at)),MPI_CHAR,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            std::fwrite(text.data(),1,text.size(),out);
        }
    }
    std::string().swap(channelText);
    if (out) std::fclose(out);
    MPI_Barrier(MPI_COMM_WORLD);
}
