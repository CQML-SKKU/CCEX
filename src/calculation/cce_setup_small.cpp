// A42: the setup of a cCCE cluster of bath dimension 2 or 4 from the small-cluster Jacobi (small_eigensolver.h), for the
// scalar path (cceSetupSmallAlone, from calCoherenceCceInto) and the A30 window (cceSetupSmall with a group's results).
// Both go through cceSetupSmall here, so a cluster's M, eigenvalues and start vectors are the same bits on either path.
// A file of its own: cceSetup and the time evolution in simulator_cce.cpp compile as before, and with CCEX_SMALL_EIG=0 the
// run is the code before A42.
#include "../../include/calculation/simulator.h"
#include "../../include/calculation/small_eigensolver.h"

void cceSetupSmall(CceWorkspace& ws, const MatrixXcd* psi, const double* const* w, const doublec* const* V){
    const int n = int(ws.Halpha.rows());
    if (w[0]) ws.up.setGiven(n,w[0],V[0]); else ws.up.compute(ws.Halpha);
    if (w[1]) ws.down.setGiven(n,w[1],V[1]); else ws.down.compute(ws.Hbeta);
    // the rest as cceSetup (simulator_cce.cpp)
    const MatrixXcd& M_up   = ws.up.eigenvectors();
    const MatrixXcd& M_down = ws.down.eigenvectors();
    MatrixXcd& M = ws.M;
    M.noalias() = M_up.adjoint()*M_down;
    ws.valUp = ws.up.eigenvalues();
    ws.valDown = ws.down.eigenvalues();
    if (psi){
        const MatrixXcd& psi0 = *psi;
        ws.rowStart.noalias() = psi0.adjoint() * M_down;
        ws.colStart.noalias() = M_up.adjoint() * psi0;
        ws.Madj = M.adjoint();
    }
}

bool cceSetupSmallAlone(CceWorkspace& ws, const MatrixXcd* psi){
    const int n = int(ws.Halpha.rows());
    if ((n != 2 && n != 4) || !smalleig::on()) return false;
    double w[2][4]; doublec V[2][16];
    const double* wOk[2] = {nullptr,nullptr};
    const doublec* VOk[2] = {nullptr,nullptr};
    long long fallbacks = 0;
    for (int s=0; s<2; ++s){
        const doublec* h = (s == 0 ? ws.Halpha : ws.Hbeta).data();
        double* wp = w[s]; doublec* Vp = V[s]; bool ok = false;
        smalleig::solve(n,1,&h,&wp,&Vp,&ok);
        if (ok){ wOk[s] = w[s]; VOk[s] = V[s]; } else ++fallbacks;
    }
    smalleig::count(2-fallbacks,fallbacks);
    cceSetupSmall(ws,psi,wOk,VOk);
    return true;
}
