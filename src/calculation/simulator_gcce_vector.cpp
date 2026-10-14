#include "../../include/calculation/simulator.h"

// A16: two pieces of the evolution=vector step (propagateStateGcce), byte-identical to what they replace.

// exp(-i H tau0) psi0 = M (phase .* (M^dag psi0)): M^dag psi0 is the same every step, so it is formed on the first call
// for the cluster (prepareGccePulses empties it) instead of every step - the same product on the same operands.
void firstFreeGcce(GcceWork* w, const Eigen::VectorXcd& psi0, const Eigen::VectorXcd& phase){
    if (w->startEigen.size() == 0) w->startEigen.noalias() = w->Madj * psi0;
    w->psiscratch = w->startEigen;
    w->psiscratch.array() *= phase.array();
    w->psi.noalias() = w->M * w->psiscratch;
}

// out = K in with K = UpulsesExpanded[ipulse] = kron(U, 1_bdim), block by block. Eigen's matrix-vector kernel forms each
// packet of the dense product as pmadd(c, alpha = 1, 0) with c = pmadd(K(:,j), x_j, ... pmadd(K(:,0), x_0, 0)); the
// entries of K off the blocks are zeros and add nothing. The same packet operations with the same entries of K in the same
// order here, at qdim^2 bdim instead of (qdim bdim)^2 products. Byte-identical only when bdim is a whole number of packets:
// otherwise the kernel takes some rows with smaller packets or scalar code, so that case (X9) is summed block by block.
void applyPulseGcce(GcceWork* w, int ipulse, const Eigen::VectorXcd& in, Eigen::VectorXcd& out){
    typedef Eigen::internal::packet_traits<doublec>::type Packet;
    const int size = Eigen::internal::packet_traits<doublec>::size;
    const MatrixXcd& K = w->UpulsesExpanded[ipulse];
    const int bdim = w->pulseBdim, qdim = bdim > 0 ? int(K.rows())/bdim : 0;
    if (bdim <= 0){
        out.noalias() = K * in;
        return;
    }
    if (bdim % size != 0){
        // X9: the same block structure, block by block in scalar order, at qdim^2 bdim products (not byte-identical to
        // the dense product: the kernel's packet/scalar rows round differently)
        out.setZero(in.size());
        for (int a=0; a<qdim; a++) for (int b=0; b<qdim; b++)
            out.segment(a*bdim,bdim).array() += K(a*bdim,b*bdim) * in.segment(b*bdim,bdim).array();
        return;
    }
    out.resize(in.size());
    const Packet zero = Eigen::internal::pset1<Packet>(doublec(0.0,0.0)), one = Eigen::internal::pset1<Packet>(doublec(1.0,0.0));
    for (int a=0; a<qdim; a++){
        for (int k=0; k<bdim; k+=size){
            Packet c = zero;
            for (int b=0; b<qdim; b++)
                c = Eigen::internal::pmadd(Eigen::internal::pset1<Packet>(K(a*bdim,b*bdim)),
                                           Eigen::internal::ploadu<Packet>(in.data()+b*bdim+k), c);
            Eigen::internal::pstoreu(out.data()+a*bdim+k, Eigen::internal::pmadd(c, one, zero));
        }
    }
}
