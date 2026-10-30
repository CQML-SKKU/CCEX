/**
 * @file utilities.cpp
 * @brief Shared helpers: error stops, logging, maths and string utilities.
 */
#include "../../include/common/utilities.h"
#include "../../include/common/memory.h"
#include <float.h>  // FLT_EPSILON
#include <cstdarg>
#include <cstdio>
#include "mpi.h"
#include <vector>
#include <algorithm>
#include <queue>
#include <map>
#include <string>
#include <cctype>
#ifdef CCEX_BLAS_DIRECT
#include "../../include/common/blas_option.h"   // the thresholds printBanner reports
#endif
#include "../../include/calculation/gpu/small_cluster_batch.h"   // smallClusterWarmDown (ccexRefuse)
#include "../../include/calculation/gpu/compute_backend.h"       // the GPU ranks (A98)

/* utils ---------------------------------------------------*/
MatrixXcd Double1dToMatrixXcd(double* val, int n){
    MatrixXcd mat = MatrixXcd::Zero(n,1);
    for (int i=0; i<n; i++){
        mat(i,0) = doublec(val[i],0.0);
    }
    return mat;
}

double* MatrixXcdToDouble1d(MatrixXcd mat){
    double* val = allocDouble1d(mat.rows());
    for (int i=0; i<mat.rows(); i++){
        val[i] = mat(i,0).real();
    }
    return val;
}

int Double_is_same_row(double* row1, double* row2) {
    for (int i = 0; i < 3; ++i) {
        if (fabs(row1[i] - row2[i]) > 1e-9) {
            return 0; // different
        }
    }
    return 1;  // same
}

/* math functions -----------------------------------------------*/

/////Calculate the parameter///////////////////////////////////
double dist(double spin1[],double spin2[]){
	double dist = 0;
	dist = sqrt(pow(spin1[0] - spin2[0],2) + pow(spin1[1] - spin2[1],2) + pow(spin1[2] - spin2[2],2));
	return dist;
}

/////polar angle from the vector_m to vector_n/////////////////
double cosTheta(double spin1[], double spin2[],double dist){ 
 //Note that spin_m[0] = gyromagnetic ratio, spin_m[0 ~ 2] --> x, y, z coordinates
    return (spin2[2] - spin1[2]) / dist;
}

double sinTheta(double spin1[], double spin2[], double dist){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    return tmp_mn / dist;
}

/////azimuthal angle from the vector_m to vector_n/////////////
double cosPhi(double spin1[], double spin2[]){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    if (tmp_mn == 0) return 1.0;
    else return (spin2[0] - spin1[0]) / tmp_mn;
}
double sinPhi(double spin1[], double spin2[]){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    if (tmp_mn == 0) return 0.0;
    else return (spin2[1] - spin1[1]) / tmp_mn;
}

/////bath-coordinate rotation//////////////////////////////////
// Tolerances. ROT_COLLINEAR_TOL is compared against |q x b| with q,b already unit,
// so it IS the sine of the angle between the two axes -- 1e-8 rejects anything
// closer than ~6e-7 degrees, where the azimuth of the new x axis is numerical noise.
// ROT_ORTHO_TOL is the post-construction self-check on R, several orders looser
// than the double-precision round-off the construction actually produces (~1e-16).
#define ROT_ZERONORM_TOL 1e-12
#define ROT_COLLINEAR_TOL 1e-8
#define ROT_ORTHO_TOL 1e-10

/**
 * @brief Normalize v into out; fatal if it is not a usable direction
 * @details Fatal if v is not three finite numbers, or has (numerically) zero length --
 *          both leave the direction undefined, and a silently wrong axis here would
 *          show up only as a wrong Hamiltonian much later.
 * @param[in]  v    the input 3-vector (a direction, any length)
 * @param[out] out  v divided by its length, a unit vector
 * @param[in]  what name of the input, used in the error message
*/
void normalizeAxisOrDie(const double v[3], double out[3], const char* what){

    for (int i=0; i<3; i++){
        if (!std::isfinite(v[i])){
            ccexRefusalPrint(stderr,"Error: %s[%d] is not a finite number (%g)\n",what,i,v[i]);
            ccexRefuse();
        }
    }

    double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (norm < ROT_ZERONORM_TOL){
        ccexRefusalPrint(stderr,"Error: %s has zero norm ( [%g, %g, %g] )\n",what,v[0],v[1],v[2]);
        ccexRefusalPrint(stderr,"An axis must be a non-zero direction vector.\n");
        ccexRefuse();
    }

    for (int i=0; i<3; i++){ out[i] = v[i]/norm; }
}

static void rotCross(const double u[3], const double v[3], double out[3]){
    out[0] = u[1]*v[2] - u[2]*v[1];
    out[1] = u[2]*v[0] - u[0]*v[2];
    out[2] = u[0]*v[1] - u[1]*v[0];
}

static double rotDot(const double u[3], const double v[3]){
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

/**
 * @brief Angle between two directions, in degrees, in [0,180]
 * @details atan2(|u x v|, u.v) rather than acos(u.v/|u||v|): acos loses almost all of
 *          its precision exactly where this is used, near 0 and 180 degrees, where its
 *          argument sits on a flat part of the cosine. atan2 stays accurate there.
 *          The result is signless but NOT folded to [0,90] -- an antiparallel pair is
 *          180 degrees apart, which is what makes a same-direction test possible.
 *          Scale-independent, so the inputs need not be normalized.
*/
double angleBetweenDeg(const double u[3], const double v[3]){

    double c[3];
    rotCross(u,v,c);
    return atan2(sqrt(rotDot(c,c)),rotDot(u,v)) * 180.0 / M_PI;
}

/**
 * @brief Build the passive transform that puts the qubit axis on the computational +z
 * @details The constructed R satisfies
 *
 *      R * normalize(qubit_axis) = [0,0,1]
 *
 *  i.e. it re-expresses a vector given in the bath file's own Cartesian frame in the
 *  frame whose z axis IS the qubit symmetry axis (a passive change of basis, not an
 *  active rotation of the system).
 *
 *      ez = q                  the new z axis is the qubit axis
 *      ex = normalize(q x z)   z = [0,0,1] of the bath file; it only fixes the AZIMUTH
 *                              of x and y about ez, which no result depends on
 *      ey = ez x ex            right-handed
 *
 *  (Until A90 the z of the second line was the input "bath_axis", default [0,0,1]; it
 *  was removed because it changed nothing but that azimuth.)
 *
 *  ex, ey, ez are the ROWS of R, so (R*v)[k] is the component of v along the k-th new
 *  axis. For qubit_axis=[1,1,1] this gives
 *
 *      ex = [ 1,-1, 0]/sqrt(2),  ey = [ 1, 1,-2]/sqrt(6),  ez = [ 1, 1, 1]/sqrt(3)
 *
 *  qubit_axis = +z gives the identity; qubit_axis = -z leaves x undefined and is refused.
 *
 * @param qubit_axis qubit symmetry axis, in the bath file frame
 * @param R          [out] 3x3 rotation matrix
 */
void buildQubitAlignedRotation(const double qubit_axis[3], double R[3][3]){

    const double b[3] = {0.0,0.0,1.0};   // the azimuth reference: the bath file's z
    double q[3];
    normalizeAxisOrDie(qubit_axis,q,"qubit_axis");

    double ex[3], ey[3], ez[3];
    ez[0] = q[0]; ez[1] = q[1]; ez[2] = q[2];

    rotCross(q,b,ex);
    double cross_norm = sqrt(rotDot(ex,ex));

    if (cross_norm < ROT_COLLINEAR_TOL){
        // q x z vanishes, so the azimuth of the new x axis is undefined. The one case
        // that still has an answer is "the qubit axis is already +z" -- then the frame
        // does not move at all and R is the identity.
        bool q_is_z = (fabs(q[0]) < ROT_COLLINEAR_TOL) && (fabs(q[1]) < ROT_COLLINEAR_TOL) && (q[2] > 0.0);

        if (q_is_z){
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){ R[i][j] = (i==j) ? 1.0 : 0.0; }
            }
            return;
        }

        ccexRefusalPrint(stderr,"Error(frame.rotation): qubit_axis is antiparallel to z ([0,0,-1])\n");
        ccexRefusalPrint(stderr,"  normalized qubit_axis : [ %g, %g, %g ]\n",q[0],q[1],q[2]);
        ccexRefusalPrint(stderr,"The computational x axis is qubit_axis x [0,0,1], which vanishes here, so the frame is undefined.\n");
        ccexRefuse();
    }

    for (int i=0; i<3; i++){ ex[i] /= cross_norm; }
    rotCross(ez,ex,ey);

    for (int j=0; j<3; j++){
        R[0][j] = ex[j];
        R[1][j] = ey[j];
        R[2][j] = ez[j];
    }

    ////////////////////////////////////////////////////////////////////////
    // Self-check : R must be a proper rotation that lands q on +z
    ////////////////////////////////////////////////////////////////////////
    double ortho_err = 0.0;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            double rtr = R[0][i]*R[0][j] + R[1][i]*R[1][j] + R[2][i]*R[2][j]; // (R^T R)_ij
            double target = (i==j) ? 1.0 : 0.0;
            if (fabs(rtr-target) > ortho_err){ ortho_err = fabs(rtr-target); }
        }
    }

    double det = R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1])
               - R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0])
               + R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);

    double zcheck_err = 0.0;
    for (int i=0; i<3; i++){
        double target = (i==2) ? 1.0 : 0.0;
        double val = R[i][0]*q[0] + R[i][1]*q[1] + R[i][2]*q[2];
        if (fabs(val-target) > zcheck_err){ zcheck_err = fabs(val-target); }
    }

    if (ortho_err > ROT_ORTHO_TOL || fabs(det-1.0) > ROT_ORTHO_TOL || zcheck_err > ROT_ORTHO_TOL){
        ccexRefusalPrint(stderr,"Error(frame.rotation): the constructed matrix is not a valid rotation\n");
        ccexRefusalPrint(stderr,"  max|R^T R - I|            = %g (tol %g)\n",ortho_err,ROT_ORTHO_TOL);
        ccexRefusalPrint(stderr,"  det(R)                    = %.15g (must be +1)\n",det);
        ccexRefusalPrint(stderr,"  max|R*qubit_axis - [001]| = %g (tol %g)\n",zcheck_err,ROT_ORTHO_TOL);
        ccexRefuse();
    }
}

/**
 * @brief Apply the rotation about a fixed point : r_new = r0 + R * (r_old - r0)
 * @param xyz [in,out] coordinate, overwritten in place
 * @param R   rotation matrix from buildQubitAlignedRotation
 * @param r0  point held fixed by the transformation (the central qubit)
 */
void rotateAboutPoint(double xyz[3], const double R[3][3], const double r0[3]){

    double d[3];
    for (int i=0; i<3; i++){ d[i] = xyz[i] - r0[i]; }

    for (int i=0; i<3; i++){
        xyz[i] = r0[i] + R[i][0]*d[0] + R[i][1]*d[1] + R[i][2]*d[2];
    }
}

/**
 * @brief Re-express a rank-2 tensor in the rotated frame : R * T * R^T
 * @details The companion of rotateAboutPoint. A tensor whose COMPONENTS are written in
 *          the bath file's Cartesian basis has to be transformed as well, or it would
 *          keep describing the old axes while the positions describe the new ones.
 *
 *          R is real, so the transpose -- not the adjoint -- is the right inverse here;
 *          they coincide numerically, but only transpose is correct in intent.
 *
 *          A tensor that was never filled in is returned untouched: BathArray leaves an
 *          unset hyperfine tensor at 0x0, and there is nothing to rotate.
 * @param T tensor in the original frame
 * @param R rotation matrix from buildQubitAlignedRotation
*/
MatrixXcd rotateTensor(const MatrixXcd T, const double R[3][3]){

    if (T.rows() != 3 || T.cols() != 3){ return T; }

    MatrixXcd Rm = MatrixXcd::Zero(3,3);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){ Rm(i,j) = doublec(R[i][j],0.0); }
    }

    return Rm * T * Rm.transpose();
}

/**
 * @brief Print the whole coordinate frame once, for the run log
 * @details Called on rank 0 only, once per run -- never per bath spin. The point is that
 *          a reader can reconstruct every frame the run used from this block alone: the
 *          two input axes, where each computational axis lies in the source frame, R and
 *          its self-checks, and the magnetic field written both ways round.
 * @param qubit_axis             qubit symmetry axis, in the bath file frame
 * @param R                      the rotation, original (bath file) frame -> computational frame
 * @param r0                     the point held fixed (the first qubit), in the original frame
 * @param bfield                 Config.bfield, i.e. the field in the COMPUTATIONAL frame
 * @param defect_axis_reference  name of the recorded reference axis, or NULL if none
*/
void reportQubitAlignedRotation(const double qubit_axis[3], const double R[3][3], const double r0[3],
                                const float bfield[3], const char* defect_axis_reference){

    double q[3];
    normalizeAxisOrDie(qubit_axis,q,"qubit_axis");

    double ortho_err = 0.0;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            double rtr = R[0][i]*R[0][j] + R[1][i]*R[1][j] + R[2][i]*R[2][j];
            double target = (i==j) ? 1.0 : 0.0;
            if (fabs(rtr-target) > ortho_err){ ortho_err = fabs(rtr-target); }
        }
    }

    double det = R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1])
               - R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0])
               + R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);

    double qz[3];
    for (int i=0; i<3; i++){ qz[i] = R[i][0]*q[0] + R[i][1]*q[1] + R[i][2]*q[2]; }

    // R^T takes a computational-frame vector back to the source frame. For a field along
    // the computational z it returns Bz * ez, i.e. Bz * normalize(qubit_axis).
    double bsrc[3];
    for (int i=0; i<3; i++){
        bsrc[i] = R[0][i]*(double)bfield[0] + R[1][i]*(double)bfield[1] + R[2][i]*(double)bfield[2];
    }

    printSubTitle("Coordinate frame rotation ...");
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","qubit axis, source frame",q[0],q[1],q[2]);
    printf("\n");
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","computational x axis, source frame",R[0][0],R[0][1],R[0][2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","computational y axis, source frame",R[1][0],R[1][1],R[1][2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","computational z axis, source frame",R[2][0],R[2][1],R[2][2]);
    printf("\n");
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 0",R[0][0],R[0][1],R[0][2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 1",R[1][0],R[1][1],R[1][2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 2",R[2][0],R[2][1],R[2][2]);
    printf("      %-32s : %.15lf\n","det(R)",det);
    printf("      %-32s : %.3e\n","max|R^T R - I|",ortho_err);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","R * normalize(qubit_axis)",qz[0],qz[1],qz[2]);
    printf("\n");
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation origin r0, source frame",r0[0],r0[1],r0[2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","B field, computational frame",
           (double)bfield[0],(double)bfield[1],(double)bfield[2]);
    printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","B field, source frame (R^T B)",bsrc[0],bsrc[1],bsrc[2]);
    printf("\n");

    if (defect_axis_reference != NULL){
        printf("      %-32s : %s\n","Defect axis reference",defect_axis_reference);
        printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","Reference axis, source frame",q[0],q[1],q[2]);
        printf("      %-32s : [ %12.8lf, %12.8lf, %12.8lf ]\n","Reference axis, calc frame",qz[0],qz[1],qz[2]);
        printf("      %-32s : recorded for future JT-axis classification,\n","Status");
        printf("      %-32s  not applied to the Defect Hamiltonian\n","");
        printf("\n");
    }

    printf("      Every qubit position and every bath-spin position gets the same rigid\n");
    printf("      transform r_new = r0 + R*(r_old - r0), about the first qubit.\n");
    printf("      Tensors with frame \"file\" (hyperfine, quadrupole, Qubit.intmap,\n");
    printf("      Defect species) and the point-dipole tensors move with them; frame \"qubit\"\n");
    printf("      stays as written -- the counts are listed below.\n");
    printf("      Read as written in the computational frame: \"bfield\" ([0,0,B] is along qubit_axis).\n");
    printf("\n");
}

///////////////////////////////////////////////////////////////
// Physics functions
///////////////////////////////////////////////////////////////

// Obtain possible substates of the spin
float* substates(float S){
    float* ms=NULL;
    int n = (int)(2*S+1);
    ms = allocFloat1d(n);
    for (int i=0; i<n; i++){
        ms[i] = S-i;
    }
    return ms;
}

bool isSubLevel(float S, float ms){

    for (int i=0; i<(int)(2*S+1); i++){
        if (ms == S-i){
            return true;
        }
    }
    return false;
}

// Obtain the spinor array
MatrixXcd getSpinor(float S, float ms){

    if (ms > S || ms < -S){
        fprintf(stderr, "Error(spinorArray): ms %g is not satisfying -S <= ms <= S (S %g)\n", ms, S);
        ccexFatal();   // a rank alone may reach this: every rank stops (D44)
    }

    int n = (int)(2*S+1);
    MatrixXcd spinor = MatrixXcd::Zero(n,1);
    //e.g. S=1
    //    i = 0, 1, 2 (ms = S-i)
    //    i = 0 -> ms = 1 -> arr = [ 1 0 0] (else = 0)
    //    i = 1 -> ms = 0 -> arr = [ 0 1 0] (else = 0)
    //    i = 2 -> ms = -1 -> arr = [ 0 0 1] (else = 0)
    bool chkms = false;
    for (int i=0; i<n; i++){
        float mstmp = S-i;
        if (mstmp == ms){
            spinor(i,0) = doublec(1.0,0.0);
            chkms = true;
        }
        else{
            spinor(i,0) = doublec(0.0,0.0);
        }
    }

    if (chkms){
        // std::cout << "spinor" << std::endl;
        // std::cout << spinor << std::endl;
        return spinor;
    }
    else{
        fprintf(stderr, "Error(spinorArray): ms %g is not a sublevel of S %g\n", ms, S);
        ccexFatal();
    }
}


MatrixXcd kron(const MatrixXcd& A, const MatrixXcd& B){
    MatrixXcd C(A.rows()*B.rows(),A.cols()*B.cols());
    for (int i=0; i<A.rows(); i++){
        for (int j=0; j<A.cols(); j++){
            C.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j)*B;
        }
    }
    return C;
}


/*! 
 * @brief Calculate the norm of the matrix
 * @details Spinor Properties : normalized set
 * @param[in] m gsl_matrix_complex
 * @return norm of the matrix
 */
double calNorm(MatrixXcd m){
    int row = m.rows();
    int col = m.cols();

    // Compare only vector
    if (col!=1){
        perror("Error(norm): the matrix is not a vector");
        exit(EXIT_FAILURE);
    }

    // Find norm^2
    doublec norm2 = (m.adjoint() * m)(0,0);
    
    // check if the norm is real number
    if (fabs(norm2.imag()) > FLT_EPSILON){
        perror("Error(norm): the norm is not a real number");
        exit(EXIT_FAILURE);
    }

    return sqrt(norm2.real());
}


int normalize(MatrixXcd* m){

    double norm = calNorm(*m);

    if (fabsf((float)norm -1.0f) <= FLT_EPSILON){
        return 0; // already normalized
    }
    else{
        *m = *m / norm;
        return 1; // give normalized matrix
    }

}

float findZbasisSubLevel(MatrixXcd spinor){
    int n = spinor.rows();
    float S = (float)(n-1)/2.0;
    float ms = 0.0;
    
    bool chkms = false;
    for (int i=0; i<n; i++){
        double spinoriReal = spinor(i,0).real();
        if (fabs(spinoriReal - 1.0) <= FLT_EPSILON){
            ms = S-i;
            if (!chkms){
                chkms = true;
            }
            else{
                fprintf(stderr, "Error(findZbasisSubLevel): ms is not unique (S %g)\n", S);
                ccexFatal();   // reached by the rank that owns a cluster: every rank stops (review 20261008-autoconv-d42 #6)
            }
        }
    }

    if (chkms){
        return ms;
    }
    else{
        fprintf(stderr, "Error(findMs): the state is not a z-basis sublevel of S %g\n", S);
        ccexFatal();
    }
}

MatrixXcd partialtrace(MatrixXcd Mij, int dimrow, int dimcol){

    // Mij : Full Matrix
    // Bkl : Block Matrix
    // Rkl : Reduced Matrix

    //        | B00 B01 B02 | 
    // Mij =  | B10 B11 B12 | 
    //        | B20 B21 B22 | 
    MatrixXcd Bkl = MatrixXcd::Zero(dimrow,dimcol);

    //        | r00 r01 r02 |
    // Rkl =  | r10 r11 r12 | = PartialTrace[Mij]
    //        | r20 r21 r22 |
    MatrixXcd Rkl = MatrixXcd::Zero(Mij.rows()/dimrow,Mij.cols()/dimcol);
    
    // Get partial trace of M after finding block matrix
    int irow = 0;
    int icol = 0;
    for (int k=0; k<Mij.rows(); k+=dimrow){
        icol=0;
        for (int l=0; l<Mij.cols(); l+=dimcol){
            // Find block matrix Bkl
            // Bkl,pq = Mij 
            // * i = (a)*dimrow + p < (a+1) * dimrow
            // * j = (b)*dimcol + q < (b+1) * dimcol
            Bkl = Mij.block(k,l,dimrow,dimcol);
            // Trace of Bkl 
            Rkl(irow,icol) = Bkl.trace();
            icol++;
        }
        irow++;
    }
    return Rkl;
}

MatrixXcd powMatrixXcdElementWise(MatrixXcd a, int n){
    // NOTE!!! 
    // here the input power is only allowing "int"
    // Because i found float n gives wrong value in this case.
    // have to use int n
    // but n<0 is okay to use
    int nrow = a.rows();
    int ncol = a.cols();
    MatrixXcd res(nrow,ncol);
 
    for (int i=0; i<nrow; i++){
        for (int j=0; j<ncol; j++){
            res(i,j) = std::pow(a(i,j),n);
        }
    }
    return res;
}

MatrixXcd mulMatrixXcdElementWise(MatrixXcd a, MatrixXcd b){
    return a.array() * b.array();
}

// a / b element by element, for the cluster correction of a density matrix:
// rho_C / rho_0. The elementwise product of gCCE assumes every element of the
// free evolution rho_0 is nonzero; an element it keeps at zero (a level the
// initial state does not touch) has no correction to apply, and the final
// product leaves it at zero. So where rho_0 is zero -- measured against the
// largest element of rho_0 at that time step, since it carries round-off --
// the correction is 1. A cluster may put a real value into such an element
// (the bath couples the level); the expansion cannot carry it, and the
// largest such value is reported through *dropped so that the run says so.
// Before 2026-09-24 these elements were round-off over zero (inf, then nan
// once an inf met a zero) or round-off over round-off (a large finite number),
// and every product they entered came out nan.
MatrixXcd divMatrixXcdElementWise(MatrixXcd a, MatrixXcd b, double* dropped){
    MatrixXcd res(a.rows(), a.cols());
    const double floor = 1e-12 * b.cwiseAbs().maxCoeff();
    for (int i=0; i<a.rows(); i++){
        for (int j=0; j<a.cols(); j++){
            if (std::abs(b(i,j)) <= floor){
                res(i,j) = doublec(1.0,0.0);
                if (dropped != nullptr && std::abs(a(i,j)) > *dropped) *dropped = std::abs(a(i,j));
            }else{
                res(i,j) = a(i,j) / b(i,j);
            }
        }
    }
    return res;
}

/* Easy print ---------------------------------------------------*/
void printInlineMatrixXcd(char* key, MatrixXcd mat){
    printf("      %-18s : [ ", key);

    if (mat.rows() * mat.cols() > 9 && mat.rows() * mat.cols() <= 400){
        int count = 0;
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                std::complex<double> z = mat(i, j);
                printf("%3.2f%+-3.2fj", z.real(), z.imag());
                if (i != mat.rows() - 1 || j != mat.cols() - 1) {
                    printf(", ");
                }
                count++;
            }
            if (count%9 == 0 && count != (mat.rows() * mat.cols())){
                printf("\n%30s"," ");
            }          
        }
    }else if (mat.rows() * mat.cols() <= 9){
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                std::complex<double> z = mat(i, j);
                printf("%3.2f%+-3.2fj", z.real(), z.imag());
                if (i != mat.rows() - 1 || j != mat.cols() - 1) {
                    printf(", ");
                }
            }
        }
    }else{
        // print the non-zero elements in the matrix
        // and print the rest as 0
        for (int i = 0; i < mat.rows(); ++i) {
            int count = 0;
            for (int j = 0; j < mat.cols(); ++j) {
                std::complex<double> z = mat(i, j);
                if (fabs(z.real()) > FLT_EPSILON || fabs(z.imag()) > FLT_EPSILON){
                    printf("(%d,%d) %3.2f%+-3.2fj", i,j,z.real(), z.imag());
                    if (i != mat.rows() - 1 || j != mat.cols() - 1) {
                        printf(", ");
                    }
                    count++;

                    if (count%9 == 0 && (i == mat.rows() - 1 && j == mat.cols() - 1)){
                        printf("\n%30s"," ");
                    }
                }
            }
        }
    }

    printf(" ]\n");
}

void printStateInDiracNot(char* key, MatrixXcd mat){
    printf("      %-18s : ", key);

    int dim = mat.rows();
    float spin = (dim-1.0)/2.0;

    for (int i=0; i <dim; i++){
        double real = mat(i,0).real();
        double imag = mat(i,0).imag();
        double abs_value = sqrt(real*real + imag*imag);

        if ( abs_value > FLT_EPSILON ){
            float ms = spin - i;
            printf("(%+g%+gj) |%2g >  ",real,imag,ms);
        }
    }
    printf("\n");
}

void printStructElementChar(char* key,char* val){
    printf("      %-18s : %s\n", key, val);
}

void printStructElementChar2d(char* key, char** val, int n){
    printf("      %-18s : [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10s", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
    if (n==0){printf("]\n");}
}

void printStructElementInt(char* key, int val){
    printf("      %-18s : %d\n", key, val);
}

void printStructElementInt1dIdx(char* key, int* val, int n){
    printf("      %-18s : [ ", key);
    for (int i = 0; i < n; i++){
        printf("%3d : %-5d",i,val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
    if (n==0){printf("]\n");}
}

void printStructElementFloat(char* key, float val){
    printf("      %-18s : %.6g\n", key, val);
}

void printStructElementFloat1d(char* key, float* val, int n){
    printf("      %-18s : [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10.2f", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
    if (n==0){printf("]\n");}
}

void printStructElementDouble(char* key, double val){
    printf("      %-18s : %.3g\n", key, val);
}

void printStructElementDouble1d(char* key, double* val, int n){
    printf("      %-18s : [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10.2f", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
    if (n==0){printf("]\n");}
}

void printStructElementBool(char* key, bool val){
    printf("      %-18s : %s\n", key, val ? "true" : "false");
}

void printLine(){
    printf("    ----------------------------------------------------------------------\n");
}

void printLineSection(){
    printf("\n    ======================================================================\n\n");
}

void printTitle(char* title){
    printf("    < %s >\n\n",title);
}

void printSubTitle(char* title){
    printf("    >> %s\n\n",title);
}

void copyBounded(char* dst, const char* src, size_t size, const char* what){
    if (strlen(src) >= size){
        ccexRefusalPrint(stderr,"Error: %s \"%s\" is longer than it can be (%zu characters)\n",what,src,size-1);
        ccexRefuse();
    }
    strcpy(dst,src);
}

void printMessage(char* title){
    printf("        %s\n",title);
}

/* Find index ---------------------------------------------------*/
int findIndexInt(int* array, int ista, int iend, int val){
    for (int i = ista; i <= iend; i++){
        if (array[i] == val){
            return i;
        }
    }
    return -1;
}

int findIndexCharFix(char array[][MAX_CHARARRAY_LENGTH], int ista, int iend, char* val){
    for (int i = ista; i <= iend; i++){
        if (strcasecmp(array[i],val) == 0){
            return i;
        }
    }
    return -1;
}

int findIndexChar(char** array, int ista, int iend, char* val){
    for (int i = ista; i <= iend; i++){
        if (strcasecmp(array[i],val) == 0){
            return i;
        }
    }
    return -1;
}

int findIndexFloat(float* array, int ista, int iend, float val){
    for (int i = ista; i <= iend; i++){
        if (array[i] == val){
            return i;
        }
    }
    return -1;
}

/* Quick sort ---------------------------------------------------*/

void QuickSort(int** d_Array, int left, int right){
    if(left<=right){
        int pivot = Partition(d_Array, left, right);
        QuickSort(d_Array, left, pivot-1);
        QuickSort(d_Array, pivot+1, right);
    }
}

int Partition(int** d_Array, int left, int right){
    int pivot = (*d_Array)[left];
    int low = left +1;
    int high = right;
    while(low<=high){
        while(low<=right && pivot >= (*d_Array)[low]){
            low++;
        }
        while(high>=(left+1) && pivot <= (*d_Array)[high]){
            high--;
        }
        if (low <= high){
            Swap(d_Array, low, high);
        }
    }
    Swap(d_Array, left, high);
    return high;
}

void Swap(int** d_Array, int a, int b){
    int temp = (*d_Array)[a];
    (*d_Array)[a] = (*d_Array)[b];
    (*d_Array)[b] = temp;
}

/* Type checker ---------------------------------------------------*/

//Check the string is the double
int isStringDouble(char *s){
    //0--> no number, 1--> is number
    size_t size = strlen(s);
    if (size == 0){return 0;} //no number
    if (size == 1){if(!isdigit(s[0])){return 0;}}

    bool CheckConti=false;
    for (int i=0;i<(int) size; i++){
        if(s[i] == '.' || s[i] =='-' || s[i] == '+' || isspace(s[i])){
            if (CheckConti==true){
                return 0; //if special string is continue, the string is not number
            }else{
                CheckConti=true;continue;
            }
        } // ok 
        if(s[i] < '0' || s[i] > '9') {return 0;} //not number
        CheckConti=false; //for check that continue special string
    }

    return 1; //this string is number
}

/* MPI ---------------------------------------------------*/
void para_range(int n1,int n2, int nprocs, int myrank, int*ista, int *iend)
{
    int iwork1, iwork2;
    iwork1 = (n2-n1+1)/nprocs;
    iwork2 = (n2-n1+1) % nprocs;
    *ista = myrank*iwork1 + n1 + min(myrank, iwork2);
    *iend = *ista + iwork1 - 1;
    if(iwork2 > myrank) *iend = *iend + 1;
}

int min(int x, int y)
{
    int v;
    if (x>=y) v = y;
    else v = x;
    return v;
}


// A31 + A45: the rows (1-based, ascending) of an order of m clusters that this rank takes. loads (nprocess entries) is
// every rank's cost so far - carried from order to order (MPI_getLocalClusters goes from the highest order down), so
// the orders together are split as one list, largest first. With a cost that is not the same for every cluster of
// the order: every cluster, the largest cost first (equal costs: the earlier in the list first), goes to the rank with
// the least load (equal: the lower rank) - longest first, so a rank holds a mix of sizes and the errors of the cost
// formula average out (A44 step 0: the contiguous A31 split left 20-27 % end idle on h-BN). Every rank computes the
// same doubles in the same order, so all agree. All costs equal, no cost, or no cluster: the equal-count contiguous
// split (para_range), so a bath of one kind of spin is distributed as before (byte-identical); loads then grow by
// each rank's count times that cost.
// A98: the split with GPU ranks (gpu_options.ranks_per_gpu). The dimensions the GPU takes are timed first, on a
// made-up cluster each (computeGpuSplitCalibrate, collective: every rank sees the same rows). Then, largest cost first, a
// cluster the GPU may take goes to the least-loaded GPU rank at its cost times the weight (computeGpuSplitWeight) or to
// the least-loaded CPU rank at its cost, whichever ends first; any other goes to the least-loaded rank of all (the GPU
// ranks too), at its cost. Equal loads: the GPU rank, then the lower rank.
static std::vector<int> rankClusterRowsGpu(int m, const std::function<double(int)>& rowCost, std::vector<double>& loads,
                                           int n, const std::function<const int*(int)>& rowOf){
    std::set<int> pending;
    for (int i=0; i<m; i++)
        if (const int d = computeGpuSplitPending(n,rowOf(i))) pending.insert(d);
    computeGpuSplitCalibrate(pending);
    // shortcut: a cost and an index per row (12 B; clusterizePairs keeps 12 B a pair), for lists far above 1e8 rows a pass per cost group
    std::vector<double> cost(m);
    for (int i=0; i<m; i++) cost[i] = rowCost(i);
    std::vector<int> order(m);
    for (int i=0; i<m; i++) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b){ return cost[a] > cost[b]; });
    const std::vector<int>& gpuRanks = computeGpuSplitRanks();
    // the least-loaded rank without a GPU (a heap: only these placements change their loads); the GPU ranks are few, scanned
    std::priority_queue<std::pair<double,int>, std::vector<std::pair<double,int>>, std::greater<std::pair<double,int>>> cpuLeast;
    std::vector<char> gpu(nprocess,0);
    for (int g : gpuRanks) gpu[g] = 1;
    for (int r=0; r<nprocess; r++) if (!gpu[r]) cpuLeast.push({loads[r], r});
    std::vector<int> rows;
    for (int i : order){
        const double w = computeGpuSplitWeight(n,rowOf(i));
        int g = gpuRanks[0];
        for (int x : gpuRanks) if (loads[x] < loads[g]) g = x;
        const int c = cpuLeast.empty() ? -1 : cpuLeast.top().second;
        // a cluster the GPU may take: on a GPU rank's card or on a CPU rank, whichever finishes it first; any other: the
        // CPU of the least-loaded rank, a GPU rank's too
        const double onG = loads[g] + cost[i]*(w > 0 ? w : 1), onC = c < 0 ? 0 : loads[c] + cost[i];
        const bool toG = c < 0 || (w > 0 ? onG <= onC : (loads[g] < loads[c] || (loads[g] == loads[c] && g < c)));
        const int r = toG ? g : c;
        if (toG) loads[g] = onG;
        else { cpuLeast.pop(); loads[c] = onC; cpuLeast.push({onC, c}); }
        if (r == rank) rows.push_back(i+1);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

std::vector<int> rankClusterRows(int m, const std::function<double(int)>& rowCost, std::vector<double>& loads,
                                 int n, const std::function<const int*(int)>& rowOf){
    std::vector<int> rows;
    loads.resize(nprocess, 0.0);
    if (rowCost && rowOf && m > 0 && computeGpuSplitOn()) return rankClusterRowsGpu(m,rowCost,loads,n,rowOf);
    bool uniform = true;
    double c0 = 0.0;
    if (rowCost && m > 0){
        // A38: no array for the test, so a uniform list (P1, 13C) costs nothing more at the peak of clusterizePairs
        c0 = rowCost(0);
        for (int i=1; i<m && uniform; i++) uniform = rowCost(i) == c0;
    }
    if (uniform){
        for (int r=0; r<nprocess; r++){
            int ista, iend;
            para_range(2, m+1, nprocess, r, &(ista) ,&(iend));
            // cluster : ista - 1 <= i < iend; none when ista - 1 >= iend
            if (r == rank) for (int i=ista-1; i<iend; i++) rows.push_back(i);
            if (iend > ista-1) loads[r] += c0*(iend-ista+1);
        }
        return rows;
    }
    std::priority_queue<std::pair<double,int>, std::vector<std::pair<double,int>>, std::greater<std::pair<double,int>>> least;
    for (int r=0; r<nprocess; r++) least.push({loads[r], r});
    auto place = [&](int i, double c){
        std::pair<double,int> top = least.top();
        least.pop();
        if (top.second == rank) rows.push_back(i+1);
        top.first += c;
        least.push(top);
    };
    // X10: few distinct costs (the cost is D^2 (D+nstep), so a mixed bath has a few): the rows grouped by cost, the
    // largest first, each group in list order - the order the stable sort below gives, so the same rows and the same
    // additions to the loads, without sorting the m rows
    constexpr size_t maxGroups = 64;
    std::map<double,std::vector<int>,std::greater<double>> groups;
    bool few = true;
    for (int i=0; i<m && few; i++){
        groups[rowCost(i)].push_back(i);
        few = groups.size() <= maxGroups;
    }
    if (few){
        for (const auto& g : groups) for (int i : g.second) place(i, g.first);
    } else {
        groups.clear();
        std::vector<double> cost(m);
        for (int i=0; i<m; i++) cost[i] = rowCost(i);
        std::vector<int> order(m);
        for (int i=0; i<m; i++) order[i] = i;
        std::stable_sort(order.begin(), order.end(), [&](int a, int b){ return cost[a] > cost[b]; });
        for (int i : order) place(i, cost[i]);
    }
    while (!least.empty()){ loads[least.top().second] = least.top().first; least.pop(); }
    std::sort(rows.begin(), rows.end());
    return rows;
}

int*** MPI_getLocalClusters(int order, int*** clusters, const std::function<double(int,const int*)>& cost,
                            int givenOrder, const std::vector<int>* givenRows,
                            const std::vector<double>* givenLoads){

    MPI_Barrier(MPI_COMM_WORLD);
    // this rank's rows of each order (rankClusterRows: the cost split, or the equal-count one)
    // A45: from the highest order down, the ranks' loads carried, so the large clusters are placed first
    std::vector<std::vector<int>> rows(order+1);
    std::vector<double> loads(nprocess, 0.0);
    for (int n=order; n>=1; n--){
        if (n == givenOrder){                       // A38: decided before the rows were dropped
            rows[n] = *givenRows;
            if (givenLoads) loads = *givenLoads;
            continue;
        }
        std::function<double(int)> rowCost;
        if (cost) rowCost = [&](int i){ return cost(n, clusters[n][i+1]); };
        rows[n] = rankClusterRows(clusters[n][0][0]-1, rowCost, loads, n, [&](int i){ return (const int*)clusters[n][i+1]; });
    }

    MPI_Barrier(MPI_COMM_WORLD);


    // make local clusters for each rank
    int*** localClusters = (int***)calloc(order+1,sizeof(int**));

    // zeroth cluster
    localClusters[0] = (int**)calloc(1,sizeof(int*));
    localClusters[0][0] = (int*)calloc(1,sizeof(int));
    localClusters[0][0][0] = clusters[0][0][0]; // = 1
    if (rank!=0){
        localClusters[0][0][0] = 0;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // n > 0 clusters
    for (int n=1; n<order+1; n++){

        int size = (int)rows[n].size() + 1;
        localClusters[n] = (int**)allocArray2d(size,n+1,sizeof(int));
        localClusters[n][0][0] = size;

        for (int i=1; i<size; i++){
            const int iroot = rows[n][i-1];
            if (clusters[n][iroot] == NULL){   // A38: a pair row this rank did not keep
                printf("rank[%d] : order %d, row %d was not kept\n",rank,n,iroot);
                perror("local cluster row missing");
                exit(1);
            }
            for (int j=0; j<n+1; j++){
                localClusters[n][i][j] = clusters[n][iroot][j];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return localClusters;
}

void MPI_freeLocalClusters(int order, int**** localClusters){
    if (localClusters == nullptr || *localClusters == nullptr) return;
    for (int n=0; n<=order; ++n){
        int** rows=(*localClusters)[n];
        if (rows != nullptr){
            const int count=n==0 ? 1 : rows[0][0];
            freeInt2d(&(*localClusters)[n],count);
        }
    }
    freeArray1d((void**)localClusters);
}

// MatrixXcd* MPI_reduceLocalResult(int nstep, MatrixXcd* local){

//     int dim = local[0].rows();

//     // Initialize the result variable   
//     MatrixXcd* result = new MatrixXcd[nstep];
//     for (int istep=0; istep<nstep; istep++){
//         result[istep] = MatrixXcd::Constant(dim,dim,doublec(1.0,0.0));
//     }
//     MPI_Barrier(MPI_COMM_WORLD);

//     // Local data reduce to root

//     if (rank==0){printf("Start MPI_Reduce ... \n");}
//     int err;
//     for (int istep=0; istep<nstep; istep++){
//         //if (rank==0){printf("Start rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         //if (rank==10){printf("Start rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         err = MPI_Reduce(local[istep].data(),result[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
//         //if (rank==0){printf("End rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         //if (rank==10){printf("End rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         MPI_Barrier(MPI_COMM_WORLD);
//     }

//     if (rank==0 && err == MPI_SUCCESS){printf("Succeed MPI_Reduce\n ... ");}
//     MPI_Barrier(MPI_COMM_WORLD);

//     if (rank!=0){
//         delete[] result;
//         return NULL;
//     }else{
//         return result;
//     }

//     //return result;
// } 

/* Print help and banner ---------------------------------------------------------*/
// A68: what this binary is, for the log: "Build : CPU | g++ 11.4.0 | Intel(R) MPI Library 2021.18 for Linux* OS | MKL on |
// preset wsl-gnu | commit 1a2b3c4d5-dirty". The compiler from its own macros, the MPI from the library (its first line, cut
// at the first comma: Open MPI's goes on with the package and the repository), the preset from the build
// (do_compile.sh, "cmake" from CMakeLists.txt), the commit from the header both builds write at build time
// (ccex_build_commit.h: "-dirty" when the tree differs from it, "unknown" outside git).
#ifndef CCEX_BUILD_PRESET
#define CCEX_BUILD_PRESET "unknown"
#endif
#if defined(__has_include)
#if __has_include("ccex_build_commit.h")
#include "ccex_build_commit.h"
#endif
#endif
#ifndef CCEX_BUILD_COMMIT
#define CCEX_BUILD_COMMIT "unknown"
#endif
#ifdef CCEX_USE_CUDA
#include "../../include/calculation/gpu/gpu_backend.h"   // gpuBuildDescription
#endif
static std::string buildLine(){
    char compiler[64];
#if defined(__INTEL_LLVM_COMPILER)
    snprintf(compiler, sizeof compiler, "icpx %d.%d.%d", __INTEL_LLVM_COMPILER/10000, __INTEL_LLVM_COMPILER/100%100, __INTEL_LLVM_COMPILER%100);
#elif defined(__INTEL_COMPILER)
    snprintf(compiler, sizeof compiler, "icpc %d.%d", __INTEL_COMPILER/100, __INTEL_COMPILER%100/10);
#elif defined(__clang__)
    snprintf(compiler, sizeof compiler, "clang %d.%d.%d", __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(__GNUC__)
    snprintf(compiler, sizeof compiler, "g++ %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#else
    snprintf(compiler, sizeof compiler, "unknown compiler");
#endif
    char library[MPI_MAX_LIBRARY_VERSION_STRING] = "unknown MPI";
    int length = 0;
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Get_library_version(library, &length);
    std::string mpi(library);
    mpi = mpi.substr(0, mpi.find_first_of("\n,"));
    while (!mpi.empty() && isspace((unsigned char)mpi.back())) mpi.pop_back();
    while (!mpi.empty() && isspace((unsigned char)mpi.front())) mpi.erase(0, 1);
#ifdef CCEX_USE_CUDA
    char device[64];
    gpuBuildDescription(device, sizeof device);
    std::string line = std::string("Build : GPU (") + device + ")";
#else
    std::string line = "Build : CPU";
#endif
#if defined(CCEX_BLAS_DIRECT) && defined(CCEX_BLAS_ACCELERATE)
    const char* mkl = "Accelerate on";
#elif defined(CCEX_BLAS_DIRECT)
    const char* mkl = "MKL on";
#else
    const char* mkl = "MKL off";
#endif
    return line + " | " + compiler + " | " + mpi + " | " + mkl + " | preset " + CCEX_BUILD_PRESET + " | commit " + CCEX_BUILD_COMMIT;
}

static FILE* quietFile = nullptr;
static int quietSaved = -1;
static void logQuietAll(){ logQuietEnd(true); }   // a run that stops while quiet: its whole output, the reason included
void logQuietBegin(bool quiet){
    if (!quiet || rank != 0 || quietFile) return;
    static bool registered = false;
    if (!registered){ atexit(logQuietAll); registered = true; }
    fflush(stdout); std::cout.flush();
    quietFile = tmpfile();
    if (!quietFile) return;   // no temporary file: everything is printed
    quietSaved = dup(STDOUT_FILENO);
    dup2(fileno(quietFile), STDOUT_FILENO);
}
void logQuietEnd(bool all){
    logAlignEnd();   // A76: an aligned region inside the quiet one (a run that stops there) ends first
    if (!quietFile) return;
    fflush(stdout); std::cout.flush();
    dup2(quietSaved, STDOUT_FILENO);
    close(quietSaved);
    quietSaved = -1;
    rewind(quietFile);
    char line[4096];
    bool keep = false;   // a line longer than the buffer comes in pieces: its continuation follows its start
    for (bool start = true; fgets(line, sizeof line, quietFile); ){
        if (start){
            const char* p = line;
            while (*p == ' ' || *p == '\t') ++p;
            // the configuration's own files and spin counts as well: Read QubitFile/GyroFile/BathFile[i]/AvaaxFile,
            // "- Number of bath spins", "Bath = f (n spins) + g (m spins)"
            keep = all || (p[0] == '[' && isalpha((unsigned char)p[1])) || strstr(line, "Warning") ||
                   strstr(line, "warning") || strstr(line, "WARNING") || strstr(line, "Note:") || strstr(line, "Error") || strstr(line, "randomly generated") ||
                   (!strncmp(p, "Read ", 5) && strstr(p, "File")) || strstr(line, "Number of bath spins") ||
                   !strncmp(p, "Bath = ", 7);
        }
        if (keep) fputs(line, stdout);
        start = strchr(line, '\n') != nullptr;
    }
    fclose(quietFile);
    quietFile = nullptr;
    fflush(stdout);
}

// A76: rank 0's standard output between logAlignBegin and logAlignEnd goes to a temporary file (as logQuietBegin), and at
// logAlignEnd its "label : value" rows (indentation 6, or 8 for a sub-item) are printed with the ':' of each table - a run
// of such rows - in one column, after the longest label of the table, and without spaces at the line ends.
static FILE* alignFile = nullptr;
static int alignSaved = -1;
static void logAlignAll(){ logAlignEnd(); }   // a run that stops inside: its output still comes out
void logAlignBegin(){
    if (rank != 0 || alignFile) return;
    static bool registered = false;
    if (!registered){ atexit(logAlignAll); registered = true; }   // registered after logQuietAll: runs before it
    fflush(stdout); std::cout.flush();
    alignFile = tmpfile();
    if (!alignFile) return;   // no temporary file: printed as it is
    alignSaved = dup(STDOUT_FILENO);
    dup2(fileno(alignFile), STDOUT_FILENO);
}
void logAlignEnd(){
    if (!alignFile) return;
    fflush(stdout); std::cout.flush();
    dup2(alignSaved, STDOUT_FILENO);
    close(alignSaved);
    alignSaved = -1;
    rewind(alignFile);
    std::vector<std::string> lines;
    {
        std::string line;
        for (int c; (c = fgetc(alignFile)) != EOF; ){
            if (c == '\n'){ lines.push_back(line); line.clear(); } else line += char(c);
        }
        if (!line.empty()) lines.push_back(line);
    }
    fclose(alignFile);
    alignFile = nullptr;
    // a row: 6 or 8 spaces, a label without ':' that does not start with a space, spaces, ':' and a space, the value
    struct Row { size_t indent; std::string label, value; };
    auto row = [](const std::string& l, Row* r){
        size_t indent = l.find_first_not_of(' ');
        if (indent != 6 && indent != 8) return false;
        size_t colon = l.find(" : ", indent);
        if (colon == std::string::npos) return false;
        std::string label = l.substr(indent, colon - indent);
        while (!label.empty() && label.back() == ' ') label.pop_back();
        if (label.empty() || label.find(':') != std::string::npos) return false;
        size_t v = l.find_first_not_of(' ', colon + 2);
        r->indent = indent; r->label = label; r->value = v == std::string::npos ? "" : l.substr(v);
        while (!r->value.empty() && r->value.back() == ' ') r->value.pop_back();
        return true;
    };
    for (size_t i = 0; i < lines.size(); ){
        std::vector<Row> table;
        Row r;
        size_t j = i;
        while (j < lines.size() && row(lines[j], &r)){ table.push_back(r); ++j; }
        if (table.empty()){
            std::string l = lines[i];
            while (!l.empty() && l.back() == ' ') l.pop_back();
            printf("%s\n", l.c_str());
            ++i;
            continue;
        }
        size_t width = 0;   // the colon column, counted from indentation 6
        for (const Row& t : table) width = std::max(width, t.indent - 6 + t.label.size());
        for (const Row& t : table)
            printf("%s%-*s :%s%s\n", std::string(t.indent, ' ').c_str(), int(width - (t.indent - 6)), t.label.c_str(),
                   t.value.empty() ? "" : " ", t.value.c_str());
        i = j;
    }
    fflush(stdout);
}

void printBanner(){

    // Print the following :
    printf("\n\n");
    printf("    ======================================================================\n\n");
    printf("                   ______   ______  _______ ___   ___\n");
    printf("                  /      | /      ||   ____|\\  \\ /  /\n");
    printf("                 |  ,----'|  ,----'|  |__    \\  V  /\n");
    printf("                 |  |     |  |     |   __|    >   <\n");
    printf("                 |  `----.|  `----.|  |____  /  .  \\\n");
    printf("                  \\______| \\______||_______|/__/ \\__\\\n");
    printf("                __   _   _   _  _           __   _  _   _  _   _ _\n");
    printf("               / _| / \\ | \\_/ || |     //  / _) | |/ / | |/ / | | |\n");
    printf("              ( (_ ( o )| \\_/ || |_   //   \\_ \\ |   <  |   <  | U |\n");
    printf("               \\__| \\_,7|_| |_||___| //    (__/ |_|\\_\\ |_|\\_\\ |___|\n");
    printf("\n");

    //Ref :                                                                              
    //https://wepplication.github.io/tools/asciiArtGen/?fontSelector=Doom&userInput=General+++CCE-X

    time_t t = time(NULL);
    struct tm *d1 = localtime(&t);
    
    printf("\n");
    printf("      Program CCEX starts on %s", asctime(d1));
    printf("      The CCEX code has compiled at '%s'\n", __DATE__);
    // A46: whether this build takes the large clusters to MKL (test/code_verification/check_numerics reads this line)
    logAlignBegin();   // A76: the two rows below in one column
#if defined(CCEX_BLAS_DIRECT) && defined(CCEX_BLAS_ACCELERATE)
    // macOS: Apple's Accelerate in MKL's place. The line keeps its words, which the checks read as "a library takes the
    // large clusters" (check_numerics' compare.py, check_adaptive's regress.py)
    printf("      MKL for large clusters : on (Apple Accelerate in its place; zgemm from dim %d, zheevd from dim %d)\n",
           ccexBlas::gemmMin, ccexBlas::eigMin);
#elif defined(CCEX_BLAS_DIRECT)
    printf("      MKL for large clusters : on (zgemm from dim %d, zheevd from dim %d)\n", ccexBlas::gemmMin, ccexBlas::eigMin);
#else
    printf("      MKL for large clusters : off\n");
#endif
    printf("      %s\n", buildLine().c_str());
    logAlignEnd();
    printf("    ======================================================================\n");
    printf("\n"); 

}

[[noreturn]] void ccexFatal(int code){
    logQuietEnd(true);   // A68: MPI_Abort runs no atexit handler
    fflush(stdout); fflush(stderr);
    int running = 0, finished = 0;
    MPI_Initialized(&running); MPI_Finalized(&finished);
    if (running && !finished) MPI_Abort(MPI_COMM_WORLD, code);
    exit(code);
}

int ccexRefusalPrint(FILE* stream, const char* format, ...){
    if (rank != 0) return 0;
    va_list args;
    va_start(args, format);
    const int n = vfprintf(stream, format, args);
    va_end(args);
    return n;
}

void ccexRefusalPerror(const char* what){
    if (rank == 0) perror(what);
}

[[noreturn]] void ccexRefuse(){
    fflush(stdout); fflush(stderr);
    int running = 0, finished = 0;
    MPI_Initialized(&running); MPI_Finalized(&finished);
    if (running && !finished){
        // review 20261014 F1: the GPU early start (smallClusterWarmUp) may still run its thread on node rank 0, page-locking
        // the shared window; exit() would destroy it joinable (std::terminate) after MPI_Finalize freed that window.
        // Collective like this function (every rank refuses here); nothing when nothing was started.
        smallClusterWarmDown();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
    exit(EXIT_FAILURE);
}
