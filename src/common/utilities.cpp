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
#ifdef CCEX_BLAS_DIRECT
#include "../../include/common/blas_option.h"   // the thresholds printBanner reports
#endif
#include "../../include/calculation/gpu/small_cluster_batch.h"   // smallClusterWarmDown (ccexRefuse)

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
 *  frame whose z axis IS the qubit symmetry axis. This is NOT the active rotation
 *  that carries bath_axis onto qubit_axis; using that one would put the field and the
 *  ZFS on the wrong axis.
 *
 *  R^T is the inverse, so what R^T maps onto qubit_axis is the computational [0,0,1] --
 *  NOT bath_axis. The two coincide only in the special case bath_axis = [0,0,1]; for any
 *  other bath_axis, R^T * bath_axis is some other vector entirely.
 *
 *      ez = q                  the new z axis is the qubit axis
 *      ex = normalize(q x b)   bath_axis only fixes the AZIMUTH of x and y about ez;
 *                              it does not change where z points
 *      ey = ez x ex            right-handed
 *
 *  ex, ey, ez are the ROWS of R, so (R*v)[k] is the component of v along the k-th new
 *  axis. For bath_axis=[0,0,1], qubit_axis=[1,1,1] this gives
 *
 *      ex = [ 1,-1, 0]/sqrt(2),  ey = [ 1, 1,-2]/sqrt(6),  ez = [ 1, 1, 1]/sqrt(3)
 *
 * @param bath_axis  reference axis of the bath file frame (azimuth reference only)
 * @param qubit_axis qubit symmetry axis, in the bath file frame
 * @param R          [out] 3x3 rotation matrix
 */
void buildQubitAlignedRotation(const double bath_axis[3], const double qubit_axis[3], double R[3][3]){

    double b[3], q[3];
    normalizeAxisOrDie(bath_axis ,b,"bath_axis");
    normalizeAxisOrDie(qubit_axis,q,"qubit_axis");

    double ex[3], ey[3], ez[3];
    ez[0] = q[0]; ez[1] = q[1]; ez[2] = q[2];

    rotCross(q,b,ex);
    double cross_norm = sqrt(rotDot(ex,ex));

    if (cross_norm < ROT_COLLINEAR_TOL){
        // q x b vanishes, so the azimuth of the new x axis is undefined. The one case
        // that still has an answer is "the qubit axis is already +z and so is the bath
        // reference" -- then the frame does not move at all and R is the identity.
        bool q_is_z = (fabs(q[0]) < ROT_COLLINEAR_TOL) && (fabs(q[1]) < ROT_COLLINEAR_TOL) && (q[2] > 0.0);
        bool b_is_z = (fabs(b[0]) < ROT_COLLINEAR_TOL) && (fabs(b[1]) < ROT_COLLINEAR_TOL) && (b[2] > 0.0);

        if (q_is_z && b_is_z){
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){ R[i][j] = (i==j) ? 1.0 : 0.0; }
            }
            return;
        }

        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): bath_axis and qubit_axis are parallel or antiparallel\n");
        ccexRefusalPrint(stderr,"  normalized bath_axis  : [ %g, %g, %g ]\n",b[0],b[1],b[2]);
        ccexRefusalPrint(stderr,"  normalized qubit_axis : [ %g, %g, %g ]\n",q[0],q[1],q[2]);
        ccexRefusalPrint(stderr,"  |qubit_axis x bath_axis| = %g\n",cross_norm);
        ccexRefusalPrint(stderr,"The cross product fixes the azimuth of the new x axis, so collinear axes\n");
        ccexRefusalPrint(stderr,"leave the frame undefined. Pick a bath_axis that is not along qubit_axis\n");
        ccexRefusalPrint(stderr,"(the only accepted collinear input is bath_axis=[0,0,1] with qubit_axis=[0,0,1],\n");
        ccexRefusalPrint(stderr," which is the identity transformation).\n");
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
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): the constructed matrix is not a valid rotation\n");
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
 * @param bath_axis              reference axis of the bath file frame, in that frame (as given to buildQubitAlignedRotation)
 * @param qubit_axis             qubit symmetry axis, in the bath file frame
 * @param R                      the rotation, original (bath file) frame -> computational frame
 * @param r0                     the point held fixed (the reference qubit), in the original frame
 * @param bfield                 Config.bfield, i.e. the field in the COMPUTATIONAL frame
 * @param defect_axis_reference  name of the recorded reference axis, or NULL if none
*/
void reportQubitAlignedRotation(const double bath_axis[3], const double qubit_axis[3],
                                const double R[3][3], const double r0[3],
                                const float bfield[3], const char* defect_axis_reference){

    double b[3], q[3];
    normalizeAxisOrDie(bath_axis ,b,"bath_axis");
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
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","bath axis, source frame" ,b[0],b[1],b[2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","qubit axis, source frame",q[0],q[1],q[2]);
    printf("\n");
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","computational x axis, source frame",R[0][0],R[0][1],R[0][2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","computational y axis, source frame",R[1][0],R[1][1],R[1][2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","computational z axis, source frame",R[2][0],R[2][1],R[2][2]);
    printf("\n");
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 0",R[0][0],R[0][1],R[0][2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 1",R[1][0],R[1][1],R[1][2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation matrix R, row 2",R[2][0],R[2][1],R[2][2]);
    printf("      %-32s: %.15lf\n","det(R)",det);
    printf("      %-32s: %.3e\n","max|R^T R - I|",ortho_err);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","R * normalize(qubit_axis)",qz[0],qz[1],qz[2]);
    printf("\n");
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","rotation origin r0, source frame",r0[0],r0[1],r0[2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","B field, computational frame",
           (double)bfield[0],(double)bfield[1],(double)bfield[2]);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","B field, source frame (R^T B)",bsrc[0],bsrc[1],bsrc[2]);
    printf("\n");

    if (defect_axis_reference != NULL){
        printf("      %-32s: %s\n","Defect axis reference",defect_axis_reference);
        printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","Reference axis, source frame",q[0],q[1],q[2]);
        printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","Reference axis, calc frame",qz[0],qz[1],qz[2]);
        printf("      %-32s: recorded for future JT-axis classification;\n","Status");
        printf("      %-32s  not applied to the Defect Hamiltonian\n","");
        printf("\n");
    }

    printf("      Every qubit position and every bath-spin position gets the same rigid\n");
    printf("      transform r_new = r0 + R*(r_old - r0), about the reference qubit.\n");
    printf("      Stored HF/QD tensors that are in the bath frame move with them, and each\n");
    printf("      Qubit.intmap tensor follows its own recorded tensor_frame / provenance --\n");
    printf("      the transformed and kept counts are listed below.\n");
    printf("      Defect entries follow their own coordinate_frame: bath-frame axes,\n");
    printf("      displacement vectors, and tensors are transformed; qubit-frame entries stay.\n");
    printf("      Read as written in the computational frame: legacy top-level \"qzfs\" and \"bfield\".\n");
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
    printf("      %-18s:   [ ", key);

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
    printf("      %-18s:   ", key);

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
    printf("      %-18s:   %-21s\n", key, val);
}

void printStructElementChar2d(char* key, char** val, int n){
    printf("      %-18s:   [ ", key);
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
    printf("      %-18s:   %-21d\n", key, val);
}

void printStructElementInt1dIdx(char* key, int* val, int n){
    printf("      %-18s:   [ ", key);
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
    printf("      %-18s:   %-21.6g\n", key, val);
}

void printStructElementFloat1d(char* key, float* val, int n){
    printf("      %-18s:   [ ", key);
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
    printf("      %-18s:   %-21.3g\n", key, val);
}

void printStructElementDouble1d(char* key, double* val, int n){
    printf("      %-18s:   [ ", key);
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
    printf("      %-18s:   %-21s\n", key, val ? "true" : "false");
}

void printLine(){
    printf("    -----------------------------------------------------------------\n");
}

void printLineSection(){
    printf("\n    ===============================================================\n\n");
}

void printTitle(char* title){
    printf("    < %s > \n\n",title);
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
std::vector<int> rankClusterRows(int m, const std::function<double(int)>& rowCost, std::vector<double>& loads){
    std::vector<int> rows;
    loads.resize(nprocess, 0.0);
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
        rows[n] = rankClusterRows(clusters[n][0][0]-1, rowCost, loads);
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
void printBanner(){

    // Print the following :
    printf("\n\n");
    printf("    ======================================================================\n\n");
    printf("                   ______   ______  _______ ___   ___ \n");
    printf("                  /      | /      ||   ____|\\  \\ /  / \n");
    printf("                 |  ,----'|  ,----'|  |__    \\  V  /  \n");
    printf("                 |  |     |  |     |   __|    >   <   \n");
    printf("                 |  `----.|  `----.|  |____  /  .  \\  \n");
    printf("                  \\______| \\______||_______|/__/ \\__\\ \n");
    printf("                __   _   _   _  _       _   _     _   _   _ _   \n");
    printf("               / _| / \\ | \\_/ || |     //  / \\   | | / \\ | | |  \n");
    printf("              ( (_ ( o )| \\_/ || |_   //  | o |n_| |( o )| U |  \n");
    printf("               \\__| \\_,7|_| |_||___| //   |_n_|\\__/  \\_/ |___|  \n");
    printf("\n");

    //Ref :                                                                              
    //https://wepplication.github.io/tools/asciiArtGen/?fontSelector=Doom&userInput=General+++CCE-X

    time_t t = time(NULL);
    struct tm *d1 = localtime(&t);
    
    printf("\n");
    printf("            Program CCEX starts on %s ", asctime(d1));
    printf("               The CCEX code has compiled at '%s' \n", __DATE__);
    // A46: whether this build takes the large clusters to MKL (test/code_verification/check_numerics reads this line)
#ifdef CCEX_BLAS_DIRECT
    printf("               MKL for large clusters : on (zgemm from dim %d, zheevd from dim %d)\n", ccexBlas::gemmMin, ccexBlas::eigMin);
#else
    printf("               MKL for large clusters : off\n");
#endif
    printf("    ======================================================================\n");
    printf("\n"); 

}

void printHelp(){

    printf("\tThis code is to simulate the many-body spin dynamics\n");
    printf("\n\n");

    printf("\tHow to use this code : CCECode [OPTION]... [FILE]...\n");
    printf("\t\t-f\t:\tUse the condition input-file\n");
    printf("\t\t-m\t:\tUse the external 'calMethod' ('single, ensemble, semi')\n");
    printf("\t\t-I\t:\tuse the external 'BathFile'\n");
    printf("\t\t-s\t:\tuse the external 'StateFile'\n");
    printf("\t\t-S\t:\tuse the external 'ExStateFile'\n");
    printf("\t\t-a\t:\tuse the external 'AvaaxFile'\n");
    printf("\t\t-B\t:\tuse the external 'B0'\n");
    printf("\t\t-o\t:\tuse the external 'Savefile' name\n");
    printf("\t\t-v\t:\tprint the information of A-&Q-tensor file\n");
    printf("\t\t-h\t:\tdisplay this help and exit\n");
    printf("\n\n");


    printf("\tAlong the 'calMethod' option, we can choose the calculation method!\n");
    printf("\t 'calMethod = single'   : use the single sample calculation method\n");
    printf("\t 'calMethod = ensemble' : use the ensemble CCE calculation method (default)\n");
    printf("\t 'calMethod = semi'     : use the semi-classical CCE calculation method\n");
    printf("\n\n");


    printf("\tTo run this code, at least you need to the three file before!\n");
    printf("\t\t1. Inputfile about nuclear configure which have spin positions and isotopes information\n");
    printf("\t\t2. Gyromagnetic ratio about isotope nuclear\n");
    printf("\t\t3. Defect position file // Write the defect position in condition file\n");
    printf("\n\n");


    printf("\tAnd then, you need another element!\n");
    printf("\t\tCommon options         : 'BathFile', 'Order', 'Pulse', 'B0', deltaT', 'nStep', 'SaveFile', HFOpt', 'QuadOpt', 'HFmediOpt', ...\n");
    printf("\t\tSingle method          : 'StateFile', 'ExStateFile', 'AvaaxFile', 'rDsrdr', 'enJTEOpt', 'nExspin'\n");
    printf("\t\t                         '&Exspin' tag related stuff (For detail See below)\n");
    printf("\t\tEnsemble method        : ''\n");
    printf("\t\tsemi-classical method  : 'IntegStep' \n");
    printf("\n\n");


    printf("\tIf you use the condition file, there are many conditions...\n");
    printf("\t\t'GyroFile'                \t: gyro magnetic ratio input-file (neccessary)\n");
    printf("\t\t'DefectFile'              \t: input-file or array related defect position (neccessary) (unit : angsrom)\n");
    printf("\t\t'BathFile'                \t: input-file of bath Configure\n");
    printf("\t\t'StateFile'               \t: spin state input file\n");
    printf("\t\t'ExStateFile'             \t: addtional spin's state input file \n");
    printf("\t\t'AvaaxFile'               \t: the numbering of available geometry of the extra spin's input file\n");
    printf("\t\t'SaveFile', 'SaveFileNoDiv', 'SaveFileWiDiv'\t: the file name to save the result of CCE\n");
    printf("\n");
}

[[noreturn]] void ccexFatal(int code){
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
