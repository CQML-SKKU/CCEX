#ifndef __CCEX_UTILITIES_UTILITIES_H_
#define __CCEX_UTILITIES_UTILITIES_H_

#include <Eigen/Dense>
#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <functional>
#include <vector>

// Stop the whole run from one rank (D38): a fatal error that only this rank meets - a file it cannot open or
// write, an eigensolver that fails - used to exit() this rank alone, and the others could wait in the next
// collective. MPI_Abort when MPI is running (the unit tests call some of this without it), exit otherwise.
[[noreturn]] void ccexFatal(int code = EXIT_FAILURE);

// An input refusal that every rank reaches alike (A55): the ranks read the same input and the same files and decide
// alike, so rank 0 alone writes the message (ccexRefusalPrint, ccexRefusalPerror) and the ranks meet, finalize and exit 1
// together (ccexRefuse). A plain exit on each rank let mpirun kill the slower one: exit 1 or 255 and one or two messages.
// Only where every rank surely decides the same: elsewhere a rank would wait at the barrier for good (ccexFatal there).
// Without MPI (before MPI_Init, the unit tests) it writes and exits.
int ccexRefusalPrint(FILE* stream, const char* format, ...) __attribute__((format(printf,2,3)));
void ccexRefusalPerror(const char* what);
[[noreturn]] void ccexRefuse();

#define _USE_MATH_DEFINES

/* MPI Global variable ------------------------------------------*/
extern int rank; /**< Rank of the process (MPI)*/
extern int nprocess; /**< Number of processes (MPI)*/
/**--------------------------------------------------------------*/

/* Print Global variable ----------------------------------------*/
extern bool verbosity; /**< Verbosity level */
/**--------------------------------------------------------------*/

/* Data Types ---------------------------------------------------*/
/**
 * @struct DoubleComplex
 * @brief Complex number with double precision
 * @note In local PC, some complex operators from complex.h are not supported.
 * 
*/
typedef struct {
    double real; /**< Real number */
    double imag; /**< Imaginary number */
} DoubleComplex;

/**
 * @struct DoubleTensor
 * @brief Tensor with double precision 
*/
typedef struct {
    double xx, xy, xz; /**< xx, xy, xz */
    double yx, yy, yz; /**< yx, yy, yz */
    double zx, zy, zz; /**< zx, zy, zz */
} DoubleTensor;

/* Eigen Library ------------------------------------------------*/
/** 
 * @typedef doublec
 * @brief Complex number with double precision in Eigen library
*/
typedef std::complex<double>doublec;

/**
 * @typedef MatrixXcd
 * @brief Matrix with complex number with double precision in Eigen library
*/
typedef Eigen::Matrix<doublec, Eigen::Dynamic, Eigen::Dynamic> MatrixXcd;

/**
 * @namespace Eigen
 * @brief Eigen namespace
*/
using namespace Eigen;


/* Math Constants -----------------------------------------------*/
// Planck constant : 1.054 571 729 x 10^-34 [J*s]
// H_BAR is the unit conversion constant (Not the Planck constant)
#define H_BAR 1.054571729

// Unit conversion
#define MHZ_TO_RADKHZ(x) (1 * (2.0 * M_PI * 1.0e+3) * x)
#define RADKHZ_TO_MHZ(x) (1 / (2.0 * M_PI * 1.0e+3) * x)
#define KHZ_TO_RADKHZ(x) (1 * (2.0 * M_PI) * x)
#define RADKHZ_TO_KHZ(x) (1 / (2.0 * M_PI) * x)

// Gyromagnetic ratio of electron spin [radkHz/G]
// Parenthesized and WITHOUT a trailing semicolon. The semicolon used to be part of the
// macro body, so it only worked in `double g = GAMMA_ELECTRON;` -- where the extra `;`
// became a null statement -- and was a syntax error anywhere an expression was wanted,
// e.g. foo(GAMMA_ELECTRON) or 2*GAMMA_ELECTRON.
#define GAMMA_ELECTRON (-17608.597050)

/* File Constants -----------------------------------------------*/
#define MAX_FILEPATH 500 /**< Maximum length of file path */

/* cce.in line constants ----------------------------------------*/
#define MAX_OPTION_ARRAY_LENGTH 500 /**< Maximum length of option array */
#define MAX_OPTION_ELEMENT_LENGTH 50 /**< Maximum length of option element */
#define MAX_CHARARRAY_LENGTH 20 /**< Maximum length of structure char array */

#define MAX_MAINTAG_SUBOPTION_NUMBER 500 /**< Maximum number of suboptions for each &MAINTAG */
#define MAX_MAINTAG_SUBOPTION_LINE_LENGTH 1000 /**< Maximum length of suboption element for each &MAINTAG */


/* utils -------------------------------------------------------*/
double* MatrixXcdToDouble1d(MatrixXcd mat);
MatrixXcd Double1dToMatrixXcd(double* val, int n);
int Double_is_same_row(double* row1, double* row2);

/* math functions -----------------------------------------------*/

// geometry
double dist(double spin1[],double spin2[]);
double cosTheta(double spin1[], double spin2[],double dist);
double sinTheta(double spin1[], double spin2[], double dist);
double cosPhi(double spin1[], double spin2[]);
double sinPhi(double spin1[], double spin2[]);

// bath-coordinate rotation
void normalizeAxisOrDie(const double v[3], double out[3], const char* what);
double angleBetweenDeg(const double u[3], const double v[3]);
void buildQubitAlignedRotation(const double bath_axis[3], const double qubit_axis[3], double R[3][3]);
void rotateAboutPoint(double xyz[3], const double R[3][3], const double r0[3]);
MatrixXcd rotateTensor(const MatrixXcd T, const double R[3][3]);
void reportQubitAlignedRotation(const double bath_axis[3], const double qubit_axis[3],
                                const double R[3][3], const double r0[3],
                                const float bfield[3], const char* defect_axis_reference);

/* Physics functions --------------------------------------------*/

// spin state control
float* substates(float S);
bool isSubLevel(float S, float ms);
MatrixXcd getSpinor(float S, float ms);
MatrixXcd kron(const MatrixXcd& a, const MatrixXcd& b);
MatrixXcd partialtrace(MatrixXcd Mij, int dimrow, int dimcol);
double calNorm(MatrixXcd m);
int normalize(MatrixXcd* m);
float findZbasisSubLevel(MatrixXcd spinor);

MatrixXcd powMatrixXcdElementWise(MatrixXcd a, int n);
MatrixXcd mulMatrixXcdElementWise(MatrixXcd a, MatrixXcd b);
MatrixXcd divMatrixXcdElementWise(MatrixXcd a, MatrixXcd b, double* dropped);
//MatrixXcd divMatrixXcdElementWise(MatrixXcd a, MatrixXcd b);

/* Easy print ---------------------------------------------------*/
void printInlineMatrixXcd(char* key, MatrixXcd mat);
void printStateInDiracNot(char* key, MatrixXcd mat);
void printStructElementChar(char* key,char* val);
void printStructElementChar2d(char* key, char** val, int n);
void printStructElementInt(char* key, int val);
void printStructElementInt1dIdx(char* key, int* val, int n);
void printStructElementFloat(char* key, float val);
void printStructElementFloat1d(char* key, float* val, int n);
void printStructElementDouble(char* key, double val);
void printStructElementDouble1d(char* key, double* val, int n);
void printStructElementBool(char* key, bool val);
void printLine();
void printLineSection();
void printTitle(char* title);
void printSubTitle(char* title);
void printMessage(char* message);

/**
 * @brief strcpy into a fixed buffer, refusing a string that does not fit
 * @details Names and paths come from files a user writes by hand. `size` is the
 *          buffer's, terminator included; `what` names the value in the error.
 */
void copyBounded(char* dst, const char* src, size_t size, const char* what);

/* Find index ---------------------------------------------------*/
int findIndexInt(int* array, int ista, int iend, int val); // find in the range of ista <= i <= iend
int findIndexCharFix(char array[][MAX_CHARARRAY_LENGTH], int ista, int iend, char* val); // find in the range of ista <= i <= iend, strcasecmp
int findIndexChar(char** array, int ista, int iend, char* val); // find in the range of ista <= i <= iend, strcasecmp
int findIndexFloat(float* array, int ista, int iend, float val); // find in the range of ista <= i <= iend

/* Type checker ---------------------------------------------------*/
int isStringDouble(char *s);

/* Quick Sort --------------------------------------------------*/
void QuickSort(int** d_Array, int left, int right);
int Partition(int** d_Array, int left, int right);
void Swap(int** d_Array, int a, int b);

/* MPI ---------------------------------------------------------*/
void para_range(int n1,int n2, int nprocs, int myrank, int*ista, int *iend);
int min(int x, int y);
// cost (A31, A45): with a cost per cluster (n, row of clusters[n]: power then n sites), the orders are split by
// rankClusterRows from the highest down with the ranks' loads carried (largest cost first, to the least-loaded rank);
// an order whose clusters all cost the same, or no cost, keeps the equal-count contiguous split (para_range). A rank's
// clusters keep the list's order.
// givenOrder (A38): that order's rows are givenRows and the loads after it givenLoads (rankClusterRows of the full
// list, made before its rows were dropped; it must be the highest order) instead of ones computed from the rows.
int*** MPI_getLocalClusters(int order, int*** clusters,
                            const std::function<double(int,const int*)>& cost = nullptr,
                            int givenOrder = 0, const std::vector<int>* givenRows = nullptr,
                            const std::vector<double>* givenLoads = nullptr);
// The rows (1-based, ascending) this rank takes of m clusters with costs rowCost(0..m-1) in list order (A45); loads
// (every rank's cost so far, resized to nprocess) is added to.
std::vector<int> rankClusterRows(int m, const std::function<double(int)>& rowCost, std::vector<double>& loads);
// Free the owned copy returned by MPI_getLocalClusters, including its headers.
// The zeroth header is allocated even on ranks whose zeroth weight is zero.
void MPI_freeLocalClusters(int order, int**** localClusters);
// MatrixXcd* MPI_reduceLocalResult(int nstep, MatrixXcd* local);

/* Print help and banner ---------------------------------------------------------*/
void printBanner();
void printHelp();

#endif // __CCEX_UTILITIES_UTILITIES_H_

