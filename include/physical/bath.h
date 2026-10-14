#ifndef __CCEX_BATH_H_
#define __CCEX_BATH_H_

#include "common/utilities.h"
#include <vector>
#include "physical/qubit.h"
#include "input/general.h"

/**
 * @struct BathSpin
 * @brief BathSpin structure contains each BathSpin's information
 * @note When you make this structure, 
 *       the following action would be done automatically.
 *       (1) state value would be allocated
 *       (2) The dimension of alpha, beta would be checked if it is 2*spin+1
 *       (3) The state would be normalized
 */
typedef struct {

    /* Basic spin information (Essential) */

    /**
     * @brief the name of spin
     * @details Read from the bath file(s) given by "bathfile" (or -I) (no default value)
    */ 
    char name[MAX_CHARARRAY_LENGTH]; 

    /**
     * @brief Spin quantum number
     * @details From the optional "gyrofile" when it lists the species, otherwise from the
     *          built-in isotope table (src/input/isotope_table.cpp); an unknown species is an error
    */
    float spin;    

    /**
     * @brief Gyromagnetic ratio (Unit : radkHz/Gauss)
     * @details From the optional "gyrofile" when it lists the species, otherwise from the
     *          built-in isotope table; the file value wins when both exist
    */
    double gyro;   

    // /**
    //  * @brief Quadrupole moment (Unit : 10^-30 [m^2])
    //  * @details Value can be obtained from &Bath -> qtensorfile (default : 0.0)
    //  *      - cf. eQ = 20.44[Q/millibarn = Q * 10^-31 m^2] , in DFT eQ = 2.044 * 10^-30 [m^2]
    //  * @note Unit in input file (cce.in) : (Unit : 10^-30 [m^2])
    //  * @todo Add documentation for cce.in : &Bath -> qtensorfile 
    //  * @todo When you read quadrupole tensor, 
    //  *       - you have to get them with unit Hatree/Bohr_radius^2 (QE) -> Think about How to do this
    //  *       - you have to get them with unit MHz (VASP) -> Think about How to do this
    //  *       - which one is better? -> Think about How to do this
    //  * @todo Now, there is no default value for eq. In future, add the information of each spin' eQ in code internally.
    // */      
    // double eq;
  
    /**
     * @brief Spin position (Unit : Angstrom)
     * @details Read from the bath file(s) given by "bathfile" (no default value)
     * @note The x,y,z position can be shifted per bath file by "bath_shift"
    */
    double xyz[3]; 


    /* Spin state (Optional)*/

    /**
     * @brief Spin sub-level m in the z-basis, a single number (used when "nstate" != 0, the sampled-state approach)
     * @details Read from "bath_state_file" when given, otherwise drawn at random ("nstate" states, "seed" fixes them)
    */
    float state; 


    /* Interaction Information (Optional) */

    /**
     * @brief The detuning value of a bath spin (Unit : radkHz)
     * @details An artificial shift of the bath spin's energy levels. There is no bath-level input key:
     *          it is 0 unless the spin is a Defect main spin, whose Defect "detuning" (MHz) is added
    */
    double detuning; 

    /**
     * @brief The mean-field effect from zz-coupling with other bath spins (Unit : radkHz)
     * @details Automatically set if &General -> "nstate" tag is non-zero (default : 0.0)
     * @todo Add documentation for cce.in : &Bath -> "disorder" tag 
    */
    double disorder; 

    /**
     * @brief The transverse (x, y) part of the same static mean field, J[:,z] m_j summed over the other bath spins
     *        (Unit : radkHz); set only with hamiltonian.common.bath_bath_mean_field "xyz", zero otherwise (2026-10-03)
    */
    double disorder_xy[2];
    bool disorder_vector;   /**< disorder_xy was set (vector_mean_field): a cluster also removes its members' x, y field */

    /**
     * @brief The hyperfine interaction between qubits and a bath spin (Unit : radkHz)
     * @details hypf[i] : the hyperfine tensor with the i-th qubit. Point dipole by default
     *          (hamiltonian.common.hyperfine "dipole": "point" without "contact"); the others read
     *          hamiltonian.common.hyperfine.file (MHz), and a Defect species' "hyperfine" sets its spins
    */
    MatrixXcd* hypf; // 2024.01.26 change naming due to &Defect.defect.hypf : hypf to intqb 보류

    /**
     * @brief The quadrupole interaction a bath spin (Unit : radkHz)
     * @details Zero by default (hamiltonian.common.quadrupole "off"); a "structure" reads
     *          hamiltonian.common.quadrupole.file (EFG), and for an
     *          electron main spin of a Defect the Defect "zfs" (MHz) is stored here
    */
    MatrixXcd quad; // 2024.01.26 change naming due to &Defect.defect.quad : quad to intself 보류

    /* &Defect */
    MatrixXcd hypf_sub; // the hyperfine interaction with bath spin (if there is dft data)
    int mainspidx; // the connected main spin index with the coupling strength of hypf_sub

    double mindist; // minimum distance between qubit set and a bath spin
} BathSpin;

/*!
 * @struct BathArray
 * @brief BathArray include all Bath spins' information
 */
typedef struct {

    /**
     * @brief The number of bath spins
     * @todo Add documentation for cce.in : &Bath -> "nspin" tag
    */
    int nspin; /**< The number of BathSpin */  

    /**
     * @brief The array of BathSpin
     * @todo Add documentation for cce.in : &Bath tag
    */
    BathSpin** bath;  

    // spin properties
    int prop_nspecies; /**< The number of species */
    char** prop_names; /**< The name of species */
    double* prop_gyros; /**< The gyromagnetic ratio of species */
    float* prop_spins; /**< The spin quantum number of species */

} BathArray;


// spin pairs
void        BathArray_connectivity(int*** cmap, float*** stmap, BathArray* ba, float rdip, float rdipcut); // connectivity between bath spins
void        makeSparsemap(int*** spmap, int** cmap, int nspin); // sparsemap

// spin interactions
MatrixXcd   BathArray_int_i_j(BathArray* ba, int i, int j);
double      BathArray_dist_i_j(BathArray* ba, int i, int j);

void        BathArray_setBathHypfs(BathArray* ba, QubitArray* qa);

void        BathArray_setBathDisorders(BathArray* ba, bool vector = false);
/**
 * @brief The same disorders (and x, y parts) as BathArray_setBathDisorders, the spins split over the MPI ranks and gathered
 *        (improvements.md A13). Every rank calls it with the same bath and states. The couplings of a rank's spins do not
 *        depend on the states; they are kept in @p cache from the first call while they fit (cache.limit bytes) and reused.
 *        Each spin's sum runs over j in the same order with the same terms, so the values are the same bits.
 */
struct DisorderCache {
    BathArray* ba = nullptr;   ///< the bath the kept couplings belong to
    int lo = 0, hi = 0;        ///< this rank's spins [lo, hi)
    bool vector = false, kept = false;
    std::vector<double> zz, xz, yz;   ///< per spin of [lo, hi), all j; xz/yz only with vector
    size_t limit = size_t(512) << 20;
};
void        BathArray_setBathDisordersRanks(BathArray* ba, bool vector, DisorderCache& cache);
// the x, y components of the static field spin jsp puts on spin isp (J[:,z] m_j), for vector_mean_field
void        BathArray_getBath_i_disorderXY_j(BathArray* ba, int isp, int jsp, double out[2]);
// ! update bath disorder if there is additional spins

double      BathArray_getBath_i_disorder_j(BathArray* ba, int isp, int jsp);

double      BathArray_getOverhaus(BathArray* ba, int iq); // overhaus field for i-th qubit
// ! update bath overhausers if there is additional spins

double      BathArray_getBath_i_overhaus_j(BathArray* ba, int isp, int iq); // overhaus field for i-th qubit

// random state generator
void        BathArray_setBathStatesRandom(BathArray* ba); // random state

// physics
int         BathArray_dim(BathArray* ba);
int         BathArray_dimBath_i(BathArray* ba, int i);
int         BathSpin_dim(BathSpin* bs);

// Hamiltonian
MatrixXcd   BathArray_ZeemanHamil(BathArray* ba, MatrixXcd** sigmas, int ib, float* bfield);
MatrixXcd   BathArray_DetuningHamil(BathArray* ba, MatrixXcd** sigmas, int ib);
MatrixXcd   BathArray_DisorderHamil(BathArray* ba, MatrixXcd** sigmas, int ib, bool rm_overlap);   // z, and x/y when set
MatrixXcd   BathArray_QuadHamil(BathArray* ba, MatrixXcd** sigmas, int ib);
MatrixXcd   BathArray_InteractionHamil(BathArray* ba, MatrixXcd** sigmas, int ib, int jb);
void        BathArray_addInteractionHamil(MatrixXcd& Hb, BathArray* ba, MatrixXcd** sigmas, int ib, int jb);
MatrixXcd** BathArray_PauliOperators(BathArray* ba);

// density matrix
MatrixXcd   BathArray_Rho0(BathArray* ba, bool isEnsemble);
MatrixXcd   BathArray_Psi0(BathArray* ba);

/* Low level functions --------------------------------------------*/

BathArray*  BathArray_init();

//report
void        BathArray_report(BathArray* ba);
void        BathArray_reportBath(BathArray* ba);
void        BathArray_reportBath_i_props(BathArray* ba, int i);
void        BathArray_reportBath_states(BathArray* ba);
void        BathArray_reportBath_detunings(BathArray* ba);
void        BathArray_reportBath_disorders(BathArray* ba);
void        BathArray_reportBath_hypf(BathArray* ba, int nqubit);
void        BathArray_reportBath_quad(BathArray* ba);
void        BathArray_reportBath_hypf_sub(BathArray* ba);
void        BathArray_reportSpinProperties(BathArray* ba);

//alloc
void        BathArray_allocBath(BathArray* ba, int nqubit);
void        BathArray_reallocBath(BathArray* ba, int nspin_old, int nspin_new, int nqubit);
void        BathArray_allocProp(BathArray* ba);
void        BathArray_reallocProp(BathArray* ba, int nspin_old, int nspin_new);
void        BathArray_allocBath_i_hypf(BathArray* ba, int i, int nqubit);

//sets
void        BathArray_setNspin(BathArray* ba, const int nspin);
void        BathArray_setBath_i(BathArray* ba, const BathSpin* bath, int i, int nqubit);
void        BathArray_setBath_i_name(BathArray* ba, const char* name, int i);
void        BathArray_setBath_i_spin(BathArray* ba, const float spin, int i);
void        BathArray_setBath_i_gyro(BathArray* ba, const double gyro, int i);
void        BathArray_setBath_i_xyz(BathArray* ba, const double* xyz, int i);
void        BathArray_setBath_i_state(BathArray* ba, const float state, int i);
void        BathArray_setBath_i_detuning(BathArray* ba, const double detuning, int i);
void        BathArray_setBath_i_disorder(BathArray* ba, const double disorder, int i);
void        BathArray_setBath_i_hypf_j(BathArray* ba, const MatrixXcd hypf, int i, int j); // j : qubit index
void        BathArray_setBath_i_quad(BathArray* ba, const MatrixXcd quad, int i);
void        BathArray_setBath_i_hypf_sub(BathArray* ba, const MatrixXcd hypf_sub, int i);
void        BathArray_setBath_i_mainspidx(BathArray* ba, const int mainspidx, int i);
void        BathArray_setProp_nspecies(BathArray* ba, const int nspecies);
void        BathArray_setProp_names_i(BathArray* ba, const char* name, const int i);
void        BathArray_setProp_gyros_i(BathArray* ba, const double gyro, const int i);
void        BathArray_setProp_spins_i(BathArray* ba, const float spin, const int i);
void        BathArray_setBath_i_mindist(BathArray* ba, const double mindist, int i);

// get
int         BathArray_getNspin(BathArray* ba);
int         BathArray_getProp_nspecies(BathArray* ba);
char**      BathArray_getProp_names(BathArray* ba);
double*     BathArray_getProp_gyros(BathArray* ba);
float*      BathArray_getProp_spins(BathArray* ba);
BathSpin*   BathArray_getBath_i(BathArray* ba, int i);
char*       BathArray_getBath_i_name(BathArray* ba, int i);
float       BathArray_getBath_i_spin(BathArray* ba, int i);
double      BathArray_getBath_i_gyro(BathArray* ba, int i);
double*     BathArray_getBath_i_xyz(BathArray* ba, int i);
float       BathArray_getBath_i_state(BathArray* ba, int i);
double      BathArray_getBath_i_detuning(BathArray* ba, int i);
double      BathArray_getBath_i_disorder(BathArray* ba, int i);
MatrixXcd   BathArray_getBath_i_hypf_j(BathArray* ba, int i, int j);
MatrixXcd   BathArray_getBath_i_quad(BathArray* ba, int i);
MatrixXcd   BathArray_getBath_i_hypf_sub(BathArray* ba, int i);
int         BathArray_getBath_i_mainspidx(BathArray* ba, int i);

// free
void        BathArray_freeAll(BathArray* ba);
void        BathArray_freeProp_names(BathArray* ba);
void        BathArray_freeProp_gyros(BathArray* ba);
void        BathArray_freeProp_spins(BathArray* ba);
void        BathArray_freeBath(BathArray* ba);
void        BathArray_freeBath_i_hypf(BathArray* ba, int i);

// BathSpin
void       BathSpin_setName(BathSpin* bs, char* name);
void       BathSpin_setName_withType(BathSpin* bs, char* name, char* type); // name : mainspin name, type : subspin
void       BathSpin_setSpin(BathSpin* bs, float spin);
void       BathSpin_setGyro(BathSpin* bs, double gyro);
void       BathSpin_setXyz(BathSpin* bs, double* xyz);
void       BathSpin_setXyz_fromRxyz(BathSpin* bs, double* xyz0, double* rxyz); // Input : xyz0 : main spin position, rxyz : sub spin relative position, Output : xyz : sub spin position
void       BathSpin_setState(BathSpin* bs, float state);
void       BathSpin_setDetuning(BathSpin* bs, double detuning);
void       BathSpin_setDisorder(BathSpin* bs, double disorder);
void       BathSpin_setHypf_i(BathSpin* bs, MatrixXcd hypf, int iq);
void       BathSpin_setQuad(BathSpin* bs, MatrixXcd quad);
void       BathSpin_setQuad_fromEFG(BathSpin* bs, MatrixXcd efg, double eq, float spin); // Input : efg : Hartree/Bohr^2, eq : 10e-30 m^2 , spin number, Output : quad : radkHz
void       BathSpin_setHypfSub(BathSpin* bs, MatrixXcd hypf_sub);
void       BathSpin_setMainspidx(BathSpin* bs, int mainspidx);

char*      BathSpin_getName(BathSpin* bs);
float      BathSpin_getSpin(BathSpin* bs);
double     BathSpin_getGyro(BathSpin* bs);
double*    BathSpin_getXyz(BathSpin* bs);
float      BathSpin_getState(BathSpin* bs);
double     BathSpin_getDetuning(BathSpin* bs);
double     BathSpin_getDisorder(BathSpin* bs);
MatrixXcd  BathSpin_getHypf_i(BathSpin* bs, int iq);
MatrixXcd  BathSpin_getQuad(BathSpin* bs);
MatrixXcd  BathSpin_getHypfSub(BathSpin* bs);
int        BathSpin_getMainspidx(BathSpin* bs);

#endif // __CCEX_BATH_H_
