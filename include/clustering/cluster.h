/**
 * @file cluster.h
 * @brief The Cluster structure: cluster settings and the clusters of the run.
 */
#ifndef __CCEX_CLUSTER_H_
#define __CCEX_CLUSTER_H_

#include "input/general.h"
#include "physical/bath.h"
#include "common/utilities.h"
#include "physical/defect.h"

/**
 * @struct Cluster
 * @brief The cluster settings of the run (method, order, rdip, ...) and the clusters made from them.
 * @details Each parameter is read from the input file or options.
*/
typedef struct {

    /**
     * @brief Clusterizing order : See details in @ref Config::order 
    */
    int order;

    /** pcce
       @brief parameters of pCCE
       @todo Hyeonsu should be write this part in detail.
    */
    int  sK; /**< Partition CCE variable pCCE(order N, sK) */
    int  max_trial;  /**< It's maximum value for kmeans algorithm opitimization. */
    int  max_iter;   /**< It's maximum value for kmeans algorithm opitimization. */

    bool kmeans_pp;   /**< If you want to use kmeans++ method, True. (Default: True) */
    bool iter_detail; /**< If you want to print detail for parition information, True. */

    /**
     * @brief Clustering algorithm : See details in @ref Config::method 
    */
    char method[MAX_CHARARRAY_LENGTH];

    /**
     * @brief The number of clusters for each order "k"
     * @details 
     * 
     *  nk is useful when you calculate CCE with order >= 3
     *  At order >=3, Creating clusters is based on rdip, a ton of cluster will be formed (at least 1e6).      
     *  This value would consider nk # clusters for each order after sorting this clusters by coupling-strength.
     * 
     * Input file format : an integer array of exactly order+1 entries, the index being the order
     *   e.g. order = 3 : "nk" : [0, 0, 0, 30000]  (keep the 30000 strongest clusters of order 3)
     *   An entry 0 means no limit for that order (all clusters within rdip).
     * 
     * Array format: 
     *   nk[order + 1] = {0,} (default, including all clusters for each order)
     *   nk[i=0] = The number of clusters at order 0 (e.i. nk[0] = 1 always)
     *   nk[i>i] = The number of clusters at order i
    */
    int* nk; 

    /**
     * @brief Include all sub-clusters of the highest order : on | off
     * @details 
     * 
     *  When calcuting the coherence in electron spin bath, 
     *  CCE with order >= 3 gives lots of diverged point in coherence fucntion.
     *  If turn on this option, 
     *  it will include all sub-clusters of the highest order 
     *  and give the more robust result.
     *  (default : on)
    */
    bool addsubclus;

    /**
     * @brief Each rank keeps only its own pair rows
     * @details
     *  localRequest: set by runOneCalculation when nothing after the split reads every pair row (no output.cluster_only, no
     *  cluster_contributions, no adaptive work request); clusterizePairs may then keep only this rank's rows.
     *  pairsLocal: it did. clusinfo[2] still has every row slot and the count, but a row is allocated only for
     *  this rank's share and the first and last four (Cluster_reportClusinfo); the other slots are NULL.
     *  pairRows (npairRows entries, ascending, 1-based): this rank's order-2 rows under the split it decided
     *  (rankClusterRows: by cost, largest first), which MPI_getLocalClusters then uses instead of the rows;
     *  pairLoads (nprocess entries): every rank's cost after that split, which order 1 starts from.
    */
    bool localRequest;
    bool pairsLocal;
    int* pairRows;
    int npairRows;
    double* pairLoads;

    /**
     * @brief The clusters for each order "k"
     * @details 
     * 
     * clusters[i][j][k]
     * i : i-th order
     * j : (if j==0) The number of clusters at order i
     *     (if j!=0) j-th clusters at order i
     * k : (if k==0) The number that you will have to 
     *                multiply(+)/divide(-) in the coherence calcualtion
     *     (if k!=0) The index of the spin in the cluster
     * 
     * e.g.
     * clusters[i=0][j=0][k=0] :
     *      The number how many you should calculate in 0-th order
     * clusters[i=1][j=0][k=0] :
     *       The number of 1st order clusters
     * clusters[i=1][j!=0][k=0] :
     *       The number how many you should calculate 
     *       the coherence for 1st order and j-th cluster
     * clusters[i=1][j!=0][k!=0] :
     *       The spin index in j-th cluster at 1st order
     * 
    */
    int*** clusinfo;

} Cluster;

/* High Level --------------------------------------------------------*/



// init
Cluster* Cluster_init();
void Cluster_clusterize(Cluster* cls, BathArray* ba, Config* config,
                       QubitArray* qa=nullptr, DefectArray* dfa=nullptr);

// A31: a cluster's cost D^2 (D + nstep), D its dimension (the qubit's too in gCCE), for MPI_getLocalClusters.
// Built from the bath as it is now: pCCE cuts the bath in clusterize, so call it after the bath is final.
std::function<double(int,const int*)> Cluster_costFunction(Config* cnf, QubitArray* qa, BathArray* ba, DefectArray* dfa);
// The D of that cost (A98: gpu_options.ranks_per_gpu sends a cluster to the GPU ranks by it)
std::function<double(int,const int*)> Cluster_dimFunction(Config* cnf, QubitArray* qa, BathArray* ba, DefectArray* dfa);

// free
void Cluster_freeAll(Cluster* cls);

// MPI
// int*** Clsuter_getLocalClusters_MPI(Cluster* cls);

// find clusters
// void Cluster_find_clusinfo_hash(Cluster* cls, BathArray* batharr, Config* config); // normal (hash) algorithm
// void Cluster_find_clusinfo_dsj(Cluster* cls, BathArray* batharr, Config* config); // disjoint algorithm
// void Cluster_find_clusinfo_itb(Cluster* cls, BathArray* batharr, Config* config); // inter-bathcluster algorithm
// void Cluster_find_clusinfo_dsjitb(Cluster* cls, BathArray* batharr, Config* config); // disjoint + inter-bathcluster algorithm

// report
void Cluster_report(Cluster* cls);
void Cluster_reportNk(Cluster* cls);
void Cluster_reportClusinfo(Cluster* cls);
void reportClusinfo(int*** clusinfo, int order);

/* Low Level --------------------------------------------------------*/

// clusinfo controller

void Cluster_allocClusinfo(Cluster* cls, int order); // alloc clusinfo[order+1] : 0 ~ order
void Cluster_freeClusinfo(Cluster* cls);

void Cluster_setClusinfo_0th(Cluster* cls); // Actually do nothing
void Cluster_setClusinfo_1th(Cluster* cls, BathArray* ba); // bathspins = clusinfo[1][j][k]
int Cluster_setClusinfo_addcluster(Cluster* cls, int order, int iter, int* cluster); // add new cluster, return the index of the cluster
void Cluster_setClusinfo_chgcluster(Cluster* cls, int order, int* cluster, int ic); // change cluster at ic-th cluster
void Cluster_setClusinfo_chgiter(Cluster* cls, int order, int iter, int ic); // change iter at ic-th cluster 

//void Cluster_setSk(Cluster* cls, int sk); //
//void Cluster_getSk(Cluster* cls); // output



// what's iter :     
// The number how many you will multiply/divide 
// in the coherence calculation for 0 th order

int Cluster_getClusinfo_ncluster(Cluster* cls, int order); // get the number of clusters at order
int* Cluster_getClusinfo_itercluster(Cluster* cls, int order, int ic); // get ic-th iter+cluster at order [0] = iter, [1~] = cluster
int Cluster_getClusinfo_iter(Cluster* cls, int order, int ic); // get ic-th iter+cluster at order [0] = iter, [1~] = cluster
int* Cluster_getClusinfo_cluster_copy(Cluster* cls, int order, int ic); // get ic-th cluster at order (copy)
int*** Cluster_getClusinfo(Cluster* cls);

// nk controller
void Cluster_allocNk(Cluster* cls);
void Cluster_freeNk(Cluster* qa);
int Cluster_getNk_order(Cluster* cls, int i); // get the number of clusters at ith-orderr at ith-order
int* Cluster_getNk(Cluster* cls);

// others controller
char* Cluster_getMethod(Cluster* cls);
int Cluster_getOrder(Cluster* cls);
bool Cluster_getAddsubclus(Cluster* cls);
int Cluster_getSk(Cluster* cls);
int Cluster_getMax_iter(Cluster* cls);
int Cluster_getMax_trial(Cluster* cls);
bool Cluster_getKmeans_pp(Cluster* cls);
bool Cluster_getIter_detail(Cluster* cls);

// set 
void Cluster_setOrder(Cluster* cls, int order);
void Cluster_setMethod(Cluster* cls, char* method);
void Cluster_setAddsubclus(Cluster* cls, bool addsubclus);
void Cluster_setNk(Cluster* cls, int* nk); // set the number of clusters at ith-order
void Cluster_setSk(Cluster* cls, int sK);
void Cluster_setMax_trial(Cluster* cls, int max_trial);
void Cluster_setMax_iter(Cluster* cls, int max_iter);
void Cluster_setKmeans_pp(Cluster* cls, bool kmeans_pp);
void Cluster_setIter_detail(Cluster* cls, bool iter_detail);

// Preflight inventory: exact raw expansion, without allocating propagation matrices.
void writeClusterInventory(Cluster*,BathArray*,QubitArray*,DefectArray*,const Config*,
                           const char* outfile,int loadedSites);

#endif // __CCEX_CLUSTER_H_
