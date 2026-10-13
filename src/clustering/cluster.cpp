#include "../../include/clustering/cluster.h"
#include "../../include/common/utilities.h"
#include "../../include/common/memory.h"
#include "../../include/clustering/cluster_hash.h"
#include "../../include/clustering/cluster_pairs.h"
#include "../../include/clustering/cluster_pcce.h"
#include "../../include/clustering/cluster_dsjitb.h"
#include "../../include/clustering/bundle/bundle_cluster.h"

/* High Level --------------------------------------------------------*/

// init
Cluster* Cluster_init(){
    Cluster* cls = (Cluster*)allocArray1d(1,sizeof(Cluster));
    Cluster_setOrder(cls, 0);
    Cluster_setAddsubclus(cls, false);
    cls->nk = NULL;
    cls->method[0] = '\0';
    cls->clusinfo = NULL;
    cls->localRequest = false;
    cls->pairsLocal = false;
    cls->pairRows = NULL;
    cls->npairRows = 0;
    cls->pairLoads = NULL;
    
    return cls;
}

std::function<double(int,const int*)> Cluster_costFunction(Config* cnf, QubitArray* qa, BathArray* ba, DefectArray* dfa){
    // A31: one D^3 eigensolve and nstep steps of D^2 work; measured on h-BN (dimensions 3..49) it follows the time per
    // cluster within about 20 %. (Moved here from main.cpp by A38, which also needs it inside clusterizePairs.)
    std::vector<double> siteDim(ba->nspin, 1.0);
    for (int i=0; i<ba->nspin; i++){
        siteDim[i] = BathSpin_dim(BathArray_getBath_i(ba,i));
        const int idf = DefectArray_findDefectIndex(dfa,BathArray_getBath_i_name(ba,i));
        for (int isp=0; idf>=0 && isp<DefectArray_getDefect_idf_naddspin(dfa,idf); isp++)
            siteDim[i] *= BathSpin_dim(DefectArray_getSubbath_i_isp(dfa,i,isp));
    }
    const double qubitDim = strcasecmp(cnf->method,"gcce") ? 1.0 : QubitArray_dim(qa);
    const double nstep = Config_getNstep(cnf);
    return [siteDim, qubitDim, nstep](int n, const int* row){
        double d = qubitDim;
        for (int j=1; j<=n; j++) d *= siteDim[row[j]];
        return d*d*(d+nstep);
    };
}

void Cluster_clusterize(Cluster* cls, BathArray* ba, Config* config,QubitArray* qa,DefectArray* dfa){

    if(!strcasecmp(config->cluster_mode,"adaptive")) {
        clusterizeAdaptive(cls,ba,qa,dfa,config);
        return;
    }

    int order = Cluster_getOrder(cls);
    int nspin = BathArray_getNspin(ba);
    char* method = Cluster_getMethod(cls);

    ///////////////////////////////////////////////////////////////////
    // Alloc. memory for storing cluster information at cls.clusinfo
    // --
    // cls.clusinfo[i]       : i = 0 ~ order  (alloc "order+1" memory)
    // cls.clusinfo[i][j]    : j = 0          (alloc "   1   " memory)
    // cls.clusinfo[i][j][k] : k = 0          (alloc "   1   " memory)
    // --
    //  >> cls.clusinfo[i][0][0] = 1          
    //  (save data as "1", 
    //   meaning that there are "1-1" number of cluster)
    ///////////////////////////////////////////////////////////////////
    Cluster_allocClusinfo(cls, order); 

    ///////////////////////////////////////////////////////////////////
    // If order is larger than 1, 
    // the following variables are required for creating clusters
    ///////////////////////////////////////////////////////////////////

    // Connectivity map : 
    // To know whether the distance between two spins within rdip
    // cmap[nspin][nspin] = 1 if connected 
    //                     = 0 else 
    int** cmap = NULL; 

    // Strength map : 
    // To know the strength between two spins 
    // stmap[nspin][nspin] = strength (radkHz)
    float** stmap = NULL; 

    // Sparse map : 
    // For cluster algorithm, cmap is shrinked if an element is 0
    // spmap[nspin][ncol], ncol = the number of connected spin + 1
    // spmap[ispin][0] = the number of connected spin + 1
    // spmap[ispin][k] = k-th connected spin index 
    int** spmap = NULL; 

    // Cutoff for second order cluster
    float rdip = Config_getRdip(config);
    float rdipcut = Config_getRdipcut(config);

    ///////////////////////////////////////////////////////////////////
    // After clusterizing, 
    // the information would stored at cls.clusinfo
    ///////////////////////////////////////////////////////////////////
    if (strcasecmp(method, "cce")==0 
            || strcasecmp(method, "gcce")==0
            || strcasecmp(method, "mecce")==0 || strcasecmp(method, "megcce")==0){

        Cluster_setClusinfo_0th(cls);    // if order == 0
        Cluster_setClusinfo_1th(cls,ba); // if order == 1

        if (order == 2 && !Cluster_getAddsubclus(cls) && !verbosity){
            // The same list without the hash and the dense maps, the pair search shared
            // between the ranks (A21); higher orders still need the maps below.
            // A38: with localRequest, only this rank's rows (and the cost split, decided from all of them first).
            clusterizePairs(cls, ba, rdip, rdipcut,
                            cls->localRequest && qa && dfa ? Cluster_costFunction(config,qa,ba,dfa) : nullptr);
        }
        else if (order > 1){             // if order >= 2
            // Get connectivity map 
            // Args Returns : cmap, stmap, spmap
            BathArray_connectivity(&cmap, &stmap, ba, rdip, rdipcut);
            makeSparsemap(&spmap, cmap, nspin);

            // Clusterize for order > 1 
            // and Solve tilde (Division) of coherence function
            // Args Returns : cls->clusinfo
            clusterizeHash(cls, nspin, spmap, stmap);

            freeInt2d(&cmap,nspin);
            freeInt2d(&spmap,nspin);
            freeFloat2d(&stmap,nspin); 

        }
    } 
    else if (strcasecmp(method, "pcce")==0){
        // Clusterize for all order 
        // and Solve tilde (Division) of coherence function
        // Args Returns : cls->clusinfo (new address)
        clusterizePcce(cls, ba, config); 
    }
    else if (  strcmp(method, "dsj")==0 
             || strcmp(method, "itb")==0
             || strcmp(method, "dsjitb")==0 ){
        // In disjoint clustering method, 
        // we will obtain clustered disjointly clustered spins
        // or interbath or both.
        BathArray_connectivity(&cmap, &stmap, ba, rdip, rdipcut);
        clusterizeDsjItb(cls, &cmap, &stmap, ba);
        freeInt2d(&cmap,nspin);
        freeInt2d(&spmap,nspin);
    }
    else{
        printf("Error: clusterize: method is not defined\n");
        exit(1);
    }

    //if (strcasecmp(method, "gcce")!=0 and strcasecmp(method, "cce")!=0){
    //    // cls -> hash
    //    HashCluster* hashClusters = NULL;
    //    convertClusinfoToHash(&hashClusters, cls);
    //    Cluster_freeClusinfo(cls);

    //    // solve tilde, hash -> cls
    //    solveTilde(&hashClusters, cls, nspin);
    //    freeHashCluster(&hashClusters, order);
    //}

    // the 0th cluster is kept for the methods with the qubit in the cluster space (gCCE, megCCE)
    if (strcasecmp(method, "gcce")!=0 && strcasecmp(method, "megcce")!=0){
        Cluster_setClusinfo_chgiter(cls, 0, 0, 0);
    }
}


// free
void Cluster_freeAll(Cluster* cls){
    if (cls == NULL) return;
    Cluster_freeClusinfo(cls);
    Cluster_freeNk(cls);
    freeArray1d((void**)&cls);
}

/* Low Level --------------------------------------------------------*/

// report
void Cluster_report(Cluster* cls){

    printTitle("Structure Cluster");

    printSubTitle("General properties");
    printStructElementInt("order",Cluster_getOrder(cls));
    printStructElementChar("method",cls->method);
    printStructElementBool("addsubclus",Cluster_getAddsubclus(cls));

    if (strcasecmp(cls->method,"pcce")==0){
        printSubTitle("Conduct pCCE method");
        printStructElementInt("max_trial",Cluster_getMax_trial(cls));
        printStructElementInt("max_iter",Cluster_getMax_iter(cls));
        printStructElementBool("kmeans_pp",Cluster_getKmeans_pp(cls));
        printStructElementBool("iter_detail",Cluster_getIter_detail(cls));
    }
    
    if (cls->nk!=NULL){
        printSubTitle("The number of clusters for each order");
        Cluster_reportNk(cls);
    }
    
    if (cls->clusinfo!=NULL){
        printSubTitle("Clustering information");
        Cluster_reportClusinfo(cls);
    }

    printf("\n");
}

void Cluster_reportNk(Cluster* cls){
    int order = Cluster_getOrder(cls);
    printStructElementInt1dIdx("nk",Cluster_getNk(cls),order+1);
}

void Cluster_reportClusinfo(Cluster* cls){

    int order = Cluster_getOrder(cls);
    printf("      %-18s:   \n\n","clusinfo");

    // 0-th order
    int iter_0th = Cluster_getClusinfo_itercluster(cls,0,0)[0];
    int ncluster_0th = Cluster_getClusinfo_ncluster(cls,0);
    printf("%10s 0-th order ... \n"," ");
    printf("%10s Cluster[0][%5d] : %3d (Iter)\n"," ",0,iter_0th);
    printf("\n%10s Total 0-Cluster # :  %d # \n\n\n"," ",ncluster_0th);

    for (int i = 1; i <= order; i++){
        printf("\n%10s %d-th order ... \n"," ",i);
        int ncluster = Cluster_getClusinfo_ncluster(cls, i);

        for (int j = 1; j<= ncluster; j++){

            if (verbosity || (j<5 || j>ncluster-4)){
                
                int* itercluster = Cluster_getClusinfo_itercluster(cls,i,j);
                int iter = itercluster[0];

                printf("%10s Cluster[%d][%5d] : %3d ["," ",i,j,iter);
                for (int k = 1; k <= i; k++){
                    printf("%-3d",itercluster[k]);
                    if (k != i){ printf(", ");} else { printf(" ]\n");}
                }
            }
            
            if (!verbosity && j==5){
                printf("%10s      : \n"," ");
            }
        }
        printf("\n%10s Total %d-Cluster # :  %d # \n\n\n"," ",i,ncluster);
    }
}

void reportClusinfo(int*** clusinfo, int order){
 
    // 0-th order
    int iter_0th = clusinfo[0][0][0];
    int ncluster_0th = 1;
    printf("%10s 0-th order ... \n"," ");
    printf("%10s Cluster[0][%5d] : %3d (Iter)\n"," ",0,iter_0th);
    printf("\n%10s Total 0-Cluster # :  %d # \n\n\n"," ",ncluster_0th);

    for (int i = 1; i <= order; i++){
        printf("\n%10s %d-th order ... \n"," ",i);
        int ncluster = clusinfo[i][0][0]-1;

        for (int j = 1; j<= ncluster; j++){

            if (verbosity || (j<5 || j>ncluster-4)){
                
                int* itercluster = &(clusinfo[i][j][0]);
                int iter = itercluster[0];

                printf("%10s Cluster[%d][%5d] : %3d ["," ",i,j,iter);
                for (int k = 1; k <= i; k++){
                    printf("%-3d",itercluster[k]);
                    if (k != i){ printf(", ");} else { printf(" ]\n");}
                }
            }
            
            if (!verbosity && j==5){
                printf("%10s      : \n"," ");
            }
        }
        printf("\n%10s Total %d-Cluster # :  %d # \n\n\n"," ",i,ncluster);
    }
}

/* Low Level --------------------------------------------------------*/

// alloc
void Cluster_allocNk(Cluster* cls){ // 0:0-th ... k:k-th order 
    int length = cls->order + 1;
    cls->nk = allocInt1d(length);
}

void Cluster_allocClusinfo(Cluster* cls, int order){
    cls->clusinfo = (int***)allocArray1d(order+1,sizeof(int**));

    for (int i=0; i<=order; i++){
        cls->clusinfo[i] = (int**)allocArray1d(1,sizeof(int*));
        cls->clusinfo[i][0] = (int*)allocArray1d(1,sizeof(int));
        cls->clusinfo[i][0][0] = 1;
    }
}

void Cluster_setClusinfo_0th(Cluster* cls){
    // Actually do nothing
    if (cls->order==0){
        if (rank==0){printf("\n\t Clustering 0th order ... \n");}
    }
}

void Cluster_setClusinfo_1th(Cluster* cls, BathArray* ba){
    //////////////////////////////////////////////
    // Simply insert bathspin index in cls.clusinfo
    // Faster than using Hash algorithm
    //////////////////////////////////////////////
    if (cls->order == 1){ 
        if (rank==0){printf("\n\t Clustering 1st order ... \n");}
        int cluster[1] = {0};
        int iter = 1;
        for (int spinidx = 0; spinidx < ba->nspin; spinidx++){
            cluster[0] = spinidx;
            Cluster_setClusinfo_addcluster(cls, 1, iter, cluster);
        }
    }
}

int Cluster_setClusinfo_addcluster(Cluster* cls, int order, int iter, int* cluster){

    if (order==0){
        printf("Error: Cluster_setClusinfo_addcluster: at 0-th order, you cannot add cluster\n");
        exit(1);
    }

    if (cls->clusinfo==NULL){
        printf("Error: Cluster_setClusinfo_addcluster: clusinfo is not allocated\n");
        exit(1);
    }

    if (iter != 1 ){
        printf("Warning: Cluster_setClusinfo_addcluster: iter is not 1 when you add new cluster ... \n");
    }

    int ncluster_old = cls->clusinfo[order][0][0];
    int ncluster_new = ncluster_old + 1;
    int ic = ncluster_new - 1; // new cluster index

    // set the number of cluster at order
    cls->clusinfo[order][0][0] = ncluster_new;

    // realloc
    reallocInt2d(&(cls->clusinfo[order]),ncluster_old,ncluster_new, order+1);

    // set new cluster
    copyInt1d(&(cls->clusinfo[order][ic][1]),cluster,order);

    // set iter at 0-th element (It is the number how many it will calculated)
    cls->clusinfo[order][ic][0] = iter;

    return ic;
}

void Cluster_setClusinfo_chgcluster(Cluster* cls, int order, int* cluster, int ic){
    if (order==0){
        printf("Error : Cluster_setClusinfo_chgcluster : at 0-th order, you cannot change cluster information\n");
        printf("For the 0th order, you can only change iter\n");
        exit(1);
    }
    copyInt1d(&(cls->clusinfo[order][ic][1]),cluster,order);
}

void Cluster_setClusinfo_chgiter(Cluster* cls, int order, int iter, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_setClusinfo_chgiter : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    cls->clusinfo[order][ic][0] = iter;
}


int Cluster_getClusinfo_ncluster(Cluster* cls, int order){
    if (order==0){
        return 1;
    }
    return cls->clusinfo[order][0][0]-1;
}

int* Cluster_getClusinfo_itercluster(Cluster* cls, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return cls->clusinfo[order][ic];
}

int Cluster_getClusinfo_iter(Cluster* cls, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return cls->clusinfo[order][ic][0];
}

int* Cluster_getClusinfo_cluster_copy(Cluster* cls, int order, int ic){
    if (order==0){
        printf("Error : Cluster_getClusinfo_cluster_copy : at 0-th order, you cannot get cluster information\n");
        exit(1);
    }
    int* cluster = allocInt1d(order);
    copyInt1d(cluster,&(cls->clusinfo[order][ic][1]),order);
    return cluster;
}

int*** Cluster_getClusinfo(Cluster* cls){
    return cls->clusinfo;
}

// get 

char* Cluster_getMethod(Cluster* cls){
    return cls->method;
}

int Cluster_getOrder(Cluster* cls){
    return cls->order;
}

bool Cluster_getAddsubclus(Cluster* cls){
    return cls->addsubclus;
}

int* Cluster_getNk(Cluster* cls){
    return cls->nk;
}

int Cluster_getNk_order(Cluster* cls, int i){
    if (i<cls->order && i>=0){
        return cls->nk[i];
    }else{
        printf("Error: Cluster_getNk_order: i is out of range\n");
        exit(1);
    }
}

int Cluster_getSk(Cluster* cls){
    return cls->sK;
}

int Cluster_getMax_iter(Cluster* cls){
    return cls->max_iter;
}

int Cluster_getMax_trial(Cluster* cls){
    return cls->max_trial;
}

bool Cluster_getKmeans_pp(Cluster* cls){
    return cls->kmeans_pp;
}

bool Cluster_getIter_detail(Cluster* cls){
    return cls->iter_detail;
}

// set 
void Cluster_setOrder(Cluster* cls, int order){
    cls->order = order;
}

void Cluster_setMethod(Cluster* cls, char* method){
    strcpy(cls->method, method);
}
void Cluster_setAddsubclus(Cluster* cls, bool addsubclus){
    cls->addsubclus = addsubclus;
}
void Cluster_setNk(Cluster* cls, int* nk){
    copyInt1d(cls->nk,nk,cls->order+1);
}

void Cluster_setSk(Cluster* cls, int sK){
    cls->sK = sK;
}                                                           

void Cluster_setMax_trial(Cluster* cls, int max_trial){
    cls->max_trial = max_trial;
}

void Cluster_setMax_iter(Cluster* cls, int max_iter){
    cls->max_iter = max_iter;
}

void Cluster_setKmeans_pp(Cluster* cls, bool kmeans_pp){
    cls->kmeans_pp= kmeans_pp;
}

void Cluster_setIter_detail(Cluster* cls, bool iter_detail){
    cls->iter_detail= iter_detail;
}

// free
void Cluster_freeNk(Cluster* cls){
    freeInt1d(&(cls->nk));
}


void Cluster_freeClusinfo(Cluster* cls){
    
    int order = cls->order;

    if (cls->clusinfo!=NULL){
        for (int i=0; i<=order; i++){

            int ncluster = Cluster_getClusinfo_ncluster(cls,i);

            if (i==0){
                freeArray2d((void***)&(cls->clusinfo[i]),1);
            }
            else{
                // Includes row zero, which stores the number of clusters.
                freeArray2d((void***)&(cls->clusinfo[i]),ncluster+1);
            }
        }
        free(cls->clusinfo);
        cls->clusinfo = NULL;
    }
    free(cls->pairRows);
    free(cls->pairLoads);
    cls->pairRows = NULL;
    cls->npairRows = 0;
    cls->pairLoads = NULL;
    cls->pairsLocal = false;
}
