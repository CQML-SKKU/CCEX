#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <climits>
#include <algorithm>
#include <vector>
#include "../../include/clustering/cluster_hash.h"
#include "../../include/clustering/cluster_pcce.h"
#include "../../include/common/memory.h"
#include "../../include/clustering/cluster.h"
#include "../../include/physical/bath.h"
#include "../../include/calculation/hamiltonian.h"
#include "../../include/clustering/pcce_partition_cache.h"


void clusterizePcce(Cluster* cls, BathArray* ba, Config* config){

    ////////////////////////////////////////////////////
    // About pCCE method :
    //  pCCE(N,K) indicates a clustering method of 
    //  the N-clusters having K number of bath spins
    ////////////////////////////////////////////////////
    //
    ////////////////////////////////////////////////////
    // pCCE related information : 
    //   sK, ncenter, pcce_nspin, rest_nspin 
    Partition_info pinfo;
    int totnspin = BathArray_getNspin(ba);

    if (cls->sK < 1 || cls->sK > totnspin){
        if (rank==0) fprintf(stderr,"Error: pcce.partition_size must be between 1 and the bath spin count (%d)\n",totnspin);
        MPI_Abort(MPI_COMM_WORLD,1);
        return;
    }
    if (cls->order < 1 || cls->order > (INT_MAX-1)/cls->sK){
        if (rank==0) fprintf(stderr,"Error: pCCE order * pcce.partition_size must be a positive integer below INT_MAX\n");
        MPI_Abort(MPI_COMM_WORLD,1);
        return;
    }

    pinfo.sK         = cls->sK;
    pinfo.ncenter    = int(totnspin / (pinfo.sK));
    pinfo.pcce_nspin = int((pinfo.ncenter) * (pinfo.sK));
    pinfo.rest_nspin = int(totnspin % (pinfo.sK));

    // Include the original ordered selection, including any excluded remainder.
    const std::string cacheIdentity = pccePartitionCacheIdentity(ba, &pinfo, cls, config);

    ////////////////////////////////////////////////////
    // Cut the # of spins (pcce bundleing)
    // Total number of spins : ba.nspin -> NK
    ////////////////////////////////////////////////////
    shrink_restspins(ba, pinfo.rest_nspin);

    ////////////////////////////////////////////////////
    // Find the center position & assigned index 
    // (Find a representative position of a cluster 
    //  to be included k spins)
    //  Arg return : best_centers, best_assigned_idx
    ////////////////////////////////////////////////////
    int ncenter              = pinfo.ncenter; 
    int pcce_nspin           = pinfo.pcce_nspin; 
    Point* best_centers      = (Point*) calloc(ncenter, sizeof(Point));
    int*   best_assigned_idx = (int*)   calloc(pcce_nspin, sizeof(int));
    pccePartitionCached(ba->bath, &pinfo, cls, best_centers, best_assigned_idx, config, cacheIdentity);
    ////////////////////////////////////////////////////
    // Clusterize for the center
    ////////////////////////////////////////////////////
    int order = Cluster_getOrder(cls);

    Cluster_setClusinfo_0th(cls);

    if (order == 1){ // For pCCE, + k=1 (order = n)
        
        if (rank==0){printf("\n\t Clustering 1st order ... \n");}
    
        int cluster[1] = {0};
        int iter = 1; // clusinfo[N][i][0]
        // center >> cls->info for o1
        for (int centeridx = 0; centeridx < ncenter; centeridx++){
            cluster[0] = centeridx; 
            Cluster_setClusinfo_addcluster(cls, order, iter, cluster);
        }
    }
    else if (order > 1){ 
    
        int** cmap    = NULL;  // connectivity map : cmap[nspin][nspin] = 1 if connected else 0
        float** stmap = NULL;  // strength map : stmap[nspin][nspin] = strength
        int** spmap   = NULL;  // sparse map : spmap[nspin][ncol] = n connected spin + 1d array of connected spins index
    
        float rdip    = Config_getRdip(config);
        float rdipcut = Config_getRdipcut(config);
    
        // center >> cls->info for on
        if (methodUsesPcceClusters(cls->method)){
            const bool spins = !strcasecmp(config->pcce_connections,"spins");
            if (spins) BathArray_connectivity_pcce_spins(&cmap, &stmap, ba, best_assigned_idx, ncenter, rdip, rdipcut);
            else       BathArray_connectivity_pcce(&cmap, &stmap, best_centers, rdip, rdipcut, ncenter);
            if (rank==0){
                long npair = 0;
                for (int a=0; a<ncenter; a++) for (int b=a+1; b<ncenter; b++) npair += cmap[a][b];
                printf("\t pCCE connections : %s (rdip=%g, rdipcut=%g): %ld of %ld partition pairs connected\n",
                       spins ? "spins" : "centers", rdip, rdipcut, npair, (long)ncenter*(ncenter-1)/2);
            }
            makeSparsemap(&spmap, cmap, ncenter);
            // main
            clusterizeHash(cls, ncenter, spmap, stmap); // clusterize + solve tilde
        }
        else{
            printf("Error: clusterize: method is not defined\n");
            exit(1);
        }
        freeInt2d(   &cmap, ncenter);
        freeInt2d(  &spmap, ncenter);
        freeFloat2d(&stmap, ncenter); 
    
    }else{
        fprintf(stderr,"Error: Cluster_clusterize: order(%d) is not defined\n",order);
        exit(1);
    }

    // ========================== //
    // Setting cls_pcce !! <= cls //
    // ========================== //
    int** cs2dArr = Cluster_setCenterIdx_spinIdx_2dArr(ncenter, cls->sK, pcce_nspin, best_assigned_idx);
    //if (rank == 0){
    //    for (int tempi=0; tempi<ncenter; tempi++){
    //        printf("cs2dArr[%d] ", tempi);
    //        for (int tempj=0; tempj<cls->sK; tempj++){
    //            printf("%d ", cs2dArr[tempi][tempj]);
    //        }
    //        printf("\n");
    //    }
    //}
    //printf("\n");
    
    int*** pcceClusInfo = convert_centerIdx_to_spinIdx(cls, cls->order, cls->sK, cs2dArr);
    freeInt2d(&cs2dArr, ncenter);
    free(best_centers);
    free(best_assigned_idx);
    
    int spinOrder = (cls->order * cls->sK);
    // nk and clusinfo must use the same order after expanding bundle indices.
    int* spinNk = allocInt1d(spinOrder+1);
    for (int i=0; i<=cls->order; ++i) spinNk[i*cls->sK] = cls->nk[i];
    Cluster_freeNk(cls);
    cls->nk = spinNk;
    Cluster_freeClusinfo(cls);
    cls->clusinfo = pcceClusInfo;
    Cluster_setOrder(cls, spinOrder);

    //printf("! ================= !\n");
    //printf("! -- cls print --!!\n\n");
    //reportClusinfo(cls->clusinfo, 4);
    //printf("\n! ================= !\n");
}

// Find the Center position & assigned index 
// arg return : center position(best_centers) and related info (pinfo, best_assigned_idx).
void simulator_cluster_partition(BathSpin** bath        , Partition_info* pinfo             , Cluster* cls,
                                 Point*     best_centers, int*            best_assigned_idx, const Config* config){
    
    clock_t start     = clock();
    char* inputfile   = NULL;
    //////////////////////////////////////////////////////////
    ////            !!  Setting MPI process  !!             //
    //////////////////////////////////////////////////////////
    int max_trial       = cls      -> max_trial;
    int trials_per_proc = max_trial / nprocess;
    if (trials_per_proc == 0){trials_per_proc = 1;}
    const bool deterministic = !strcasecmp(config->pcce_partition,"deterministic");

    //int start_trial     = rank      * trials_per_proc;
    //int end_trial       = (rank == nprocess - 1) ? max_trial : start_trial + trials_per_proc;

    //for (int i=0; i<nprocess; i++){
    //    if (rank == i){
    //        printf(" || Rank: %d\n", rank);
    //        printf(" || trials_per_proc: %d\n", trials_per_proc);
    //        //printf("start_trial: %d\n", start_trial);
    //        //printf("end_trial: %d\n", end_trial);
    //        print_BD("=", 55);printf("\n");
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    //}
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    ////        !!  Setting the information array !!        //
    //////////////////////////////////////////////////////////
    int pcce_nspin = pinfo -> pcce_nspin;  // ncenter * sK
    int ncenter    = pinfo -> ncenter;     // nspin / sK (truncated)

    Point* local_best_spins      = (Point*) calloc(pcce_nspin, sizeof(Point));
    Point* local_best_centers    = (Point*) calloc(ncenter,    sizeof(Point));
    int* local_best_assigned_idx = (int*)   calloc(pcce_nspin, sizeof(int));
    int* local_best_shuffle_idx  = (int*)   calloc(pcce_nspin, sizeof(int));

    // Reuse O(ncenter) assignment scratch instead of storing all N*K distances.
    PcceTrialWorkspace trial_workspace(pcce_nspin,ncenter);
    // Converged, cycle, limit-only trials; logical, executed, skipped iterations.
    unsigned long long local_search[6]={}, global_search[6]={};
    
    //////////////////////////////////////////////////////////
    ////                     !!  Run !!                     //
    //////////////////////////////////////////////////////////
    double local_best_inertia =  INFINITY;
    double local_best_sil     = -INFINITY;
    int    local_best_trial   = INT_MAX;
    int local_trials = 0;
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    // Global trials are independent of the number of ranks. Some ranks can
    // legitimately have no work; all still participate in the winner selection.
    for (int64_t trial64=deterministic ? rank : 0;
         trial64<(deterministic ? max_trial : trials_per_proc);
         trial64+=deterministic ? nprocess : 1){
        const int trial=int(trial64);
        PcceRandom random(config->seed,trial);
        PcceRandom* rng=deterministic ? &random : nullptr;
        ++local_trials;


        Point* local_spins          = (Point*) calloc(pcce_nspin, sizeof(Point));
        Point* local_centers        = (Point*) calloc(ncenter,    sizeof(Point));
        int* local_assigned_idx     = (int*)   calloc(pcce_nspin, sizeof(int));
        int* local_shuffle_idx      = (int*)   calloc(pcce_nspin, sizeof(int));

        copy_bath_to_point_type(bath, local_spins, pcce_nspin);

        for (int i=0; i<pcce_nspin; i++){ 
            local_shuffle_idx[i] = i;
        }

        // Shuffle "spins" & initialize the "local_centers". // 
        shuffle_spins(local_spins, local_shuffle_idx, pcce_nspin, rng);
        choose_initial_method(cls->kmeans_pp, local_spins, local_centers, pcce_nspin, ncenter, rng);

        const PcceTrialResult result=pcceOptimizeTrial(local_spins,local_centers,
            local_assigned_idx,pcce_nspin,ncenter,cls->max_iter,trial_workspace);
        const int status=result.converged ? 0 : result.cycle_period ? 1 : 2;
        ++local_search[status];
        local_search[3]+=result.iterations;
        local_search[4]+=result.executed;
        local_search[5]+=result.skipped;
        if(cls->iter_detail && rank==0 && local_trials<=config->pcce_detail_limit)
            printf("[partition-iteration] rank=%d trial=%d status=%s iterations=%d executed=%d skipped=%d period=%d\n",
                   rank,trial,status==0 ? "converged" : status==1 ? "cycle" : "limit",
                   result.iterations,result.executed,result.skipped,result.cycle_period);

        // Find the local optimization centers configuration. // 
        double local_inertia = calculateInertia(local_spins, local_centers, local_assigned_idx, pcce_nspin, ncenter);
        if (local_inertia < local_best_inertia){

            double local_sil = silhouetteCoefficient(local_spins, local_assigned_idx, local_centers, pcce_nspin, ncenter);
            //printf("\n\n\nLOCAL & BEST: %.2lf & %.2lf\n\n\n", local_sil, local_best_sil);
            local_best_trial   = trial;
            local_best_inertia = local_inertia;
            local_best_sil     = local_sil;

            copy_best_centers(local_centers, local_best_centers, ncenter);
            copy_best_assigned(local_assigned_idx, local_best_assigned_idx, pcce_nspin);
            copy_best_idxList(local_shuffle_idx, local_best_shuffle_idx, pcce_nspin);
            copy_best_spins(local_spins, local_best_spins, pcce_nspin);
        }
        //else if (local_inertia < local_best_inertia) 
        else if (local_inertia == local_best_inertia) {
            double local_sil = silhouetteCoefficient(local_spins, local_assigned_idx, local_centers, pcce_nspin, ncenter);
            if (local_sil > local_best_sil){
                local_best_trial   = trial;
                local_best_inertia = local_inertia;
                local_best_sil     = local_sil;

                copy_best_centers(local_centers, local_best_centers, ncenter);
                copy_best_assigned(local_assigned_idx, local_best_assigned_idx, pcce_nspin);
                copy_best_spins(local_spins, local_best_spins, pcce_nspin);
                copy_best_idxList(local_shuffle_idx, local_best_shuffle_idx, pcce_nspin);
                //printf("\n\n\nLOCAL & BEST: %.2lf & %.2lf\n\n\n", local_sil, local_best_sil);
            }
        }

        ////////////////////////////////////////////////////////
        //clock_t trial_end = clock();
        //for (int i=0; i < nprocess; i++){
        //    if (rank==i){
        //        print_optimize_kmeans_info(trial, iteration, local_inertia, local_best_inertia, local_best_sil, rank);
        //        printf("\n||---------------- Optimization Time ----------------||\n\n");
        //        printf("                 "); print_time(trial_start, trial_end);
        //        print_BD("=", 55); printf("\n");
        //        fflush(stdout);
        //    }
        //MPI_Barrier(MPI_COMM_WORLD);
        //}
        ////////////////////////////////////////////////////////
        free(local_spins);
        free(local_centers);
        free(local_assigned_idx);
        free(local_shuffle_idx);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //for (int i=0; i<nprocess; i++){
    //    if (rank == i){
    //        printf("\n||-------------------  Rank : %d | Best Local -------------------||\n\n", rank);
    //        print_BD("=", 55);
    //        print_centers(ncenter, local_best_centers, rank);
    //        print_assign_idx(local_best_assigned_idx, pcce_nspin, rank);
    //        print_spins(local_best_spins, pcce_nspin, rank);
    //        print_BD("=", 55);
    //        fflush(stdout);
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    //}

    ////////////////////////////////////////////////////////
    //  !!  Setting the "All & Global best" varibles  !!  //
    ////////////////////////////////////////////////////////
    double*  all_best_inertia = (double*)calloc(nprocess, sizeof(double));
    double*  all_best_sil     = (double*)calloc(nprocess, sizeof(double));
    int*     all_best_trial   = (int*)   calloc(nprocess, sizeof(int));

    int    best_rank           =  0;
    int    global_best_trial   = INT_MAX;
    double global_best_inertia =  INFINITY;
    double global_best_sil     = -INFINITY;

    ////////////////////////////////////////////////////////
    //   Gathering local best parameters in each procs.   //
    ////////////////////////////////////////////////////////
    MPI_Gather(&local_best_inertia, 1, MPI_DOUBLE, all_best_inertia, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_best_sil    , 1, MPI_DOUBLE, all_best_sil    , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_best_trial  , 1, MPI_INT   , all_best_trial  , 1, MPI_INT   , 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////
    //           !!  Find the best rank index !!          //
    ////////////////////////////////////////////////////////
    if (rank == 0) {
        for (int i = 0; i < nprocess; i++){
            if (pcce_better_partition(all_best_inertia[i], all_best_sil[i],
                                      global_best_inertia, global_best_sil,deterministic) ||
                (deterministic && all_best_inertia[i]==global_best_inertia &&
                 all_best_sil[i]==global_best_sil && all_best_trial[i]<global_best_trial)) {

                best_rank = i;
                update_global_best(i, all_best_inertia, all_best_sil, all_best_trial, &global_best_inertia, &global_best_sil, &global_best_trial);
            }
        
        }

        //printf("\n"); print_BD("=", 55); printf("\n");
        //printf("Best rank: %d, Iner.: %.2lf, Sil.: %.2lf, Trial: %d\n", best_rank, global_best_inertia, global_best_sil, global_best_trial);
        //printf("\n"); print_BD("=", 55); printf("\n");
    }
    free(all_best_inertia);
    free(all_best_sil);
    free(all_best_trial);

    int total_trials=0;
    MPI_Reduce(&local_trials,&total_trials,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(local_search,global_search,6,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank==0){
        if (!std::isfinite(global_best_inertia) || !std::isfinite(global_best_sil)){
            fprintf(stderr,"Error: pCCE partition search produced no finite candidate\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        if(cls->iter_detail && local_trials>config->pcce_detail_limit)
            printf("[partition-detail] rank=0 suppressed=%d limit=%d; global totals follow\n",
                   local_trials-config->pcce_detail_limit,config->pcce_detail_limit);
        printf("[partition] mode=%s total_trials=%d best_trial=%d best_rank=%d\n",
               config->pcce_partition,total_trials,global_best_trial,best_rank);
        printf("[partition-search] converged=%llu cycles=%llu limit=%llu iterations=%llu executed=%llu skipped=%llu\n",
               global_search[0],global_search[1],global_search[2],global_search[3],global_search[4],global_search[5]);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    ////////////////////////////////////////////////////////
    //     !!  Send & Receive the best rank resuts !!     //
    ////////////////////////////////////////////////////////
    // Chnage the indecies for "local_spins => ba->bath"
    Point* global_best_centers              = local_best_centers;

    Point* global_best_inverse_spins        = (Point*)malloc(pcce_nspin*sizeof(Point));
    int*   global_best_inverse_assigned_idx = (int*)calloc(pcce_nspin, sizeof(int));
    int*   global_best_inverse_shuffle_idx  = (int*)calloc(pcce_nspin, sizeof(int)); 
    changeIdx_inverse_shuffle_idx(local_best_shuffle_idx, global_best_inverse_shuffle_idx, pcce_nspin);
    
    if (rank == best_rank){
        for (int i=0; i<pcce_nspin; i++){
            int inverse_idx                               = local_best_shuffle_idx[i];
            global_best_inverse_assigned_idx[inverse_idx] = local_best_assigned_idx[i];
            global_best_inverse_spins[inverse_idx]        = local_best_spins[i];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(global_best_inverse_spins       , pcce_nspin* 3, MPI_DOUBLE, best_rank,  MPI_COMM_WORLD);
    MPI_Bcast(global_best_centers             , ncenter  * 3 , MPI_DOUBLE, best_rank,  MPI_COMM_WORLD);
    MPI_Bcast(global_best_inverse_assigned_idx, pcce_nspin   , MPI_INT   , best_rank,  MPI_COMM_WORLD);

    for (int i=0; i<ncenter; i++){
        best_centers[i].x = global_best_centers[i].x; 
        best_centers[i].y = global_best_centers[i].y; 
        best_centers[i].z = global_best_centers[i].z; 
    }
    for (int i=0; i<pcce_nspin; i++){
        best_assigned_idx[i] = global_best_inverse_assigned_idx[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //if (rank == 0){
    //    printf("\n||-------------------------------------------------------||\n");
    //    printf("\n||----------------! FINISH OPTIMIZATION !----------------||\n");
    //    printf("\n||-------------------------------------------------------||\n\n");
    //    printf("\n||-----------------! Center Information !----------------||\n\n");
    //    print_centers(ncenter, best_centers, rank);
    //    printf("\n||--------------! Assigned Index for Center !------------||\n");
    //    print_assign_idx(best_assigned_idx, pcce_nspin, rank);
    //    print_BD("=", 55);
    //}
    //MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    ////      !!  Close the dynamical allocation !!         //
    //////////////////////////////////////////////////////////
    FILE *c_savefile; 
    FILE *s_savefile; 
    
    //if (rank == 0){
    //    //write_centers(&c_savefile, &(ba->centers), ncenter);
    //    //write_spins(&s_savefile, &global_best_inverse_spins, &(ba->assigned_idx), pcce_nspin);
    //    //write_centers(&c_savefile, &best_centers, ncenter);
    //    //write_spins(&s_savefile, &global_best_inverse_spins, &best_assigned_idx, pcce_nspin);
    //    print_best_centers(&best_centers, ncenter);
    //    print_best_spins(&global_best_inverse_spins, &best_assigned_idx, pcce_nspin);
    //}

    free(global_best_centers);
    free(global_best_inverse_spins);
    free(global_best_inverse_assigned_idx);
    free(global_best_inverse_shuffle_idx);

    free(local_best_spins);
    free(local_best_assigned_idx);
    free(local_best_shuffle_idx);

    ////////////////////////////////////////////////////////
    //                   !!  Finish  !!                   //
    ////////////////////////////////////////////////////////
    //clock_t end = clock();

    //if (rank==0) {
    //    printf("\n||--------------- Total Execution Time ---------------||\n\n");
    //    printf("                 "); print_time(start, end);
    //    print_BD("=", 55); printf("\n");
    //}   
}

void shrink_restspins(BathArray* ba, int rest_nspin){

    if (rest_nspin < 0 || rest_nspin > ba->nspin){
        fprintf(stderr,"Error: invalid pCCE remainder spin count\n");
        MPI_Abort(MPI_COMM_WORLD,1);
        return;
    }
    
    // Destroy the spins that fall off the end BEFORE shrinking the pointer array. The old
    // code just shrank the array and dropped them: they were already leaked, and now that
    // BathSpin is a properly constructed C++ object their MatrixXcd members would leak
    // with them. createBathArray deep-copies every spin it takes, so nothing else aliases
    // these -- they are ours to delete.
    int old_nspin = ba->nspin;
    ba->nspin -= rest_nspin;
    for (int i=ba->nspin; i<old_nspin; i++){
        BathArray_freeBath_i_hypf(ba,i);
        delete ba->bath[i];
    }
    BathSpin** temp = (BathSpin**)realloc(ba->bath, ba->nspin * sizeof(BathSpin*));
    // NB: the old check tested ba->bath, not temp -- it could never fire.
    if (temp == NULL && ba->nspin>0){
        perror("Failed to realloc ba->bath\n");
        exit(EXIT_FAILURE);
    }
    ba->bath = temp;
}

///////////////////////////////////////////////////
void update_global_best(int index, double *inertia, double *sil, 
                        int *trial, double *global_best_inertia, double *global_best_sil, int *global_best_trial) {

    *global_best_inertia = inertia[index];
    *global_best_sil     = sil[index];
    *global_best_trial   = trial[index];
}

bool pcce_better_partition(double inertia, double silhouette,
                           double best_inertia, double best_silhouette,bool deterministic){
    return inertia < best_inertia ||
           (inertia == best_inertia && (deterministic ? silhouette > best_silhouette : silhouette < best_silhouette));
}

///////////////////////////////////////////////
//       Util. & Setting the spin number     // 
///////////////////////////////////////////////
double calc_dist(double* arr1, double* arr2){

    double dx = arr1[0] - arr2[0];
    double dy = arr1[1] - arr2[1];
    double dz = arr1[2] - arr2[2];

    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

void print_pcce_info(int bath_nspin, int nqubit, int ncenter, int pcce_nspin, int rest_nspin){
    print_BD("-", 30);
    printf("Total Bath Spin: %d\n", bath_nspin);
    printf("Qubit #: %d\n", nqubit);
    printf("ncenter: %d\n", ncenter);
    printf("pcce_nspin: %d\n", pcce_nspin);
    printf("rest_nspin: %d\n", rest_nspin);
    print_BD("-", 30);
}

///////////////////////////////////////////////
//       Constrained k-means function.       // 
///////////////////////////////////////////////
int if_contain(int* arr, int arr_size, int value) {
    for (int k=0; k < arr_size;  k++){
        if (arr[k] == value) {
            return 1;
        }
    }
    return 0;           
}

// Inititialize the centers
void initializeCentroids(Point* spins, Point* centers, int ncenter, int pcce_nspin, PcceRandom* rng) {
    
    // Make and initialize the Random number index list // 
    int count = 0;
    int* RandomList = (int*) calloc(ncenter, sizeof(int));
    for (int i=0; i<ncenter; i++){
        RandomList[i] = (pcce_nspin+100);
    }
    
    // Random vales => srand //
    // srand(time(NULL)+i*i);
    for (int count=0; count <ncenter;) {
        int random_index = rng ? rng->index(pcce_nspin) : rand() % pcce_nspin;
        //printf("%d", random_index);
    
        if (!if_contain(RandomList, count, random_index)) {
            RandomList[count] = random_index;
            centers[count].x = spins[random_index].x;
            centers[count].y = spins[random_index].y;
            centers[count].z = spins[random_index].z;
            count++;
        }
    }
    //printf("\n\n Finish making the random number index list!!\n\n");
    free(RandomList);
}

void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int pcce_nspin, int ncenter, PcceRandom* rng){
    //printf("||------------------- Start Iter. -------------------||\n\n");
    if (kmeans_pp == true){
        //printf(" - Centroids Initialization: K-Means++ \n\n");
        //printf(" || Centroids Initialization: K-Means++ ||\n\n");
        //print_BD("=", 55);
        kMeansPlusPlus(spins, centers, pcce_nspin, ncenter, rng);
    }
    else {
        //printf(" - Centroids Initialization: Random \n");
        //printf(" || Centroids Initialization: Random ||\n\n");
        //print_BD("=", 55);
        initializeCentroids(spins, centers, ncenter, pcce_nspin, rng);
    }
    //printf("\n");
}

void kMeansPlusPlus(Point *spins, Point *centers, int pcce_nspin, int ncenter, PcceRandom* rng) {

    int index = rng ? rng->index(pcce_nspin) : rand() % pcce_nspin;
    
    centers[0].x = spins[index].x;
    centers[0].y = spins[index].y;
    centers[0].z = spins[index].z;
    //centers[0] = spins[index];  // 첫 번째 중심점 선택

    // Cache the distance to the nearest chosen center. Updating with only the
    // newest center reduces initialization from O(N*K*K) to O(N*K). Keep the
    // same sqrt/square, summation order and rand() calls as the legacy method.
    std::vector<double> nearest(pcce_nspin, INFINITY);
    std::vector<double> dists(pcce_nspin);

    for (int i = 1; i < ncenter; i++) {
        double sum = 0.0;
        // 각 포인트에 대한 최소 거리를 제곱하여 계산
        for (int j = 0; j < pcce_nspin; j++) {
            double d = distance(spins[j], centers[i-1]);
            if (d < nearest[j]) nearest[j] = d;
            d = nearest[j];
            dists[j] = d * d;
            sum += dists[j];
        }

        // 확률적으로 다음 중심점 선택
        double r = (rng ? rng->unit() : (double)rand() / RAND_MAX) * sum;
        sum = 0;
        for (int j = 0; j < pcce_nspin; j++) {
            sum += dists[j];
            if (sum >= r) {
                //centers[i] = spins[j];
                centers[i].x = spins[j].x;
                centers[i].y = spins[j].y;
                centers[i].z = spins[j].z;
                break;
            }
        }
    }
}
void shuffle_spins(Point *array, int* shuffle_idxs, int pcce_nspin, PcceRandom* rng) {
    // 배열의 마지막 요소부터 시작하여 첫 번째 요소 이전까지 반복
    for (int i = pcce_nspin - 1; i > 0; i--) {
        // 0부터 i까지의 무작위 인덱스 선택
        int j = rng ? rng->index(i+1) : rand() % (i + 1);

        Point temp = array[i];
        array[i]   = array[j];
        array[j]   = temp;

        int temp_idx    = shuffle_idxs[i];
        shuffle_idxs[i] = shuffle_idxs[j];
        shuffle_idxs[j] = temp_idx;

    }
}

void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        int j = shuffle_idx[i];
        inverse_shuffle_idx[j] = i;
    }
}

void assignCentroids(int* assigned_cluster_idx, double** distances, int pcce_nspin, int ncenter) {
    int spinsPerCenter = pcce_nspin / ncenter;
    int *centroidCounts = (int *)malloc(ncenter * sizeof(int));
    for (int i = 0; i < ncenter; i++) {
        centroidCounts[i] = 0;
    }

    for (int i = 0; i < pcce_nspin; i++) {
        int minIndex = -1;
        double minDistance = HUGE_VAL;
        for (int j = 0; j < ncenter; j++) {
            if (distances[i][j] < minDistance && centroidCounts[j] < spinsPerCenter) {
                minDistance = distances[i][j];
                minIndex = j;
            }
        }
        if (minIndex == -1) {
            // It may happen that no centroid can take more points, distribute remaining points
            for (int j = 0; j < ncenter; j++) {
                if (centroidCounts[j] < spinsPerCenter) {
                    minIndex = j;
                    break;
                }
            }
        }
        assigned_cluster_idx[i] = minIndex;
        centroidCounts[minIndex]++;
    }
    free(centroidCounts);
}

void calc_dist_to_centers(Point* spins, Point* centers, double** distances, int pcce_nspin, int ncenter) {
    for (int i = 0; i < pcce_nspin; i++) {
        for (int j = 0; j < ncenter; j++) {
            double dx = spins[i].x - centers[j].x;
            double dy = spins[i].y - centers[j].y;
            double dz = spins[i].z - centers[j].z;
            //dists[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            distances[i][j] = (dx * dx + dy * dy + dz * dz);
            //printf("|| DISTANCE[%d][%d] =%lf || \n", i, j, (dx * dx + dy * dy + dz * dz));
        }
    }
}

PcceTrialResult pcceOptimizeTrial(Point* spins,Point* centers,int* assigned,
                                  int nspin,int ncenter,int max_iter,PcceTrialWorkspace& workspace) {
    PcceTrialResult result;
    workspace.cycle.reset();
    while(result.iterations<max_iter) {
        assignCentroidsGrid(assigned,spins,centers,nspin,ncenter,workspace);
        result.converged=updateCentroids(spins,assigned,centers,nspin,ncenter);
        ++result.iterations; ++result.executed;
        // A real convergence takes precedence over cycle detection. In a proven
        // nonconvergent cycle every phase has already passed this same test.
        if(result.converged) break;
        if(!result.cycle_period) {
            result.cycle_period=int(workspace.cycle.observe(centers,assigned));
            if(result.cycle_period) {
                const int remaining=max_iter-result.iterations;
                result.skipped=(remaining/result.cycle_period)*result.cycle_period;
                result.iterations+=result.skipped;
                // Keep the CURRENT full state. Skipping only whole periods lands
                // at precisely the old state; execute the remaining partial period.
                // No RNG is consumed inside the iteration being skipped.
            }
        }
    }
    return result;
}

void assignCentroidsFromPoints(int* assigned, const Point* spins, const Point* centers,
                               int pcce_nspin, int ncenter, int* counts, int* active) {
    const int capacity = pcce_nspin / ncenter;
    for (int c=0; c<ncenter; ++c) { counts[c]=0; active[c]=c; }
    int available=ncenter;
    for (int i=0; i<pcce_nspin; ++i) {
        // Keep ascending center IDs: strict < must select the same first center
        // in a distance tie as the original full scan. A full center cannot win.
        int best=0;
        double minimum=HUGE_VAL;
        for (int a=0; a<available; ++a) {
            const int c=active[a];
            const double dx=spins[i].x-centers[c].x;
            const double dy=spins[i].y-centers[c].y;
            const double dz=spins[i].z-centers[c].z;
            const double squared=dx*dx+dy*dy+dz*dz;
            if (squared<minimum) { minimum=squared; best=a; }
        }
        // best=0 also preserves the old first-nonfull fallback when all
        // distances are infinite/NaN. Production N is exactly capacity*K.
        const int center=active[best];
        assigned[i]=center;
        if (++counts[center]==capacity) {
            --available;
            std::move(active+best+1,active+available+1,active+best);
        }
    }
}

// assignCentroidsFromPoints gives each spin, in order, the open center with the smallest
// squared distance, the smallest center ID among equal distances. Here the open centers
// sit in grid cells; rings of cells around the spin are searched until no unsearched
// cell can hold a center at that distance or less, and the squared distance of every
// center looked at is the same expression as there. Pruning only uses a lower bound
// with a margin far above the rounding of the cell geometry and of that expression.
void assignCentroidsGrid(int* assigned, const Point* spins, const Point* centers,
                         int pcce_nspin, int ncenter, PcceTrialWorkspace& w) {
    const int capacity = pcce_nspin / ncenter;
    double lo[3]={INFINITY,INFINITY,INFINITY}, hi[3]={-INFINITY,-INFINITY,-INFINITY};
    // Few centers, a remainder, or coordinates whose squares could overflow or that are
    // not numbers: the reference (its all-infinite and remainder behaviour included).
    bool usable = ncenter>=64 && capacity*ncenter==pcce_nspin;
    auto widen=[&](const Point& p) {
        const double v[3]={p.x,p.y,p.z};
        for (int a=0; a<3; ++a) {
            if (!(std::fabs(v[a])<=1e100)) usable=false;
            lo[a]=std::min(lo[a],v[a]); hi[a]=std::max(hi[a],v[a]);
        }
    };
    for (int i=0; i<pcce_nspin && usable; ++i) widen(spins[i]);
    for (int c=0; c<ncenter && usable; ++c) widen(centers[c]);
    double ext=0, maxabs=0;
    if (usable) for (int a=0; a<3; ++a) {
        ext=std::max(ext,hi[a]-lo[a]);
        maxabs=std::max(maxabs,std::max(std::fabs(lo[a]),std::fabs(hi[a])));
    }
    const int side=std::max(1,int(std::ceil(std::cbrt(double(ncenter)))));
    const double h=usable ? ext/side : 0;
    // a cell width that is not positive and finite (ext 0, or a subnormal ext/side rounded to 0: P2-01 of review
    // 20261012-233058-452278) has no grid; checked before any conversion to int
    if (!usable || !(h>0 && std::isfinite(h))) {
        assignCentroidsFromPoints(assigned,spins,centers,pcce_nspin,ncenter,w.counts.data(),w.active.data());
        return;
    }
    const double slack=1e-9*(maxabs+ext);
    int g[3];
    for (int a=0; a<3; ++a) g[a]=std::min(side,int((hi[a]-lo[a])/h)+1);
    auto cellIndex=[&](double v,int a) {
        int ix=int((v-lo[a])/h);
        return ix<0 ? 0 : ix>=g[a] ? g[a]-1 : ix;
    };
    const int ncell=g[0]*g[1]*g[2];
    w.cell_start.assign(ncell+1,0); w.cell_open.assign(ncell,0);
    w.cell_items.resize(ncenter); w.center_cell.resize(ncenter); w.center_slot.resize(ncenter);
    w.open.resize(ncenter); w.open_slot.resize(ncenter);
    int* counts=w.counts.data();
    for (int c=0; c<ncenter; ++c) {
        const int cell=(cellIndex(centers[c].z,2)*g[1]+cellIndex(centers[c].y,1))*g[0]+cellIndex(centers[c].x,0);
        w.center_cell[c]=cell; ++w.cell_start[cell+1];
        counts[c]=0; w.open[c]=c; w.open_slot[c]=c;
    }
    for (int cell=0; cell<ncell; ++cell) w.cell_start[cell+1]+=w.cell_start[cell];
    for (int c=0; c<ncenter; ++c) {
        const int cell=w.center_cell[c], slot=w.cell_start[cell]+w.cell_open[cell]++;
        w.cell_items[slot]=c; w.center_slot[c]=slot;
    }
    int nopen=ncenter;
    for (int i=0; i<pcce_nspin; ++i) {
        const double v[3]={spins[i].x,spins[i].y,spins[i].z};
        const int s[3]={cellIndex(v[0],0),cellIndex(v[1],1),cellIndex(v[2],2)};
        int best=-1;
        double minimum=INFINITY;
        auto consider=[&](int c) {
            const double dx=spins[i].x-centers[c].x;
            const double dy=spins[i].y-centers[c].y;
            const double dz=spins[i].z-centers[c].z;
            const double squared=dx*dx+dy*dy+dz*dz;
            if (squared<minimum || (squared==minimum && c<best)) { minimum=squared; best=c; }
        };
        long visited=0;
        for (int r=0;; ++r) {
            const int z0=std::max(0,s[2]-r), z1=std::min(g[2]-1,s[2]+r);
            const int y0=std::max(0,s[1]-r), y1=std::min(g[1]-1,s[1]+r);
            const int x0=std::max(0,s[0]-r), x1=std::min(g[0]-1,s[0]+r);
            for (int z=z0; z<=z1; ++z) for (int y=y0; y<=y1; ++y) {
                const bool face=std::abs(z-s[2])==r || std::abs(y-s[1])==r;
                const int step=face ? 1 : 2*r;
                for (int x=face ? x0 : s[0]-r; x<=x1; x+=step) {
                    if (x<x0) continue;
                    const int cell=(z*g[1]+y)*g[0]+x;
                    ++visited;
                    for (int k=w.cell_start[cell], e=k+w.cell_open[cell]; k<e; ++k) consider(w.cell_items[k]);
                }
            }
            // Every unsearched cell lies beyond one face of the searched block.
            double bound=INFINITY;
            for (int a=0; a<3; ++a) {
                if (s[a]-r>0)      bound=std::min(bound,v[a]-(lo[a]+(s[a]-r)*h));
                if (s[a]+r<g[a]-1) bound=std::min(bound,lo[a]+(s[a]+r+1)*h-v[a]);
            }
            if (bound==INFINITY) break;
            bound-=slack;
            if (best>=0 && bound>0 && bound*bound*(1-1e-9)-1e-290>minimum) break;
            if (visited>nopen) {
                for (int k=0; k<nopen; ++k) consider(w.open[k]);
                break;
            }
        }
        assigned[i]=best;
        if (++counts[best]==capacity) {
            const int cell=w.center_cell[best], last=w.cell_start[cell]+(--w.cell_open[cell]);
            const int slot=w.center_slot[best], moved=w.cell_items[last];
            w.cell_items[slot]=moved; w.center_slot[moved]=slot;
            w.cell_items[last]=best; w.center_slot[best]=last;
            const int j=w.open_slot[best], tail=w.open[--nopen];
            w.open[j]=tail; w.open_slot[tail]=j;
        }
    }
}

bool updateCentroids(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) {

    bool converged = true;
    int *counts = (int *)malloc(ncenter * sizeof(int));
    Point *sums = (Point *)malloc(ncenter * sizeof(Point));
    
    // Initialize //
    for (int i = 0; i < ncenter; i++) {
        sums[i].x = sums[i].y = sums[i].z = 0;
        counts[i] = 0;
    }

    for (int i = 0; i < pcce_nspin; i++) {
        int cid = assignments[i];
        sums[cid].x += spins[i].x;
        sums[cid].y += spins[i].y;
        sums[cid].z += spins[i].z;
        counts[cid]++;
    }

    double TOLERANCE = 0.0001;
    for (int i = 0; i < ncenter; i++) {
        if (counts[i] > 0) {
            double newX = sums[i].x / counts[i];
            double newY = sums[i].y / counts[i];
            double newZ = sums[i].z / counts[i];
            // Check if centroid moved significantly
            if (fabs(newX - centers[i].x) > TOLERANCE || fabs(newY - centers[i].y) > TOLERANCE || fabs(newZ - centers[i].z) > TOLERANCE){
                converged = false;
            }
            centers[i].x = newX;
            centers[i].y = newY;
            centers[i].z = newZ;
        }
    }
    free(counts);
    free(sums);
    return converged;
}

double calculateInertia(Point *spins, Point *centers, int *assignments, int pcce_nspin, int ncenter) {
    double inertia = 0.0;
    for (int i = 0; i < pcce_nspin; i++) {
        int cid = assignments[i];
        double dx = spins[i].x - centers[cid].x;
        double dy = spins[i].y - centers[cid].y;
        double dz = spins[i].z - centers[cid].z;
        inertia += dx * dx + dy * dy + dz * dz;
    }
    return inertia;
}

double silhouetteCoefficient(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) {
    // A single bundle has no neighboring bundle; singleton bundles contribute 0.
    if (pcce_nspin < 2 || ncenter < 2) return 0.0;
    std::vector<int> counts(ncenter, 0);
    for (int i=0; i<pcce_nspin; ++i) ++counts[assignments[i]];
    std::vector<double> sums(ncenter);
    double total = 0.0;
    for (int i=0; i<pcce_nspin; ++i){
        const int own = assignments[i];
        if (counts[own] <= 1) continue;
        std::fill(sums.begin(), sums.end(), 0.0);
        // Accumulate all bundle distances in one pass, preserving each bundle's
        // original spin order. This removes the O(N*N*K) nested bundle scan.
        for (int j=0; j<pcce_nspin; ++j)
            sums[assignments[j]] += distance(spins[i], spins[j]);
        const double a = sums[own] / (counts[own]-1);
        double b = INFINITY;
        for (int c=0; c<ncenter; ++c)
            if (c != own && counts[c] > 0)
                b = std::min(b, sums[c]/counts[c]);
        const double scale = std::max(a,b);
        if (std::isfinite(b) && scale > 0) total += (b-a)/scale;
    }
    return total/pcce_nspin;
}

double distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/////////////////////////////////////////////////////////////////////////////
// Function: Copy now array --> best array //
void copy_best_centers(Point* centers, Point* best_centers, int ncenter){
    for (int i=0; i<ncenter; i++){
        best_centers[i] = centers[i];
    }
}

void copy_best_assigned(int* assigned, int* best_assigned, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_assigned[i] = assigned[i];
    }
}

void copy_best_idxList(int* idxList, int* best_idxList, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_idxList[i] = idxList[i];
    }
}

void copy_best_spins(Point* spins, Point* best_spins, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_spins[i] = spins[i];
    }
}

void copy_bath_to_point_type(BathSpin** bath, Point* spins, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        spins[i].x = bath[i]->xyz[0];
        spins[i].y = bath[i]->xyz[1];
        spins[i].z = bath[i]->xyz[2];
    }
}


///////////////////////////////////////////////
//            Print & Write Function         // 
///////////////////////////////////////////////
void print_BD(char* str, int num){
    printf("\n");
    for (int b=0; b<num; b++){
        fputs(str, stdout);
    }
    printf("\n");
}

///////////////////////////////////
void print_centers(int ncenter, Point* centers, int rank){
    print_BD("-", 55);
    printf("  - Rank: %d\n\n", rank);
    for (int i=0; i < ncenter; i++){
        printf(" Point %5d: %10.3lf %10.3lf %10.3lf\n", i, centers[i].x, centers[i].y, centers[i].z);
    }
}

///////////////////////////////////
void print_kmeans_converge_info(int trial, int iteration, int ncenter, Point** centers){
    print_BD("-", 55);
    printf(" || Iter. %d - (%d) ||\n", trial, iteration);
    for (int i=0; i < ncenter; i++){
        printf(" Centroid %d: %10.5lf %10.5lf %10.5lf\n", i, (*centers)[i].x, (*centers)[i].y, (*centers)[i].z);
    }
}

///////////////////////////////////
void print_optimize_kmeans_info(int trial, int iteration, double inertia, double best_inertia, double best_sil, int rank){

    //print_BD("=", 55); 
    printf("\n||------------- Optimization Infomation -------------||\n\n");
    //printf(" ||  Iteration  infomation  ||\n\n");
    printf("    - Rank: %d \n\n", rank);
    printf("    - %d trials & %d iterations.\n\n", trial, iteration);
    printf("    - Iner. & Best Iner. : %.2lf & %.2lf\n\n", inertia, best_inertia);
    printf("    - Best Sil.  : %.2lf\n\n", best_sil);
    //print_BD("=", 55); 
}

///////////////////////////////////
void print_assign_idx(int *assign_idx, int pcce_nspin, int rank) {
    printf("\n || Point assignments to centers ||\n\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i = 0; i < pcce_nspin ; i++) {
        printf(" Point %5d -> Center %4d\n", i, assign_idx[i]);
    }
    printf("\n\n");
}

void print_spins(Point *spins, int pcce_nspin, int rank){
    //printf("shuffle!!\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i=0; i<pcce_nspin; i++){
        printf("Point %5d - %10.5lf %10.5lf %10.5lf\n", i, spins[i].x, spins[i].y, spins[i].z);
    }
}

void print_time(clock_t start, clock_t end){
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    int hours      = (int) cpu_time_used / 3600;
    int minutes    = ((int) cpu_time_used % 3600) / 60;
    double seconds = cpu_time_used - (hours * 3600) - (minutes*60);
    printf("    %d : %d : %.3f \n\n", hours, minutes, seconds);
}

///////////////////////////////////
void write_centers(FILE **c_savefile, Point **centers, int ncenter){

    *c_savefile = fopen("./centers.txt", "w");
    //*c_savefile = fopen("./out/centers.txt", "w");
    if (*c_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*c_savefile, " x   y   z \n");
    for (int p=0; p<ncenter; p++){
        fprintf(*c_savefile, " %10.5lf %10.5lf %10.5lf\n", (*centers)[p].x, (*centers)[p].y, (*centers)[p].z);
    }
    fclose(*c_savefile);
}

///////////////////////////////////
void write_spins(FILE **s_savefile, Point **spins, int** assigned_idx, int pcce_nspin){

    *s_savefile = fopen("./spins.txt", "w");
    if (*s_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*s_savefile, " x   y   z   center Idx\n");
    for (int p=0; p<pcce_nspin; p++){
        fprintf(*s_savefile, " %10.5lf %10.5lf %10.5lf %10d\n", (*spins)[p].x, (*spins)[p].y, (*spins)[p].z, (*assigned_idx)[p]);
    }
    fclose(*s_savefile);
}

///////////////////////////////////
void print_best_spins(Point **spins, int** assigned_idx, int pcce_nspin){

    printf(" ===============\n");
    printf(" ==== Spins ====\n");
    printf(" ===============\n");
    printf(" x   y   z   center Idx\n");
    for (int p=0; p<pcce_nspin; p++){
        printf(" %10.5lf %10.5lf %10.5lf %10d\n", (*spins)[p].x, (*spins)[p].y, (*spins)[p].z, (*assigned_idx)[p]);
    }
    printf(" ===============\n");
} 


void print_best_centers(Point **centers, int ncenter){

    printf(" ===========\n");
    printf(" = Centers =\n");
    printf(" ===========\n");
    printf(" x   y   z \n");
    for (int p=0; p<ncenter; p++){
        printf(" %10.5lf %10.5lf %10.5lf\n", (*centers)[p].x, (*centers)[p].y, (*centers)[p].z);
    }
    printf(" ===========\n");
}
///////////////////////////////////
void BathArray_connectivity_pcce(int*** cmap, float*** stmap, Point* centers, float rdip, float rdipcut, int ncenter){
     
    //Point* centers = ba->centers;
    // cmap  [ncenter][ncenter] = 1 if connected else 0
    // stmap [ncenter][ncenter] = strength of the connection

    // Connectivity Map and Strength Map
    *cmap  =   allocInt2d(ncenter, ncenter);
    *stmap = allocFloat2d(ncenter, ncenter);

    // Find connectivity Map and strength Map
    double r = 0.0;
    float strength = 0.0;
    MatrixXcd tensor;

    for (int sp1=0; sp1<ncenter; sp1++){
        for (int sp2=sp1+1; sp2<ncenter; sp2++){
            
            Point c1 = centers[sp1];  double cxyz1[3] = {c1.x, c1.y, c1.z};
            Point c2 = centers[sp2];  double cxyz2[3] = {c2.x, c2.y, c2.z};
        
            r        = dist(cxyz1, cxyz2);
            tensor   = calPointDipoleTensor(cxyz1, cxyz2, 10000, 10000);
            strength = tensor(2,2).real();

            //printf("r = %10.5lf\n", r);
            if ( r < rdip && r > rdipcut){
                (*cmap)[sp1][sp2] = 1;
                (*cmap)[sp2][sp1] = 1;
            }

            (*stmap)[sp1][sp2] = fabs(strength);
            (*stmap)[sp2][sp1] = fabs(strength);
        }
    }
}

///////////////////////////////////
// pcce_connections = spins: two partitions are connected when any spin of the
// one and any spin of the other are within (rdipcut, rdip), as the bundles of
// the adaptive mode are (bundle_cluster.cpp). The strength is the sum of |T_zz|
// over every cross pair of the two partitions, from the same pair tensor that
// the standard CCE connectivity uses.
void BathArray_connectivity_pcce_spins(int*** cmap, float*** stmap, BathArray* ba, int* assigned, int ncenter, float rdip, float rdipcut){

    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        if (assigned[i] < 0 || assigned[i] >= ncenter){
            fprintf(stderr,"Error: invalid pCCE center assignment (spin %d -> %d of %d)\n", i, assigned[i], ncenter);
            MPI_Abort(MPI_COMM_WORLD,1);
            return;
        }
    }
    *cmap  =   allocInt2d(ncenter, ncenter);
    *stmap = allocFloat2d(ncenter, ncenter);

    for (int i=0; i<nspin; i++){
        for (int j=i+1; j<nspin; j++){
            int a = assigned[i], b = assigned[j];
            if (a == b) continue;
            double r = BathArray_dist_i_j(ba,i,j);
            if ( r < rdip && r > rdipcut){
                (*cmap)[a][b] = 1;
                (*cmap)[b][a] = 1;
            }
            float strength = fabs(BathArray_int_i_j(ba,i,j)(2,2).real());
            (*stmap)[a][b] += strength;
            (*stmap)[b][a] += strength;
        }
    }
}

///////////////////////////////////
int** Cluster_setCenterIdx_spinIdx_2dArr(int ncenter, int sK, int pcce_nspin, int* best_assigned_idx){

    int** cs2dArr = allocInt2d(ncenter, sK);
    int*  assignNumArr = allocInt1d(ncenter);

    for (int i=0; i<pcce_nspin; i++){
        int centerIdx = best_assigned_idx[i];
        if (centerIdx < 0 || centerIdx >= ncenter || assignNumArr[centerIdx] >= sK){
            fprintf(stderr,"Error: invalid pCCE center assignment or bundle capacity\n");
            MPI_Abort(MPI_COMM_WORLD,1);
            return NULL;
        }
        cs2dArr[centerIdx][assignNumArr[centerIdx]++] = i;
    } 
    //for (int i=0; i<ncenter; i++){
    //    for (int j=0; j<sK; j++){
    //        printf("cs2dArr[%d][%d] = %d\n", i, j, cs2dArr[i][j]);
    //    }
    //}
    //printf("! ======================= !\n\n");
    free(assignNumArr);
    return cs2dArr;
}

void print_clusinfo_sorder_i(int*** pcceClusInfo, int sorder, int n_cInfo){
    for (int n_cinfo=0; n_cinfo<n_cInfo; n_cinfo++){
        printf("pcceClusInfo[%d][%d] = %p\n", sorder, n_cinfo, (void*)pcceClusInfo[sorder][n_cinfo]);
    }
    printf("=======\n");
}

void print_clusinfo_sorder_i_j(int*** pcceClusInfo, int sorder, int n_cinfo, int cIdx){
    //printf("=======\n");
    //printf("sorder == %d\n", sorder);
    //printf("sorder / sK = 0\n");
    //printf("=======\n");
    for (int cidx=0; cidx<cIdx; cidx++){
        printf("pcceClusInfo[%d][%d][%d] = %d\n", sorder, n_cinfo, cidx, pcceClusInfo[sorder][n_cinfo][cidx]);
    }
}

int*** convert_centerIdx_to_spinIdx(Cluster* cls, int cOrder, int sK, int** cs2dArr){
    //int*** pcceClusInfo;

    int sOrder              = cOrder * sK; // (sOrder = spin Order)  &  (cOrder = center Order <== cce.input file) //
    int*** pcceClusInfo     = (int***)allocArray1d(sOrder+1, sizeof(int**)); 

    // When cOrder = 0 // 
    pcceClusInfo[0]       = (int**)allocArray1d(1, sizeof(int*));
    pcceClusInfo[0][0]    = (int* )allocArray1d(1, sizeof(int));
    pcceClusInfo[0][0][0] = 1;
    
    // When cOrder > 0 // 
    for (int sorder = 1; sorder <= sOrder; sorder++){

        if (sorder % sK==0){
            int corder                 = int(sorder/sK);
            int n_cInfo                = cls->clusinfo[corder][0][0] - 1;               // n_clnfo: The number of "center clusinfo"'s n-th order cluster //
            pcceClusInfo[sorder]       = (int**)allocArray1d(n_cInfo+1, sizeof(int*));
            //print_clusinfo_sorder_i(pcceClusInfo, sorder, n_cInfo+1);

            if (n_cInfo != 0) {
                pcceClusInfo[sorder][0]    = (int* )allocArray1d(1, sizeof(int));
                pcceClusInfo[sorder][0][0] = (n_cInfo+1);
                //printf("pcceClusInfo[%d][0][0] = %d\n", sorder, pcceClusInfo[sorder][0][0]);

                for (int n_cinfo=1; n_cinfo<(n_cInfo+1); n_cinfo++){
                    //printf("!=========================!\n");
                    //printf("sorder & n_cinfo & n_cInfo: %d & %d & %d \n\n", sorder, n_cinfo, n_cInfo);
                    pcceClusInfo[sorder][n_cinfo] = (int* )allocArray1d((sorder+1), sizeof(int));

                    int bd = 1;
                    pcceClusInfo[sorder][n_cinfo][0] = cls->clusinfo[corder][n_cinfo][0];
                    //printf("cls      -> clustinfo[%d][%d][0] = %d\n"  , corder, n_cinfo, cls      -> clusinfo[corder][n_cinfo][0]);
                    //printf("pcceClusInfo[%d][%d][0] = %d\n\n", sorder, n_cinfo, pcceClusInfo[sorder][n_cinfo][0]);
                    for (int cidx=1; cidx<corder+1; cidx++){
                        int centerIdx = cls->clusinfo[corder][n_cinfo][cidx];
                        //printf("corder, n_cinfo, cidx: %d, %d, %d\n", corder, n_cinfo, cidx);
                        for (int sk=0; sk<sK; sk++){
                            pcceClusInfo[sorder][n_cinfo][bd++] = cs2dArr[centerIdx][sk];
                        }
                    }
                    //print_clusinfo_sorder_i_j(pcceClusInfo, sorder, n_cinfo, sOrder+1);
                    //printf("\n!=========================!\n");
                    //printf("AFTER SETTING nsij!!\n");
                    //printf("---------------------\n");
                }
            }

            else {
                pcceClusInfo[sorder][0]      = (int* )allocArray1d(1, sizeof(int));
                pcceClusInfo[sorder][0][0]   = (n_cInfo+1);
            }
        }

        else {
            pcceClusInfo[sorder]       = (int**)allocArray1d(1,sizeof(int*));
            pcceClusInfo[sorder][0]    = (int*)allocArray1d(1,sizeof(int));
            pcceClusInfo[sorder][0][0] = 1;

            //printf("=======\n");
            //printf("sorder == %d\n", sorder);
            //printf("sorder !~ sK != 0\n");
            //printf("pcceClusInfo[sorder][0][0] = %d\n", pcceClusInfo[sorder][0][0]);
            //printf("=======\n");
        }
    }
    //printf("!! FINSIH: Function \"Convert center indices => spin indices. \" !!\n");
    return pcceClusInfo;
}
