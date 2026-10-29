/**
 * @file cluster_pairs.cpp
 * @brief The order-2 cluster list of CCE/gCCE made without the hash, shared out between the ranks.
 */
#include "../../include/clustering/cluster_pairs.h"
#include "../../include/common/memory.h"

#include <math.h>
#include <mpi.h>
#include <vector>

// Why this exists (improvements.md A21): at order 2 every rank built the same list through
// BathArray_connectivity (two dense nspin x nspin maps), makeSparsemap and clusterizeHash
// (one hash entry, id string and pair per cluster, then a linked-list merge sort of all pairs).
// For the production P1 bath (6,035 spins, 4,503,988 pairs) that was 24.9 s of a 131 s run on
// 64 ranks, and 64 copies of the same random-access work also slowed each other down.
// The list made here is the same, row for row: the same pairs, the same order, the same counts.

namespace {

// by_strength (cluster_hash.cpp) on the two strengths: descending, 0 on a tie.
int byStrength(float a, float b){
    if (a == b) return 0;
    return a < b ? 1 : -1;
}

// HASH_SORT (uthash's HASH_SRT) as an array sort: the same bottom-up merge - runs of 1, 2, 4, ...
// from the front, the left run taken while the comparison is <= 0, stop after the pass with one
// merge. Doing the same merges with the same comparisons gives the same order whatever the
// strengths are (std::stable_sort would agree only for a consistent ordering, i.e. without NaN).
void sortLikeHash(std::vector<int>& order, const std::vector<float>& strength){
    const size_t n = order.size();
    if (n == 0) return;
    std::vector<int> merged(n);
    for (size_t width = 1;; width *= 2){
        size_t nmerges = 0;
        for (size_t lo = 0; lo < n; lo += 2*width){
            nmerges++;
            size_t p = lo, pend = (lo+width < n) ? lo+width : n;
            size_t q = pend, qend = (pend+width < n) ? pend+width : n;
            size_t out = lo;
            while (p < pend || q < qend){
                if (p == pend)                                              merged[out++] = order[q++];
                else if (q == qend)                                         merged[out++] = order[p++];
                else if (byStrength(strength[order[p]],strength[order[q]]) <= 0) merged[out++] = order[p++];
                else                                                        merged[out++] = order[q++];
            }
        }
        order.swap(merged);
        if (nmerges <= 1) break;
    }
}

// updateNk's report (cluster_hash.cpp), printed at the same three places.
void reportNk(const int* nk, int order){
    if (rank==0 && verbosity){
        printf("        update \"Nk\"\n        ");
        for (int i=0; i<=order; i++){
            printf("%s%d:%d",i ? "  " : "",i,nk[i]);
        }
        printf("\n\n");
    }
}

} // namespace

void clusterizePairs(Cluster* cls, BathArray* ba, float rdip, float rdipcut,
                     const std::function<double(int,const int*)>& cost){

    const int order = Cluster_getOrder(cls);
    const int nspin = BathArray_getNspin(ba);
    int* nks = Cluster_getNk(cls);

    // makeHashClusterO1 + updateNk: one single per spin.
    if (rank==0){printf("\n      CCE1 ...\n");}
    if (nspin > 0) nks[1] = nspin;
    reportNk(nks,order);

    // The pairs in the order makeHashClusterO2 inserted them: rows ascending, then columns above
    // the row ascending. Each rank searches a block of rows (about the same number of candidate
    // pairs per rank) and the blocks are put back together in rank order, so the result is the
    // whole list in that same order on every rank. Connection and strength are the expressions of
    // BathArray_connectivity; the tensor is now made only for connected pairs, the only ones whose
    // strength is read at order 2.
    if (rank==0){printf("\n      CCE2 ...\n");}
    const long long candidates = (long long)nspin*(nspin-1)/2;
    int first = nspin, last = nspin;
    {
        long long seen = 0;
        const long long from = candidates*rank/nprocess, to = candidates*(rank+1)/nprocess;
        for (int row = 0; row < nspin; row++){
            if (first == nspin && seen >= from) first = row;
            if (seen >= to){ last = row; break; }
            seen += nspin-1-row;
        }
        if (first > last) first = last;
    }
    std::vector<int> localRow, localCol;
    std::vector<float> localStrength;
    for (int sp1 = first; sp1 < last; sp1++){
        for (int sp2 = sp1+1; sp2 < nspin; sp2++){
            double r = BathArray_dist_i_j(ba,sp1,sp2);
            if ( r < rdip && r > rdipcut){
                float strength = BathArray_int_i_j(ba,sp1,sp2)(2,2).real();
                localRow.push_back(sp1);
                localCol.push_back(sp2);
                localStrength.push_back(fabs(strength));
            }
        }
    }

    int localCount = (int)localRow.size();
    std::vector<int> counts(nprocess), displs(nprocess);
    MPI_Allgather(&localCount,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD);
    long long total = 0;
    for (int i = 0; i < nprocess; i++){
        displs[i] = (int)total;
        total += counts[i];
    }
    // The cluster list counts its rows in int; so did the hash path.
    if (total >= 2147483647LL){
        if (rank==0) fprintf(stderr,"Error: clusterize: %lld pairs within rdip do not fit the cluster list\n",total);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    const int npair = (int)total;
    std::vector<int> pairRow(npair), pairCol(npair);
    std::vector<float> strength(npair);
    MPI_Allgatherv(localRow.data(),localCount,MPI_INT,pairRow.data(),counts.data(),displs.data(),MPI_INT,MPI_COMM_WORLD);
    MPI_Allgatherv(localCol.data(),localCount,MPI_INT,pairCol.data(),counts.data(),displs.data(),MPI_INT,MPI_COMM_WORLD);
    MPI_Allgatherv(localStrength.data(),localCount,MPI_FLOAT,strength.data(),counts.data(),displs.data(),MPI_FLOAT,MPI_COMM_WORLD);
    // not read again; their storage freed now, not at the end of the function (peak memory; P2-02 of review 20261012-233058-452278)
    std::vector<int>().swap(localRow); std::vector<int>().swap(localCol); std::vector<float>().swap(localStrength);

    std::vector<int> sorted(npair);
    for (int i = 0; i < npair; i++) sorted[i] = i;
    sortLikeHash(sorted,strength);
    if (npair > 0) nks[2] = npair;   // updateNk sets it only when there is a pair
    reportNk(nks,order);
    reportNk(nks,order);             // after makeHashClusterOn, which has nothing to do at order 2

    // A38: with a cost, this rank allocates only the rows it will calculate - the rows of
    // MPI_getLocalClusters, decided here from every pair with the same cost function, in the same
    // order (rankClusterRows), so the split is the one it would make from the full list - and the first
    // and last four, which rank 0 prints (Cluster_reportClusinfo). The other slots stay NULL.
    // A45: the rows are no longer one range (largest cost first); a cost that differs keeps 12 B a pair while it decides.
    std::vector<float>().swap(strength);   // not read again; freed before the rows (peak memory)
    std::vector<int> held;
    if (cost){
        int row[3] = {1,0,0};
        std::vector<double> loads(nprocess, 0.0);   // order 2 is the highest: MPI_getLocalClusters splits it first
        held = rankClusterRows(npair,[&](int i){
            row[1] = pairRow[sorted[i]];
            row[2] = pairCol[sorted[i]];
            return cost(2,row);
        },loads,2,[&](int i){
            row[1] = pairRow[sorted[i]];
            row[2] = pairCol[sorted[i]];
            return (const int*)row;
        });
        free(cls->pairRows);
        cls->npairRows = (int)held.size();
        cls->pairRows = held.empty() ? NULL : (int*)malloc(held.size()*sizeof(int));
        for (size_t r = 0; r < held.size(); r++) cls->pairRows[r] = held[r];
        free(cls->pairLoads);
        cls->pairLoads = (double*)malloc(nprocess*sizeof(double));
        for (int r = 0; r < nprocess; r++) cls->pairLoads[r] = loads[r];
        cls->pairsLocal = true;
    }
    size_t nextHeld = 0;                  // held is ascending, the rows below are walked in order

    // solveTilde at order 2: every pair is a row with count 1, in the sorted order; a single's
    // count is 1 minus the pairs it is in, and the singles with count 0 are left out, the rest
    // kept in spin order.
    cls->clusinfo[2] = (int**)realloc(cls->clusinfo[2],(npair+1)*sizeof(int*));
    std::vector<int> single(nspin,1);
    for (int i = 0; i < npair; i++){
        const int p = sorted[i];
        const int index = i+1;
        int* row = NULL;
        const bool mine = nextHeld < held.size() && held[nextHeld] == index;
        if (mine) nextHeld++;
        if (!cost || mine || index < 5 || index > npair-4){
            row = (int*)calloc(3,sizeof(int));
            row[0] = 1;
            row[1] = pairRow[p];
            row[2] = pairCol[p];
        }
        cls->clusinfo[2][index] = row;
        single[pairRow[p]]--;
        single[pairCol[p]]--;
    }
    cls->clusinfo[2][0][0] = npair+1;

    int nsingle = 0;
    for (int i = 0; i < nspin; i++) if (single[i] != 0) nsingle++;
    cls->clusinfo[1] = (int**)realloc(cls->clusinfo[1],(nsingle+1)*sizeof(int*));
    int k = 1;
    for (int i = 0; i < nspin; i++){
        if (single[i] == 0) continue;
        int* row = (int*)calloc(2,sizeof(int));
        row[0] = single[i];
        row[1] = i;
        cls->clusinfo[1][k++] = row;
    }
    cls->clusinfo[1][0][0] = nsingle+1;
}
