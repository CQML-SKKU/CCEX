// rdip ladder (include/auto_convergence/rdip_ladder.h; (internal notes)).
#include "../../include/auto_convergence/rdip_ladder.h"
#include "../../include/clustering/pcce_partition_cache.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <strings.h>

namespace rdipladder {


// Every rank reaches the same refusal (the checks read the input and the clusters every rank builds alike), so the ranks
// meet and stop with exit 1 instead of MPI_Abort: an abort let mpirun drop the message now and then (check_rdip_ladder
// bad_saveall failed twice on 2026-10-08 with no message in the log).
[[noreturn]] void refuse(const char* why){
    ccexRefusalPrint(stderr,"Error: rdip_ladder %s\n",why);
    ccexRefuse();   // the same stop as every input refusal (A55, include/common/utilities.h)
}

// pCCE: the partition this run's clusterize made or read (pcce_partition_cache.cpp), after the remainder spins were cut.
Units::Units(Config* c,Cluster* cls,BathArray* bath):cut(Config_getRdipcut(c)),ba(bath){
    if(c->nrdip_ladder==0 || strcasecmp(c->method,"pcce"))return;
    const PccePartition& p=pcceLastPartition();
    const int nspin=BathArray_getNspin(ba);
    if(int(p.assigned.size())!=nspin || p.centers.empty())refuse("found no pCCE partition of this bath");
    size=cls->sK;centres=strcasecmp(c->pcce_connections,"spins")!=0;
    of=p.assigned;centre=p.centers;
    members.assign(centre.size(),{});
    for(int i=0;i<nspin;++i){
        if(of[i]<0 || of[i]>=int(members.size()))refuse("found a spin outside the pCCE partitions");
        members[of[i]].push_back(i);
    }
    for(auto& m:members)if(int(m.size())!=size)refuse("found a pCCE partition whose size is not sK");
}

// The run's own connectivity (BathArray_connectivity, BathArray_connectivity_pcce, BathArray_connectivity_pcce_spins):
// connected when rdipcut < d < rdip, d the spin distance, the centre distance, or any cross spin distance.
double Units::distance(int a,int b) const{
    constexpr double never=std::numeric_limits<double>::infinity();
    if(of.empty()){const double d=BathArray_dist_i_j(ba,a,b);return d>cut ? d : never;}
    if(centres){
        double ca[3]={centre[a].x,centre[a].y,centre[a].z},cb[3]={centre[b].x,centre[b].y,centre[b].z};
        const double d=dist(ca,cb);return d>cut ? d : never;
    }
    double m=never;
    for(int i:members[a])for(int j:members[b]){
        const double d=i<j ? BathArray_dist_i_j(ba,i,j) : BathArray_dist_i_j(ba,j,i);
        if(d>cut && d<m)m=d;
    }
    return m;
}

// After clusterize, on every rank: refuse what the ladder does not cover and re-add the dropped single units at power 0.
void prepare(Config* c,Cluster* cls,BathArray* ba,Output* op){
    if(c->nrdip_ladder==0)return;
    const bool pcce=!strcasecmp(c->method,"pcce");
    if((strcasecmp(c->method,"cce") && strcasecmp(c->method,"gcce") && !pcce) || strcasecmp(c->quantity,"coherence"))
        refuse("needs method CCE, gCCE or pCCE and quantity coherence");
    if(Cluster_getOrder(cls)!=(pcce ? 2*cls->sK : 2))refuse("needs order 2");
    if(strcasecmp(c->cluster_mode,"standard"))refuse("needs cluster_mode standard");
    if(c->cluster_contributions)refuse("cannot be combined with cluster_contributions");
    if(strcasecmp(Output_getSavemode(op),"coherence"))refuse("needs savemode coherence");
    const float rdip=Config_getRdip(c),cut=Config_getRdipcut(c);
    for(int k=0;k<c->nrdip_ladder;++k)
        if(!(c->rdip_ladder[k]>cut && c->rdip_ladder[k]<rdip))refuse("rungs must lie between rdipcut and rdip");
    const Units units(c,cls,ba);
    const int n=units.count(),size=units.size;
    long pairs=0;                                   // the run's pairs must be every pair within (rdipcut, rdip): no Nk cut
    for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)if(units.distance(a,b)<rdip)++pairs;
    if(pairs!=long(cls->clusinfo[2*size][0][0])-1)refuse("needs every pair within rdip as a cluster (no Nk cut)");
    int**& rows=cls->clusinfo[size];
    const int old=rows[0][0];
    std::vector<char> present(n,0);
    for(int i=1;i<old;++i)present[units.unit(rows[i][1])]=1;
    std::vector<int> missing;
    for(int u=0;u<n;++u)if(!present[u])missing.push_back(u);
    if(missing.empty())return;
    int** grown=allocInt2d(old+missing.size(),size+1);
    for(int i=1;i<old;++i)for(int s=0;s<=size;++s)grown[i][s]=rows[i][s];
    grown[0][0]=old+int(missing.size());
    for(size_t m=0;m<missing.size();++m){                    // spins ascending, as the run's own single units
        grown[old+m][0]=0;
        for(int s=0;s<size;++s)grown[old+m][1+s]=units.of.empty() ? missing[m] : units.members[missing[m]][s];
    }
    freeInt2d(&rows,old);rows=grown;
    if(rank==0)printf("[rdip-ladder] %d rungs below rdip %g; %zu single %s re-added at power 0\n",c->nrdip_ladder,double(rdip),
                      missing.size(),pcce ? "partitions" : "spins");
}



// Rank 0: write rung k as <outfile>_rdip<r_k>_state<s>_{wiDiv,noDiv}.
void save(Output* op,float rung,MatrixXcd* wD,MatrixXcd* nD,int nstep,double deltat,int istate){
    const std::string original=Output_getOutfile(op);
    char suffix[64];snprintf(suffix,sizeof suffix,"_rdip%g",double(rung));
    strcat(Output_getOutfile(op),suffix);
    Output_save(op,wD,nD,nstep,deltat,istate);
    strcpy(Output_getOutfile(op),original.c_str());
}

}  // namespace rdipladder
