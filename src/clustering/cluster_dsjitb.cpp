/**
 * @file cluster_dsjitb.cpp
 * @brief The disjoint (dsj) and inter-bath-cluster (itb) clustering variants.
 */
#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "../../include/clustering/cluster.h"
#include "../../include/common/memory.h"
#include "../../include/clustering/cluster_dsjitb.h"

extern int rank;

void clusterizeDsjItb(Cluster* cls, int*** cmap, float*** stmap, BathArray* ba){

    MPI_Barrier(MPI_COMM_WORLD);

    int order = Cluster_getOrder(cls);
    int nspin = BathArray_getNspin(ba);
    char* method = Cluster_getMethod(cls);

    if (rank == 0){
        printSubTitle("Make a disjoint clusters");
    }


    //make 2d Duo13C with coupling strength
    // NOTE! in this case, the index is larger than "1" actual index in BathArray->bath
    //Duo13CwithCS[0] = [ nspin + 1 , 0.0, 0.0 ]
    //Duo13CwithCS[1] = [ PairIdx[0] + 1 , PairIdx[1] + 1, CouplingStr ]
    double **Duo13CwithCS = NULL;
    makeDuo13CwithCS(&Duo13CwithCS, nspin, cmap, stmap); //

    //printf("\n ! Duo13CwithCS \n"); 
    //for(int j=0;j<(int)Duo13CwithCS[0][0];j++){
    //    printf("[%d] :",j); 
    //    printf("%lf %lf %lf",Duo13CwithCS[j][0],Duo13CwithCS[j][1],Duo13CwithCS[j][2] );
    //    printf("\n");
    //}
    //printf("\n");

    //sort Duo13CwithCS based on 2th column (last column)
    QuickSort_Double_2d(&(Duo13CwithCS),3,1,int((Duo13CwithCS[0][0])-1));

    //printf("\n ! Duo13CwithCS \n"); 
    //for(int j=0;j<(int)Duo13CwithCS[0][0];j++){
    //    printf("[%3d] :",j); 
    //    printf("%3d - 1    %3d - 1  %10.4lf",int(Duo13CwithCS[j][0]),int(Duo13CwithCS[j][1]),Duo13CwithCS[j][2] );
    //    printf("\n");
    //}
    //printf("\n");

    // interbath pairs
    int** InterBathPairs = allocInt2d(1,2);
    InterBathPairs[0][0] = 1; 

    // find disjoint Groups
    // Groups (group# + 1 , order + 1): 
    // [0][0] : the number of groups + 1
    // [i][0] : the number of spins in the n-th group + 1 = N[G_n] + 1
    // (i!=0)
    int line_groups = 1;
    int** Groups = allocInt2d(line_groups,order+1);
    Groups[line_groups-1][0] = line_groups;

    // Check if the spin is in the group
    // spinGroupChecker : 
    //  [0] : the number of nuclei in the n-th group = N[G_n]
    //  [1] : group number n 
    int **GroupChecker = allocInt2d((int)(nspin + 1),2);
    GroupChecker[0][0] = (int)(nspin + 1);
    
    for (int i = 1; i < (int)Duo13CwithCS[0][0]; i++){

        int spin1 = (int)Duo13CwithCS[i][0];
        int spin2 = (int)Duo13CwithCS[i][1];

        int GroupChecker_nSpin_spin1 = GroupChecker[spin1][0];
        int GroupChecker_Group_spin1 = GroupChecker[spin1][1];

        int GroupChecker_nSpin_spin2 = GroupChecker[spin2][0];
        int GroupChecker_Group_spin2 = GroupChecker[spin2][1];


		//printf("Spin : %d, %d\n",spin1,spin2);
		//printf("GroupChecker_Group : %d, %d\n",GroupChecker_Group_spin1,GroupChecker_Group_spin2);
		//printf("GroupChecker_nSpin : %d, %d\n",GroupChecker_nSpin_spin1,GroupChecker_nSpin_spin2);

        ////////////////////////////////////////////////////////////
        //printf("Groups : \n");
        //for(int igroup=0; igroup<Groups[0][0]; igroup++){
        //    printf("[%4d]",igroup);
        //    for (int jgroup=0; jgroup<order+1;jgroup++){
        //        printf("%4d ",Groups[igroup][jgroup]);
        //    }
        //    printf("\n");
        //}                
        //////////////////////////////////////////////////////////

        // Case1 - i, j not in group 
        if (GroupChecker_Group_spin1 == 0 && GroupChecker_Group_spin2 == 0){
            // make new group
            int* group;
            jointSpinSpin(&(group),order+1,spin1,spin2);
			//printf("make Group : %d %d %d\n",group[0],group[1],group[2]);
            //add group to Groups
            int newGroupNumber;
            newGroupNumber = addNewGroup(&(Groups),group, order);

            // update GroupChecker
            updateGroupChecker(&(GroupChecker),group,newGroupNumber);

			//printf("- %d, %d (case1)\n",newGroupNumber,newGroupNumber);
        }
        // Case2-1 - i in group , j not in group  and Ì†µÌ≤©[Ì†µÌ¥æ(i)] < order 
        else if ( GroupChecker_Group_spin1 != 0 \
               && GroupChecker_Group_spin2 == 0 \
               && GroupChecker_nSpin_spin1 < order){

            // add spin2 to group
            jointGroupSpin(&(Groups[GroupChecker_Group_spin1]),spin2);

            // update GroupChecker
            updateGroupChecker(&(GroupChecker),Groups[GroupChecker_Group_spin1],GroupChecker_Group_spin1);
               
			//printf("- %d, %d (case2-1)\n",GroupChecker_Group_spin1,GroupChecker_Group_spin1);
        }
        // Case2-2 - i in group , j not in group and Ì†µÌ≤©[Ì†µÌ¥æ(i)] >= order 
        else if ( GroupChecker_Group_spin1 != 0 \
               && GroupChecker_Group_spin2 == 0 \
               && GroupChecker_nSpin_spin1 >= order){

            int* group = (int*)calloc(order+1,sizeof(int)); 
			group[0] = 2;
			group[1] = spin2;

            int newGroupNumber;
			newGroupNumber = addNewGroup(&Groups,group,order);
            updateGroupChecker(&(GroupChecker),group,newGroupNumber);
       
            addInterBathPair(&(InterBathPairs),spin1,spin2);
			//printf("- %d, %d (case2-2)\n",GroupChecker_Group_spin1,newGroupNumber);
        }
        // Case3-1 - i not in group , j in group  and Ì†µÌ≤©[Ì†µÌ¥æ(j)] < order 
        else if ( GroupChecker_Group_spin1 == 0 \
               && GroupChecker_Group_spin2 != 0 \
               && GroupChecker_nSpin_spin2 < order){

            // add spin1 to group
            jointGroupSpin(&(Groups[GroupChecker_Group_spin2]),spin1);

            // update GroupChecker
            updateGroupChecker(&(GroupChecker),Groups[GroupChecker_Group_spin2],GroupChecker_Group_spin2);
			//printf("- %d, %d (case3-1)\n",GroupChecker_Group_spin2,GroupChecker_Group_spin2);
               
        }
        // Case3-2 - i not in group , j in group  and Ì†µÌ≤©[Ì†µÌ¥æ(j)]  >= order 
        else if ( GroupChecker_Group_spin1 == 0 \
               && GroupChecker_Group_spin2 != 0 \
               && GroupChecker_nSpin_spin2 >= order){

            int* group = (int*)calloc(order+1,sizeof(int)); 
			group[0] = 2;
			group[1] = spin1;
			int newGroupNumber = 0;
			newGroupNumber = addNewGroup(&Groups,group,order);
            updateGroupChecker(&(GroupChecker),group,newGroupNumber);
            addInterBathPair(&(InterBathPairs),spin1,spin2);
			//printf("- %d, %d (case3-2)\n",newGroupNumber,GroupChecker_Group_spin2);
        }
        // Case4-1 - I,j in different group and Ì†µÌ≤©[Ì†µÌ¥æ(i)] + Ì†µÌ≤©[Ì†µÌ¥æ(j)] <= order 
        else if ( GroupChecker_Group_spin1 != 0 \
               && GroupChecker_Group_spin2 != 0 \
               && GroupChecker_Group_spin1 != GroupChecker_Group_spin2 \
               && GroupChecker_nSpin_spin1 + GroupChecker_nSpin_spin2 <= order){

            // we'll add the spins into the group that the other spin include
            // the target group that add other spins would have smaller group number that other spins' group
            // e.i. updating group should be smaller than removing group
            int updatingGroup = 0;
            int removingGroup = 0;
            if (GroupChecker_Group_spin1 < GroupChecker_Group_spin2){
                updatingGroup = GroupChecker_Group_spin1;
                removingGroup = GroupChecker_Group_spin2;
            }
            else{
                updatingGroup = GroupChecker_Group_spin2;
                removingGroup = GroupChecker_Group_spin1;
            }

            // joint group
            jointGroupGroup(&(Groups[updatingGroup]),&(Groups[removingGroup]));
            // remove group
            int newrow = 0;
            newrow = Delete_Int_2dArr(&(Groups),Groups[0][0],order+1,removingGroup);
            Groups[0][0] = newrow;
            
            // update GroupChecker for newly jointed group
            updateGroupChecker(&(GroupChecker),Groups[updatingGroup],updatingGroup);

            // update GroupChecker from deleted group to last group
            for (int j = removingGroup; j < Groups[0][0]; j++){
                updateGroupChecker(&(GroupChecker),Groups[j],j);
            }
			//printf("- Updated group : %d, Deleted group %d (case4-1)\n",updatingGroup,removingGroup);
        }
        // Case4-2 - I,j in different group and Ì†µÌ≤©[Ì†µÌ¥æ(i)] + Ì†µÌ≤©[Ì†µÌ¥æ(j)] > order
        else if ( GroupChecker_Group_spin1 != 0 \
               && GroupChecker_Group_spin2 != 0 \
               && GroupChecker_Group_spin1 != GroupChecker_Group_spin2 \
               && GroupChecker_nSpin_spin1 + GroupChecker_nSpin_spin2 > order){
			//printf("- Not made a group (case4-2)\n");
            addInterBathPair(&(InterBathPairs),spin1,spin2);
        }
		// Case 5 - i,j in the same group
        else if ( GroupChecker_Group_spin1 != 0 \
               && GroupChecker_Group_spin2 != 0 \
               && GroupChecker_Group_spin1 == GroupChecker_Group_spin2){
			;
			//printf("- Same group (case5)\n");
        }
        else{
            printf("Error in createDisjointClusters\n");
            exit(1);
        }
    }

	// not paired spin (not in duo.)
	for (int i=1; i<GroupChecker[0][0]; i++){
		if (GroupChecker[i][1] == 0){
            int* group = (int*)calloc(order+1,sizeof(int)); 
			group[0] = 2;
			group[1] = i;
			int newGroupNumber = 0;
			newGroupNumber = addNewGroup(&Groups,group,order);
            updateGroupChecker(&(GroupChecker),group,newGroupNumber);
		}
    }

    // -------// dsj algo done
    /////////////////////////////////////////////////////////////////////////////////
    //
	//printf("----- Groups ----- \n");
    //printf("%d \n",(int)ID->Mono13C[0][0]);
    int* nSpinChk = (int*)calloc((int)(nspin + 1),sizeof(int));
    nSpinChk[0] = int(nspin + 1);

    int nSpinInGroup=0;
    for (int i=1; i<Groups[0][0]; i++){
        for (int j=1; j<Groups[i][0]; j++){
            //printf("%5d ",Groups[i][j]);
            nSpinChk[Groups[i][j]]++;
        }
		//printf("\n");
        nSpinInGroup=nSpinInGroup+Groups[i][0]-1;
    }
	//printf("\nTotal Spin # in Group : %d \n",nSpinInGroup);

	if (nSpinInGroup != nspin){
		printf("Error, spin# in disjoint group(%d) != total spin#(%d)\n",nSpinInGroup,(int)nspin);
		for (int i=1; i<nSpinChk[0]; i++){
			if (nSpinChk[i]!=1){
				printf("Error spin : spin \"%d - 1\" - %d#\n",i,nSpinChk[i]);
			}
		}
		exit(1);
	}
	free(nSpinChk);
	///printf("\n");
    ///
    //printf("nSpinChk : \n");
    //for (int i=0; i<nSpinChk[0]; i++){
    //    printf("[%d] %d\n",i,nSpinChk[i]);

    //}


        
    /////////////////////////////////////////////////////////////////////////////////
    // Add clustered spins in clusinfo
    /////////////////////////////////////////////////////////////////////////////////
    //
	// cce1 -> all [0] = 0
	// because all spins are grouped, so initial value should be 0
    //Cluster_reportClusinfo(cls);
    //
    reallocInt2d(&(cls->clusinfo[1]),(int)cls->clusinfo[1][0][0],nspin+1,2); 
    cls->clusinfo[1][0][0] = nspin + 1; 
	for (int i=1; i<cls->clusinfo[1][0][0]; i++){
		cls->clusinfo[1][i][0]  = 0;	
		cls->clusinfo[1][i][1]  = i-1;	
	}

    //Cluster_reportClusinfo(cls);
    // if itb + dsj or dsj, then do 
    if (strcasecmp(method,"dsjitb")==0
            || strcasecmp(method,"dsj")==0){


        // update cls->clusinfo[order] with Groups
        for (int i=1; i<Groups[0][0]; i++){
            int nSpin = Groups[i][0] - 1;
            int clusIdx = 0;

            if (nSpin !=1){
                // clsinfo alloc.
                reallocInt2d(&(cls->clusinfo[nSpin]),cls->clusinfo[nSpin][0][0],cls->clusinfo[nSpin][0][0]+1,nSpin+1);

                // Give values
                cls->clusinfo[nSpin][0][0] = cls->clusinfo[nSpin][0][0] + 1; //Total #
                clusIdx = cls->clusinfo[nSpin][0][0] - 1; 
                cls->clusinfo[nSpin][clusIdx][0] = 1;

                for (int j=1; j<Groups[i][0]; j++){
                    int spin = Groups[i][j] - 1;   // actual index need "-1"
                	cls->clusinfo[nSpin][clusIdx][j] = spin;
                }
            }else{
                clusIdx = Groups[i][1];
                int spin = Groups[i][1] - 1;   // actual index need "-1"
                cls->clusinfo[nSpin][clusIdx][0] = 1;
                cls->clusinfo[nSpin][clusIdx][1] = spin;
            }

            QuickSort(&(cls->clusinfo[nSpin][clusIdx]),1,nSpin);
        }
        // Write the clusters
	    //for (int i=1; i<order+1; i++){

//      //      sort2DArray(&(CCE->ClusInfo[i]),CCE->ClusInfo[i][0][0],i+1,1);
	    //	char orderStr[5];
	    //	char fname[500] = "Clusters_disjoint_CCE";
	    //	sprintf(orderStr,"%d",i);
	    //	strcat(fname,orderStr);
	    //	//writeClusters(CCE->ClusInfo[i],i,fname);
	    //}
    }

    // if itb + dsj or itb, then
    if (strcasecmp(method,"dsjitb")==0
            || strcasecmp(method,"itb")==0){

        // write InterBathPairs
	    //char orderStr[5];
	    //char fname[500] = "Clusters_InterBathPair";
	    //writeClusters(InterBathPairs,1,fname);

        // update cls->clusinfo[order] with InterBathPairs
        for (int i=1; i<InterBathPairs[0][0]; i++){

            int spin1 = InterBathPairs[i][0] - 1; // actual index need "-1"
            int spin2 = InterBathPairs[i][1] - 1; // actual index need "-1"
	    
	    	// L^~_ij
            reallocInt2d(&(cls->clusinfo[2]),cls->clusinfo[2][0][0],cls->clusinfo[2][0][0]+1,3);
            cls->clusinfo[2][0][0] = cls->clusinfo[2][0][0] + 1;

            int clusIdx = cls->clusinfo[2][0][0] - 1;
            cls->clusinfo[2][clusIdx][0] = 1;
            cls->clusinfo[2][clusIdx][1] = spin1;
            cls->clusinfo[2][clusIdx][2] = spin2;

	    	// L~_ij = Lij / Li, Lj. Row s+1 of clusinfo[1] is spin s (row 0 holds the count); the rows were indexed by
	    	// s until 2026-11-05, which took the power from the spin before and, for spin 0, shortened the list
            cls->clusinfo[1][spin1+1][0]--;
            cls->clusinfo[1][spin2+1][0]--;
        }
    }
	
    shrinkClusterO1(cls);    


	//write clusters after adding interBath clusters
	//for (int i=1; i<order+1; i++){
//  //      sort2DArray(&(CCE->ClusInfo[i]),CCE->ClusInfo[i][0][0],i+1,1);
	//	char orderStr[5];
	//	char fname[500] = "Clusters_final_CCE";
	//	sprintf(orderStr,"%d",i);
	//	strcat(fname,orderStr);
	//	writeClusters(CCE->ClusInfo[i],i,fname);
	//}

    //printf("\n<CCE->ClusInfo> \n"); 
    //for(int i=0;i<order+1;i++){
    //    for(int j=0;j<CCE->ClusInfo[i][0][0];j++){
    //        printf("CCE->ClusInfo[%d][%d] :",i,j); 
    //        for(int k=0;k<i+1;k++){
    //            printf("%d ",CCE->ClusInfo[i][j][k]);
    //        }
    //        printf("\n");
    //    }
    //    printf("\n");
    //}

    freeDouble2d(&Duo13CwithCS, int(Duo13CwithCS[0][0]));
	freeInt2d(&InterBathPairs,InterBathPairs[0][0]);
	freeInt2d(&Groups,Groups[0][0]);
	freeInt2d(&GroupChecker,GroupChecker[0][0]);


    MPI_Barrier(MPI_COMM_WORLD);

}

//double calCouplingStrength(double* spin1, double* spin2){
//
//    double r = dist(spin1,spin2);
//    double st12 = sinTheta(spin1,spin2,r);
//    double ct12 = cosTheta(spin1,spin2,r);
//    double sp12 = sinPhi(spin1,spin2);
//    double cp12 = cosPhi(spin1,spin2);
//
//    double Axx = 1 - 3 * pow(st12,2) * pow(cp12,2);
//    double Ayy = 1 - 3 * pow(st12,2) * pow(sp12,2);
//    double Azz = 1 - 3 * pow(ct12,2);
//
//    return pow(Axx,2) + pow(Ayy,2) + pow(Azz,2);
//
//}

void makeDuo13CwithCS(double*** Duo13CwithCS, int nspin, int*** cmap, float*** stmap ){

    // copy Duo13C to (*Duo13CwithCS) with coupling strength for last column
    //
    int idx_duo13C = 1;
    (*Duo13CwithCS) =  allocDouble2d(idx_duo13C,3);
    (*Duo13CwithCS)[0][0] = idx_duo13C;
    (*Duo13CwithCS)[0][1] = 0.0;
    (*Duo13CwithCS)[0][2] = 0.0;


    for (int i = 0; i < nspin; i++){
        for (int j = i+1; j < nspin; j++){
            if ((*cmap)[i][j] != 0){
                reallocDouble2d(Duo13CwithCS,idx_duo13C,idx_duo13C+1,3);
                (*Duo13CwithCS)[idx_duo13C][0] = (double)i+1; // Index is larger than 1 compared to actual index
                (*Duo13CwithCS)[idx_duo13C][1] = (double)j+1; // Index is larger than 1 compared to actual index
                (*Duo13CwithCS)[idx_duo13C][2] = (double)((*stmap)[i][j]);
                idx_duo13C++;
            }
        }
    }

    (*Duo13CwithCS)[0][0] = idx_duo13C;
    (*Duo13CwithCS)[0][1] = 0.0;
    (*Duo13CwithCS)[0][2] = 0.0;

    //if (idx_duo13C != nspin + 1){
    //    printf("Error! idx_duo13C(%d) != nspin + 1 (%d + 1) \n", idx_duo13C,nspin);
    //}

}

void jointSpinSpin(int** group, int len, int spin1, int spin2){

    (*group) = allocInt1d(len);

    (*group)[0] = 3;
    (*group)[1] = spin1;
    (*group)[2] = spin2;
}

void jointGroupSpin(int** group, int spin){

    int current_nspin = (*group)[0] - 1;
    int new_nspin = current_nspin + 1;

    (*group)[0] = new_nspin + 1;
    (*group)[new_nspin] = spin;
    
}

void jointGroupGroup(int** group1, int** group2){

    int current_nspin_group1 = (*group1)[0] - 1;
    int current_nspin_group2 = (*group2)[0] - 1;

    int new_nspin_group1 = current_nspin_group1 + current_nspin_group2 ;

    (*group1)[0] = new_nspin_group1 + 1;
//
    for (int i = 1; i < current_nspin_group2 + 1; i++){
        (*group1)[current_nspin_group1 + i] = (*group2)[i];
        (*group2)[i] = 0; //Changed
    }

    (*group2)[0] = 1;
}

void updateGroupChecker(int*** GroupChecker, int* newGroup, int newGroupNumber){

    int nSpin = newGroup[0] - 1;
    for (int i=1; i< newGroup[0]; i++){
        int spin = newGroup[i];
        (*GroupChecker)[spin][0] = nSpin;
        (*GroupChecker)[spin][1] = newGroupNumber;
    }

}

int Delete_Int_2dArr(int*** arr, int row, int col, int delrow){

    memmove(&((*arr)[delrow]), &((*arr)[delrow+1]), sizeof(int*)*(row-delrow-1));
    row--;

    *arr = (int**)realloc(*arr, sizeof(int*)*row);

    return row;
}

int addNewGroup(int*** Groups, int* group, int order){
    // Return group Number of new group

    int current_nGroup = (*Groups)[0][0] - 1;
    int new_nGroup = current_nGroup + 1;

    reallocInt2d(Groups,current_nGroup+1, new_nGroup+1,order+1);

    (*Groups)[0][0] = new_nGroup + 1;
    (*Groups)[new_nGroup] = group;

    return new_nGroup;
}

void addInterBathPair(int*** Pairs, int spin1, int spin2){

    int current_nPair = (*Pairs)[0][0] - 1;
    int new_nPair = current_nPair + 1;

    reallocInt2d(Pairs,current_nPair+1, new_nPair+1,2);

    (*Pairs)[0][0] = new_nPair + 1;
    (*Pairs)[new_nPair][0] = spin1;
    (*Pairs)[new_nPair][1] = spin2;

}

//void writeClusters(int** arr, int order, const char* fname){
//
//	FILE *fp2 = NULL;
//    fp2 = fopen(fname,"w");
//
//    int nCluster = arr[0][0];
//    for (int j=0; j<nCluster; j++){
//        for (int k = 0; k<order+1; k++){
//            fprintf(fp2,"%d ",arr[j][k]);
//        }
//        fprintf(fp2,"\n");
//    }
//    fclose(fp2);
//
//}

//void writeCoherenceEachCluster(double _Complex* arr, int len, const char* fname){
//
//	FILE *fp2 = NULL;
//    fp2 = fopen(fname,"w");
//
//    for (int j=0; j<len; j++){
//        fprintf(fp2,"%30.20lf \t %30.20lf\n",creal(arr[j]),cimag(arr[j]));  
//    }
//    fclose(fp2);
//}

//void writeBathSpins(InputData* ID, Gyro* Gy){
//
//    //write bath
//
//    printf("\n\t-----Writing Bath spins-----\n");
//
//    FILE *fp = NULL;
//    char filename[100];
//    sprintf(filename,"./BathSpins.csv");
//    fp = fopen(filename,"w");
//
//    for (int i=0; i<(int)ID->Mono13C[0][0]; i++){
////	        printf("%15.5lf\t%15.5lf\t%15.5lf\t%15.5lf\n"
////			,ID->Mono13C[i][0],ID->Mono13C[i][1],ID->Mono13C[i][2]
////	   	    ,ID->Mono13C[i][3]);
//
//
//		if (i!=0){
//	        fprintf(fp,"%15.5lf\t%15.5lf\t%15.5lf\t%15.5lf\n"
//			,ID->Mono13C[i][0],ID->Mono13C[i][1],ID->Mono13C[i][2]
//	   	    ,Gy->Atom_gyro[(int)ID->Mono13C[i][3]]);
//		}
//		else{
//	        fprintf(fp,"%10.5lf\t%10.5lf\t%10.5lf%10.5lf\n"	
//	        ,ID->Mono13C[i][0]
//			,ID->Mono13C[i][1]
//			,ID->Mono13C[i][2]
//    	    ,ID->Mono13C[i][3]);
//
//		}
//    }
//    fclose(fp);
//}


void sort2DArray(int*** arr,int ROWS, int COLS, int col){
    for (int i = 0; i < ROWS - 1; i++) {
        for (int j = 0; j < ROWS - i - 1; j++) {
            if ((*arr)[j][col] > (*arr)[j + 1][col]) {
                // ÌäπÏ†ï Ïó¥ÏùÑ Í∏∞Ï§ÄÏúºÎ°ú ÌñâÏùÑ ÍµêÌôò
                for (int k = 0; k < COLS; k++) {
                    int temp = (*arr)[j][k];
                    (*arr)[j][k] = (*arr)[j + 1][k];
                    (*arr)[j + 1][k] = temp;
                }
            }
            else if ((*arr)[j][col] == (*arr)[j + 1][col]) {
                if ((*arr)[j][col+1] > (*arr)[j + 1][col+1]) {
                    // ÌäπÏ†ï Ïó¥ÏùÑ Í∏∞Ï§ÄÏúºÎ°ú ÌñâÏùÑ ÍµêÌôò
                    for (int k = 0; k < COLS; k++) {
                        int temp = (*arr)[j][k];
                        (*arr)[j][k] = (*arr)[j + 1][k];
                        (*arr)[j + 1][k] = temp;
                    }
                }
            }

        }
    }
}


void QuickSort_Double_2d(double ***arr, int col, int left, int right) {
    int i = left, j = right;
    double pivot = (*arr)[(left + right) / 2][col - 1]; // pivotÏùÄ Î∞∞Ïó¥Ïùò Í∞ÄÏö¥Îç∞ Í∞íÏúºÎ°ú ÏÑ§Ï†ï

    // partition
    while (i <= j) {
        while ((*arr)[i][col - 1] > pivot) i++;
        while ((*arr)[j][col - 1] < pivot) j--;
        if (i <= j) {
            // swap
            double *temp = (*arr)[i];
            (*arr)[i] = (*arr)[j];
            (*arr)[j] = temp;
            i++;
            j--;
        }
    }

    // recursion
    if (left < j) QuickSort_Double_2d(arr, col, left, j);
    if (i < right) QuickSort_Double_2d(arr, col, i, right);
}

void shrinkClusterO1(Cluster* cls){

    int i=1;

    int write_j = 1;
    int origin_line = cls->clusinfo[i][0][0];
    for (int j = 1; j < origin_line; j++) {
        int* current = cls->clusinfo[i][j];
        if (current != NULL && current[0] != 0) {
            // copy to write_j
            if (write_j != j) {
                cls->clusinfo[i][write_j] = cls->clusinfo[i][j];
            }
            write_j++;
        } else {
            // free
            free(cls->clusinfo[i][j]);
        }
    }

    cls->clusinfo[i][0][0] = write_j;
    cls->clusinfo[i] = (int**)realloc(cls->clusinfo[i], write_j * sizeof(int*));
}
