#include "../../include/common/output.h"
#include "../../include/common/memory.h"
#include "../../include/common/utilities.h"
#include <string>

// A result file that could not be written in full has to stop the run: a full disk, a quota
// or a file-size limit used to leave a truncated file behind and end with exit status 0.
static void closeOutput(FILE* fp, const char* path){
    const bool failed = ferror(fp) != 0;
    if (fclose(fp) != 0 || failed){
        fprintf(stderr, "Error: cannot write file %s\n", path);
        ccexFatal();   // rank 0 alone writes: the others must not wait at the next barrier (D38)
    }
}

Output* Output_init(){
    Output* op = (Output*)allocArray1d(1,sizeof(Output));
    op->savemode[0] = '\0';
    op->outfile = NULL;
    op->save_nodiv = true;
    return op;
}

void Output_save(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, double deltat, int istate){

    char* outfile = Output_getOutfile(op);

    // istate to string
    // "_state" and any int fit in 32; the old buffer of 10 overflowed from _state1000 on (review 78d583 R1)
    char istate_str[32];
    if (istate==0){
        istate_str[0] = '\0';
    }
    else{
        snprintf(istate_str, sizeof(istate_str), "_state%d", istate);
    }

    int dim = result_wD[0].rows();

    char* outfile_wD = allocChar1d(MAX_FILEPATH);
    char* outfile_nD = allocChar1d(MAX_FILEPATH);

    strcpy(outfile_wD, outfile);
    strcat(outfile_wD, istate_str);
    strcat(outfile_wD, "_wiDiv");

    strcpy(outfile_nD, outfile);
    strcat(outfile_nD, istate_str);
    strcat(outfile_nD, "_noDiv");

    FILE* fp_wD = fopen(outfile_wD, "w");
    FILE* fp_nD = op->save_nodiv ? fopen(outfile_nD, "w") : NULL;   // savemode.nodiv false: no _noDiv file (2026-10-08)
    
    if (fp_wD == NULL || (op->save_nodiv && fp_nD == NULL)){
        printf("Error: cannot open file %s\n", outfile_wD);
        printf("Error: cannot open file %s\n", outfile_nD);
        ccexFatal(1);
    }

    // write format :
    // for all nstate , write time dm00 dm01 dm02 ... dmnn (write only)
    for (int i=0; i<nstep; i++){
        int c = 0;
        fprintf(fp_wD, "%15g\t", i*deltat);
        if (fp_nD) fprintf(fp_nD, "%15g\t", i*deltat);
        for (int j=0; j<dim; j++){
            for (int k=0; k<dim; k++){

                fprintf(fp_wD, "%+15.10lf", result_wD[i](j,k).real());
                if (fp_nD) fprintf(fp_nD, "%+15.10lf", result_nD[i](j,k).real());

                if (std::isnan(result_wD[i](j,k).imag())){
                    fprintf(fp_wD, "%+-4lfj%10s\t", result_wD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_wD, "%+-13.10lfj\t", result_wD[i](j,k).imag());
                }

                if (std::isnan(result_nD[i](j,k).imag())){
                    if (fp_nD) fprintf(fp_nD, "%+-4lfj%10s\t", result_nD[i](j,k).imag()," ");
                }else{
                    if (fp_nD) fprintf(fp_nD, "%+-13.10lfj\t", result_nD[i](j,k).imag());
                }

                if ((c%3==2) && (j*k!=(dim-1)*(dim-1))){
                    fprintf(fp_wD, "\n");
                    fprintf(fp_wD, "%15s\t"," ");

                    if (fp_nD) fprintf(fp_nD, "\n");
                    if (fp_nD) fprintf(fp_nD, "%15s\t"," ");
                }

                c++;
            }
        }
        fprintf(fp_wD, "\n");
        if (fp_nD) fprintf(fp_nD, "\n");
    }

    closeOutput(fp_wD, outfile_wD);
    if (fp_nD) closeOutput(fp_nD, outfile_nD);

    freeChar1d(&outfile_wD);
    freeChar1d(&outfile_nD);
    
    return;
}

void Output_save_all(Output* op, MatrixXcd* result, int nstep, double deltat, int* cluster, int nspin, int istate){

    //////////////////////////////
    // Save the original one
    char* outfile_original = allocChar1d(MAX_FILEPATH);
    strcpy(outfile_original, Output_getOutfile(op));
    //////////////////////////////

    //////////////////////////////
    // Change outfile name

    // cluster to string (built to its length: a fixed 100-character label overflowed for large clusters with long spin
    // indices; hard-coding check 2026-10-09), and the whole name has to fit the outfile buffers with its suffixes
    std::string cluster_str;
    for (int i=0; i<nspin; i++){
        cluster_str += "_" + std::to_string(cluster[i]);
    }
    if (strlen(op->outfile) + cluster_str.size() + 64 >= MAX_FILEPATH){
        fprintf(stderr, "Error: savemode cluster_coherence: the file name of cluster %s is longer than %d characters; use a shorter outfile\n",
                cluster_str.c_str(), MAX_FILEPATH - 1);
        ccexFatal(1);
    }

    // add string outfilehead, cluster, istate to outfile
    strcat(op->outfile, cluster_str.c_str());
    
    Output_save(op, result, result, nstep, deltat, istate);
    //////////////////////////////

    //////////////////////////////
    // Restore outfile name
    strcpy(op->outfile, outfile_original);
    //////////////////////////////

    freeChar1d(&outfile_original);

    return;
}

Output Output_fullOutput(const Output* op, char* name, const int* cluster, int n, int istate){
    std::string full(op->outfile);
    for (int i=0; i<n; i++){
        full += "_" + std::to_string(cluster[i]);
    }
    full += "_full";
    if (istate!=0) full += "_state" + std::to_string(istate);
    if (full.size() + 16 >= MAX_FILEPATH){   // Output_save appends _wiDiv / _noDiv
        fprintf(stderr, "Error: savemode cluster_density_matrix: the file name %s is longer than %d characters; use a shorter outfile\n",
                full.c_str(), MAX_FILEPATH - 17);
        ccexFatal(1);
    }
    strcpy(name, full.c_str());
    Output out = *op;
    out.outfile = name;
    return out;
}

void Output_save_info(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, double deltat, int istate, double Azx, double Azz, char* bathfile){

    char* savemode = Output_getSavemode(op);

    if (strcasecmp(savemode, "coherence_tagged") != 0){
        printf("Error: savemode is not coherence_tagged\n");
        ccexFatal(1);
    }

    char* outfile = Output_getOutfile(op);

    // istate to string
    // "_state" and any int fit in 32; the old buffer of 10 overflowed from _state1000 on (review 78d583 R1)
    char istate_str[32];
    if (istate==0){
        istate_str[0] = '\0';
    }
    else{
        snprintf(istate_str, sizeof(istate_str), "_state%d", istate);
    }


    int dim = result_wD[0].rows();

    char* outfile_wD = allocChar1d(MAX_FILEPATH);
    char* outfile_nD = allocChar1d(MAX_FILEPATH);

    strcpy(outfile_wD, outfile);
    strcat(outfile_wD, istate_str);
    strcat(outfile_wD, "_wiDiv");

    strcpy(outfile_nD, outfile);
    strcat(outfile_nD, istate_str);
    strcat(outfile_nD, "_noDiv");

    FILE* fp_wD = fopen(outfile_wD, "w");
    FILE* fp_nD = op->save_nodiv ? fopen(outfile_nD, "w") : NULL;   // savemode.nodiv false: no _noDiv file (2026-10-08)
    
    if (fp_wD == NULL || (op->save_nodiv && fp_nD == NULL)){
        printf("Error: cannot open file %s\n", outfile_wD);
        printf("Error: cannot open file %s\n", outfile_nD);
        ccexFatal(1);
    }

    for (int i=0; i<nstep; i++){
        int c = 0;
        fprintf(fp_wD, "%15g\t", i*deltat);
        if (fp_nD) fprintf(fp_nD, "%15g\t", i*deltat);
        for (int j=0; j<dim; j++){
            for (int k=0; k<dim; k++){

                fprintf(fp_wD, "%+15.10lf", result_wD[i](j,k).real());
                if (fp_nD) fprintf(fp_nD, "%+15.10lf", result_nD[i](j,k).real());

                if (std::isnan(result_wD[i](j,k).imag())){
                    fprintf(fp_wD, "%+-4lfj%10s\t", result_wD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_wD, "%+-13.10lfj\t", result_wD[i](j,k).imag());
                }

                if (std::isnan(result_nD[i](j,k).imag())){
                    if (fp_nD) fprintf(fp_nD, "%+-4lfj%10s\t", result_nD[i](j,k).imag()," ");
                }else{
                    if (fp_nD) fprintf(fp_nD, "%+-13.10lfj\t", result_nD[i](j,k).imag());
                }

                if ((c%3==2) && (j*k!=(dim-1)*(dim-1))){
                    fprintf(fp_wD, "\n");
                    fprintf(fp_wD, "%15s\t"," ");

                    if (fp_nD) fprintf(fp_nD, "\n");
                    if (fp_nD) fprintf(fp_nD, "%15s\t"," ");
                }

                c++;
            }
        }
        // write Azx, Azz, bathfile
        fprintf(fp_wD, "%+.10lf\t%+.10lf\t%s\n", Azx, Azz, bathfile);
        if (fp_nD) fprintf(fp_nD, "%+.10lf\t%+.10lf\t%s\n", Azx, Azz, bathfile);
    }

    closeOutput(fp_wD, outfile_wD);
    if (fp_nD) closeOutput(fp_nD, outfile_nD);
    
    freeChar1d(&outfile_wD);
    freeChar1d(&outfile_nD);
    return;
}

void Output_setSavemode(Output* op, char* savemode){

    // the names before 2026-10-13 are refused with the name that replaced them (docs/user/InputFile_Output.md)
    const char* oldNames[] = {"normal", "info", "all", "allfull"};
    const char* names[] = {"coherence", "coherence_tagged", "cluster_coherence", "cluster_density_matrix"};
    for (int k=0; k<4; ++k)
        if (strcasecmp(savemode, oldNames[k]) == 0){
            ccexRefusalPrint(stderr, "Error: savemode \"%s\" is now \"%s\"\n", oldNames[k], names[k]);
            ccexRefuse();   // read from the input on every rank (A55)
        }
    int k = 0;
    while (k<4 && strcasecmp(savemode, names[k])) ++k;
    if (k == 4) {
        ccexRefusalPrint(stderr, "Error: savemode \"%s\" is not one of coherence, coherence_tagged, cluster_coherence, "
                        "cluster_density_matrix\n", savemode);
        ccexRefuse();
    }
    strcpy(op->savemode, savemode);   // one of the names above (any case), which fits
}

char* Output_getSavemode(Output* op){
    return op->savemode;
}

void Output_allocOutfile(Output* op){
    op->outfile = allocChar1d(MAX_FILEPATH);
}

void Output_freeOutfile(Output* op){
    freeArray1d((void**)&(op->outfile));
}

void Output_setOutfile(Output* op, char* outfile){
    // Output_save and Output_save_all append a cluster label (under 100
    // characters), "_state<n>" and "_wiDiv" to this inside buffers of the same
    // size, so the name itself has to leave them room.
    copyBounded(op->outfile,outfile,MAX_FILEPATH-128,"outfile");
}

char* Output_getOutfile(Output* op){
    return op->outfile;
}

void Output_freeAll(Output* op){
    Output_freeOutfile(op);
    freeArray1d((void**)&op);
}

void Output_report(Output* op){

    printTitle("Structure Output");

    printStructElementChar("savemode", op->savemode);
    printStructElementChar("outfile", op->outfile);
    printStructElementBool("savemode.nodiv", op->save_nodiv);

}
