#include "../../include/input/reader.h"
#include "../../include/input/isotope_table.h"
#include "../../include/common/memory.h"
#include "../../include/common/utilities.h"
#include "../../include/clustering/cluster.h"
#include "../../include/calculation/hamiltonian.h"

static int addIsotopeSpecies(BathArray* ba, const char* name);   // below readGyrofile

void readQubitfile(QubitArray* qa, Config* cnf){

    char* fname = Config_getQubitfile(cnf);
    double* xyz = QubitArray_getQubit_i_xyz(qa,0); // only add the info to 0-th qubit

    // file open
    FILE* data;

    data    =   fopen(fname,"r");

    // A qubitfile the input names has to be read: a typo in its path used to fall through to "already set" and run
    // the qubit at (0,0,0) (review 4dd5fd R1). Only without the key (an empty name) does the Qubit block stand alone.
    if (data == NULL && fname != NULL && fname[0] != '\0'){
        ccexRefusalPrint(stderr,"Error(readQubitFile): cannot open the qubit file '%s': %s\n",fname,strerror(errno));
        ccexRefuse();
    }
    if (data == NULL && QubitArray_getNqubit(qa) == 0){
        ccexRefusalPerror("Error(readQubitFile): Cannot open the qubit file\n");
        ccexRefusalPrint(stdout,"Current WRONG qubitfile name : %s\n",fname);
        ccexRefusalPrint(stdout,"Check if the qubit file exists\n");
        ccexRefuse();
    }else if (data == NULL && QubitArray_getNqubit(qa) != 0){
        if (rank==0){printMessage("QubitArray already set with option ...\n");}
        return;
    }else{
        if (rank==0){
            char message[MAX_FILEPATH];
            snprintf(message,sizeof message,"Read QubitFile : %s ... ",fname);
            printMessage(message);
        }
        
    }

    // read xyz from QubitFile
    size_t max = 1000;
    char lline[max];
    int count = 0;
    while(!feof(data)){
        fgets(lline, max, data);
        count = sscanf(lline, "%lf %lf %lf\n", &(xyz[0]), &(xyz[1]), &(xyz[2]));
        if (count == 3){
            fclose(data);

            if (rank==0){
                printStructElementDouble1d("  Qubit[0].xyz ", xyz, 3);
                printf("\n");
            }
            return;
        }
    }
    fclose(data);

    // if no data
    ccexRefusalPerror("Error(readQubitFile): Cannot read the qubit file\n");
    ccexRefusalPrint(stdout,"Current WRONG qubitfile name : %s\n",fname);
    ccexRefusalPrint(stdout,"Check if your qubit file has 'xyz' values in a right form : x y z \n");
    ccexRefuse();
}

/**
 * @brief Stops when two qubits share a name, since reference_qubit resolves by exact name
*/
static void checkQubitNamesUnique(QubitArray* qa){

    int nqubit = QubitArray_getNqubit(qa);
    for (int i=0; i<nqubit; i++){
        for (int j=i+1; j<nqubit; j++){
            if (strcmp(QubitArray_getQubit_i_name(qa,i),QubitArray_getQubit_i_name(qa,j)) == 0){
                ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): qubit[%d] and qubit[%d] share the name \"%s\"\n",
                        i,j,QubitArray_getQubit_i_name(qa,i));
                ccexRefusalPrint(stderr,"reference_qubit resolves by exact name, so the names must be unique.\n");
                ccexRefuse();
            }
        }
    }
}

/**
 * @brief Index of the qubit the rotation is taken about
 * @details "reference_qubit" names it; without the key the first qubit is used, which
 *          is unambiguous for as long as nqubit is restricted to 1. The lookup is by
 *          exact name, so duplicate names are refused first -- otherwise a name would
 *          silently resolve to whichever qubit happened to be listed earlier.
 *
 *          Written to survive the nqubit=1 restriction being lifted: nothing here
 *          assumes there is only one qubit.
*/
static int resolveReferenceQubit(QubitArray* qa, Config* cnf){

    int nqubit = QubitArray_getNqubit(qa);
    checkQubitNamesUnique(qa);

    char* refname = Config_getRot_reference_qubit(cnf);
    if (refname[0] == '\0'){ return 0; } // not given : the first qubit

    int iref = QubitArray_getQubitIdx_fromName(qa,refname);
    if (iref < 0){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): reference_qubit \"%s\" is not one of the qubits\n",refname);
        ccexRefusalPrint(stderr,"Defined qubit names are :");
        for (int i=0; i<nqubit; i++){ ccexRefusalPrint(stderr," \"%s\"",QubitArray_getQubit_i_name(qa,i)); }
        ccexRefusalPrint(stderr,"\n");
        ccexRefusalPrint(stderr,"A single qubit given through \"Qubit.file\" is named \"q0\" by default.\n");
        ccexRefuse();
    }
    return iref;
}

/**
 * @brief Check the coordinate_frame_rotation input against the QubitArray
 * @details PURE VALIDATION. It opens nothing, reads no file and changes no state; it
 *          only compares what the input file and the command line have already produced.
 *          Called from main once every option source is in but before readQubitfile,
 *          readBathfiles, readHftensorfile and readQdtensorfile -- so an unsupported
 *          combination stops the run before a single tensor file is opened.
 *
 *          The checks live here rather than in cJSON_readOptionConfig because they need
 *          BOTH structures: the Config knows the rotation and the readmodes, the
 *          QubitArray knows how many qubits there are and what they are called, and
 *          those are parsed by two different functions.
 *
 *          The multi-qubit restriction on external HF and QD tensor files is NOT new
 *          here. readHftensorfile and readQdtensorfile have always refused nqubit > 1,
 *          because one tensor file is written around one defect and applying it to
 *          several qubits has no defined meaning. Those guards stay where they are; this
 *          one only reports the same thing earlier and in terms of the rotation.
 *
 *          A no-op unless the rotation is enabled, so nothing about an existing input
 *          changes.
*/
void validateCoordinateFrameRotationInputs(QubitArray* qa, Config* cnf){

    if (!Config_getRot_enabled(cnf)){ return; }

    int nqubit = QubitArray_getNqubit(qa);
    if (nqubit < 1){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): no qubit is defined (nqubit = %d)\n",nqubit);
        ccexRefuse();
    }

    checkQubitNamesUnique(qa);

    ////////////////////////////////////////////////////////////////////////
    // Single qubit : everything below is optional, exactly as before
    ////////////////////////////////////////////////////////////////////////
    if (nqubit == 1){
        // reference_qubit and qubit_position_frame may both be omitted -- there is one
        // qubit to rotate about and one frame it can be in. An explicit reference_qubit
        // is still resolved (and still has to match) when the rotation runs.
        return;
    }

    ////////////////////////////////////////////////////////////////////////
    // Several qubits : nothing may be left to a default
    ////////////////////////////////////////////////////////////////////////
    if (Config_getRot_reference_qubit(cnf)[0] == '\0'){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): nqubit = %d, so reference_qubit is required\n",nqubit);
        ccexRefusalPrint(stderr,"It names the qubit the rotation is taken about, and with more than one qubit\n");
        ccexRefusalPrint(stderr,"there is no sensible default. Defined qubit names are :");
        for (int i=0; i<nqubit; i++){ ccexRefusalPrint(stderr," \"%s\"",QubitArray_getQubit_i_name(qa,i)); }
        ccexRefusalPrint(stderr,"\n");
        ccexRefuse();
    }

    // Resolves by exact name and is fatal if it does not match; the index is discarded
    // here because this function only validates.
    (void)resolveReferenceQubit(qa,cnf);

    if (Config_getRot_qubit_position_frame(cnf)[0] == '\0'){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): nqubit = %d, so qubit_position_frame is required\n",nqubit);
        ccexRefusalPrint(stderr,"State the frame every Qubit.xyz is written in. The only supported value is\n");
        ccexRefusalPrint(stderr,"\"bath\": the same Cartesian basis and the same origin as the bath files, so that\n");
        ccexRefusalPrint(stderr,"one rotation moves every qubit and every bath spin together.\n");
        ccexRefuse();
    }

    ////////////////////////////////////////////////////////////////////////
    // External tensor files : a pre-existing restriction, caught earlier
    ////////////////////////////////////////////////////////////////////////
    if (Config_getHf_readmode(cnf) != 0){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): a hyperfine tensor file (hamiltonian.common.hyperfine) with %d qubits is not supported\n",
                nqubit);
        ccexRefusalPrint(stderr,"readHftensorfile already refuses more than one qubit with a hyperfine tensor file:\n");
        ccexRefusalPrint(stderr,"the file is written around a single defect, so which qubit it belongs to is not\n");
        ccexRefusalPrint(stderr,"defined. Use hyperfine dipole \"point\" with contact false, or run one qubit at a time.\n");
        ccexRefuse();
    }
    if (Config_getQd_readmode(cnf) != 0){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): a quadrupole tensor file (hamiltonian.common.quadrupole) with %d qubits is not supported\n",
                nqubit);
        ccexRefusalPrint(stderr,"readQdtensorfile already refuses more than one qubit with a quadrupole tensor file,\n");
        ccexRefusalPrint(stderr,"for the same reason. Use quadrupole \"off\", or run one qubit at a time.\n");
        ccexRefuse();
    }

    ////////////////////////////////////////////////////////////////////////
    // Every explicit intmap tensor has to say which frame it is in
    ////////////////////////////////////////////////////////////////////////
    for (int i=0; i<nqubit; i++){
        for (int j=i; j<nqubit; j++){
            if (QubitArray_getIntmap_i_j_frame(qa,i,j) != INTMAP_EXPLICIT_UNSPECIFIED){ continue; }
            ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): the intmap tensor between \"%s\" and \"%s\"\n",
                    QubitArray_getQubit_i_name(qa,i),QubitArray_getQubit_i_name(qa,j));
            ccexRefusalPrint(stderr,"has no \"tensor_frame\", and with nqubit = %d it cannot be guessed.\n",nqubit);
            ccexRefusalPrint(stderr,"  \"bath\"  : components are in the original bath Cartesian basis, and are\n");
            ccexRefusalPrint(stderr,"            transformed as R*T*R^T\n");
            ccexRefusalPrint(stderr,"  \"qubit\" : components are already in the common computational basis --\n");
            ccexRefusalPrint(stderr,"            reference-qubit-aligned, NOT an individual local frame -- and\n");
            ccexRefusalPrint(stderr,"            are left as written\n");
            ccexRefuse();
        }
    }
}

void readBathfiles(BathArray* ba, QubitArray* qa, Config* cnf){

    char message[MAX_FILEPATH];

    int     nqubit      = QubitArray_getNqubit(qa);
    int     nbathfiles  = Config_getNbathfiles(cnf);
    double  rbath       = Config_getRbath(cnf);
    double  rbathcut    = Config_getRbathcut(cnf);
    
    // bath spin properties
    int     nspecies    = BathArray_getProp_nspecies(ba);
    char**  names       = BathArray_getProp_names(ba);
    float*  spins       = BathArray_getProp_spins(ba);
    double* gyros       = BathArray_getProp_gyros(ba);

    // Read bath files
    int nspin_prev = 0;
    int nspin = 0;

    for (int i=0; i<nbathfiles; i++){

        FILE*   data;
        char*   fname       = Config_getBathfiles_i(cnf,i);
        double* bathadjust  = Config_getBathadjust_i(cnf,i);

        // file open        
        data = fopen(fname,"r");
        if(data==NULL){
            ccexRefusalPerror("Error(readBathFile): Cannot open the bath file\n");
            ccexRefusalPrint(stdout,"Current WRONG bathfile name : %s\n",fname);
            ccexRefusalPrint(stdout,"Check if the bath file exists\n");
            ccexRefuse();
        }

        // read xyz from BathFile
        double xyz[3] = {0.0,};
        char name[100];

        // Read first line (the number of spins)
        int fline = 1; // the first line is the number of spins
        fscanf(data, "%*f %*f %*f %*s\n");
        //
        //double nline_d;
        //fscanf(data, "%lf %*lf %*lf %*s\n", &nline_d);
        //int nspins_inBathFile = ( (int)nline_d ) - 1;
        //int lineIdx = 0;
        //Config_alloc_flines(cnf, nspins_inBathFile); // alloc nline
        //printf("nline = %d\n", nspins_inBathFile);
        //
        while(!feof(data)){
        
            int count = fscanf(data, "%lf %lf %lf %99s\n", &xyz[0], &xyz[1], &xyz[2], name);
            xyz[0] += bathadjust[0];
            xyz[1] += bathadjust[1];
            xyz[2] += bathadjust[2];
            fline++; // the line is started from 2

            if (count == 4){
                double r = QubitArray_mindist(xyz, qa);
                if ((r <= rbath) && (r >= rbathcut)){

                    //////////////////////////////
                    // increase the number of spins
                    nspin++;
                    BathArray_setNspin(ba, nspin);                
                    // allocate the BathArray->Bath
                    if (nspin == 1){
                        BathArray_allocBath(ba, nqubit); // alloc bath spins
                        Config_alloc_flines(cnf, nspin); // alloc fline
                    }else{
                        BathArray_reallocBath(ba, nspin-1, nspin, nqubit);
                        Config_realloc_flines(cnf, nspin-1, nspin);
                    }

                    // set the bath
                    BathArray_setBath_i_name(ba, name, nspin-1);
                    BathArray_setBath_i_xyz(ba, xyz, nspin-1);
                    BathArray_setBath_i_mindist(ba, r, nspin-1);

                    // set fline
                    Config_set_nflines(cnf, nspin);
                    Config_set_flines_i(cnf, fline, nspin-1);

                    // set the bath spin properties
                    int ispeceis = findIndexChar(names,0,nspecies-1,name);

                    if (ispeceis == -1){                 // not in the gyrofile: the built-in isotope table
                        ispeceis = addIsotopeSpecies(ba, name);
                        if (ispeceis == -1){
                            ccexRefusalPrint(stderr,"Error(readBathFile): the species \"%s\" is neither in the gyrofile nor in the built-in "
                                    "isotope table (names like 13C, 14N, 11B); give it in the gyrofile (name spin gyro)\n",name);
                            ccexRefusalPrint(stderr,"The bath spin is : %s %lf %lf %lf \n",name,xyz[0],xyz[1],xyz[2]);
                            ccexRefuse();
                        }
                        nspecies = BathArray_getProp_nspecies(ba);    // the property arrays were reallocated
                        names    = BathArray_getProp_names(ba);
                        spins    = BathArray_getProp_spins(ba);
                        gyros    = BathArray_getProp_gyros(ba);
                    }
                    BathArray_setBath_i_spin(ba, spins[ispeceis], nspin-1);
                    BathArray_setBath_i_gyro(ba, gyros[ispeceis], nspin-1);
                    // overlap..?
                } //else{
                  //  Config_set_flines_i(cnf, 0, lineIdx);
                //}
            } else{
                ccexRefusalPerror("Error(readBathFile): BathFile format is NOT x y z name");
                ccexRefuse();
            }
            //lineIdx++;
        }
        fclose(data);

        ////////////////////////////////////////////////////////////////////////
        // Print
        if(rank==0){
            printLine();
            snprintf(message,sizeof message,"Read BathFile[%d] : %s\n",i,fname); printMessage(message);
            snprintf(message,sizeof message,"- Number of bath spins : %d (nspin = %d)",nspin-nspin_prev,nspin); printMessage(message);
            for (int isp=nspin_prev; isp<nspin; isp++){
                if (verbosity || (isp<nspin_prev + 3 || isp>nspin-3)){ 
                    // printf("        ( fline %5d )",Config_get_flines_i(cnf,isp));
                    BathArray_reportBath_i_props(ba, isp);
                }
                if (!verbosity && isp==nspin_prev + 3){
                    snprintf(message,sizeof message,"   : \n"); printMessage(message);
                }
            }
            printLine();

            nspin_prev = BathArray_getNspin(ba);
            printf("\n");
        }
        ////////////////////////////////////////////////////////////////////////
    }
    
    if (Config_get_nflines(cnf) != nspin){
        ccexRefusalPrint(stderr,"Error(readBathFile): The number of bath spins is not consistent\n");
        ccexRefusalPrint(stderr,"The number of bath spins from the file : %d\n",nspin);
        ccexRefusalPrint(stderr,"The number of bath spins from the file line : %d\n",Config_get_nflines(cnf));
        ccexRefuse();
    }

    // qsort(ba->bath, ba->nspin, sizeof(BathSpin*), compare_dist);
    //BathArray_report(ba);
    //for (int i=0; i < 1858; i++){
    //    printf("%d - ", i);
    //    printf("Config_set_flines_i(cnf, %d): ", i);
    //    std::cout << Config_get_flines_i(cnf, i) << std::endl;
    //}
}


/**
 * @brief Print the whole coordinate frame, including every qubit, once per run
 * @details Wraps the pure-geometry report in utilities.cpp and adds what only the
 *          QubitArray can say: which qubit the rotation was taken about, where every
 *          qubit sits in each frame, and how the interaction tensors were handled.
 * @param qa   the qubits, already transformed to the computational frame
 * @param cnf  the run's settings (axes, frames, field)
 * @param R    the rotation, original frame -> computational frame
 * @param r0   the point held fixed, in the original frame
 * @param iref index of the reference qubit the rotation was taken about
 * @param qsrc qubit positions BEFORE the rotation; the QubitArray already holds the
 *             transformed ones by the time this runs
 * @param defect_axis_reference name of the recorded reference axis, or NULL if none
 * @param nint_rotated number of qubit intmap tensors that were transformed
 * @param nint_kept    number of qubit intmap tensors left as they were
 * @param nprov        count of intmap tensors by IntmapProvenance (five kinds)
*/
static void reportCoordinateFrame(QubitArray* qa, Config* cnf, const double R[3][3],
                                  const double r0[3], int iref, double** qsrc,
                                  const char* defect_axis_reference,
                                  int nint_rotated, int nint_kept, const int nprov[5]){

    int nqubit = QubitArray_getNqubit(qa);

    reportQubitAlignedRotation(Config_getRot_bath_axis(cnf),Config_getRot_qubit_axis(cnf),
                               R,r0,Config_getBfield(cnf),defect_axis_reference);

    char* posframe = Config_getRot_qubit_position_frame(cnf);
    printf("      %-32s: %s\n","Qubit position frame",posframe[0] ? posframe : "bath (default)");
    printf("      %-32s  all Qubit.xyz values use the same Cartesian basis and\n","Meaning:");
    printf("      %-32s  origin as the bath files\n","");
    printf("      %-32s  the physical axis of reference_qubit in the source frame\n","Qubit axis:");
    printf("      %-32s  common reference-qubit-aligned computational frame,\n","tensor_frame = \"qubit\":");
    printf("      %-32s  not an individual local frame\n","");
    printf("\n");
    printf("      %-32s: \"%s\"  (index %d of %d)\n","Reference qubit",
           QubitArray_getQubit_i_name(qa,iref),iref,nqubit);
    printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","r0, source frame",r0[0],r0[1],r0[2]);
    printf("\n");

    for (int iq=0; iq<nqubit; iq++){
        double* now = QubitArray_getQubit_i_xyz(qa,iq);
        printf("      qubit[%d] \"%s\"%s\n",iq,QubitArray_getQubit_i_name(qa,iq),
               (iq==iref) ? "   (reference)" : "");
        printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","  source frame",
               qsrc[iq][0],qsrc[iq][1],qsrc[iq][2]);
        printf("      %-32s: [ %12.8lf, %12.8lf, %12.8lf ]\n","  computational frame",
               now[0],now[1],now[2]);
    }
    printf("\n");

    printf("      %-32s: %d transformed as R*T*R^T, %d left as written\n",
           "Qubit intmap tensors",nint_rotated,nint_kept);
    printf("      %-32s: %d default-zero, %d auto (source geometry),\n",
           "  by provenance",nprov[INTMAP_DEFAULT_ZERO],nprov[INTMAP_AUTO_SOURCE_GEOMETRY]);
    printf("      %-32s  %d explicit bath, %d explicit qubit, %d explicit unspecified\n","",
           nprov[INTMAP_EXPLICIT_BATH],nprov[INTMAP_EXPLICIT_QUBIT],nprov[INTMAP_EXPLICIT_UNSPECIFIED]);
    printf("\n");
}

/**
 * @brief Move the bath into the qubit-aligned frame
 * @details One pipeline step of its own, run once the bath positions and every stored
 *          tensor are in place. Doing it here rather than inside readBathfiles keeps
 *          the frame change to a single moment: bath POSITIONS are in the source frame
 *          until this runs and in the qubit-aligned frame after it, and the tensors
 *          that are still in the bath frame are transformed in the same pass, so
 *          positions and the tensors beside them never disagree.
 *
 *          Tensor COMPONENTS are not all in one frame beforehand: hf_tensor_frame and
 *          qd_tensor_frame say which basis the tensor FILES use, and a point-dipole
 *          tensor CCEX computed itself is in the bath frame whatever they say. This
 *          step normalizes them, transforming the bath-frame ones and leaving the
 *          qubit-frame ones alone. hypf_bathframe is what tells the two apart.
 *
 *          Defect entries carry their own coordinate_frame. Their axis and rxyzs are
 *          vectors, while hypf/efg/zfs are rank-2 tensors; bath-frame entries receive
 *          the same R here before Defect construction. Scalar detuning does not move.
 *
 *          The tensor-file lookups happen before this, against un-moved coordinates,
 *          because the position columns of those files are in the source bath frame.
 *
 *          The magnetic field is parsed independently of this transformation. After all
 *          CLI overrides are applied, Config_validateBfieldAlignment checks that it lies
 *          along the computational z axis.
 *
 *          rbath / rbathcut selection happens earlier, in the source frame, which is
 *          equivalent: the rotation is taken about the reference qubit and an
 *          orthogonal transformation about that point leaves every distance to it
 *          unchanged. mindist is recomputed here anyway, so the value stored on each
 *          spin describes the coordinates actually stored next to it.
 *
 *          A no-op unless coordinate_frame_rotation is enabled.
*/
void applyCoordinateFrameRotation(BathArray* ba, QubitArray* qa, DefectArray* dfa, Config* cnf,
                                  HypfProvenance hypf_bathframe){

    if (!Config_getRot_enabled(cnf)){ return; }

    double R[3][3];
    buildQubitAlignedRotation(Config_getRot_bath_axis(cnf),Config_getRot_qubit_axis(cnf),R);

    int iref = resolveReferenceQubit(qa,cnf);
    int nqubit = QubitArray_getNqubit(qa);

    // COPY the origin before anything moves: r0 aliases Qubit[iref].xyz otherwise, and
    // the first rotateAboutPoint would rewrite the very point it rotates about.
    double r0[3];
    copyDouble1d(r0,QubitArray_getQubit_i_xyz(qa,iref),3);

    // Kept for the report, which shows every qubit in both frames.
    double** qsrc = allocDouble2d(nqubit,3);
    for (int iq=0; iq<nqubit; iq++){
        copyDouble1d(qsrc[iq],QubitArray_getQubit_i_xyz(qa,iq),3);
    }

    // Qubits first, bath second. mindist is the distance to the NEAREST qubit, so it can
    // only be recomputed once every qubit has moved. The reference qubit has r_old = r0
    // and therefore does not move at all.
    for (int iq=0; iq<nqubit; iq++){
        rotateAboutPoint(QubitArray_getQubit_i_xyz(qa,iq),R,r0);
    }


    // Which stored tensors are in the bath frame and therefore have to move with it.
    //
    // For the hyperfine tensors this is decided PER SPIN, not per run, because a single
    // run can hold both kinds. hf_tensor_frame describes only what came out of the
    // tensor file; anything readHftensorfile computed itself -- a point-dipole tensor,
    // whether as the whole answer (readmode 0), as the anisotropic part next to an
    // isotropic Fermi contact (readmode 1), or as the fallback for a spin the file did
    // not cover (readmode 2/3) -- is in the source frame regardless. hypf_bathframe
    // carries that distinction over from the reader.
    //
    // The point-dipole tensor is covariant, so R*T*R^T is the same thing as recomputing
    // it from the rotated positions; test_rotation_unit.cpp checks that identity.
    //
    // The quadrupole side stays a single flag: readQdtensorfile has no computed-here
    // fallback to distinguish, and qd_readmode = 0 leaves every quad at zero.
    bool hf_file_is_bathframe = (strcasecmp(Config_getRot_hf_tensor_frame(cnf),"bath") == 0);
    bool rotate_qd = (Config_getQd_readmode(cnf) == 0)
                  || (strcasecmp(Config_getRot_qd_tensor_frame(cnf),"bath") == 0);

    if (hypf_bathframe == NULL){
        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): the hyperfine provenance is missing.\n");
        ccexRefusalPrint(stderr,"readHftensorfile must run, and hand its provenance over, before this step.\n");
        ccexRefuse();
    }

    ////////////////////////////////////////////////////////////////////////
    // Qubit-qubit interaction tensors
    ////////////////////////////////////////////////////////////////////////
    // Decided per entry from its recorded provenance, never from what kind of tensor it
    // is: a self/ZFS tensor and a pair tensor are both 3x3 and both can be given in
    // either frame. A qubit whose crystallographic orientation differs from the
    // reference qubit is expressed by ITS OWN ZFS tensor, transformed by the same global
    // R as everything else; off-diagonal elements surviving that is normal and correct.
    int nint_rotated = 0, nint_kept = 0;
    int nprov[5] = {0,0,0,0,0};

    for (int i=0; i<nqubit; i++){
        for (int j=i; j<nqubit; j++){

            IntmapProvenance frame = QubitArray_getIntmap_i_j_frame(qa,i,j);
            if ((int)frame >= 0 && (int)frame < 5){ nprov[(int)frame]++; }

            bool rotate_this;
            switch (frame){
                case INTMAP_DEFAULT_ZERO:         rotate_this = false; break; // zero either way
                case INTMAP_AUTO_SOURCE_GEOMETRY: rotate_this = true;  break;
                case INTMAP_EXPLICIT_BATH:        rotate_this = true;  break;
                case INTMAP_EXPLICIT_QUBIT:       rotate_this = false; break;
                case INTMAP_EXPLICIT_UNSPECIFIED:
                    // validateCoordinateFrameRotationInputs makes this unreachable for
                    // nqubit > 1. At nqubit == 1 it is the legacy "qzfs" case, which was
                    // never rotated and must keep not being rotated.
                    if (nqubit > 1){
                        ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): intmap[%d][%d] has no tensor_frame\n",i,j);
                        ccexRefusalPrint(stderr,"and nqubit = %d. This should have been caught before any file was read.\n",nqubit);
                        ccexRefuse();
                    }
                    rotate_this = false;
                    break;
                default:
                    ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): intmap[%d][%d] has an unusable\n",i,j);
                    ccexRefusalPrint(stderr,"provenance (state %d).\n",(int)frame);
                    ccexRefuse();
            }

            if (rotate_this){
                QubitArray_setIntmap_i_j(qa,rotateTensor(QubitArray_getIntmap_i_j(qa,i,j),R),i,j);
                nint_rotated++;
            }else{
                nint_kept++;
            }
        }
    }

    int ndefect_rotated = DefectArray_applyCoordinateFrameRotation(dfa,R);

    if (rank==0){
        const char* defaxis =
            (Config_getDefect_axis_reference(cnf) == DEFECT_AXIS_REFERENCE_QUBIT_AXIS) ? "qubit_axis" : NULL;
        reportCoordinateFrame(qa,cnf,R,r0,iref,qsrc,defaxis,nint_rotated,nint_kept,nprov);
        printf("      %-28s: %d transformed from bath to qubit frame\n",
               "Defect templates",ndefect_rotated);
        for (int idf=0; idf<DefectArray_getNdefect(dfa); idf++){
            if (DefectArray_getDefect_idf_axis_set(dfa,idf)){
                double* axis = DefectArray_getDefect_idf_axis(dfa,idf);
                printf("      Defect[%s].axis, calc frame : [ %12.8lf, %12.8lf, %12.8lf ]\n",
                       DefectArray_getDefect_idf_dfname(dfa,idf),axis[0],axis[1],axis[2]);
            }
        }
        printf("\n");
    }

    freeDouble2d(&qsrc,nqubit);

    int nspin  = BathArray_getNspin(ba);
    int nhf_rotated = 0, nhf_kept = 0;

    for (int i=0; i<nspin; i++){

        double* xyz = BathArray_getBath_i_xyz(ba,i);
        rotateAboutPoint(xyz,R,r0);
        BathArray_setBath_i_mindist(ba,QubitArray_mindist(xyz,qa),i);

        for (int j=0; j<nqubit; j++){

            bool computed_here;
            switch (hypf_bathframe[i][j]){
                case HYPF_SOURCE_GEOMETRY: computed_here = true;  break;
                case HYPF_FILE_TENSOR:     computed_here = false; break;
                default:
                    // Including HYPF_UNSET. Guessing either way is a silently wrong
                    // Hamiltonian -- a missed rotation or a doubled one.
                    ccexRefusalPrint(stderr,"Error(coordinate_frame_rotation): hypf[%d][%d] has no usable\n",i,j);
                    ccexRefusalPrint(stderr,"provenance (state %d). Expected %d (from the tensor file) or %d\n",
                            hypf_bathframe[i][j],HYPF_FILE_TENSOR,HYPF_SOURCE_GEOMETRY);
                    ccexRefusalPrint(stderr,"(computed from the source geometry).\n");
                    ccexRefuse();
            }

            if (hf_file_is_bathframe || computed_here){
                BathArray_setBath_i_hypf_j(ba,rotateTensor(BathArray_getBath_i_hypf_j(ba,i,j),R),i,j);
                nhf_rotated++;
            }else{
                nhf_kept++;
            }
        }

        if (rotate_qd){
            BathArray_setBath_i_quad(ba,rotateTensor(BathArray_getBath_i_quad(ba,i),R),i);
        }
    }

    if (rank==0){
        printf("      %-28s: %d transformed as R*T*R^T, %d left as written\n",
               "hyperfine tensors",nhf_rotated,nhf_kept);
        printf("      %-28s: %s\n","quadrupole tensors",
               rotate_qd ? "transformed as R*T*R^T" : "left as written (qd_tensor_frame = qubit)");
        printf("\n");

        // readHftensorfile and readQdtensorfile each print their tensors when they finish,
        // but that is BEFORE this step, so those values are the ones in the source frame.
        // Print them again here, so the log ends on the tensors that actually reach the
        // Hamiltonian rather than on the ones that were read.
        printSubTitle("Tensors after the coordinate frame rotation");
        BathArray_reportBath_hypf(ba,nqubit);
        BathArray_reportBath_quad(ba);
    }
}

void setBathStates(BathArray* ba, Config* cnf, int i){

    /**
     * @param ba : BathArray
     * @param cnf : Config
     * @param i : index of the state file ( 1 <= i <= nstate )
     * @details 
     *  If i is zero, this calculation is for the ensemble approach
     *  Therefore, the bath states shouldn't be set
    */

    // Check the number of bath states if i is zero
    if (i == 0){
        // err
        ccexRefusalPrint(stderr,"Error(setBathState): The index of the state file is zero\n");
        ccexRefusalPrint(stderr,"This calculation is for the ensemble approach\n");
        ccexRefusalPrint(stderr,"Therefore, the bath states shouldn't be set\n");
        ccexRefuse();
    }

    // For print
    char message[MAX_FILEPATH];

    // Set filename as fname_i
    char fname[MAX_FILEPATH];
    snprintf(fname,sizeof fname,"%s%d",Config_getStatefile(cnf),i);

    // Read the state file
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);
    if(data==NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                snprintf(message,sizeof message,"Warning(readStateFile): Cannot open the state file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                snprintf(message,sizeof message,"Warning(readStateFile): The number of bathfiles is larger than 1"); printMessage(message);
            }
            snprintf(message,sizeof message,"The bath state is randomly generated\n"); printMessage(message);
        }
        BathArray_setBathStatesRandom(ba);
        return;

    }else{

        // Read the number of lines
        int idx = 0;
        int count = fscanf(data, "%*s\n");
        int fline = 1; // the first line is the length of the file
        int nspin = BathArray_getNspin(ba);

        while(!feof(data)){

            //Read the bath state from the file
            //We will get the bath state at the same line of bathfile that we get the bath spins
            fline++;
            // Read the bath state
            float state = 0.0;
            count = fscanf(data, "%f\n", &state);

            if (fline == Config_get_flines_i(cnf,idx)){
                if (count == 1){
                    BathArray_setBath_i_state(ba, state, idx);
                    idx++;
                }else{
                    ccexRefusalPrint(stderr,"Error(readStateFile): StateFile format is NOT `float` \n");
                    ccexRefuse();
                }
            }
            if (idx == nspin){
                if (rank==0){
                    snprintf(message,sizeof message,"Read StateFile : %s\n",fname); printMessage(message);
                }
                break;
            }
        }
        fclose(data);
        if (idx != nspin){
            if (rank==0){
                fprintf(stderr,"Error(readStateFile): The number of bath spins is not consistent\n");
                fprintf(stderr,"The number of bath spins from bathfile : %d\n",idx);
                fprintf(stderr,"The number of bath spins from the file line : %d\n",nspin);
            }
            ccexFatal();   // every rank reads the file; rank 0 alone used to exit (review b467fe R2)
        }
    }
    return;
}

int compare_dist(const void *a, const void *b){
    BathSpin* s1 = *(BathSpin**)a;
    BathSpin* s2 = *(BathSpin**)b;
    if (s1->mindist < s2->mindist) return -1;
    if (s1->mindist > s2->mindist) return  1;
    return 0;
}

void setDefectPaxes(DefectArray* dfa, BathArray* ba, Config* cnf){

    // For print
    char message[MAX_FILEPATH];

    // Read the state file
    char* fname = Config_getAvaaxfile(cnf);
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);

    if(data==NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                snprintf(message,sizeof message,"Warning(setPaxes): Cannot open the avaax file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                snprintf(message,sizeof message,"Warning(setPaxes): The number of bathfiles is larger than 1\n"); printMessage(message);
            }
            snprintf(message,sizeof message,"The paxes is randomly generated\n"); printMessage(message);
        }
        DefectArray_setPaxesRandom(dfa, ba);
        return;
    }

    // Read the number of lines
    int nspin = Config_get_nflines(cnf);
    int idx = 0;
    int count = fscanf(data, "%*s\n");
    int fline = 1; // the first line is the length of the file

    while(!feof(data)){

        // Read the bath axes from the file
        // We will get the bath axes at the same line of bathfile that we get the bath spins
        fline++;
        // Read the bath axes
        float paxestmp = 0;
        int paxes = 0;
        count = fscanf(data, "%f\n", &paxestmp);
        paxes = int(paxestmp);
        if (fline == Config_get_flines_i(cnf,idx)){
            if (count == 1){
                //////////////////////////////////////////////
                char* dfname = BathArray_getBath_i_name(ba,idx);
                int idf = DefectArray_findDefectIndex(dfa,dfname);
                if (idf == -1){
                    // err if the defect name is not found
                    ccexRefusalPrint(stderr,"Error(setPaxes): Cannot find the defect name in the defect array\n");
                    ccexRefusalPrint(stderr,"The bath spin[%d] is : %s \n",idx,dfname);
                    ccexRefuse();
                }
                // stored as the input's navaax + 1 (index 0 is the main spin, option.cpp): the axes are 1 .. navaax-1;
                // `paxes > navaax` let navaax itself through and the tables were read one past their end (2026-10-09)
                int navaax = DefectArray_getDefect_idf_navaax(dfa,idf);
                if (paxes > navaax - 1 || paxes < 1){
                    // err if the principal axis is out of range
                    ccexRefusalPrint(stderr,"Error(setPaxes): The principal axis is out of range\n");
                    ccexRefusalPrint(stderr,"The bath spin[%d] is : %s \n",idx,dfname);
                    ccexRefusalPrint(stderr,"The number of principal axis : %d \n",navaax - 1);
                    ccexRefusalPrint(stderr,"The principal axis from the file : %d \n",paxes);
                    ccexRefuse();
                }
                //////////////////////////////////////////////
                DefectArray_setPaxes_i(dfa, idx, paxes);
                idx++;
            }else{
                ccexRefusalPrint(stderr,"Error(setPaxes): AvaaxFile format is NOT `int` \n");
                ccexRefuse();
            }
        }
        //Check if the number of bath spins is consistent
        if (Config_get_nflines(cnf) == idx){
            if (rank==0){
                snprintf(message,sizeof message,"Read AvaaxFile : %s\n",fname); printMessage(message);
            }
            break;
        }
    }
    fclose(data);

    // Check if the number of bath spins is consistent
    if (Config_get_nflines(cnf) != idx){
        if (rank==0){
            ccexRefusalPrint(stderr,"Error(setPaxes): The number of bath spins is not consistent\n");
            ccexRefusalPrint(stderr,"The number of bath spins from bathfile : %d\n",idx);
            ccexRefusalPrint(stderr,"The number of bath spins from the file line : %d\n",Config_get_nflines(cnf));
        }
        ccexRefuse();   // every rank counts the same file; rank 0 alone used to exit (A55)
    }
    return;
}

void setSubbathStates(DefectArray* dfa, BathArray* ba, Config* cnf, int i){

    /**
     * @param ba : BathArray
     * @param dfa : DefectArray
     * @param cnf : Config
     * @param i : index of the state file ( 1 <= i <= nstate )
     * @details 
     *  If i is zero, this calculation is for the ensemble approach
     *  Therefore, the bath states shouldn't be set
    */

    // Check the number of bath states if i is zero
    if (i == 0){
        // err
        ccexRefusalPrint(stderr,"Error(setBathState): The index of the state file is zero\n");
        ccexRefusalPrint(stderr,"This calculation is for the ensemble approach\n");
        ccexRefusalPrint(stderr,"Therefore, the bath states shouldn't be set\n");
        ccexRefuse();
    }

    // For print
    char message[MAX_FILEPATH];

    // File
    char fname[MAX_FILEPATH];
    snprintf(fname,sizeof fname,"%s%d",Config_getExstatefile(cnf),i);

    // Read the state file
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);
    if (data == NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                
                snprintf(message,sizeof message,"Warning(setSubbathStates): Cannot open the exstate file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                snprintf(message,sizeof message,"Warning(setSubbathStates): The number of bathfiles is larger than 1"); printMessage(message);
            }
            snprintf(message,sizeof message,"The subbath states are randomly generated\n"); printMessage(message);
        }
        DefectArray_setSubbathStatesRandom(dfa, ba);
        return;
    }else{
            
        // Read the number of lines
        int nspin = BathArray_getNspin(ba);
        int idx = 0;
        int count = fscanf(data, "%*s\n");
        int fline = 1; // the first line is the length of the file

        while(!feof(data)){

            // Read the bath state from the file
            // We will get the bath state at the same line of bathfile that we get the bath spins
            fline++;
            // Read the bath state
            const int max = 1000;
            char lline[max];
            count = fscanf(data, "%999s\n", lline);

            if (fline == Config_get_flines_i(cnf,idx)){
                int naddspin = DefectArray_getNaddspins_i(dfa,idx);

                if (naddspin != 0){

                    char* ptr = strtok(lline," ");
                    int isp=0;
                     
                    while (ptr !=NULL){
                        //state add & check
                        float ms = float(atof(ptr));
                        BathSpin* bs = DefectArray_getSubbath_i_isp(dfa,idx,isp);
                        float  S  = BathSpin_getSpin(bs);
                        BathSpin_setState(bs,ms); // Check if the state is valid in this function
                        ptr = strtok(NULL," ");
                        isp++;
                    }

                    if (isp != naddspin){
                        ccexRefusalPrint(stderr,"Error(setSubbathStates): The number of subbath spins is not consistent\n");
                        ccexRefusalPrint(stderr,"The number of subbath spins from the file line : %d\n",isp);
                        ccexRefusalPrint(stderr,"The number of subbath spins from the defect array : %d\n",naddspin);
                        ccexRefuse();
                    }
                }
                idx++;
            }
            if (idx == nspin){
                if (rank==0){
                    snprintf(message,sizeof message,"Read ExstateFile : %s\n",fname); printMessage(message);
                }
                break;
            }
        }
        fclose(data);
        if (idx != nspin){
            if (rank==0){
                fprintf(stderr,"Error(setSubbathStates): The number of bath spins is not consistent\n");
                fprintf(stderr,"The number of bath spins from bathfile : %d\n",idx);
                fprintf(stderr,"The number of bath spins from the file line : %d\n",nspin);
            }
            ccexRefuse();   // every rank reads the file; rank 0 alone used to exit (review 20261014 F2)
        }
    }
    return;
}

// A bath species the gyrofile does not name, from the built-in isotope table: appended to the species properties.
// Returns its index, or -1 when the table does not have it.
static int addIsotopeSpecies(BathArray* ba, const char* name){
    const IsotopeEntry* iso = findIsotope(name);
    if (iso == NULL){ return -1; }
    int nspecies = BathArray_getProp_nspecies(ba) + 1;
    BathArray_setProp_nspecies(ba, nspecies);
    BathArray_reallocProp(ba, nspecies-1, nspecies);
    BathArray_setProp_names_i(ba, name, nspecies-1);
    BathArray_setProp_spins_i(ba, iso->spin, nspecies-1);
    BathArray_setProp_gyros_i(ba, iso->gyro, nspecies-1);
    if (rank==0){
        printf("        %s from the built-in isotope table: spin %.1f, gyro %.8f rad/ms/G\n", name, iso->spin, iso->gyro);
    }
    return nspecies-1;
}

void readGyrofile(BathArray* ba, Config* cnf){
    
    char* fname = Config_getGyrofile(cnf);

    // No gyrofile: every bath species comes from the built-in isotope table (readBathfiles, addIsotopeSpecies)
    if (fname == NULL || fname[0] == '\0'){
        if(rank==0){printf("        No gyrofile: bath species are taken from the built-in isotope table (EasySpin, %d isotopes)\n\n",
                           isotopeTableSize());}
        return;
    }

    FILE* data;
    data = fopen(fname,"r");

    if(data==NULL){
        ccexRefusalPrint(stderr,"Error(readGyroFile): Cannot open the gyro file\n"); 
        ccexRefusalPrint(stderr,"Current WRONG gyrofile name : %s\n",fname);
        ccexRefusalPrint(stderr,"Check if the gyro file exists\n");
        ccexRefuse();
    }

    // read gyrofile set spin properties
    int nspecies = 0;
    while(!feof(data)){
        // Wider than a name may be: BathArray_setProp_names_i refuses one that
        // is too long, where reading it into its own size overran the stack.
        char    name[100] = "\0";
        float   spin = 0.0;
        double  gyro = 0.0;
        int     count = fscanf(data, "%99s %f %lf\n", name, &spin, &gyro);
        if (count == 3){
            nspecies++;
            BathArray_setProp_nspecies(ba, nspecies);
            if (nspecies == 0){ 
                BathArray_allocProp(ba); }
            else{
                BathArray_reallocProp(ba, nspecies-1, nspecies);
            }
            BathArray_setProp_names_i(ba, name, nspecies-1);
            BathArray_setProp_spins_i(ba, spin, nspecies-1);
            BathArray_setProp_gyros_i(ba, gyro, nspecies-1);
            // the gyrofile wins over the built-in table; a clear disagreement is reported
            const IsotopeEntry* iso = findIsotope(name);
            if (rank==0 && iso != NULL && (iso->spin != spin || fabs(gyro - iso->gyro) > 1e-3 * fabs(iso->gyro))){
                printf("        Note: %s in the gyrofile (spin %.1f, gyro %.6f) differs from the built-in isotope table "
                       "(spin %.1f, gyro %.6f); the gyrofile is used\n", name, spin, gyro, iso->spin, iso->gyro);
            }
        }
        else{
            ccexRefusalPrint(stderr,"Error(readGyroFile): GyroFile format is NOT name spin gyro\n");
            ccexRefuse();
        }
    }
    fclose(data);
    if(rank==0){printf("        Read GyroFile : %s\n",fname);}
    if(rank==0){BathArray_reportSpinProperties(ba);}
    if(rank==0){printf("\n");}
}

void readHftensorfile(BathArray* ba, QubitArray* qa, Config* cnf, HypfProvenance* hypf_bathframe){

    // Required bath spin properties 
    int nspecies = BathArray_getProp_nspecies(ba);
    char** names = BathArray_getProp_names(ba);

    // hf-related parameters
    int hf_readmode = Config_getHf_readmode(cnf);
    char* hf_tensorfile = Config_getHf_tensorfile(cnf);
    int hf_ignore_oor = Config_getHf_ignore_oor(cnf);

    // the defect's spin is the qubit's (the input key DefectTotSpin had to equal it until A47); CorrTotSpin, used for
    // old-format VASP files only, was removed with it and is 0 here, as it was by default
    double DefectTotSpin = QubitArray_getQubit_i_spin(qa,0);
    double CorrTotSpin = 0.0;
    double SpinFactor = 0.0;

    // Provenance of every tensor this function is about to store. Both states occur in
    // the same run: a spin outside the tensor file's range falls back to a point-dipole
    // tensor even when every other spin was matched.
    //
    // Everything starts UNSET, deliberately. allocInt2d zeroes, and 0 is HYPF_FILE_TENSOR
    // -- a real state -- so a storage path that forgot to record its provenance would
    // otherwise pass silently as "came from the file" and skip a rotation it needed.
    int nspin_prov  = BathArray_getNspin(ba);
    int nqubit_prov = QubitArray_getNqubit(qa);
    *hypf_bathframe = allocInt2d(nspin_prov,nqubit_prov);
    for (int i=0; i<nspin_prov; i++){
        for (int j=0; j<nqubit_prov; j++){ (*hypf_bathframe)[i][j] = HYPF_UNSET; }
    }

    if (rank==0){
        printSubTitle("Read the Hyperfine interaction from DFT inputfile...");
    }

    //the information in Atensor file
    if (hf_readmode==0){
        if (rank==0){
            printf("      %-18s:   %4d ( Point-dipole tensor )\n\n", "HF Readmode ", hf_readmode);
        }
        BathArray_setBathHypfs(ba,qa);

        // Every one of them was built from the source geometry.
        for (int i=0; i<nspin_prov; i++){
            for (int j=0; j<nqubit_prov; j++){ (*hypf_bathframe)[i][j] = HYPF_SOURCE_GEOMETRY; }
        }
    }
    else if (hf_readmode==1 || hf_readmode==2 || hf_readmode==3){

        if (rank==0){

            if (hf_readmode==1){
                printf("      %-18s:   %d ( Fermi-contact term + point-dipole tensor )\n", "HF Readmode ", hf_readmode);
            }
            else if (hf_readmode==2){
                printf("      %-18s:   %d ( DFT dipolar tensor )\n", "HF Readmode ", hf_readmode);
            }
            else if (hf_readmode==3){
                printf("      %-18s:   %d ( Fermi-contact term + DFT dipole tensor )\n", "HF Readmode ", hf_readmode);
            }
            printStructElementChar("HF Tensor file ", hf_tensorfile);
            printStructElementInt("HF ignoring range ", hf_ignore_oor);
            printf("\n");

        }

        //////////////////////////////////////////////////////////////////
        // Tensor information
        //////////////////////////////////////////////////////////////////

        // main information
        double*  A_Gfactor;     //g-factor___ ;
        double** A_Etc;         //etc____ ; Axx, Axy, Axz,..., Ayz, Azz
        double** AtensorArray;  //x,y,z,Axx,Axy,Axz,...,Ayz,Azz in wiDefect file

        // Boundary condition (vertex version)
        double** A_vertex;      //vertex about A-tensor
        double** A_center;      //center of plane about A-tensor
        double** A_normal;      //normal vector of plane about A-tensor
        char* vertex_condi[8] = {"v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"};

        // Boundary condition (xyz version)
        double  A_MinDif[3] = {0.0,0.0,0.0};    // MinDif[A]
        double  A_MaxDif[3] = {0.0,0.0,0.0};    // MaxDif[A]

        // Check if the BD is used with vertex or xyz version
        bool isVertex = false;

        ////////////////////////////////////////////////////////////
        // Read the tensor files
        ////////////////////////////////////////////////////////////
    
        int version = 0;
        double fileS = 0.0;
        if (READ_Tensor_provenance(hf_tensorfile, DefectTotSpin, &version, &fileS)){
            // format v3: the values are the hyperfine of the spin-S defect; DefectTotSpin and CorrTotSpin do not scale them
            SpinFactor = 1.0;
            if (rank==0){
                printf("\n\t\tTensor file format v3 (provenance line): the hyperfine of the spin-%g defect, from %s\n",
                       fileS, version == 1 ? "VASP" : "Quantum ESPRESSO");
                printf("\t\t-->SpinFactor = 1 (read as written); DefectTotSpin (%g) is only checked against S,"
                       " CorrTotSpin (%g) is not used\n\n", DefectTotSpin, CorrTotSpin);
            }
        }else{
            // an older file: the rule before format v3, unchanged
            version = READ_Tensor_ver(hf_tensorfile,&SpinFactor ,DefectTotSpin, &CorrTotSpin);

            if (rank==0){
                printHfInfo_version(version, DefectTotSpin, CorrTotSpin, SpinFactor);
                printf("\t\t!!! WARNING: this tensor file has no provenance line (a file before format v3).\n"
                       "\t\t!!! Its factor %g comes from the first line alone: VASP -> x 1 (or CorrTotSpin/DefectTotSpin),\n"
                       "\t\t!!! anything else -> x 0.5/DefectTotSpin. A SpinTensorReader v1 file (OUTCAR x 2) with a\n"
                       "\t\t!!! VASP line is read 2x. Remake the file with example/tools/SpinTensorReader (format v3).\n\n",
                       SpinFactor);
                fprintf(stderr, "Warning: hyperfine tensor file %s has no provenance line (format before v3); "
                        "read x %g by the first-line rule - remake it with SpinTensorReader\n", hf_tensorfile, SpinFactor);
            }
        }

        //the information to checking BD (in A-&Q-tensor)
        isVertex = \
        READ_BD_vertex(hf_tensorfile, &(A_vertex), &(A_center), &(A_normal), (vertex_condi));

        if (!isVertex){
            bool Check_min = READ_Tensor_array(hf_tensorfile,&(A_MinDif),"MinDif[A]",3);
            bool Check_max = READ_Tensor_array(hf_tensorfile,&(A_MaxDif),"MaxDif[A]",3);
            if ((Check_min==false) || (Check_max==false)){
                ccexRefusalPrint(stdout,"\tThe Boundary Condition (vertex and Range) can not read from tensor file (%s)\n",hf_tensorfile);
                ccexRefusalPrint(stdout,"\tPlz, Check the Tensor file !!!\n");
                ccexRefuse();
            }
        }

        if (rank==0 && verbosity){
            printHfInfo_BD(A_vertex, A_center, A_normal, A_MinDif, A_MaxDif, isVertex);
        }

        // Read gfactor and etc
        READ_Tensor_const(hf_tensorfile,names,nspecies,&(A_Gfactor),"g-factor___"); // length = nspecies+1
        READ_Tensor_etc(hf_tensorfile,names,nspecies,&(A_Etc),"etc____",1); // length = nspecies+1

        if (rank==0 && verbosity){
            printHfInfo_etc(A_Etc, A_Gfactor, names, nspecies, 0);
        }

        // Read the tensors
        READ_Tensor(hf_tensorfile,&(AtensorArray),13); // length = AtensorArray[0][0]

        if (rank==0 && verbosity){
            printHfInfo_tensor(AtensorArray, 0);
        }
        
        ////////////////////////////////////////////////////////////
        // Set Hyperfine tensor to bath array
        ////////////////////////////////////////////////////////////
        int nspin = BathArray_getNspin(ba);
        int nqubit = QubitArray_getNqubit(qa);

        if (nqubit > 1){
            ccexRefusalPrint(stderr,"Error(readHftensorfile): The number of qubit is more than 1\n");
            ccexRefusalPrint(stderr,"In this case, DONOT use Atensorfile\n");
            ccexRefusalPrint(stderr,"Instead of that, put the tensor directly.\n");
            ccexRefusalPrint(stderr,"Recommanded keyword is Bath : {}, use cce.in as json file\n");
            ccexRefuse();
        }

        for (int iqubit=0; iqubit<nqubit; iqubit++){
            for (int ispin=0; ispin<nspin; ispin++){

                // resultant paramter
                double fc = 0.0;
                MatrixXcd Adip = MatrixXcd::Zero(3,3);
                MatrixXcd Atot = MatrixXcd::Zero(3,3);

                // qubit information
                double qgyro = QubitArray_getQubit_i_gyro(qa,iqubit);
                double* qxyz = QubitArray_getQubit_i_xyz(qa,iqubit);


                // bath spin information
                double spgyro = BathArray_getBath_i_gyro(ba,ispin);
                double* spxyz = BathArray_getBath_i_xyz(ba,ispin);
                char* spname = BathArray_getBath_i_name(ba,ispin);

                // g-factor of the bath spin
                int speciesIdx = findIndexChar(names,0,nspecies-1,spname);
                double gfactor = A_Gfactor[speciesIdx+1]; // length = nspecies+1

                // relative distance between qubit and bath spin
                double rxyz[3];
                rxyz[0] = spxyz[0] - qxyz[0];
                rxyz[1] = spxyz[1] - qxyz[1];
                rxyz[2] = spxyz[2] - qxyz[2];
                
                // Check if the bath spin is in the maximum, minimum range
                bool isExist = false;
                double err = 0.001; // tolerance for the range

                if (isVertex){
                    isExist = CheckBD_vertex(rxyz, A_vertex, A_center, A_normal, err);
                }else{
                    isExist = CheckBD_Range(rxyz, A_MinDif, A_MaxDif, err);
                }

                // Check if the bath spin is matched with the relative distance of the file
                bool isFound = false;
                int tensorIdx = 0;

                if (isExist){
                    //AtesnorArray[][0~2]  : atomic position
                    isFound = FIND_AtomPosi(rxyz, AtensorArray, &tensorIdx);
                }

                // read the tensor information
                if (isExist && isFound){
                    //AtesnorArray[][3]    : fermi contact term
                    //AtesnorArray[][4~12] : Atensor value
                    fc = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][3] * gfactor * SpinFactor);
                    Adip(0,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][4] * gfactor * SpinFactor);
                    Adip(0,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][5] * gfactor * SpinFactor);
                    Adip(0,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][6] * gfactor * SpinFactor);
                    Adip(1,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][7] * gfactor * SpinFactor);
                    Adip(1,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][8] * gfactor * SpinFactor);
                    Adip(1,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][9] * gfactor * SpinFactor);
                    Adip(2,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][10] * gfactor * SpinFactor);
                    Adip(2,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][11] * gfactor * SpinFactor);
                    Adip(2,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][12] * gfactor * SpinFactor);

                    if (hf_readmode == 1){
                        // fermi contact term + Point-dipole approximation
                        // NOTE the DFT tensor read just above is DISCARDED here: Adip is
                        // overwritten by a point-dipole tensor built from the source
                        // geometry, and only the isotropic fc survives from the file. So
                        // what gets stored is a bath-frame tensor no matter what
                        // hf_tensor_frame says -- fc is a scalar and has no frame.
                        Adip = calPointDipoleTensor(qxyz, spxyz, qgyro, spgyro);
                        Atot = Adip + fc * MatrixXcd::Identity(3,3);
                        (*hypf_bathframe)[ispin][iqubit] = HYPF_SOURCE_GEOMETRY;
                    }else if (hf_readmode == 2){
                        // DFT dipolar tensor : straight from the file
                        Atot = Adip;
                        (*hypf_bathframe)[ispin][iqubit] = HYPF_FILE_TENSOR;
                    }else if (hf_readmode == 3){
                        // fermi contact term + DFT dipolar tensor : the anisotropic part
                        // is the file's, and the isotropic part has no frame
                        Atot = Adip + fc * MatrixXcd::Identity(3,3);
                        (*hypf_bathframe)[ispin][iqubit] = HYPF_FILE_TENSOR;
                    }
                    BathArray_setBath_i_hypf_j(ba, Atot, ispin, iqubit);

                }else if (isExist && !isFound){
                    if (hf_ignore_oor==0){
                        // You CANNOT ignore the atoms that deviated from the range.
                        ccexRefusalPrint(stderr,"Error(readHftensorfile): The nuclear spin is within the range, but doesn't exist in A-file!!\n");
                        ccexRefusalPrint(stderr,"A_MinDif : %lf %lf %lf\n",A_MinDif[0],A_MinDif[1],A_MinDif[2]);
                        ccexRefusalPrint(stderr,"A_MaxDif : %lf %lf %lf\n",A_MaxDif[0],A_MaxDif[1],A_MaxDif[2]);
                        ccexRefusalPrint(stderr,"difXYZ : %lf %lf %lf\n", rxyz[0],rxyz[1],rxyz[2]);
                        ccexRefusalPrint(stderr,"Nuclear spin : %lf %lf %lf\n",spxyz[0],spxyz[1],spxyz[2]);
                        ccexRefusalPrint(stderr,"In the A-file --> %d-line is the same!!", tensorIdx);
                        ccexRefuse();
                    }else{
                        // You CAN ignore the atoms that exist in the range.
                        Adip = calPointDipoleTensor(qxyz, spxyz, qgyro, spgyro);
                        BathArray_setBath_i_hypf_j(ba, Adip, ispin, iqubit);
                        (*hypf_bathframe)[ispin][iqubit] = HYPF_SOURCE_GEOMETRY;
                    }
                }else{
                    // If the atom in bath file doesn't exist in tensor file, 
                    // then, the hyperfine interaction would be simply dipolar interaction
                    Adip = calPointDipoleTensor(qxyz, spxyz, qgyro, spgyro);
                    BathArray_setBath_i_hypf_j(ba, Adip, ispin, iqubit);
                    (*hypf_bathframe)[ispin][iqubit] = HYPF_SOURCE_GEOMETRY;
                }
            }
        }

        ////////////////////////////////////////////////////////////
        // AInfo memory free
        ////////////////////////////////////////////////////////////

        int length_AtensorArray = (int)AtensorArray[0][0];
        int length_others = nspecies+1;
        
        freeDouble1d(&A_Gfactor);
        freeDouble2d(&A_Etc,length_others);
        freeDouble2d(&AtensorArray, length_AtensorArray);

        if (isVertex){
            freeDouble2d(&A_vertex, 8);
            freeDouble2d(&A_center, 8);
            freeDouble2d(&A_normal, 8);
        }

    }

    ////////////////////////////////////////////////////////////
    // Every tensor stored above must have said where it came from
    ////////////////////////////////////////////////////////////
    // An internal invariant, not an input check: Config_setHf_readmode already rejects
    // anything outside 0-3, so no input should be able to reach this. It stays because
    // the cost of being wrong is a silently mis-rotated Hamiltonian -- if a future
    // storage path forgets to record its provenance, this is what says so.
    int nunset = 0, first_i = -1, first_j = -1;
    for (int i=0; i<nspin_prov; i++){
        for (int j=0; j<nqubit_prov; j++){
            if ((*hypf_bathframe)[i][j] == HYPF_UNSET){
                if (nunset == 0){ first_i = i; first_j = j; }
                nunset++;
            }
        }
    }
    if (nunset > 0){
        ccexRefusalPrint(stderr,"Error(readHftensorfile): %d of %d hyperfine tensors were stored without\n",
                nunset,nspin_prov*nqubit_prov);
        ccexRefusalPrint(stderr,"recording which frame they are in; the first is hypf[%d][%d].\n",first_i,first_j);
        ccexRefusalPrint(stderr,"hf_readmode = %d -- the supported values are 0, 1, 2 and 3.\n",
                Config_getHf_readmode(cnf));
        ccexRefuse();
    }

    if (rank==0){
        int nqubit = QubitArray_getNqubit(qa);
        BathArray_reportBath_hypf(ba, nqubit); //report updated hyperfine tensors
    }
}

/**
 * @brief readHftensorfile for callers that do not need the provenance
 * @details The provenance only matters to applyCoordinateFrameRotation. This overload
 *          keeps the older three-argument call working by owning the array itself.
*/
void readHftensorfile(BathArray* ba, QubitArray* qa, Config* cnf){

    HypfProvenance hypf_bathframe = NULL;
    readHftensorfile(ba,qa,cnf,&hypf_bathframe);
    freeInt2d(&hypf_bathframe,BathArray_getNspin(ba));
}

void readQdtensorfile(BathArray* ba, QubitArray* qa, Config* cnf){

    // Required bath spin properties 
    int nspecies = BathArray_getProp_nspecies(ba);
    char** names = BathArray_getProp_names(ba);

    // qd-related parameters
    int qd_readmode = Config_getQd_readmode(cnf);
    char* qd_tensorfile = Config_getQd_tensorfile(cnf);
    char* qd_tensorfile_woqubit = Config_getQd_tensorfile_woqubit(cnf);
    double* qd_cellpara = Config_getQd_cellpara(cnf);

    if (rank==0){
        printSubTitle("Read the Quadrupole interaction from DFT inputfile...");
    }

    //the information in Qtensor file
    if (qd_readmode == 0){
        if (rank==0){
            printf("      %-18s:   %4d ( No Quadrupole interaction )\n\n", "Qd Readmode ", qd_readmode);
            printf("\n");
        }
    }
    else if (qd_readmode == 1){
        if (rank==0){
            ccexRefusalPrint(stdout,"      %-18s:   %d ** Removed option ( Quadrupole Experimental constant P_para and P_perp (MHz) \n", "Qd Readmode ", qd_readmode);
        }
        ccexRefuse();
    } 
    else if (qd_readmode == 2 || qd_readmode == 3 || qd_readmode == 4){ 
         if (rank==0){
            if (qd_readmode==2){
                printf("      %-18s:   %d ( Quadrupole DFT File (QuadFileWiDefect) )\n", "Qd Readmode ", qd_readmode);
                printStructElementChar("Qd Tensor file ", qd_tensorfile);
                printf("\n");
            }
            else if (qd_readmode==3){
                printf("      %-18s:   %d \n", "Qd Readmode ", qd_readmode);
                printf("\t  Cell Parameter with perpendicular Direction(X,Y,Z) (QuadCellPara)               : %lf %lf %lf\n",qd_cellpara[0],qd_cellpara[1],qd_cellpara[2]);
                printf("\t  Quadrupole DFT File for strained structure (QuadFileWiDefect)                   : %s\n",qd_tensorfile);
                printf("\t  Quadrupole DFT File for strained structure without defect (QuadFileWoDefect)    : %s\n",qd_tensorfile_woqubit);
            }
            else if (qd_readmode==4){
                printf("      %-18s:   %d \n", "Qd Readmode ", qd_readmode);
                printf("\t  Periodic Cell Parameter with perpendicular Direction(X,Y,Z) (QuadCellPara)      : {%lf, %lf, %lf}\n",qd_cellpara[0],qd_cellpara[1],qd_cellpara[2]);
                printf("\t  Quadrupole DFT File for insert Defect strcuture (QuadFileWiDefect)              : %s\n",qd_tensorfile);
                printf("\t  Quadrupole DFT File for strained structure without defect (QuadFileWoDefect)    : %s\n",qd_tensorfile_woqubit);
            }
            printf("\n");
        }      

        //////////////////////////////////////////////////////////////////
        // Tensor information
        //////////////////////////////////////////////////////////////////

        // main information
        //the inforamtion in Quadrupole file
        double*  Q_eQ = NULL;  //eQ___ ;
        double** Q_Etc = NULL; //etc____ ; Qxx,Qxy,Qxz,...,Qyz,Qzz

        double** Q_vertexwiDef = NULL;     //vertext about Q-tensor with Defect
        double** Q_centerwiDef = NULL;     //center of plane about Q-tensor with Defect
        double** Q_normalwiDef = NULL;     //normal vector of plane about Q-tensor with Defect
        double  Q_MinDifwiDef[3]={0.0,}; //MinDif[A]
        double  Q_MaxDifwiDef[3]={0.0,}; //MaxDif[A]
        double** QtensorArrayWiDef = NULL; //x,y,z,Qxx,Qxy,Qxz,...,Qyz,Qzz in wiDefect file

        double** Q_vertexwoDef = NULL;     //vertext about Q-tensor without Defect
        double** Q_centerwoDef = NULL;     //center of plane about Q-tensor without Defect
        double** Q_normalwoDef = NULL;     //normal vector of plane about Q-tensor without Defect
        double  Q_MinDifwoDef[3]={0.0,}; //MinDif[A]
        double  Q_MaxDifwoDef[3]={0.0,}; //MaxDif[A]
        double** QtensorArrayWoDef = NULL; //x,y,z,Qxx,Qxy,Qxz,...,Qyz,Qzz in woDefect file

        char* vertex_condi[8] = {"v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"};

        // Check if the BD is used with vertex or xyz version
        bool isVertex_wiDef = false;
        bool isVertex_woDef = false;

        ////////////////////////////////////////////////////////////
        // Read the tensor files
        ////////////////////////////////////////////////////////////
    
        //////////////////////////////////////////////////////
        // woDef
        //
        //the information to checking BD (in A-&Q-tensor)
        isVertex_wiDef = \
        READ_BD_vertex(qd_tensorfile, &(Q_vertexwiDef), &(Q_centerwiDef), &(Q_normalwiDef), (vertex_condi));

        if (!isVertex_wiDef){
            bool Check_min = READ_Tensor_array(qd_tensorfile,&(Q_MinDifwiDef),"MinDif[A]",3);
            bool Check_max = READ_Tensor_array(qd_tensorfile,&(Q_MaxDifwiDef),"MaxDif[A]",3);
            if ((Check_min==false) || (Check_max==false)){
                ccexRefusalPrint(stdout,"\tThe Boundary Condition (vertex and Range) can not read from tensor file (%s)\n",qd_tensorfile);
                ccexRefusalPrint(stdout,"\tPlz, Check the Tensor file !!!\n");
                ccexRefuse();
            }
        }

        if (rank==0 && verbosity){
            printHfInfo_BD(Q_vertexwiDef, Q_centerwiDef, Q_normalwiDef, Q_MinDifwiDef, Q_MaxDifwiDef, isVertex_wiDef);
        }

        // Read eQ, etc
        READ_Tensor_const(qd_tensorfile,names,nspecies,&(Q_eQ),"eQ___"); // length = nspecies+1
        READ_Tensor_etc(qd_tensorfile,names,nspecies,&(Q_Etc),"etc____",9); // length = nspecies+1
        if (rank==0 && verbosity){
            printHfInfo_etc(Q_Etc, Q_eQ, names, nspecies, 1);
        }

        // Read the tensors
        READ_Tensor(qd_tensorfile,&(QtensorArrayWiDef),12); // length = AtensorArray[0][0]
        if (rank==0 && verbosity){
            printHfInfo_tensor(QtensorArrayWiDef, 1); // mode=1 : wiDef
        }

        //////////////////////////////////////////////////////
        // woDef
        if (qd_readmode == 3 || qd_readmode == 4){
            //the information to checking BD for woqubit 
            isVertex_woDef = \
            READ_BD_vertex(qd_tensorfile_woqubit, &(Q_vertexwoDef), &(Q_centerwoDef), &(Q_normalwoDef), (vertex_condi));
    
            if (!isVertex_woDef){
                bool Check_min = READ_Tensor_array(qd_tensorfile_woqubit,&(Q_MinDifwoDef),"MinDif[A]",3);
                bool Check_max = READ_Tensor_array(qd_tensorfile_woqubit,&(Q_MaxDifwoDef),"MaxDif[A]",3);
                if ((Check_min==false) || (Check_max==false)){
                    ccexRefusalPrint(stdout,"\tThe Boundary Condition (vertex and Range) can not read from tensor file (%s)\n",qd_tensorfile_woqubit);
                    ccexRefusalPrint(stdout,"\tPlz, Check the Tensor file !!!\n");
                    ccexRefuse();
                }
            }
    
            if (rank==0 && verbosity){
                printHfInfo_BD(Q_vertexwoDef, Q_centerwoDef, Q_normalwoDef, Q_MinDifwoDef, Q_MaxDifwoDef, isVertex_woDef);
            }

            // Read the tensors
            READ_Tensor(qd_tensorfile_woqubit,&(QtensorArrayWoDef),12); // length = AtensorArray[0][0]
            if (rank==0 && verbosity){
                printHfInfo_tensor(QtensorArrayWoDef, 2); // mode=2 : woDef
            }
        }

        ////////////////////////////////////////////////////////////
        // Set EFG tensor to bath array
        ////////////////////////////////////////////////////////////
        int nspin = BathArray_getNspin(ba);
        int nqubit = QubitArray_getNqubit(qa);

        if (nqubit > 1){
            ccexRefusalPrint(stderr,"Error(readQdtensorfile): The number of qubit is more than 1\n");
            ccexRefusalPrint(stderr,"In this case, DONOT use Qtensorfile\n");
            ccexRefusalPrint(stderr,"Instead of that, put the tensor directly.\n");
            ccexRefusalPrint(stderr,"Recommanded keyword is Bath : {}, use cce.in as json file\n");
            ccexRefuse();
        }

        for (int iqubit=0; iqubit<nqubit; iqubit++){
            for (int ispin=0; ispin<nspin; ispin++){

                // resultant paramter
                MatrixXcd efg = MatrixXcd::Zero(3,3);

                // qubit information
                double qgyro = QubitArray_getQubit_i_gyro(qa,iqubit);
                double* qxyz = QubitArray_getQubit_i_xyz(qa,iqubit);


                // bath spin information
                double spgyro = BathArray_getBath_i_gyro(ba,ispin);
                float spspin = BathArray_getBath_i_spin(ba,ispin);
                double* spxyz = BathArray_getBath_i_xyz(ba,ispin);
                char* spname = BathArray_getBath_i_name(ba,ispin);

                // eQ of the bath spin
                int speciesIdx = findIndexChar(names,0,nspecies-1,spname);
                double eQ = Q_eQ[speciesIdx+1]; // length = nspecies+1
                double* efg_etc = Q_Etc[speciesIdx+1];

                // relative distance between qubit and bath spin
                double rxyz[3];
                rxyz[0] = spxyz[0] - qxyz[0];
                rxyz[1] = spxyz[1] - qxyz[1];
                rxyz[2] = spxyz[2] - qxyz[2];

                // Required parameters
                bool isExist_wiDef = false;
                double err = 0.001; // tolerance for the range
                bool isFound_wiDef = false;
                int tensorIdx_wiDef = 0;

                bool isExist_woDef = false;
                bool isFound_woDef = false;
                int tensorIdx_woDef = 0;

                //////////////////////////////////////////////////////
                // Find matched tensor
                switch (qd_readmode){
                    
                    case 2:
                        // Check if the bath spin is in the maximum, minimum range
                        if (isVertex_wiDef){
                            isExist_wiDef = CheckBD_vertex(rxyz, Q_vertexwiDef, Q_centerwiDef, Q_normalwiDef, err);
                        }else{
                            isExist_wiDef = CheckBD_Range(rxyz, Q_MinDifwiDef, Q_MaxDifwiDef, err);
                        }

                        // Check if the bath spin is matched with the relative distance of the file
                        if (isExist_wiDef){
                            //AtesnorArray[][0~2]  : atomic position
                            isFound_wiDef = FIND_AtomPosi(rxyz, QtensorArrayWiDef, &tensorIdx_wiDef);
                        }


                        // read the tensor information
                        if (isExist_wiDef){
                            if (isFound_wiDef){
                                //QtesnorArray[][3-11] : Qtensor value
                                efg(0,0) = QtensorArrayWiDef[tensorIdx_wiDef][3] ;
                                efg(0,1) = QtensorArrayWiDef[tensorIdx_wiDef][4] ;
                                efg(0,2) = QtensorArrayWiDef[tensorIdx_wiDef][5] ;
                                efg(1,0) = QtensorArrayWiDef[tensorIdx_wiDef][6] ;
                                efg(1,1) = QtensorArrayWiDef[tensorIdx_wiDef][7] ;
                                efg(1,2) = QtensorArrayWiDef[tensorIdx_wiDef][8] ;
                                efg(2,0) = QtensorArrayWiDef[tensorIdx_wiDef][9] ;
                                efg(2,1) = QtensorArrayWiDef[tensorIdx_wiDef][10] ;
                                efg(2,2) = QtensorArrayWiDef[tensorIdx_wiDef][11] ;
                            }else{
                                ccexRefusalPrint(stdout,"error in EFG(Opt:%d), the atom is within the range but doesn't exist in defect EFG file \n",qd_readmode);
                                ccexRefuse();
                            }
                        }else{
                            // Out of range
                            efg(0,0) = efg_etc[0];     
                            efg(0,1) = efg_etc[1];     
                            efg(0,2) = efg_etc[2];
                            efg(1,0) = efg_etc[3];     
                            efg(1,1) = efg_etc[4];     
                            efg(1,2) = efg_etc[5];
                            efg(2,0) = efg_etc[6];     
                            efg(2,1) = efg_etc[7];     
                            efg(2,2) = efg_etc[8];
                        }
                        break;
                    case 3:
                        // Check if the bath spin is in the maximum, minimum range
                        if (isVertex_wiDef){
                            isExist_wiDef = CheckBD_vertex(rxyz, Q_vertexwiDef, Q_centerwiDef, Q_normalwiDef, err);
                        }else{
                            isExist_wiDef = CheckBD_Range(rxyz, Q_MinDifwiDef, Q_MaxDifwiDef, err);
                        }

                        // Check if the bath spin is matched with the relative distance of the file
                        if (isExist_wiDef){
                            //AtesnorArray[][0~2]  : atomic position
                            isFound_wiDef = FIND_AtomPosi(rxyz, QtensorArrayWiDef, &tensorIdx_wiDef);
                        }

                        // read the tensor information
                        if (isExist_wiDef){
                            if (isFound_wiDef){
                                //QtesnorArray[][3-11] : Qtensor value
                                efg(0,0) = QtensorArrayWiDef[tensorIdx_wiDef][3] ;
                                efg(0,1) = QtensorArrayWiDef[tensorIdx_wiDef][4] ;
                                efg(0,2) = QtensorArrayWiDef[tensorIdx_wiDef][5] ;
                                efg(1,0) = QtensorArrayWiDef[tensorIdx_wiDef][6] ;
                                efg(1,1) = QtensorArrayWiDef[tensorIdx_wiDef][7] ;
                                efg(1,2) = QtensorArrayWiDef[tensorIdx_wiDef][8] ;
                                efg(2,0) = QtensorArrayWiDef[tensorIdx_wiDef][9] ;
                                efg(2,1) = QtensorArrayWiDef[tensorIdx_wiDef][10] ;
                                efg(2,2) = QtensorArrayWiDef[tensorIdx_wiDef][11] ;
                            }else{
                                ccexRefusalPrint(stdout,"error1 in EFG(Opt:%d), the atom is withIn the range of defect cell(curve or bubble),\n",qd_readmode);
                                ccexRefusalPrint(stdout,"                           but is not the same atoms in the defect cell(curve or bubble)!!\n");
                                ccexRefuse();
                            }
                        }else{
                            if ((qd_cellpara[0] == 0.0) && (qd_cellpara[1]==0.0) && (qd_cellpara[2]==0.0)){
                                printf("Something is wrong to Hdata structure!!\n");
                                printf("qd_cellpara : %lf %lf %lf\n",qd_cellpara[0],qd_cellpara[1],qd_cellpara[2]);
                            }
                            if(qd_cellpara[0]!=0.0){rxyz[0]=ReDefinediff(rxyz[0],Q_MinDifwoDef[0],Q_MaxDifwoDef[0],qd_cellpara[0]);}
                            if(qd_cellpara[1]!=0.0){rxyz[1]=ReDefinediff(rxyz[1],Q_MinDifwoDef[1],Q_MaxDifwoDef[1],qd_cellpara[1]);}
                            if(qd_cellpara[2]!=0.0){rxyz[2]=ReDefinediff(rxyz[2],Q_MinDifwoDef[2],Q_MaxDifwoDef[2],qd_cellpara[2]);}

                            // Check if the bath spin is in the maximum, minimum range
                            if (isVertex_wiDef){
                                isExist_wiDef = CheckBD_vertex(rxyz, Q_vertexwiDef, Q_centerwiDef, Q_normalwiDef, err);
                            }else{
                                isExist_wiDef = CheckBD_Range(rxyz, Q_MinDifwiDef, Q_MaxDifwiDef, err);
                            }
                            
                            // Check if the bath spin is matched with the relative distance of the file
                            if (isExist_wiDef){
                                //AtesnorArray[][0~2]  : atomic position
                                isFound_wiDef = FIND_AtomPosi(rxyz, QtensorArrayWiDef, &tensorIdx_wiDef);
                            }
                            
                            // read the tensor information
                            if (isExist_wiDef){
                                if (isFound_wiDef){
                                    //QtesnorArray[][3-11] : Qtensor value
                                    efg(0,0) = QtensorArrayWiDef[tensorIdx_wiDef][3] ;
                                    efg(0,1) = QtensorArrayWiDef[tensorIdx_wiDef][4] ;
                                    efg(0,2) = QtensorArrayWiDef[tensorIdx_wiDef][5] ;
                                    efg(1,0) = QtensorArrayWiDef[tensorIdx_wiDef][6] ;
                                    efg(1,1) = QtensorArrayWiDef[tensorIdx_wiDef][7] ;
                                    efg(1,2) = QtensorArrayWiDef[tensorIdx_wiDef][8] ;
                                    efg(2,0) = QtensorArrayWiDef[tensorIdx_wiDef][9] ;
                                    efg(2,1) = QtensorArrayWiDef[tensorIdx_wiDef][10] ;
                                    efg(2,2) = QtensorArrayWiDef[tensorIdx_wiDef][11] ;
                                }else{
                                    ccexRefusalPrint(stdout,"error2 in EFG(Opt:%d), the atom is withOut the range of defect cell,\n",qd_readmode);
                                    ccexRefusalPrint(stdout,"                           and is withIn the range of No defect cell(curve or bubble),\n");
                                    ccexRefusalPrint(stdout,"                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                                    ccexRefuse();

                                }
                            }else{
                                // Out of range
                                efg(0,0) = efg_etc[0];     
                                efg(0,1) = efg_etc[1];     
                                efg(0,2) = efg_etc[2];
                                efg(1,0) = efg_etc[3];     
                                efg(1,1) = efg_etc[4];     
                                efg(1,2) = efg_etc[5];
                                efg(2,0) = efg_etc[6];     
                                efg(2,1) = efg_etc[7];     
                                efg(2,2) = efg_etc[8];
                            }
                        }
                        break;
                    case 4:
                        // Check if the bath spin is in the maximum, minimum range
                        if (isVertex_wiDef){
                            isExist_wiDef = CheckBD_vertex(rxyz, Q_vertexwiDef, Q_centerwiDef, Q_normalwiDef, err);
                        }else{
                            isExist_wiDef = CheckBD_Range(rxyz, Q_MinDifwiDef, Q_MaxDifwiDef, err);
                        }

                        // Check if the bath spin is matched with the relative distance of the file
                        if (isExist_wiDef){
                            //AtesnorArray[][0~2]  : atomic position
                            isFound_wiDef = FIND_AtomPosi(rxyz, QtensorArrayWiDef, &tensorIdx_wiDef);
                        }

                        // read the tensor information
                        if (isExist_wiDef){
                            if (isFound_wiDef){
                                //QtesnorArray[][3-11] : Qtensor value
                                efg(0,0) = QtensorArrayWiDef[tensorIdx_wiDef][3] ;
                                efg(0,1) = QtensorArrayWiDef[tensorIdx_wiDef][4] ;
                                efg(0,2) = QtensorArrayWiDef[tensorIdx_wiDef][5] ;
                                efg(1,0) = QtensorArrayWiDef[tensorIdx_wiDef][6] ;
                                efg(1,1) = QtensorArrayWiDef[tensorIdx_wiDef][7] ;
                                efg(1,2) = QtensorArrayWiDef[tensorIdx_wiDef][8] ;
                                efg(2,0) = QtensorArrayWiDef[tensorIdx_wiDef][9] ;
                                efg(2,1) = QtensorArrayWiDef[tensorIdx_wiDef][10] ;
                                efg(2,2) = QtensorArrayWiDef[tensorIdx_wiDef][11] ;
                            }else{
                                if (isVertex_woDef){
                                    isExist_woDef = CheckBD_vertex(rxyz, Q_vertexwoDef, Q_centerwoDef, Q_normalwoDef, err);
                                }else{
                                    isExist_woDef = CheckBD_Range(rxyz, Q_MinDifwoDef, Q_MaxDifwoDef, err);
                                }
        
                                // Check if the bath spin is matched with the relative distance of the file
                                if (isExist_woDef){
                                    //AtesnorArray[][0~2]  : atomic position
                                    isFound_woDef = FIND_AtomPosi(rxyz, QtensorArrayWoDef, &tensorIdx_woDef);
                                }
        
                                if (isExist_woDef){
                                    if (isFound_woDef){
                                        //QtesnorArray[][3-11] : Qtensor value
                                        efg(0,0) = QtensorArrayWoDef[tensorIdx_woDef][3] ;
                                        efg(0,1) = QtensorArrayWoDef[tensorIdx_woDef][4] ;
                                        efg(0,2) = QtensorArrayWoDef[tensorIdx_woDef][5] ;
                                        efg(1,0) = QtensorArrayWoDef[tensorIdx_woDef][6] ;
                                        efg(1,1) = QtensorArrayWoDef[tensorIdx_woDef][7] ;
                                        efg(1,2) = QtensorArrayWoDef[tensorIdx_woDef][8] ;
                                        efg(2,0) = QtensorArrayWoDef[tensorIdx_woDef][9] ;
                                        efg(2,1) = QtensorArrayWoDef[tensorIdx_woDef][10] ;
                                        efg(2,2) = QtensorArrayWoDef[tensorIdx_woDef][11] ;
                                    }else{
                                        ccexRefusalPrint(stdout,"error1 in EFG(Opt:%d), the atom is withIn the range of special defect cell,\n",qd_readmode);
                                        ccexRefusalPrint(stdout,"                           and is withIn the range of No defect cell(curve or bubble),\n");
                                        ccexRefusalPrint(stdout,"                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                                        ccexRefuse();
                                    }
                                }else{
                                    if ((qd_cellpara[0] == 0.0) && (qd_cellpara[1]==0.0) && (qd_cellpara[2]==0.0)){
                                        printf("Something is wrong to Hdata structure!!\n");
                                        printf("qd_cellpara : %lf %lf %lf\n",qd_cellpara[0],qd_cellpara[1],qd_cellpara[2]);
                                    }
                                    if(qd_cellpara[0]!=0.0){rxyz[0]=ReDefinediff(rxyz[0],Q_MinDifwoDef[0],Q_MaxDifwoDef[0],qd_cellpara[0]);}
                                    if(qd_cellpara[1]!=0.0){rxyz[1]=ReDefinediff(rxyz[1],Q_MinDifwoDef[1],Q_MaxDifwoDef[1],qd_cellpara[1]);}
                                    if(qd_cellpara[2]!=0.0){rxyz[2]=ReDefinediff(rxyz[2],Q_MinDifwoDef[2],Q_MaxDifwoDef[2],qd_cellpara[2]);}
        
                                    // Check if the bath spin is in the maximum, minimum range
                                    if (isVertex_woDef){
                                        isExist_woDef = CheckBD_vertex(rxyz, Q_vertexwoDef, Q_centerwoDef, Q_normalwoDef, err);
                                    }else{
                                        isExist_woDef = CheckBD_Range(rxyz, Q_MinDifwoDef, Q_MaxDifwoDef, err);
                                    }
                                    
                                    // Check if the bath spin is matched with the relative distance of the file
                                    if (isExist_woDef){
                                        //AtesnorArray[][0~2]  : atomic position
                                        isFound_woDef = FIND_AtomPosi(rxyz, QtensorArrayWoDef, &tensorIdx_woDef);
                                    }
                                    
                                    // read the tensor information
                                    if (isExist_woDef){
                                        if (isFound_woDef){
                                            //QtesnorArray[][3-11] : Qtensor value
                                            efg(0,0) = QtensorArrayWoDef[tensorIdx_woDef][3] ;
                                            efg(0,1) = QtensorArrayWoDef[tensorIdx_woDef][4] ;
                                            efg(0,2) = QtensorArrayWoDef[tensorIdx_woDef][5] ;
                                            efg(1,0) = QtensorArrayWoDef[tensorIdx_woDef][6] ;
                                            efg(1,1) = QtensorArrayWoDef[tensorIdx_woDef][7] ;
                                            efg(1,2) = QtensorArrayWoDef[tensorIdx_woDef][8] ;
                                            efg(2,0) = QtensorArrayWoDef[tensorIdx_woDef][9] ;
                                            efg(2,1) = QtensorArrayWoDef[tensorIdx_woDef][10] ;
                                            efg(2,2) = QtensorArrayWoDef[tensorIdx_woDef][11] ;
                                        }else{
                                            ccexRefusalPrint(stdout,"error2 in EFG(Opt:%d), the atom is withIn the range of special defect cell,\n",qd_readmode);
                                            ccexRefusalPrint(stdout,"                           and is withOut the range of No defect cell(curve or bubble),\n");
                                            ccexRefusalPrint(stdout,"                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                                            ccexRefuse();
        
                                        }
                                    }else{
                                        // Out of range
                                        efg(0,0) = efg_etc[0];     
                                        efg(0,1) = efg_etc[1];     
                                        efg(0,2) = efg_etc[2];
                                        efg(1,0) = efg_etc[3];     
                                        efg(1,1) = efg_etc[4];     
                                        efg(1,2) = efg_etc[5];
                                        efg(2,0) = efg_etc[6];     
                                        efg(2,1) = efg_etc[7];     
                                        efg(2,2) = efg_etc[8];
                                    }
                                }
                            }
                        //consider the bubble cell Q-tensor along x,y,z-axis (2D & 3D)
                        }else{
                            // Check if the bath spin is in the maximum, minimum range
                            if (isVertex_woDef){
                                isExist_woDef = CheckBD_vertex(rxyz, Q_vertexwoDef, Q_centerwoDef, Q_normalwoDef, err);
                            }else{
                                isExist_woDef = CheckBD_Range(rxyz, Q_MinDifwoDef, Q_MaxDifwoDef, err);
                            }
                            
                            // Check if the bath spin is matched with the relative distance of the file
                            if (isExist_woDef){
                                //AtesnorArray[][0~2]  : atomic position
                                isFound_woDef = FIND_AtomPosi(rxyz, QtensorArrayWoDef, &tensorIdx_woDef);
                            }

                            if (isExist_woDef){
                                if (isFound_woDef){
                                    //QtesnorArray[][3-11] : Qtensor value
                                    efg(0,0) = QtensorArrayWoDef[tensorIdx_woDef][3] ;
                                    efg(0,1) = QtensorArrayWoDef[tensorIdx_woDef][4] ;
                                    efg(0,2) = QtensorArrayWoDef[tensorIdx_woDef][5] ;
                                    efg(1,0) = QtensorArrayWoDef[tensorIdx_woDef][6] ;
                                    efg(1,1) = QtensorArrayWoDef[tensorIdx_woDef][7] ;
                                    efg(1,2) = QtensorArrayWoDef[tensorIdx_woDef][8] ;
                                    efg(2,0) = QtensorArrayWoDef[tensorIdx_woDef][9] ;
                                    efg(2,1) = QtensorArrayWoDef[tensorIdx_woDef][10] ;
                                    efg(2,2) = QtensorArrayWoDef[tensorIdx_woDef][11] ;
                                }else{
                                    ccexRefusalPrint(stdout,"error3 in EFG(Opt:%d), the atom is withOut the range of special defect cell,\n",qd_readmode);
                                    ccexRefusalPrint(stdout,"                           and is withIn the range of No defect cell(curve or bubble),\n");
                                    ccexRefusalPrint(stdout,"                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                                    ccexRefuse();
                                }
                            }else{
                            
                                if ((qd_cellpara[0] == 0.0) && (qd_cellpara[1]==0.0) && (qd_cellpara[2]==0.0)){
                                    printf("Something is wrong to Hdata structure!!\n");
                                    printf("qd_cellpara : %lf %lf %lf\n",qd_cellpara[0],qd_cellpara[1],qd_cellpara[2]);
                                }
                                if(qd_cellpara[0]!=0.0){rxyz[0]=ReDefinediff(rxyz[0],Q_MinDifwoDef[0],Q_MaxDifwoDef[0],qd_cellpara[0]);}
                                if(qd_cellpara[1]!=0.0){rxyz[1]=ReDefinediff(rxyz[1],Q_MinDifwoDef[1],Q_MaxDifwoDef[1],qd_cellpara[1]);}
                                if(qd_cellpara[2]!=0.0){rxyz[2]=ReDefinediff(rxyz[2],Q_MinDifwoDef[2],Q_MaxDifwoDef[2],qd_cellpara[2]);}

                                // Check if the bath spin is in the maximum, minimum range
                                if (isVertex_woDef){
                                    isExist_woDef = CheckBD_vertex(rxyz, Q_vertexwoDef, Q_centerwoDef, Q_normalwoDef, err);
                                }else{
                                    isExist_woDef = CheckBD_Range(rxyz, Q_MinDifwoDef, Q_MaxDifwoDef, err);
                                }
                            
                            // Check if the bath spin is matched with the relative distance of the file
                                if (isExist_woDef){
                                    //AtesnorArray[][0~2]  : atomic position
                                    isFound_woDef = FIND_AtomPosi(rxyz, QtensorArrayWoDef, &tensorIdx_woDef);
                                }
                                
                                // read the tensor information
                                if (isExist_woDef){
                                    if (isFound_woDef){
                                        //QtesnorArray[][3-11] : Qtensor value
                                        efg(0,0) = QtensorArrayWoDef[tensorIdx_woDef][3] ;
                                        efg(0,1) = QtensorArrayWoDef[tensorIdx_woDef][4] ;
                                        efg(0,2) = QtensorArrayWoDef[tensorIdx_woDef][5] ;
                                        efg(1,0) = QtensorArrayWoDef[tensorIdx_woDef][6] ;
                                        efg(1,1) = QtensorArrayWoDef[tensorIdx_woDef][7] ;
                                        efg(1,2) = QtensorArrayWoDef[tensorIdx_woDef][8] ;
                                        efg(2,0) = QtensorArrayWoDef[tensorIdx_woDef][9] ;
                                        efg(2,1) = QtensorArrayWoDef[tensorIdx_woDef][10] ;
                                        efg(2,2) = QtensorArrayWoDef[tensorIdx_woDef][11] ;
                                    }else{
                                        ccexRefusalPrint(stdout,"error4 in EFG(Opt:%d), the atom is withOut the range of special defect cell,\n",qd_readmode);
                                        ccexRefusalPrint(stdout,"                           and is withOut the range of No defect cell(curve or bubble),\n");
                                        ccexRefusalPrint(stdout,"                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                                        ccexRefuse();
                                    }
                                }else{
                                    // Out of range
                                    efg(0,0) = efg_etc[0];     
                                    efg(0,1) = efg_etc[1];     
                                    efg(0,2) = efg_etc[2];
                                    efg(1,0) = efg_etc[3];     
                                    efg(1,1) = efg_etc[4];     
                                    efg(1,2) = efg_etc[5];
                                    efg(2,0) = efg_etc[6];     
                                    efg(2,1) = efg_etc[7];     
                                    efg(2,2) = efg_etc[8];
                                }
                            }
                        }

                        break;
                }
                //////////////////////////////////////////////////////
                BathSpin* bs = BathArray_getBath_i(ba,ispin);
                BathSpin_setQuad_fromEFG(bs,efg,eQ,spspin);
            }
        }

        ////////////////////////////////////////////////////////////
        // AInfo memory free
        ////////////////////////////////////////////////////////////

        int length_QtensorArrayWiDef = (int)QtensorArrayWiDef[0][0];
        int length_others = nspecies+1;

        freeDouble1d(&Q_eQ);
        freeDouble2d(&Q_Etc,length_others);
        freeDouble2d(&QtensorArrayWiDef, length_QtensorArrayWiDef);

        int length_QtensorArrayWoDef = 0; 
        if (QtensorArrayWoDef != NULL){
            length_QtensorArrayWoDef = (int)QtensorArrayWoDef[0][0];
            freeDouble2d(&QtensorArrayWoDef, length_QtensorArrayWoDef);

             if (isVertex_woDef){
                 freeDouble2d(&Q_vertexwoDef, 8);
                 freeDouble2d(&Q_centerwoDef, 8);
                 freeDouble2d(&Q_normalwoDef, 8);
             }
        }

        if (isVertex_wiDef){
            freeDouble2d(&Q_vertexwiDef, 8);
            freeDouble2d(&Q_centerwiDef, 8);
            freeDouble2d(&Q_normalwiDef, 8);
        }
    }

    if (rank==0){
        int nqubit = QubitArray_getNqubit(qa);
        BathArray_reportBath_quad(ba); //report updated hyperfine tensors
    }
}

///////////////////////////////////////////////////////////////////
//This script is related with the judging the range and atom
//
//In A-&Q-tensor, we should judge the atom is in the range or not.
//And also we match the nuclear spin and tensor data!!
///////////////////////////////////////////////////////////////////


bool READ_BD_vertex(const char* inputfile, double*** vertex, double*** center, double*** normal, char** vertex_condi){

    //char* vertex_condi = {"v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"};
    
    //Allocate the vertex
    double tmp[3] = {0.0,}; bool Vertex_OK = false;
    (*vertex) = allocDouble2d(8,3);

    //Find the vertex from inputfile
    for (int i = 0; i <8; i++){
        Vertex_OK = READ_Tensor_array(inputfile, &(tmp), vertex_condi[i],3);

        //Fail to find the vertex
        if (Vertex_OK == false){
            for (int j=0; j<8; j++){free((*vertex)[j]);}
            free(*vertex); vertex = NULL;
            return false;
        }

        //copy the data into vertex
        for (int j=0; j<3; j++){
            (*vertex)[i][j] = tmp[j];
        }
        tmp[0] = 0.0; tmp[1] = 0.0; tmp[2] = 0.0;
    }
    
    //Sort the vertex along value
    // vertex is like as,
    //  0 : (-x,-y,-z)    1 : (+x,-y,-z)
    //  2 : (-x,+y,-z)    3 : (+x,+y,-z)
    //  4 : (-x,-y,+z)    5 : (+x,-y,+z)
    //  6 : (-x,+y,+z)    7 : (+x,+y,+z)
    //
    //  QuickSort_2d(&(*vertex),2,0,7);
    //  QuickSort_2d(&(*vertex),1,0,7);
    //  QuickSort_2d(&(*vertex),0,0,7);


    //Create the Plane information
    (*center) = allocDouble2d(8,3);
    (*normal) = allocDouble2d(8,3);
    CreatePlaneInfo(*vertex, center, normal);
    
    return true;
}



////////////////////
//related Hyperfine
////////////////////

// The tensor readers below split a line at any blank. Until 2026-10-13 they split at spaces only, and the tool writes
// "{:10.5f}" fields between tabs: a value of 10 characters or more (-123.45678, 1234.56789) ran into the field before
// it and the row was dropped or read shifted ((internal notes) "텐서 도구").
static const char* TENSOR_BLANKS = " \t\r\n";

//Read the tensor(3x3) value in tensor file
void READ_Tensor(const char* inputfile, double*** Tensor,int numCol){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    //allocation Tensor
    (*Tensor) = allocDouble2d(1,numCol);
    (*Tensor)[0][0]=1;

    double test[numCol];
    memset(test, 0, numCol*sizeof(test[0]));
    //read and save the tensor value
    const int max=1000;
    char lline[max];
    while(!feof(data)){
        fgets(lline,max,data);

        char *ptr = strtok(lline, TENSOR_BLANKS);
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            double tempNum=atof(ptr);
            if (tempNum != 0){
                test[count]=tempNum;
                count++;
            }else{
                if (isStringDouble(ptr) != 0){test[count]=tempNum;count++;}
            }
            ptr = strtok(NULL, TENSOR_BLANKS);
            if (count==numCol){break;}
        }

        //if the number of value is same, write the value in Tensor array
        if (count==numCol){
            int NowRow=(int)(*Tensor)[0][0];
            reallocDouble2d(Tensor, NowRow, NowRow+1,numCol);
            for (int i=0; i<numCol;i++){
                (*Tensor)[NowRow][i]=test[i];
            }
            (*Tensor)[0][0]=NowRow+1;
            
            ////for test the Tensor
            //for (int i=0;i<numCol;i++){printf("%lf ",(*Tensor)[NowRow-1][i]);}
            //printf("\n");
        }
        //initalize the test array
        memset(test, 0, numCol*sizeof(test[0])); 

    }
    rewind(data);
    fclose(data);

    if ((*Tensor)[0][0]==1){
        printf("\tThe Tensor tensor can not read from tensor file (%s)\n",inputfile);
        printf("\tPlz, Check the Tensor file !!!\n");
    }
}

//Read the etc value in tensor file
void READ_Tensor_etc(const char* inputfile, char** names, int nspecies, double*** Tensor,char* condition, int numCol){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    //allocate the memory of Array
    (*Tensor) = allocDouble2d(nspecies+1,numCol);
    (*Tensor)[0][0]=nspecies+1;

    //read and save the tensor value
    const int max=1000;
    char lline[max];
    bool CheckIs=true;
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckCondi=false;
        bool CheckAtom=false;
        char *ptr = strtok(lline, TENSOR_BLANKS);
        int count=0;
        int col=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            //find the conditions
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckCondi=true;
            }
            if (CheckAtom == false){
                for (int i=1; i<nspecies+1;i++){
                    char *Atom = names[i-1];
                    if (strncmp(ptr, Atom, strlen(Atom))==0){
                        CheckAtom=true;
                        count=i;
                        break;
                    }
                }
            }
            // write the value in array
            if (CheckCondi && CheckAtom){
                double tempNum=atof(ptr);
                if (tempNum != 0){
                    (*Tensor)[count][col]=tempNum;col++;;
                }else{
                    if (isStringDouble(ptr) !=0){
                        (*Tensor)[count][col]=tempNum;col++;
                    }
                }
            }
            ptr = strtok(NULL, TENSOR_BLANKS);
            if (numCol==col){CheckIs=false;break;}
        }
    }
    ////print array for test
    //for (int i=1; i<Gy->Gyro_size;i++){
    //    printf("Tensor[%d] : ",i);
    //    for (int j=0; j<numCol; j++){
    //        printf(" %lf ",(*Tensor)[i][j]);
    //    }
    //    printf("\n");
    //}

    rewind(data);
    fclose(data);

    //Check the conditions is in the input file
    if (CheckIs){
        ccexRefusalPrint(stdout,"\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
        ccexRefusalPrint(stdout,"\tPlz, Check the Tensor file !!!\n");
        ccexRefuse();
    }
}

//Read the constant value in tensor file (g-factor, eQ)
void READ_Tensor_const(const char* inputfile, char** names, int nspecies, double** Array,char* condition){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    //allocate the memory of Array
    *Array=(double*)calloc((nspecies+1),sizeof(double));
    (*Array)[0] = nspecies+1;
    //read and save the tensor value
    const int max=1000;
    char lline[max];
    bool CheckIs=true;
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckCondi=false;
        bool CheckAtom=false;
        char *ptr = strtok(lline, TENSOR_BLANKS);
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            //find the conditions
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckCondi=true;
            }
            if (CheckAtom==false){
                for (int i=1; i<nspecies+1;i++){
                    char *Atom = names[i-1];
                    if (strncmp(ptr, Atom, strlen(Atom))==0){
                        CheckAtom=true;
                        count=i;
                        break;
                    }
                }
            }
            // write the value in array
            if (CheckCondi && CheckAtom){
                double tempNum=atof(ptr);
                CheckIs=false;
                if (tempNum != 0){
                    (*Array)[count]=tempNum;break;
                }else{
                    if (isStringDouble(ptr) !=0){(*Array)[count]=tempNum;break;}
                }
            }
            ptr = strtok(NULL, TENSOR_BLANKS);
        }
    }
    //print array for test
    //for (int i=1; i<Gy->Gyro_size;i++){printf("Array : %lf\n",(*Array)[i]);}

    rewind(data);
    fclose(data);

    //Check the conditions is in the input file
    if (CheckIs){
        ccexRefusalPrint(stdout,"\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
        ccexRefusalPrint(stdout,"\tPlz, Check the Tensor file !!!\n");
        ccexRefuse();
    }
}

//Read the condition's value in tensor file (MinDif, MaxDif)
bool READ_Tensor_array(const char* inputfile, double(*Array)[3],char* condition,int numCol){
    
    //check file exit
    FILE* data;
    data = fopen(inputfile, "r");
    if(data==NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    //read and save the tensor value
    const int max=1000;
    char lline[max];
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckOption=false;
        char *ptr = strtok(lline, TENSOR_BLANKS);
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckOption=true;
            }
            if (CheckOption){
                double tempNum=atof(ptr);
                if (tempNum != 0){
                    (*Array)[count]=tempNum;
                    count++;
                }else{
                    if (isStringDouble(ptr) !=0){(*Array)[count]=tempNum;count++;}
                }
            }
            ptr = strtok(NULL, TENSOR_BLANKS);
            if (count==numCol){rewind(data);fclose(data);return true;} //success to get data
        }
    }
    rewind(data);
    fclose(data);
    
    //fail to get data
    return false; 

    //printf("\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
    //printf("\tPlz, Check the Tensor file !!!\n");
    //exit(0);
}

// The provenance line of tensor file format v3 (2026-10-13, (internal notes) "텐서 파일 형식 정리").
// Before it the factor came from the first line alone (VASP -> x 1, anything else -> x 0.5/DefectTotSpin), and
// nothing in the file said which of SpinTensorReader's v1 (OUTCAR x 2) or v2 (OUTCAR) it held: the h-BN tutorial
// file was v1 data with a VASP line, read 2x. A v3 file says what its numbers are, and CCEX reads them x 1.
[[noreturn]] static void refuseTensorFile(const char* inputfile, const char* why){
    if (rank == 0){
        ccexRefusalPrint(stderr, "Error: hyperfine tensor file %s: %s\n", inputfile, why);
    }
    ccexRefuse();
}

bool READ_Tensor_provenance(const char* inputfile, double DefectTotSpin, int* version, double* fileS){
    FILE* data = fopen(inputfile, "r");
    if (data == NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    const char* prefix = "# CCEX-tensor";
    const char* v3 = "# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ";
    char first[1000] = "";
    char line[1000];
    char prov[1000] = "";
    int lineno = 0, provline = 0;
    while (fgets(line, sizeof line, data) != NULL){
        lineno++;
        if (lineno == 1){ strncpy(first, line, sizeof first - 1); }
        if (strncmp(line, prefix, strlen(prefix)) == 0){
            if (provline != 0){ fclose(data); refuseTensorFile(inputfile, "more than one '# CCEX-tensor' line"); }
            provline = lineno;
            strncpy(prov, line, sizeof prov - 1);
        }
    }
    fclose(data);
    if (provline == 0){ return false; }

    char why[1400];
    // strip the line end
    prov[strcspn(prov, "\r\n")] = '\0';
    if (provline != 2){
        snprintf(why, sizeof why, "the '# CCEX-tensor' line is line %d; it must be line 2, under the VASP or "
                 "QUANTUMESPRESSO line", provline);
        refuseTensorFile(inputfile, why);
    }
    if (strncmp(prov, v3, strlen(v3)) != 0){
        snprintf(why, sizeof why, "'%s' is not a tensor format this version reads (expected \"%s<VASP OUTCAR|QE "
                 "hyperfine.out>, S = <S>\")", prov, v3);
        refuseTensorFile(inputfile, why);
    }
    // source = <...>, S = <number>
    const char* src = prov + strlen(v3);
    const char* sep = strstr(src, ", S = ");
    if (sep == NULL){
        snprintf(why, sizeof why, "the provenance line '%s' has no \", S = <S>\"", prov);
        refuseTensorFile(inputfile, why);
    }
    std::string source(src, sep - src);
    const char* snum = sep + strlen(", S = ");
    char* end = NULL;
    double S = strtod(snum, &end);
    while (end != NULL && *end != '\0' && isspace((unsigned char)*end)){ end++; }
    if (end == snum || end == NULL || *end != '\0' || !(S > 0) || std::fabs(2*S - std::round(2*S)) > 1e-9){
        snprintf(why, sizeof why, "the provenance line '%s': S must be a positive half-integer", prov);
        refuseTensorFile(inputfile, why);
    }

    // the tag of line 1, first word
    char tag[1000] = "";
    sscanf(first, "%999s", tag);
    int fromSource;
    if (source == "VASP OUTCAR"){ fromSource = 1; }
    else if (source == "QE hyperfine.out"){ fromSource = 0; }
    else{
        snprintf(why, sizeof why, "the provenance line names source '%s'; known: 'VASP OUTCAR', 'QE hyperfine.out'",
                 source.c_str());
        refuseTensorFile(inputfile, why);
    }
    int fromTag = -1;
    if (strcmp(tag, "VASP") == 0){ fromTag = 1; }
    else if (strcmp(tag, "QUANTUMESPRESSO") == 0){ fromTag = 0; }
    if (fromTag != fromSource){
        snprintf(why, sizeof why, "line 1 is '%s' but the provenance line says source = %s (expected line 1 %s); "
                 "the file contradicts itself, remake it with SpinTensorReader", tag, source.c_str(),
                 fromSource == 1 ? "VASP" : "QUANTUMESPRESSO");
        refuseTensorFile(inputfile, why);
    }
    if (std::fabs(S - DefectTotSpin) > 1e-9){
        snprintf(why, sizeof why, "the file holds the hyperfine of a spin S = %g defect (provenance line) but "
                 "the qubit spin is %g; the file and the input describe different defects", S, DefectTotSpin);
        refuseTensorFile(inputfile, why);
    }
    *version = fromSource;
    *fileS = S;
    return true;
}

//Read the version tag of A-tensor of Q-tensor file (Q.E or VASP)
int READ_Tensor_ver(const char* inputfile,double* SpinFactor, double DefectTotSpin, double* CorrTotSpin){
    
    //check file exit
    FILE* data;
    data = fopen(inputfile, "r");
    if(data==NULL){ccexRefusalPrint(stdout,"\tThere is no interaction-file, %s!!\n",inputfile);ccexRefuse();}

    //version tag name
    char QE[] = "QUANTUMESPRESSO";  // will return 0; (with no-find)
    char VASP[] = "VASP";           // will return 1;
    
    //read and save the tensor value
    const int max=1000;     char lline[max];
    bool CheckOption=false;
    while(!feof(data)){
        fgets(lline,max,data);
        char *ptr = strtok(lline," ");
        int count=0;

        while(ptr !=NULL){
            //QE check
            if (strncmp(ptr,QE,strlen(QE))==0){
                CheckOption=true;
                if (rank==0){
                    printf("\t\t  Successfully, The version tag is read from Tensorfile (Q.E)\n");
                }
                break; // Q.E. version
            }
            
            //VASP check
            if (strncmp(ptr,VASP,strlen(VASP))==0){
                CheckOption=true;
                if (rank==0){
                    printf("\t\t  Successfully, The version tag is read from Tensorfile (VASP)\n");
                }
                //just use total spin written in A-tensor file (CorrTotSpin, which scaled it, was removed in A47)
                (*SpinFactor) = 1;
                if (rank==0){
                    printf("\t\t  You will use the raw data of A-tensor file\n");
                }
                rewind(data);fclose(data);
                return 1;   // VASP version
            }
            
            ptr = strtok(NULL, " ");
        }
    }
    rewind(data);
    fclose(data);

    //if you can search the version tag, 
    //The 'Version Tag' will be recongize as 'Q.E.' version
    if (CheckOption == false && rank == 0){
        printf("\t\t  The version tage is not finded from Tensorfile, Now we will use 'Q.E.'\n");
    }
    
    //calculate the spin factor with total spin & correlation
    (*CorrTotSpin) = 0.5; //automatically, set S=1/2 in Q.E.
    (*SpinFactor) = (*CorrTotSpin) / DefectTotSpin;
    if (DefectTotSpin != 1 && rank == 0){
        printf("\t\t  DefectTotSpin is different from default value(1)!!\n");
        printf("\t\t  To check the detail option, set '-v' options\n");
    }
    return 0;   // Q.E.

}

//print information from inputfile

//Atensor
void printHfInfo_version(int version, double DefectTotSpin, double CorrTotSpin, double SpinFactor){

    printf("\n");

    //print A-tensor version (Q.E. or VASP)
    if (version == 0){ //Q.E.
        printf("\t\tYour Hyperfine file version is 'Q.E.'\n");
        printf("\t\tCentral Tot spin & Correlation spin = %lf & %lf\n",DefectTotSpin, CorrTotSpin);
        printf("\t\t-->SpinFactor = (Correlation Spin) / (Central Tot spin) = %lf\n",SpinFactor);

    }else if (version == 1){    //VASP
        printf("\t\tYour Hyperfine file version is 'VASP'\n");
        printf("\t\tNow doesn't consider the SpinFactor(%lf) due to using raw data of A-file!!\n",SpinFactor);

    }else{
        printf("\tThe version option(%d) is weird!!\n",version);
        printf("\tNeed to check the reading version from A-tensor File");
        ccexFatal(1);   // rank 0 alone prints this (readHftensorfile), so no barrier
    }
    printf("\n");
}

void printHfInfo_BD(double** vertex, double** center, double** normal, double MinDif[3], double MaxDif[3], bool Usingvertex){
    if (Usingvertex == true){
        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),vertex[i][0], vertex[i][1], vertex[i][2]);
        }
        printf("\n");

        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),center[i][0], center[i][1], center[i][2]);
        }
        printf("\n");

        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),normal[i][0], normal[i][1], normal[i][2]);
        }
        printf("\n");


    }else{
        printf("\t[Min Range] (angstrom) : %+10lf  %+10lf  %+10lf\n",MinDif[0], MinDif[1],MinDif[2]);
        printf("\t[Max Range] (angstrom) : %+10lf  %+10lf  %+10lf\n",MaxDif[0], MaxDif[1],MaxDif[2]);
        printf("\n");
    }
}

void printHfInfo_etc(double** Etc, double* Gfactor, char** names, int nspecies, int mode){

    char string1[200];
    char string2[200];

    int etc_length=0;
    switch(mode){
        case 0:
            strcpy(string1,"\t[g-factor] (no-unit) :\n");
            strcpy(string2,"\t[A etc] (MHz) :\n\tIso  \t\tA-tensor\n");
            etc_length=1;
            break;
        case 1:
            strcpy(string1,"\t[eQ] (Q/milibarn*10^-1) :\n");
            strcpy(string2,"\t[Q etc] (MHz) :\n\tIso   \t\tEFG-tensor\n");
            etc_length=9;
            break;
    }
    
    printf("%s",string1);
    for (int i =1; i<nspecies+1; i++){
        printf("\t%s  =  %12lf\n",names[i-1],Gfactor[i]);
    }
    printf("\n");

    //print etc
    printf("%s",string2);
    for (int i =1; i<nspecies+1; i++){
        printf("\t%s  : ",names[i-1]);
        for (int j=0; j<etc_length; j++){
            printf(" %12lf ",Etc[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printHfInfo_tensor(double** tensorArray, int mode){

    char string1[200];
    
    switch(mode){
        case 0:
            printf("\t[A-tensor] (MHz) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\tFermi\t\t\t\tA-tensor\n");
            break;
        case 1:
            printf("\t[Q-tensor wi defect] (angstrom) & (Hatree/Bohr_radius^2) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\t\tEFG-tensor\n");
            break;

        case 2:
            printf("\t[Q-tensor wo defect] (angstrom) & (Hatree/Bohr_radius^2) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\t\tEFG-tensor\n");
    }   

    for (int i=1; i< int(tensorArray[0][0]); i++){
        //num of ion
        printf("\t%3d  \t",i);
        //ion position
        for (int j =0; j<3; j++){
            printf("%+12lf ",tensorArray[i][j]);
        }
        printf("    ");
        // Hyperfine files have a Fermi-contact column at index 3; quadrupole
        // files do not. Both formats then contain exactly nine tensor entries.
        int tensor_begin = 3;
        if (mode == 0){
            printf("%12lf ",tensorArray[i][3]);
            printf("    ");
            tensor_begin = 4;
        }
        for (int j=tensor_begin; j<tensor_begin+9; j++){

            printf("%12lf ",tensorArray[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

////
//For checking the atom and tensor's posi
bool CheckPosition(const double Posi[3],double* refPosi){

    bool Is_there=false;
    double err=0.05;

    if(fabs(Posi[0] - refPosi[0]) <= err){
        if(fabs(Posi[1] - refPosi[1]) <= err){
            if(fabs(Posi[2] - refPosi[2]) <= err){
                Is_there=true;   
            }
        }
    }
    return Is_there;
}
bool FIND_AtomPosi(const double refPosi[3], double** TensorValue, int* num){

    //TensorValue ==> " x, y, z, ~~~~~"
    
    int totcol=(int)TensorValue[0][0];
    for(int i=1;i<totcol;i++){
        if(CheckPosition(refPosi,TensorValue[i])){
            *num=i;
            return true;
        }
    }
    return false;
}

//Check the BD using maxRange & minRange
bool CheckBD_Range(const double Posi[3],const double minRange[3],const double maxRange[3],double err){

    bool Is_there=false;

    if(Posi[0] >= minRange[0]-err && Posi[0] <= maxRange[0]+err){
        if(Posi[1] >= minRange[1]-err && Posi[1] <= maxRange[1]+err){
            if(Posi[2] >= minRange[2]-err && Posi[2] <= maxRange[2]+err){
                Is_there=true;
            }
        }
    }
         
    return Is_there;
}

//Check the BD using 8-vertex
bool CheckBD_vertex(const double Posi[3],double** vertex, double** center, double** normal, double err){

    bool Is_there = false;

    //Set the BD to use cutting the atoms
    //
    //Max positive & negative value
    double OuterPosiBD[3] = {0,};    double OuterNegaBD[3] = {0,}; 
    
    //find the max & min value of vertex
    for (int i=0; i<8; i++){
        for (int j=0; j<3; j++){
            if (vertex[i][j] > OuterPosiBD[j]) {OuterPosiBD[j] = vertex[i][j];}
            if (vertex[i][j] < OuterNegaBD[j]) {OuterNegaBD[j] = vertex[i][j];}
        }
    }

    //1. Check that the atom is out of Boundary
    bool tmp_there = false;
    if(Posi[0] >= OuterNegaBD[0]-err && Posi[0] <= OuterPosiBD[0]+err){
        if(Posi[1] >= OuterNegaBD[1]-err && Posi[1] <= OuterPosiBD[1]+err){
            if(Posi[2] >= OuterNegaBD[2]-err && Posi[2] <= OuterPosiBD[2]+err){
                tmp_there=true;
            }
        }
    }
    if (tmp_there == false){return false;} //the atom is out of Outer BD

    //2. Check that the atom is in the unique BD (related with Hexagonal shape, ...)
    Is_there = CheckingUniqeBD(Posi, center, normal);

    //  if (Is_there == false){
    //      printf("\tCaution!, The Atom is between OuterBD and InnerBD!!\n");
    //      printf("\tIf the shape of A- & Q-tensor is orthgonal, Something is werid!!\n");
    //      printf("\tCheck plz!\n\n");
    //  }else{
    //      printf("\tSucess the checking! Posi : %lf %lf %lf\n",Posi[0],Posi[1],Posi[2]);
    //  }

    return Is_there;
}

//
bool CheckingUniqeBD(const double Posi[3], double** center, double** normal){
    //Ref : https://github.com/fgomez03/hexacheck/blob/main/hexacheck.py
    //
    //This function find the atom is in the BD or not
    //To find that, we use the normal vector and center of plane

    //It checks if all dot products are non-positive
    //return 'true' if they are and -1 otherwise.
    double c =0; double vect[3] = {0.0,};
    for (int i=0; i<6; i++){
        for (int j=0; j<3; j++){
            vect[j] = Posi[j] - center[i][j];
        }

        c = vector_dot(vect, normal[i]);
        
        if (c>0){return false;} //point is outside
    }
    return true; //the point is inside
}

//Create the normal vector and center position
void CreatePlaneInfo(double** vertex, double*** center, double*** normal){
    //Ref : https://github.com/fgomez03/hexacheck/blob/main/hexacheck.py
    //
    //The vertex of plane is like as,
    //  0 : (-x,-y,-z)    1 : (+x,-y,-z)
    //  2 : (-x,+y,-z)    3 : (+x,+y,-z)
    //  4 : (-x,-y,+z)    5 : (+x,-y,+z)
    //  6 : (-x,+y,+z)    7 : (+x,+y,+z)
    //    
    //            z
    //            |
    //       .6------7  y
    //     .' |    .'|.'
    //    4---+--5'  |    
    //    |   |  | --|----->x 
    //    |  .2--+---3
    //    |.'    | .'
    //    0------1'

    //Create the center
    //  0 : (0,1,2,3), 1 : (0,2,4,6), 2 : (0,1,4,5)
    //  3 : (1,3,5,7), 4 : (2,3,6,7), 5 : (4,5,6,7)
    CalCenter(&((*center)[0]), vertex[0], vertex[1], vertex[2], vertex[3]);
    CalCenter(&((*center)[1]), vertex[0], vertex[2], vertex[4], vertex[6]);
    CalCenter(&((*center)[2]), vertex[0], vertex[1], vertex[4], vertex[5]);
    CalCenter(&((*center)[3]), vertex[1], vertex[3], vertex[5], vertex[3]);
    CalCenter(&((*center)[4]), vertex[2], vertex[3], vertex[6], vertex[7]);
    CalCenter(&((*center)[5]), vertex[4], vertex[5], vertex[6], vertex[7]);

    //Create the normal
    // 0 : (2-0) x (1-0), 1 : (4-0) x (2-0), 2 : (1-0) x (4-0)
    // 3 : (3-1) x (5-1), 4 : (2-3) x (7-3), 5 : (5-4) x (6-4)
    double u[3] = {0.0,}; double v[3] = {0.0,};

    vector_diff(&(u),vertex[2],vertex[0]); vector_diff(&v,vertex[1],vertex[0]); 
    CalNormal(&((*normal)[0]), u, v);
    vector_diff(&u,vertex[4],vertex[0]); vector_diff(&v,vertex[2],vertex[0]); 
    CalNormal(&((*normal)[1]), u, v);
    vector_diff(&u,vertex[1],vertex[0]); vector_diff(&v,vertex[4],vertex[0]); 
    CalNormal(&((*normal)[2]), u, v);
    vector_diff(&u,vertex[3],vertex[1]); vector_diff(&v,vertex[5],vertex[1]); 
    CalNormal(&((*normal)[3]), u, v);
    vector_diff(&u,vertex[2],vertex[3]); vector_diff(&v,vertex[7],vertex[3]); 
    CalNormal(&((*normal)[4]), u, v);
    vector_diff(&u,vertex[5],vertex[4]); vector_diff(&v,vertex[6],vertex[4]); 
    CalNormal(&((*normal)[5]), u, v);

}

////
//For Q-tensor without defect, re-define the difXYZ
double ReDefinediff(double difXYZ,const  double minRange,const double maxRange,const double Copy_Length){

    if(maxRange-minRange > Copy_Length){
        ccexRefusalPrint(stdout,"Some problem in Copy_Length\n");
        ccexRefusalPrint(stdout,"maxRange-minRange : %lf\nCopy_Length : %lf\n",maxRange-minRange,Copy_Length);
        ccexRefuse();
    }

    double err=0.05;
    int check_exit=0;

    while(difXYZ > maxRange+err || difXYZ < minRange-err){
        int pre_exit=check_exit;

        if(difXYZ > maxRange+err){
            difXYZ=difXYZ-Copy_Length;
            check_exit++;
            continue;
        }
        if(difXYZ < minRange-err){
            difXYZ=difXYZ+Copy_Length;
            check_exit--;
            continue;
        }
        
        if(pre_exit == check_exit ){
            ccexRefusalPrint(stdout,"error, spin standard position is weird\n");
            ccexRefusalPrint(stdout,"Need to check the Copy_Length & curve atom position\n");
            ccexRefuse();
        }

    }
    double newdifXYZ=difXYZ;

    return newdifXYZ;
}

////
//calculate the center position using 4-vertex
void CalCenter(double** center, double p1[], double p2[], double p3[], double p4[]){
    for (int i =0; i<3; i++){
        (*center)[i] = (p1[i] + p2[i] + p3[i] + p4[i])/4;
    }
}
//calculate the normal vector using 2-vector
void CalNormal(double** normal, double u[], double v[]){
    vector_cross(&(*normal), u, v, true);
}

////
//Vector cross & dot
void vector_diff(double (*result)[3], double p1[], double p2[]){
    //making vector 
    //result = p1 - p2
    for (int i=0; i<3; i++){
        (*result)[i] = p1[i] - p2[i];
    }
}

double vector_dot(double u[], double v[]){
    double result = 0.0;
    for (int i=0; i<3; i++){
        result += u[i]*v[i];
    }
    return result;
}
void vector_cross(double** cross_P, double u[], double v[], bool norm){
    //vector cross product
    //cross_P = u x v
    (*cross_P)[0] = u[1] * v[2] - u[2] * v[1];
    (*cross_P)[1] = u[2] * v[0] - u[0] * v[2];
    (*cross_P)[2] = u[0] * v[1] - u[1] * v[0];

    //nomalization
    if (norm == true){
        double val=pow((*cross_P)[0],2) + pow((*cross_P)[1],2) + pow((*cross_P)[2],2);
        (*cross_P)[0] /= val;
        (*cross_P)[1] /= val;
        (*cross_P)[2] /= val;
    }
}
