#include "../../include/physical/pulse.h"
#include "../../include/common/memory.h"
#include <iostream>
#include <fstream>
#include <string>

#define MAX_SEQUENCE 100

/* Low-level functions ------------------------------------------------------*/

Pulse* Pulse_init(){
    Pulse* pulse = (Pulse*)allocArray1d(1, sizeof(Pulse));
    pulse->npulse = 0;
    //pulse->pulse_time     = 0.0;
    pulse->pulseiter = false;
    pulse->pulsename[0] = '\0';
    pulse->detuning_factor = 1.0;

    pulse->pulse_axes       = NULL;
    pulse->pulse_angles     = NULL;
    pulse->sequence_indices = NULL;
    pulse->sequence         = NULL;

    // DEER bath pulses (pulse.bath): none until the block asks for them
    pulse->bnpulse = 0;
    pulse->bpulse_defect[0] = '\0';
    pulse->bspin = 0.5; pulse->balphams = 0.5; pulse->bbetams = -0.5;
    pulse->n_bpulse_shift = 0;
    pulse->bpulse_energy_shifts = NULL;
    pulse->bsequence_fracs      = NULL;
    pulse->bpulse_axes          = NULL;
    pulse->bpulse_angles        = NULL;
    pulse->tot_bpulse_axes      = NULL;
    pulse->tot_bpulse_angles    = NULL;
    return pulse;
}

void Pulse_setNpulse(Pulse* pulse, int npulse){
    pulse->npulse = npulse;
}

void Pulse_setPulseTime(Pulse* pulse, double pulse_time){
    pulse->pulse_time = pulse_time;
}

void Pulse_setPulseDetuningFactor(Pulse* pulse, double detuning_factor){
    pulse->detuning_factor = detuning_factor;
}

void Pulse_setPulsename(Pulse* pulse, char* pulsename){
    copyBounded(pulse->pulsename,pulsename,sizeof pulse->pulsename,"pulsename");
}

void Pulse_setPulseiter(Pulse* pulse, bool pulseiter){
    pulse->pulseiter = pulseiter;
}

void assign_sequence_indices(Pulse* pulse) {
    for (int i = 0; i < (pulse->npulse)+1; ++i) {
        int found = 0;
        for (int j = 0; j < i; ++j) {
            if ((pulse->sequence)[i][2] == (pulse->sequence)[j][2]) {
                pulse->sequence_indices[i] = pulse->sequence_indices[j];
                found = 1;
                break;
            }
        }
        if (!found) {
            (pulse->sequence_indices)[i] = i;
        }
    }
}

int Pulse_getNpulse(Pulse* pulse){
    return pulse->npulse;
}

double Pulse_getPulseTime(Pulse* pulse){
    return pulse->pulse_time;
}

double Pulse_getPulseDetuningFactor(Pulse* pulse){
    return pulse->detuning_factor;
}

char* Pulse_getPulsename(Pulse* pulse){
    return pulse->pulsename;
}

bool Pulse_getPulseiter(Pulse* pulse){
    return pulse->pulseiter;
}

double** Pulse_getSequence(Pulse* pulse){
    return pulse->sequence;
}

char* Pulse_getPulseAxes(Pulse* pulse){
    return pulse->pulse_axes;
}

double* Pulse_getPulseAngles(Pulse* pulse){
    return pulse->pulse_angles;
}

int* Pulse_getSequenceIndices(Pulse* pulse){
    return pulse->sequence_indices;
}

// DEER bath pulses (pulse.bath), and the merged qubit+bath sequence once Pulse_resolveBathSequence has installed it
static void reportBathPulses(Pulse* pulse){
    printSubTitle((char*)"Bath pulses (pulse.bath)");
    printStructElementInt((char*)"bath pulse count", pulse->bnpulse);
    printStructElementChar((char*)"bath spin driven", pulse->bpulse_defect);
    printStructElementDouble((char*)"bath spin S", pulse->bspin);
    printStructElementDouble((char*)"bath alphams", pulse->balphams);
    printStructElementDouble((char*)"bath betams", pulse->bbetams);
    for (int k=0; k<pulse->n_bpulse_shift; k++)
        printf("%27s   pump tone [%d] : %.9g MHz\n"," ",k,RADKHZ_TO_MHZ(pulse->bpulse_energy_shifts[k]));
    if (pulse->tot_bpulse_axes == NULL) return;
    printSubTitle((char*)"Merged qubit + bath sequence (each pulse at the end of its interval)");
    printLine();
    printf("%14s   %9s   %9s   %9s   %9s\n", "Diff", "Q.Axis", "Q.Angle", "B.Axis", "B.Angle");
    for (int i=0; i<(pulse->npulse)+1; i++)
        printf("%14.10f   %9c   %9.2f   %9c   %9.2f\n", pulse->sequence[i][2], pulse->pulse_axes[i],
               pulse->pulse_angles[i], pulse->tot_bpulse_axes[i], pulse->tot_bpulse_angles[i]);
    printLine();
}

void Pulse_report(Pulse* pulse){

    int      npulse    = Pulse_getNpulse(pulse);
    char*    pulsename = Pulse_getPulsename(pulse);
    bool     pulseiter = Pulse_getPulseiter(pulse);
    double** sequence  = Pulse_getSequence(pulse);

    printTitle("Structure Pulse");

    printStructElementChar("pulsename", Pulse_getPulsename(pulse));
    printStructElementInt("npulse", Pulse_getNpulse(pulse));
    printStructElementDouble("detuning_factor", Pulse_getPulseDetuningFactor(pulse));
    printStructElementDouble("pulse_time", Pulse_getPulseTime(pulse));
    
    printStructElementBool("pulseiter", Pulse_getPulseiter(pulse));
    printf("%27s * true  - pulse is applied to each qubit\n"," ");
    printf("%27s * false - pulse is applied to qubit-array\n"," ");
    
    if (sequence != NULL){
        printSubTitle("Sequence");
        printLine();
        //printf("%6s               %8s %8s %8s %8s \n"," ","start","end","diff\n");
        //for (int i=0; i<npulse+1; i++){
        //    printf("%6s sequence[%d] : %8.3lf %8.3lf %8.3lf\n"
        //    ," ",i,sequence[i][0],sequence[i][1],sequence[i][2]);
        //}
        printf("%9s   %9s   %9s   %9s   %9s   %9s\n", "Start", "End", "Diff", "Axis", "Angle", "Index");
        for (int i=0; i<(pulse->npulse)+1; i++){
            printf("%9.3f   ", pulse->sequence[i][0]);
            printf("%9.3f   ", pulse->sequence[i][1]);
            printf("%9.3f   ", pulse->sequence[i][2]);
            printf("%9c   "  , pulse->pulse_axes[i]);
            printf("%9.2f   ", pulse->pulse_angles[i]);
            printf("%9d\n", pulse->sequence_indices[i]);
        }
        printLine();
    }
    if (pulse->bnpulse > 0) reportBathPulses(pulse);   // DEER (pulse.bath); nothing printed without bath pulses
}

//////////////////////////
//      Allocation      //
//////////////////////////
void Pulse_allocSequence(Pulse* pulse){
    pulse->sequence = allocDouble2d((pulse->npulse)+1, 3);
}

void Pulse_allocAxes(Pulse* pulse){
    pulse->pulse_axes = allocChar1d((pulse->npulse)+1);
}

void Pulse_allocAngles(Pulse* pulse){
    pulse->pulse_angles = allocDouble1d((pulse->npulse)+1);
}

void Pulse_allocSequenceIndices(Pulse* pulse){
    pulse->sequence_indices = allocInt1d((pulse->npulse)+1);
}


//////////////////////////
//         Free         //
//////////////////////////
void Pulse_freeSequence(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    freeDouble2d(&(pulse->sequence), npulse+1);
}

void Pulse_freeAxes(Pulse* pulse){
    freeChar1d(&(pulse->pulse_axes));
}

void Pulse_freeAngles(Pulse* pulse){
    freeDouble1d(&(pulse->pulse_angles));
}

void Pulse_freeSequenceIndices(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    freeInt1d(&(pulse->sequence_indices));
}

void Pulse_freeAll(Pulse* pulse){
    //Pulse_freeTauFractions(pulse);
    Pulse_freeSequence(pulse);
    Pulse_freeAxes(pulse);
    Pulse_freeAngles(pulse);
    Pulse_freeSequenceIndices(pulse);
    // DEER bath pulses (calloc'd by the pulse.bath reader and Pulse_resolveBathSequence)
    freeArray1d((void**)&pulse->bpulse_energy_shifts);
    freeArray1d((void**)&pulse->bsequence_fracs);
    freeArray1d((void**)&pulse->bpulse_axes);
    freeArray1d((void**)&pulse->bpulse_angles);
    freeArray1d((void**)&pulse->tot_bpulse_axes);
    freeArray1d((void**)&pulse->tot_bpulse_angles);
    freeArray1d((void**)&pulse);
}

void allocateDefaultSequence(Pulse* pulse) {
    if (strcasecmp(pulse->pulsename, "WAHUHA") == 0) {
        Pulse_setPulsename(pulse,"WAHUHA");
        generateSequenceWAHUHA(pulse);
    } else{
        if (pulse->npulse == 0 || strcasecmp(pulse->pulsename, "Ramsey") == 0) {
            Pulse_setPulsename(pulse,"Ramsey");
            generateSequenceRamsey(pulse);
        }
        else if (pulse->npulse == 1 || strcasecmp(pulse->pulsename, "HahnEcho") == 0) {
            Pulse_setPulsename(pulse,"HahnEcho");
            generateSequenceHE(pulse);
        }
        else if (pulse->npulse > 1 || strcasecmp(pulse->pulsename, "CPMG") == 0) {
            Pulse_setPulsename(pulse,"CPMG");
            generateSequenceCPMG(pulse);
        }
    }
}

void generateSequenceRamsey(Pulse* pulse){
    
    (pulse->sequence)[0][0]      = 0.0; 
    (pulse->sequence)[0][1]      = 1.0; 
    (pulse->sequence)[0][2]      = 1.0; 
    (pulse->pulse_axes)[0]       = 'I';
    (pulse->pulse_angles)[0]     = 0.0; 
    (pulse->sequence_indices)[0] =   0;

}

void generateSequenceHE(Pulse* pulse){
    (pulse->sequence)[0][0]    =   0.0; 
    (pulse->sequence)[0][1]    =   0.5; 
    (pulse->sequence)[0][2]    =   0.5; 
    (pulse->pulse_axes)[0]     = 'X';
    (pulse->pulse_angles)[0]     = 180.0; 
    (pulse->sequence_indices)[0] =     0;

    (pulse->sequence)[1][0]    =   0.5; 
    (pulse->sequence)[1][1]    =   1.0; 
    (pulse->sequence)[1][2]    =   0.5; 
    (pulse->pulse_axes)[1]     = 'I';
    (pulse->pulse_angles)[1]     =   0; 
    (pulse->sequence_indices)[1] =   0;
}

void generateSequenceWAHUHA(Pulse* pulse){

    (pulse->sequence)[0][0]      =       0.0; 
    (pulse->sequence)[0][1]      =   1.0/6.0; 
    (pulse->sequence)[0][2]      =   1.0/6.0; 
    (pulse->pulse_axes)[0]       = 'X';
    (pulse->pulse_angles)[0]     =      90.0; 
    (pulse->sequence_indices)[0] =         0;

    (pulse->sequence)[1][0]      =   1.0/6.0; 
    (pulse->sequence)[1][1]      =   2.0/6.0; 
    (pulse->sequence)[1][2]      =   1.0/6.0; 
    (pulse->pulse_axes)[1]       = 'Y';
    (pulse->pulse_angles)[1]     =     270.0; 
    (pulse->sequence_indices)[1] =         0;

    (pulse->sequence)[2][0]      =   2.0/6.0; 
    (pulse->sequence)[2][1]      =   4.0/6.0; 
    (pulse->sequence)[2][2]      =   2.0/6.0; 
    (pulse->pulse_axes)[2]       = 'Y';
    (pulse->pulse_angles)[2]     =      90.0; 
    (pulse->sequence_indices)[2] =         2;

    (pulse->sequence)[3][0]      =   4.0/6.0; 
    (pulse->sequence)[3][1]      =   5.0/6.0; 
    (pulse->sequence)[3][2]      =   1.0/6.0; 
    (pulse->pulse_axes)[3]       = 'X';
    (pulse->pulse_angles)[3]     =     270.0; 
    (pulse->sequence_indices)[3] =         0;

    (pulse->sequence)[4][0]      =   5.0/6.0; 
    (pulse->sequence)[4][1]      =       1.0; 
    (pulse->sequence)[4][2]      =   1.0/6.0; 
    (pulse->pulse_axes)[4]       = 'I';
    (pulse->pulse_angles)[4]     =       0.0; 
    (pulse->sequence_indices)[4] =         0;
}

void generateSequenceCPMG(Pulse* pulse){

    // n : pulse
    // c : sequence i
    // sequence (2c-1)/2n * T

    int npulse = pulse -> npulse;
    double start = 0.0;
    double end = 1.0;

    for (int i=0; i<(npulse+1); i++){

        if (i==0 && i!=npulse){
            start = 0.0;
            end = double(2*(i+1)-1)/double(2*npulse); 
        }
        else if (i!=0 && i==npulse){
            start = double(2*i-1)/double(2*npulse);
            end = 1.0;
        }
        else if (i==0 && i==npulse){
            start = 0.0;
            end = 1.0;
        }
        else{
            start = double(2*i-1)/double(2*npulse);
            end = double(2*(i+1)-1)/double(2*npulse); 
        }

        (pulse->sequence)[i][0]   = start; 
        (pulse->sequence)[i][1]   = end;
        (pulse->sequence)[i][2]   = end - start; 
        (pulse->pulse_axes)[i]    = 'X';
        (pulse->pulse_angles)[i]  =       180.0; 

        // Find 3-rd value
        if (i==0 || i==npulse){
            pulse->sequence_indices[i] = 0; 
        }
        else{
            pulse->sequence_indices[i] = 1; 
        }
    }
}


//void generateSequenceEqual(double*** sequence, int npulse){
void generateSequenceEqual(Pulse* pulse){

    // n : pulse
    // c : sequence i
    // sequence c/n * T
    int npulse = pulse->npulse;
    double start = 0.0;
    double end = 1.0;
    double diff = 1.0/double(npulse+1);

    for (int i=0; i<(npulse+1); i++){

        (pulse->sequence)[i][0] = diff * i; 
        (pulse->sequence)[i][1] = diff * (i+1);
        (pulse->sequence)[i][2] = diff;

        (pulse->sequence_indices)[i]= 0; 

    }
}

