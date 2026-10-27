/**
 * @file pulse.h
 * @brief The Pulse structure: the pulse sequence on the qubit and the DEER bath pulses.
 */
#ifndef __CCEX_PULSE_H_
#define __CCEX_PULSE_H_

#include "../../include/common/utilities.h"

/**
 * @struct Pulse
 * @brief The pulse sequence of the run: the pulses on the qubit (and on bath spins for DEER) and their timing.
 */
// The default half width of the window around a DEER pump tone (pulse.bath.bandwidth, A67): the fixed 10 rad*kHz of before
// (adding_MECCE src/simulator_gcce.cpp:546-553, an arbitrary value; = 10/(2 pi 1000) MHz, about 0.0015915 MHz). The physical
// value is to be decided with the DEER author.
constexpr double kBathPulseBandwidthRadkHz = 10.0;

typedef struct {

    int npulse; /**< The number of pulse (default : -)*/
    bool pulseiter; /**< The pulse iteration (default : false)*/
    /**<
     * About npulse (Details)
     * npulse = 0 : Ramsey (default)
     * npulse = 1 : HahnEcho (default)
     * npulse > 1 : CPMG (default)
     * npulse = 4 & pulsename = "WAHUHA" : WAHUHA (optional)
    */
    char pulsename[100]; /**< The pulsename (default : "None")*/
    /**<
     * About pulsename (Details)
     * pulsename = HahnEcho | Ramsey | CPMG | WAHUHA | Manual
     * above pulsename option without Maunal  doesn't need a sequence 
     * "Manual" option must need a sequence
     * "Manual" with rotations other than X or Y 180 degrees needs the "gCCE" method;
     * cCCE (and pCCE) calculate every pulse as an ideal pi flip and refuse anything else
    */

    /**
     *  @brief : Related with pulse sequence.
     */
    double** sequence; /**< pulse sequence (default : NULL)*/
    double*  pulse_angles; 
    /**<
     * puulse angles: Array of rotation angle of pulse.
     * Unit: degree, in the input file and in this array; converted to radian where applied.
     *   len(sequence) = len(pulse_angles)
     */
    char*    pulse_axes;
    /**<
     * pulse axes: Array of rotation axis of pulse.
     * you can choose the axis of rotation angle. 
     * Input file parameters: "X", "Y", "Z", "I"(Identity) 
     *   len(sequence) = len(pulse_axes)
     */
    int*     sequence_indices;
    /**<
     *  sequence_indices: The Index that have the same difference value
     *                        if there is no the same difference value 
     *                        then give the current index
     *                        ( This would reduce the calculational time cost ) 
     *   len(sequence) = len(sequence_indices)
     */
    double pulse_time;
    double detuning_factor;

    //
    /**<
     * About sequence (Details)
     * 
     * If you want to adjust the pulse sequence, that you have to set the pulse name : "Manual", 
     * and you have to give the "sequence" variable.
     *
     * -------------------------------------------------------------------
     * ex)          !                                            !(t_free)
     * Pulse Index  0    1    2    3    4    5    6    7    8    |
     * Pulse (8#)   |____|____|____|____|____|____|____|____|____|
     *              |tau1|tau2|           .......           |tau9|
     *              |seq1|seq2|           .......           |seq9|
     * it's okay for tau1 != tau2.
     * -------------------------------------------------------------------
     *
     * Sequence format:
     *  "pulsename" : "Manual",
     *  "sequence"  : [["frac_1",     "axis_1", float(angle_1) (deg)],
     *                 ["frac_2",     "axis_2", float(angle_2) (deg)],
     *                 ["frac_3",     "axis_3", float(angle_3) (deg)],
     *                                     ...
     *                 ["frac_n",          "I",                  0.0]]
     *
     *     where Summation_n(frac_n) must be float(1.0).
     *  And the rotation axis and angle of last sequence(sequence[-1]) is 0
     * 
     * ex1) Ramsey sequence 
     *    inputfile tag:
     *       - "npulse": 0
     *           or
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1.0", "I", 0.0]]
     *
     *
     * ex2) Hahn-Echo sequence 
     *    inputfile tag:
     *       - "npulse" : 1
     *           or
     *       - "pulsename" : "Manual",
     *       - "sequence"  : [["0.5", "X", 180.0],
     *                        ["0.5", "I",   0.0]]
     *
     * ex3) CPMG sequence 
     *    inputfile tag:
     *       - "npulse" > 1
     *           or
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1/{1tau}", "X", 180.0],
     *                       ["1/{2tau}", "X", 180.0],
     *                       ["1/{2tau}", "X", 180.0],
     *                                  ...
     *                       ["1/{2tau}", "X", 180.0],
     *                       ["1/{1tau}", "I",   0.0]]
     *
     * ex3) WAHUHA sequence 
     *    inputfile tag:
     *       - "npulse"   : 4
     *       - "pulsename": "WAHUHA",
     *           or
     *       - "npulse"   : 4
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1/6", "X",  90.0],
     *                       ["1/6", "Y", 270.0],
     *                       ["2/6", "Y",  90.0],
     *                       ["1/6", "X", 270.0],
     *                       ["1/6", "I",   0.0]]
     *
     * | ------------------------------------------------------------- |
     * | ----------------- !        WARNING       ! ------------------ | 
     * | ------------------------------------------------------------- |
     * | - In CCEX CODE same variable name(sequence), but differnt!  - |
     *
     * In code, Array format:
     *  sequence is started with 0.0 and end with 1.0
     *  sequence[npulse+1][3] (defulat : NULL)
     *  sequence[ipulse][0] = Fraction of Previous pulse
     *  sequence[ipulse][1] = Fraction of current pulse
     *  sequence[ipulse][2] = Difference of Fractions "[i][1] - [i][0]"
     *
     * e.g. sequence parameter : Ramsey
     *      sequence parameter :
     *          sequence[0]      = {0.0 , 1.0 , 1.0}
     *          pulse_angles     = {0.0}
     *          pulse_axes       = {"I"}
     *          sequence_indices = {0}
     *
     * e.g. sequence parameter : Hahn-echo
     *          sequence[0]      = {"0.0", "0.5" , 0.5 }
     *          sequence[1]      = {"0.5", "1.0" , 0.5 }
     *          pulse_angles     = {180.0}
     *          pulse_axes       = {"X"}
     *          sequence_indices = {0, 0}
     *
     * e.g. sequence parameter : CPMG
     *         sequence[0]     = 0.0 , 0.25 , 0.25 
     *         sequence[1]     = 0.25 , 0.75 , 0.5 
     *         sequence[2]     = 0.75 , 1.0 , 0.25 
     *         pulse_angles    = {180.0,  180.0}
     *         pulse_axes      = {"X", "X"}
     *         sequence_indices = {0, 1, 0}
     *
     * e.g. sequence parameter : WAHUHA (generateSequenceWAHUHA)
     *         sequence[0]     = 0   , 1/6 , 1/6
     *         sequence[1]     = 1/6 , 2/6 , 1/6
     *         sequence[2]     = 2/6 , 4/6 , 2/6
     *         sequence[3]     = 4/6 , 5/6 , 1/6
     *         sequence[4]     = 5/6 , 1   , 1/6
     *         pulse_angles    = {90.0, 270.0, 90.0, 270.0, 0.0}
     *         pulse_axes      = {"X", "Y", "Y", "X", "I"}
     *         sequence_indices = {0, 0, 2, 0, 0}
    */

    // DEER bath pulses (bnpulse = 0: none), from the input block pulse.bath. As given in the input:
    int     bnpulse;             // number of bath pulses
    char    bpulse_defect[50];   // bath spin name driven (case-insensitive)
    double  bspin, balphams, bbetams;
    int     n_bpulse_shift;      // number of pump tones
    double* bpulse_energy_shifts;// pump tones [rad*kHz] (input MHz)
    double  bpulse_bandwidth;    // half width of the window around a tone [rad*kHz] (input MHz, pulse.bath.bandwidth; A67;
                                 // default kBathPulseBandwidthRadkHz)
    double* bsequence_fracs;     // bnpulse+1 interval fractions
    char*   bpulse_axes;         // bnpulse+1
    double* bpulse_angles;       // bnpulse+1 [deg]
    // After Pulse_resolveBathSequence: aligned with the merged sequence (npulse+1 rows)
    char*   tot_bpulse_axes;
    double* tot_bpulse_angles;   // [deg]
    // Continuous drive of bath spins (pulse.bath_drive, calculation/bath_drive.h); NULL: none
    struct BathDrive* bdrive;
} Pulse;

Pulse* Pulse_init();

void Pulse_report(Pulse* pulse);

void Pulse_setNpulse(Pulse* pulse, int npulse);
void Pulse_setPulseTime(Pulse* pulse, double pulse_time);
void Pulse_setPulseDetuningFactor(Pulse* pulse, double detuning_factor);
void Pulse_setPulsename(Pulse* pulse, char* pulsename);
void assign_sequence_indices(Pulse* pulse);
void Pulse_setPulseiter(Pulse* pulse, bool pulseiter);

int Pulse_getNpulse(Pulse* pulse);
double Pulse_getPulseTime(Pulse* pulse);
double Pulse_getPulseDetuningFactor(Pulse* pulse);
char* Pulse_getPulsename(Pulse* pulse);
char* Pulse_getPulseAxes(Pulse* pulse);
double* Pulse_getPulseAngles(Pulse* pulse);
int* Pulse_getSequenceIndices(Pulse* pulse);
double** Pulse_getSequence(Pulse* pulse);
bool Pulse_getPulseiter(Pulse* pulse);

// allocation
void Pulse_allocSequence(Pulse* pulse);
void Pulse_freeSequence(Pulse* pulse);
void Pulse_freeAll(Pulse* pulse);

void Pulse_allocAxes(Pulse* pulse);
void Pulse_allocAngles(Pulse* pulse);
void Pulse_allocSequenceIndices(Pulse* pulse);
void Pulse_freeAxes(Pulse* pulse);
void Pulse_freeAngles(Pulse* pulse);
void Pulse_freeSequenceIndices(Pulse* pulse);

// generate pulse sequence
void generateSequenceRamsey(Pulse* pulse);
void generateSequenceHE(Pulse* sequence);
void generateSequenceWAHUHA(Pulse* sequence);
void generateSequenceCPMG(Pulse* pulse);
void allocateDefaultSequence(Pulse* pulse);


#endif // __CCEX_PULSE_H_
