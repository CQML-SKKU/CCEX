#ifndef __CCEX_OPTION_H_
#define __CCEX_OPTION_H_

#include "input/general.h"
#include "physical/qubit.h"
#include "physical/defect.h"
#include "clustering/cluster.h"
#include "physical/pulse.h"
#include "common/output.h"
#include "input/json.h"
#include "common/utilities.h"


void cJSON_readOptionConfig(Config* cnf, char* fccein); // general.h
void cJSON_readOptionQubitArray(QubitArray* qa, char* fccein); //qubit.h
void cJSON_readOptionCluster(Cluster* clus, char* fccein); // cluster.h
// void cJSON_readOptionBathArray(BathArray* ba, char* fccein); // bath.h
void cJSON_readOptionPulse(Pulse* pulse, char* fccein); // pulse.h
void cJSON_readOptionOutput(Output* output, char* fccein); // output.h
void cJSON_readSavemode(cJSON* item, Output* op, bool readNodiv); // savemode as a string or {"mode", "nodiv"}; readNodiv false takes the mode alone
void cJSON_readOptionDefectArray(DefectArray* dfa, char* fccein); // defect.h
void cJSON_applyMethodOption(Config* cnf, Cluster* cls, const char* method); // -m: cCCE, gCCE, meCCE or megCCE, the input's cluster_mode kept (A47)
void cJSON_setMethodOverride(const char* method); // the last -m (or NULL), read in place of the file's method before the readers judge it with cluster_mode

char* cJSON_ReadFccein(char* fccein);
void cJSON_checkInputKeys(char* fccein);   // after every cJSON_readOption*: a likely typo of a key stops the run, other unread keys are reported
char* cJSON_ReadFilePath(cJSON* root, char* key, bool _default, char* default_value);
char* cJSON_ReadString(cJSON* root, char* key, bool _default, char* default_value);

int cJSON_ReadInt(cJSON* root, char* key, bool _default, int default_value);
float cJSON_ReadFloat(cJSON* root, char* key, bool _default, float default_value);
double cJSON_ReadDouble(cJSON* root, char* key, bool _default, double default_value);
bool cJSON_ReadBool(cJSON* root, char* key, bool _default, bool default_value);

char** cJSON_ReadString1d(cJSON* root, char* key, bool _default, char** default_value, int isze); // memory free is needed
char** cJSON_ReadFilePath1d(int* length, cJSON* root, char* key, bool _default, char** default_value); // memory free is needed
int* cJSON_ReadInt1d(cJSON* root, char* key, bool _default, int* default_value, int size); // memory free is needed
float* cJSON_ReadFloat1d(cJSON* root, char* key, bool _default, float* default_value, int size); // memory free is needed
double* cJSON_ReadDouble1d(cJSON* root, char* key, bool _default, double* default_value, int size); // memory free is needed
double** cJSON_ReadDouble2d(cJSON* root, char* key, bool _default, double** default_value, int row, int col); // memory free is needed

MatrixXcd cJSON_ReadTensor(cJSON* root, char* key, bool _default, MatrixXcd default_value);

void cJSON_ReadDefectInfo_IntCharDoubleArray(cJSON* root, char* key, int valuecount, double**** array, char** types, int navaax, int naddspin);
void cJSON_ReadDefectInfo_IntCharMatrixXcd2d(cJSON* root, char* key, int valuecount, MatrixXcd*** array, char** types, int navaax, int naddspin);
void cJSON_ReadDefectInfo_IntCharMatrixXcd1d(cJSON* root, char* key, int valuecount, MatrixXcd** array, int navaax);
void cJSON_ReadDefectInfo_IntCharDouble(cJSON* root, char* key, double** array, int navaax);

//void cJSON_PrintNotFount(

// void checkComparable(Config* cnf, QubitArray* qa, BathArray* ba, Cluster* clus, Pulse* pulse, Output* output);
// DefectTotalSpin vs qa->qubit[i]->spin

#endif // __CCEX_OPTION_H_
