#include "../../include/input/mpi_seed.h"
#include "../../include/input/option.h"
#include "../../include/input/json.h"
#include "../../include/common/memory.h"
#include "../../include/calculation/superoperator.h"
#include "../../include/auto_convergence/rdip_ladder.h"

#include <cmath>
#include <limits>
#include <fstream>
#include <unistd.h>
#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <string>
#include <vector>

// Every key the readers below look up, with the block it was looked up in, so that
// cJSON_checkInputKeys can tell a key nobody reads from one that was read (D24, 2026-10-02). A block
// is named by its path from the top, lower-cased (cJSON compares names without case), array elements
// as "[]": "qubit.qubit[]" for each entry of Qubit.qubit. Every parse here records the path of every
// object it made, so a lookup knows its block whichever of the six readers made the tree.
static std::set<std::string>& lookedUpKeys(){ static std::set<std::string> keys; return keys; }
static std::map<const cJSON*, std::string>& objectPaths(){ static std::map<const cJSON*, std::string> paths; return paths; }
// the same paths as written (case kept), for the messages: "Error: pulse.count is not found" (A47, 2026-10-13)
static std::map<const cJSON*, std::string>& shownPaths(){ static std::map<const cJSON*, std::string> paths; return paths; }
static std::string lowerKey(const char* key){
    std::string out(key ? key : "");
    for (char& c : out) c = (char)std::tolower((unsigned char)c);
    return out;
}
static std::string joinKey(const std::string& path, const std::string& key){ return path.empty() ? key : path + "." + key; }
static void recordPaths(const cJSON* node, const std::string& path, const std::string& shown){
    objectPaths()[node] = path;
    shownPaths()[node] = shown;
    for (const cJSON* c = node->child; c; c = c->next){
        const bool member = cJSON_IsObject(node) && c->string;
        const std::string childPath = member ? joinKey(path, lowerKey(c->string)) : path + "[]";
        const std::string childShown = member ? joinKey(shown, c->string) : shown + "[]";
        if (cJSON_IsObject(c) || cJSON_IsArray(c)) recordPaths(c, childPath, childShown);
    }
}
// key inside object, as the input writes the path to it ("Qubit.file"; a top-level key alone)
static std::string shownKey(const cJSON* object, const char* key){
    const auto where = shownPaths().find(object);
    return where == shownPaths().end() ? std::string(key ? key : "") : joinKey(where->second, key ? key : "");
}
static cJSON* trackedParse(const char* text){
    cJSON* root = cJSON_Parse(text);
    if (root) recordPaths(root, "", "");
    return root;
}
static cJSON* trackedGetObjectItem(const cJSON* object, const char* key){
    if (key){
        const auto where = objectPaths().find(object);
        lookedUpKeys().insert(joinKey(where == objectPaths().end() ? "?" : where->second, lowerKey(key)));
    }
    return cJSON_GetObjectItem(object, key);
}
#define cJSON_Parse trackedParse
#define cJSON_GetObjectItem trackedGetObjectItem

/**
 * @brief Read a 3-vector axis, which must be exactly three JSON numbers
 * @details A string, null or boolean element (valuedouble 0) is refused rather than read,
 *          so "qubit_axis": [1,"1",1] cannot quietly become [1,0,1].
 * @param section      the enclosing JSON object
 * @param section_name name of the enclosing section, for the error message
 * @param key          the key looked up in section
 * @param[in,out] out  left at its default value if the key is absent and not required
 * @param required     whether a missing key is an error
*/
static void readAxis3(cJSON* section, const char* section_name, const char* key,
                     double out[3], bool required){

    cJSON* item = cJSON_GetObjectItem(section,key);
    if (item == NULL){
        if (required){
            ccexRefusalPrint(stderr,"Error: %s.%s is required\n",section_name,key);
            ccexRefuse();
        }
        return; // absent and optional : keep the default
    }

    if (!cJSON_IsArray(item) || cJSON_GetArraySize(item) != 3){
        ccexRefusalPrint(stderr,"Error: %s.%s must be an array of exactly 3 numbers\n",section_name,key);
        ccexRefuse();
    }

    int i = 0;
    cJSON* element;
    cJSON_ArrayForEach(element,item){
        if (!cJSON_IsNumber(element)){
            ccexRefusalPrint(stderr,"Error: %s.%s[%d] is not a number\n",section_name,key,i);
            ccexRefusalPrint(stderr,"Write the axis as plain numbers, e.g. [0, 0, 1] -- not strings, null or booleans.\n");
            ccexRefuse();
        }
        out[i] = element->valuedouble;
        i++;
    }
}

// A JSON boolean and nothing else. cJSON_ReadBool falls back to its default for any
// other type, which would turn "enabled": "true" into a silent no-op.
static bool readStrictBool(cJSON* section, const char* section_name, const char* key, bool default_value){

    cJSON* item = cJSON_GetObjectItem(section,key);
    if (item == NULL){ return default_value; }
    if (!cJSON_IsBool(item)){
        ccexRefusalPrint(stderr,"Error: %s.%s must be a boolean (true or false)\n",section_name,key);
        ccexRefuse();
    }
    return (bool)item->valueint;
}

static void refuseWrongType(cJSON* item, const char* key, bool right, const char* what);
static double finiteNumber(cJSON* item, const char* key, int i, int j, bool asFloat);

////////////////////////////////////////////////////////////////////////////////////////////////
// A47 (2026-10-13, (internal notes)): the input keys were renamed and grouped into blocks. An old name
// is refused, not read, and the message names the new one (as savemode's old values are): a run that took a default for a
// key it no longer knows would look like a finished one. The converter rewrites an old input.
static const char* const kConvertHint =
    "Old inputs: python3 example/tools/InputKeys/convert_input_keys.py FILE... rewrites them with the new keys";

// A lookup the unread-key check does not record (cJSON_checkInputKeys): the refusals look for names nobody reads, and a
// recorded old name would be offered as "did you mean" for a typo.
static cJSON* untrackedItem(const cJSON* object, const char* key){
    if (!object || !cJSON_IsObject(object)) return NULL;
    for (cJSON* c = object->child; c; c = c->next)
        if (c->string && !strcasecmp(c->string,key)) return c;
    return NULL;
}

static void refuseOldInputKeys(cJSON* root){
    struct Renamed { const char* block; const char* key; const char* now; };
    static const Renamed renamed[] = {
        {"", "qubitfile", "is now \"Qubit.file\" (\"Qubit\": {\"file\": ...})"},
        {"", "qspin", "is now \"Qubit.spin\" (with Qubit.file)"},
        {"", "qzfs", "is now \"Qubit.zfs\" (with Qubit.file)"},
        {"", "qalphams", "is now \"Qubit.alphams\" (with Qubit.file)"},
        {"", "qbetams", "is now \"Qubit.betams\" (with Qubit.file)"},
        {"", "rbathcut", "is now the inner radius of \"rbath\": [inner, outer]"},
        {"", "rdipcut", "is now the inner radius of \"rdip\": [inner, outer]"},
        {"", "bathadjust", "is now \"bath_shift\""},
        {"", "statefile", "is now \"bath_state_file\""},
        {"", "exstatefile", "is now \"Defect.state_file\""},
        {"", "avaaxfile", "is now \"Defect.axis_file\""},
        {"", "addsubclus", "is now \"add_subclusters\" (true or false)"},
        {"", "nk", "is now \"max_clusters_per_order\""},
        {"", "vector_mean_field", "is now \"hamiltonian.common.bath_bath_mean_field\" (\"xyz\" for true, \"z\" for false)"},
        {"", "hf_readmode", "is now \"hamiltonian.common.hyperfine\" {\"dipole\": \"point\"|\"dft\", \"contact\": true|false} "
                            "(0 point/false, 1 point/true, 2 dft/false, 3 dft/true)"},
        {"", "hf_tensorfile", "is now \"hamiltonian.common.hyperfine.file\""},
        {"", "hf_ignore_oor", "is now \"hamiltonian.common.hyperfine.missing_row\" (0 \"stop\", 1 \"point_dipole\")"},
        {"", "qd_readmode", "is now \"hamiltonian.common.quadrupole\" (0 \"off\"; 2, 3, 4 {\"structure\": \"plain\", "
                            "\"strained\", \"bubble\"})"},
        {"", "qd_tensorfile", "is now \"hamiltonian.common.quadrupole.file\""},
        {"", "qd_tensorfile_woqubit", "is now \"hamiltonian.common.quadrupole.file_without_defect\""},
        {"", "qd_cellpara", "is now \"hamiltonian.common.quadrupole.cell\""},
        {"", "knight", "is now \"hamiltonian.cCCE.knight\""},
        {"", "hfmedi", "is now \"hamiltonian.cCCE.hfmedi\""},
        {"", "DefectTotSpin", "was removed: the qubit spin is used (also for the S of a v3 tensor file)"},
        {"", "CorrTotSpin", "was removed: it served old-format VASP tensor files only; remake such a file with "
                            "example/tools/SpinTensorReader (format v3)"},
        {"", "sK", "is now \"pcce.partition_size\""},
        {"", "max_trial", "is now \"pcce.trials\""},
        {"", "max_iter", "is now \"pcce.max_iterations\""},
        {"", "kmeans_pp", "is now \"pcce.kmeans_pp\""},
        {"", "iter_detail", "is now \"pcce.print_trials\""},
        {"", "pcce_detail_limit", "is now \"pcce.print_limit\""},
        {"", "pcce_partition", "is now \"pcce.partition\""},
        {"", "pcce_connections", "is now \"pcce.connections\""},
        {"", "pcce_cache_mode", "is now \"pcce.cache.mode\""},
        {"", "pcce_cache_file", "is now \"pcce.cache.file\""},
        {"", "bundle_max_dim", "is now \"bundle.max_dim\""},
        {"", "bundle_merge_cutoff", "is now \"bundle.merge_cutoff\""},
        {"", "bundle_connections", "is now \"bundle.connections\""},
        {"", "bundle_connection_budget", "is now \"bundle.connection_budget\""},
        {"", "expansion_max_clusters", "is now \"bundle.max_clusters\""},
        {"", "expansion_max_work", "is now \"bundle.max_work\""},
        {"", "adaptive_work_file", "is now \"adaptive.work_file\""},
        {"", "npulse", "is now \"pulse.count\""},
        {"", "pulsename", "is now \"pulse.name\""},
        {"", "sequence", "is now \"pulse.sequence\""},
        {"", "pulseiter", "is now \"pulse.per_qubit\""},
        {"", "pulse_time", "is now \"pulse.duration\""},
        {"", "detuning_factor", "is now \"pulse.detuning_factor\""},
        // DEER bath pulses: adding_MECCE's flat keys (src/option.cpp:677-730 there) are the block pulse.bath here
        {"", "bnpulse", "is now \"pulse.bath.count\""},
        {"", "bsequence", "is now \"pulse.bath.sequence\""},
        {"", "bpulse_defect", "is now \"pulse.bath.defect\""},
        {"", "bpulse_energy_shift", "is now \"pulse.bath.energy_shift\""},
        {"", "bspin", "is now \"pulse.bath.spin\""},
        {"", "balphams", "is now \"pulse.bath.alphams\""},
        {"", "bbetams", "is now \"pulse.bath.betams\""},
        // the master-equation dissipators: adding_MECCE's top-level key is a key of the method's hamiltonian section here
        {"", "jump_operators", "is now \"hamiltonian.meCCE.jump_operators\" (or hamiltonian.megCCE.jump_operators)"},
        {"", "gcce_refinement_max_sweeps", "is now \"gcce_refinement.max_sweeps\""},
        {"", "cluster_contribution_max_mb", "is now \"cluster_contributions.max_mb\""},
        {"Qubit", "nqubit", "was removed: the number of qubits is the length of Qubit.qubit"},
        {"Qubit", "overhaus", "is now \"hamiltonian.common.qubit_bath_mean_field\""},
        {"Qubit", "alphaidx", "is now \"Qubit.alpha_state\": {\"index\": n}"},
        {"Qubit", "betaidx", "is now \"Qubit.beta_state\": {\"index\": n}"},
        {"Qubit", "psia", "is now \"Qubit.alpha_state\": {\"vector\": [...]}"},
        {"Qubit", "psib", "is now \"Qubit.beta_state\": {\"vector\": [...]}"},
    };
    static const Renamed species[] = {
        {"", "dfname", "is now \"name\""},
        {"", "apprx", "was removed (it changed nothing)"},
        {"", "naddspin", "is now \"n_nuclear_spins\""},
        {"", "navaax", "is now \"n_configurations\""},
        {"", "eqs", "is now \"quadrupole_moments\""},
        {"", "rxyzs", "is now \"nuclear_positions\""},
        {"", "hypf", "is now \"hyperfine\""},
    };
    std::vector<std::string> found;
    for (const Renamed& r : renamed){
        cJSON* block = *r.block ? untrackedItem(root,r.block) : root;
        if (untrackedItem(block,r.key))
            found.push_back(std::string("input key \"") + (*r.block ? std::string(r.block) + "." : "") + r.key + "\" " + r.now);
    }
    if (cJSON* method = untrackedItem(root,"method")){
        const char* v = cJSON_IsString(method) ? method->valuestring : "";
        if (!strcasecmp(v,"cce")) found.push_back(std::string("method \"") + v + "\" is now \"cCCE\"");
        for (const char* shape : {"pcce", "dsj", "itb", "dsjitb"})
            if (!strcasecmp(v,shape))
                found.push_back(std::string("method \"") + v + "\" is now \"method\": \"cCCE\" with \"cluster_mode\": \""
                                + shape + "\"");
    }
    if (cJSON* mode = untrackedItem(root,"cluster_mode"))
        if (cJSON_IsString(mode) && !strcasecmp(mode->valuestring,"adaptive"))
            found.push_back("cluster_mode \"adaptive\" is now \"bundle\" (the bundle partition; the \"adaptive\" block "
                            "asks for the refinement)");
    if (cJSON* item = untrackedItem(root,"gcce_refinement"))
        if (cJSON_IsBool(item)) found.push_back("\"gcce_refinement\": true|false is now \"gcce_refinement\": {\"enabled\": true|false}");
    if (cJSON* item = untrackedItem(root,"cluster_contributions"))
        if (cJSON_IsBool(item))
            found.push_back("\"cluster_contributions\": true|false is now \"cluster_contributions\": {\"enabled\": true|false}");
    if (cJSON* defect = untrackedItem(root,"Defect")){
        cJSON* list = defect;
        if (cJSON_IsArray(defect))
            found.push_back("\"Defect\": [...] is now \"Defect\": {\"species\": [...]} (with \"state_file\" and \"axis_file\" "
                            "beside the list)");
        else list = untrackedItem(defect,"species");
        int i = 0;
        for (cJSON* entry = cJSON_IsArray(list) ? list->child : NULL; entry; entry = entry->next, ++i)
            for (const Renamed& r : species)
                if (untrackedItem(entry,r.key))
                    found.push_back("input key \"Defect.species[" + std::to_string(i) + "]." + r.key + "\" " + r.now);
    }
    if (found.empty()) return;
    if (rank==0){
        for (const std::string& f : found) ccexRefusalPrint(stderr,"Error: %s\n",f.c_str());
        ccexRefusalPrint(stderr,"%s (docs/user/InputFile.md has the table of old and new names)\n",kConvertHint);
    }
    ccexRefuse();
}

// The calculation method and the cluster shape (A47): "method" is cCCE, gCCE, meCCE or megCCE (the master-equation
// methods, adding_MECCE), "cluster_mode" standard, pcce, bundle, dsj, itb or dsjitb. Inside, the code keeps the names it
// compares against - CCE, gCCE, pCCE, dsj, itb, dsjitb, meCCE, megCCE, pgCCE, pmeCCE, pmegCCE, and cluster_mode adaptive for
// bundle - so nothing past the readers changed (and the -P inventory writes what it wrote).
// cluster_mode pcce runs with cCCE, gCCE, meCCE and megCCE: inside pCCE, pgCCE, pmeCCE and pmegCCE (methodUsesPcceClusters
// and methodKernel, general.h); the adding_MECCE branch's method names pgCCE, pmeCCE and pmegCCE are refused with the new
// form. gCCE also runs with bundle; dsj, itb and dsjitb run with cCCE only (gCCE with them is not implemented yet;
// improvements.md A47). noise is a reserved method name (A49): accepted by the reader, refused as not implemented yet.
// A method name: cCCE, gCCE, meCCE or megCCE (true for gCCE); reserved, old and unknown names are refused.
static bool checkedMethodName(const char* method, const char* where){
    const bool gcce = !strcasecmp(method,"gcce");
    const bool mecce = !strcasecmp(method,"mecce"), megcce = !strcasecmp(method,"megcce");
    if (!strcasecmp(method,"noise")){
        ccexRefusalPrint(stderr,"Error: %s \"%s\" is not implemented yet (a reserved name; the methods that run are cCCE, gCCE, meCCE "
                       "and megCCE)\n",where,method);
        ccexRefuse();
    }
    if (!gcce && !mecce && !megcce && strcasecmp(method,"ccce")){
        static const char* branch[][2] = {{"pgcce","gCCE"}, {"pmecce","meCCE"}, {"pmegcce","megCCE"}};   // adding_MECCE
        const char* now = nullptr;
        for (const auto& b : branch) if (!strcasecmp(method,b[0])) now = b[1];
        if (now)
            ccexRefusalPrint(stderr,"Error: %s \"%s\" is now method \"%s\" with \"cluster_mode\": \"pcce\" in the input\n",where,method,now);
        else if (!strcasecmp(method,"cce"))
            ccexRefusalPrint(stderr,"Error: %s \"%s\" is now \"cCCE\"\n",where,method);
        else if (!strcasecmp(method,"pcce") || !strcasecmp(method,"dsj") || !strcasecmp(method,"itb") || !strcasecmp(method,"dsjitb"))
            ccexRefusalPrint(stderr,"Error: %s \"%s\" is now method \"cCCE\" with \"cluster_mode\": \"%s\" in the input\n",where,method,
                    lowerKey(method).c_str());
        else
            ccexRefusalPrint(stderr,"Error: %s must be cCCE, gCCE, meCCE or megCCE (noise is reserved, not implemented yet), not \"%s\"\n",
                    where,method);
        ccexRefuse();
    }
    return gcce;
}

static void internalMethod(const char* method, const char* shape, const char* where, char* methodOut, char* modeOut){
    const bool gcce = checkedMethodName(method,where);
    const bool mecce = !strcasecmp(method,"mecce"), megcce = !strcasecmp(method,"megcce");
    static const char* shapes[] = {"standard", "pcce", "bundle", "dsj", "itb", "dsjitb"};
    int k = 0;
    while (k < 6 && strcasecmp(shape,shapes[k])) ++k;
    if (k == 6){
        if (!strcasecmp(shape,"adaptive"))
            ccexRefusalPrint(stderr,"Error: cluster_mode \"adaptive\" is now \"bundle\"\n");
        else
            ccexRefusalPrint(stderr,"Error: cluster_mode must be standard, pcce, bundle, dsj, itb or dsjitb, not \"%s\"\n",shape);
        ccexRefuse();
    }
    if ((mecce || megcce) && k != 0 && k != 1){
        ccexRefusalPrint(stderr,"Error: method %s with cluster_mode \"%s\" is not implemented (cluster_mode standard or pcce)\n",
                mecce ? "meCCE" : "megCCE", shapes[k]);
        ccexRefuse();
    }
    if (gcce && k != 0 && k != 1 && k != 2){
        ccexRefusalPrint(stderr,"Error: method gCCE with cluster_mode \"%s\" is not implemented yet; cluster_mode %s runs with method cCCE\n",
                shapes[k],shapes[k]);
        ccexRefuse();
    }
    static const char* internal[] = {"CCE", "pCCE", "CCE", "dsj", "itb", "dsjitb"};
    if (k == 1 && (gcce || mecce || megcce)) strcpy(methodOut,gcce ? "pgCCE" : mecce ? "pmeCCE" : "pmegCCE");
    else strcpy(methodOut,gcce ? "gCCE" : mecce ? "meCCE" : megcce ? "megCCE" : internal[k]);
    strcpy(modeOut,k == 2 ? "adaptive" : "standard");
}

// The last -m of the command line, or NULL (cJSON_setMethodOverride). It replaces the file's method before the
// combination with cluster_mode is judged, so a file with gCCE + pcce runs under -m cCCE (review 20261013-180049 I02).
static const char* methodOverride = NULL;

static void readMethod(cJSON* root, char* methodOut, char* modeOut){
    cJSON* mode = cJSON_GetObjectItem(root,"cluster_mode");
    if (mode && (!cJSON_IsString(mode) || !mode->valuestring)){
        ccexRefusalPrint(stderr,"Error: cluster_mode must be a string (standard, pcce, bundle, dsj, itb, dsjitb)\n");
        ccexRefuse();
    }
    const char* method = cJSON_ReadString(root,(char*)"method",true,(char*)"cCCE");
    if (methodOverride){
        checkedMethodName(method,"method");   // the file's own value is still a method name
        method = methodOverride;
    }
    internalMethod(method,mode ? mode->valuestring : "standard",methodOverride ? "-m" : "method",methodOut,modeOut);
}

void cJSON_setMethodOverride(const char* method){
    methodOverride = method;
}

void cJSON_applyMethodOption(Config* cnf, Cluster* cls, const char* method){
    // -m changes the method and keeps the input's cluster_mode
    const char* shape = methodUsesPcceClusters(cnf->method) ? "pcce" : "standard";
    for (const char* s : {"dsj", "itb", "dsjitb"}) if (!strcasecmp(cnf->method,s)) shape = s;
    if (!strcasecmp(cnf->cluster_mode,"adaptive")) shape = "bundle";
    char internal[MAX_CHARARRAY_LENGTH], mode[16];
    internalMethod(method,shape,"-m",internal,mode);
    Config_setMethod(cnf,internal);
    Cluster_setMethod(cls,internal);
}

// An input block ("pulse", "hamiltonian", ...): absent is NULL, anything but an object is refused.
static cJSON* readBlock(cJSON* parent, const char* key){
    cJSON* item = parent ? cJSON_GetObjectItem(parent,key) : NULL;
    if (item && !cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: %s must be an object\n",shownKey(parent,key).c_str());
        ccexRefuse();
    }
    return item;
}

// rbath and rdip: a number (the outer radius) or [inner, outer] (A47; the inner radius was rbathcut, rdipcut). The outer
// item comes back, the inner one in *inner (NULL for a number).
static cJSON* radiusItem(cJSON* root, const char* key, cJSON** inner){
    cJSON* item = cJSON_GetObjectItem(root,key);
    *inner = NULL;
    if (!cJSON_IsArray(item)) return item;
    if (cJSON_GetArraySize(item) != 2){
        ccexRefusalPrint(stderr,"Error: %s must be a number or [inner, outer]\n",key);
        ccexRefuse();
    }
    *inner = cJSON_GetArrayItem(item,0);
    return cJSON_GetArrayItem(item,1);
}

// cJSON_ReadFloat of an item that is there
static float readFloatItem(cJSON* item, const char* label){
    refuseWrongType(item,label,cJSON_IsNumber(item),"a number");
    return finiteNumber(item,label,-1,-1,true);
}

static DefectCoordinateFrame readDefectCoordinateFrame(cJSON* defect, const char* dfname,
                                                       bool* was_explicit){

    cJSON* item = cJSON_GetObjectItem(defect,"coordinate_frame");
    *was_explicit = (item != NULL);

    // Backward compatibility: legacy Defect tensors were consumed as already being in
    // the computational frame, because there was no Defect frame conversion at all.
    if (item == NULL){ return DEFECT_COORDINATE_FRAME_QUBIT; }

    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].coordinate_frame must be \"bath\" or \"qubit\"\n",dfname);
        ccexRefuse();
    }
    if (strcasecmp(item->valuestring,"bath") == 0){ return DEFECT_COORDINATE_FRAME_BATH; }
    if (strcasecmp(item->valuestring,"qubit") == 0){ return DEFECT_COORDINATE_FRAME_QUBIT; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].coordinate_frame = \"%s\" is not supported\n",
            dfname,item->valuestring);
    ccexRefusalPrint(stderr,"Use \"bath\" for crystal/bath-file components or \"qubit\" for computational-frame components.\n");
    ccexRefuse();
}

static double readDefectFrequencyScaleToMHz(cJSON* object, const char* dfname,
                                            const char* field_name){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; } // documented default
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.unit must be a string\n",dfname,field_name);
        ccexRefuse();
    }

    if (strcasecmp(item->valuestring,"Hz")  == 0){ return 1.0e-6; }
    if (strcasecmp(item->valuestring,"kHz") == 0){ return 1.0e-3; }
    if (strcasecmp(item->valuestring,"MHz") == 0){ return 1.0; }
    if (strcasecmp(item->valuestring,"GHz") == 0){ return 1.0e3; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.unit = \"%s\" is not supported\n",
            dfname,field_name,item->valuestring);
    ccexRefusalPrint(stderr,"Supported frequency units are Hz, kHz, MHz and GHz.\n");
    ccexRefuse();
}

static double readIntmapFrequencyScaleToKHz(cJSON* intmap, const char* qubit1_name,
                                            const char* qubit2_name){

    cJSON* item = cJSON_GetObjectItem(intmap,"unit");
    if (item == NULL){ return 1.0; } // legacy Qubit.intmap tensors are in kHz
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].unit must be a string\n",
                qubit1_name,qubit2_name);
        ccexRefuse();
    }

    if (strcasecmp(item->valuestring,"Hz")  == 0){ return 1.0e-3; }
    if (strcasecmp(item->valuestring,"kHz") == 0){ return 1.0; }
    if (strcasecmp(item->valuestring,"MHz") == 0){ return 1.0e3; }
    if (strcasecmp(item->valuestring,"GHz") == 0){ return 1.0e6; }

    ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].unit = \"%s\" is not supported\n",
            qubit1_name,qubit2_name,item->valuestring);
    ccexRefusalPrint(stderr,"Supported frequency units are Hz, kHz, MHz and GHz.\n");
    ccexRefuse();
}

static double readFiniteNumber(cJSON* item, const char* what){
    if (!cJSON_IsNumber(item) || !std::isfinite(item->valuedouble)){
        ccexRefusalPrint(stderr,"Error: %s must be a finite JSON number\n",what);
        ccexRefuse();
    }
    return item->valuedouble;
}

static MatrixXcd readDefectTensor3x3(cJSON* item, const char* what, double scale){

    if (!cJSON_IsArray(item)){
        ccexRefusalPrint(stderr,"Error: %s must be a flat 9-number or nested 3x3 array\n",what);
        ccexRefuse();
    }

    MatrixXcd tensor = MatrixXcd::Zero(3,3);
    int outer_size = cJSON_GetArraySize(item);

    if (outer_size == 9){
        for (int k=0; k<9; k++){
            char element_label[220];
            snprintf(element_label,sizeof(element_label),"%s[%d]",what,k);
            tensor(k/3,k%3) = doublec(scale * readFiniteNumber(cJSON_GetArrayItem(item,k),element_label),0.0);
        }
        return tensor;
    }

    if (outer_size == 3){
        for (int i=0; i<3; i++){
            cJSON* row = cJSON_GetArrayItem(item,i);
            if (!cJSON_IsArray(row) || cJSON_GetArraySize(row) != 3){
                ccexRefusalPrint(stderr,"Error: %s row %d must contain exactly 3 numbers\n",what,i);
                ccexRefuse();
            }
            for (int j=0; j<3; j++){
                char element_label[220];
                snprintf(element_label,sizeof(element_label),"%s[%d][%d]",what,i,j);
                tensor(i,j) = doublec(scale * readFiniteNumber(cJSON_GetArrayItem(row,j),element_label),0.0);
            }
        }
        return tensor;
    }

    ccexRefusalPrint(stderr,"Error: %s must be a flat 9-number or nested 3x3 array\n",what);
    ccexRefuse();
}

static void readDefectZfs(cJSON* defect, DefectArray* dfa, int idf, int nconfig,
                          bool coordinate_frame_explicit){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"zfs");

    if (item == NULL){
        cJSON_ReadDefectInfo_IntCharMatrixXcd1d(defect,(char*)"zfs",9,&df->zfs,nconfig+1);
        DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,DEFECT_ZFS_NONE);
        return;
    }

    if (cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharMatrixXcd1d(defect,(char*)"zfs",9,&df->zfs,nconfig+1);
        DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,DEFECT_ZFS_INDEXED_LEGACY);
        return;
    }

    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs must be an object or the legacy indexed array\n",df->dfname);
        ccexRefuse();
    }
    if (!coordinate_frame_explicit){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s] uses an object-form zfs but has no coordinate_frame\n",df->dfname);
        ccexRefusalPrint(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the tensor basis is explicit.\n");
        ccexRefuse();
    }
    if (nconfig < 1){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s] has an object-form zfs but n_configurations = %d\n",df->dfname,nconfig);
        ccexRefusalPrint(stderr,"At least one indexed configuration is required.\n");
        ccexRefuse();
    }

    cJSON* D_item = cJSON_GetObjectItem(item,"D");
    cJSON* E_item = cJSON_GetObjectItem(item,"E");
    cJSON* tensor_item = cJSON_GetObjectItem(item,"tensor");
    if ((D_item != NULL) == (tensor_item != NULL)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs must contain exactly one of D or tensor\n",df->dfname);
        ccexRefuse();
    }
    if (tensor_item != NULL && E_item != NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs.E cannot be combined with zfs.tensor\n",df->dfname);
        ccexRefuse();
    }

    double scale_to_mhz = readDefectFrequencyScaleToMHz(item,df->dfname,"zfs");
    MatrixXcd zfs = MatrixXcd::Zero(3,3);
    DefectZfsInputMode mode;

    if (D_item != NULL){
        if (!df->axis_set){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs uses D/E form but Defect.species[%s].axis is missing\n",
                    df->dfname,df->dfname);
            ccexRefusalPrint(stderr,"The symmetry axis is required to construct the Cartesian ZFS tensor.\n");
            ccexRefuse();
        }

        char what[160];
        snprintf(what,sizeof(what),"Defect.species[%s].zfs.D",df->dfname);
        double D_mhz = scale_to_mhz * readFiniteNumber(D_item,what);
        double E_mhz = 0.0;
        if (E_item != NULL){
            snprintf(what,sizeof(what),"Defect.species[%s].zfs.E",df->dfname);
            E_mhz = scale_to_mhz * readFiniteNumber(E_item,what);
        }
        if (fabs(E_mhz) > 1.0e-12){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs.E is non-zero (%g MHz)\n",df->dfname,E_mhz);
            ccexRefusalPrint(stderr,"An axis fixes z but not the transverse x/y directions needed by E. Use zfs.tensor instead.\n");
            ccexRefuse();
        }

        // H_ZFS = S . [D (n n^T - I/3)] . S. For n=[1,1,1]/sqrt(3),
        // D=2870 MHz gives the familiar +/-956.6667 MHz off-diagonal tensor.
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                double value = D_mhz * (df->axis[i]*df->axis[j] - ((i == j) ? 1.0/3.0 : 0.0));
                zfs(i,j) = doublec(value,0.0);
            }
        }
        mode = DEFECT_ZFS_SHARED_DE;
    }else{
        char tensor_label[160];
        snprintf(tensor_label,sizeof(tensor_label),"Defect.species[%s].zfs.tensor",df->dfname);
        zfs = readDefectTensor3x3(tensor_item,tensor_label,scale_to_mhz);
        mode = DEFECT_ZFS_SHARED_TENSOR;
    }

    df->zfs[0] = MatrixXcd::Zero(3,3);
    for (int iax=1; iax<=nconfig; iax++){ df->zfs[iax] = zfs; }
    DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,mode);
}

static void readDefectDetuning(cJSON* defect, DefectArray* dfa, int idf, int nconfig){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"detuning");

    // The absent and legacy-array paths preserve the original MHz behavior exactly.
    if (item == NULL || cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharDouble(defect,(char*)"detuning",&df->detuning,nconfig+1);
        return;
    }

    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning must be an object or the legacy indexed array\n",
                df->dfname);
        ccexRefuse();
    }

    cJSON* values = cJSON_GetObjectItem(item,"values");
    if (!cJSON_IsArray(values)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values must be an indexed array\n",df->dfname);
        ccexRefusalPrint(stderr,"Example: \"detuning\": {\"values\": [[1,\"e\",-2.14]], \"unit\": \"MHz\"}\n");
        ccexRefuse();
    }

    for (int iax=0; iax<=nconfig; iax++){ df->detuning[iax] = 0.0; }
    double scale_to_mhz = readDefectFrequencyScaleToMHz(item,df->dfname,"detuning");

    int nentry = cJSON_GetArraySize(values);
    for (int i=0; i<nentry; i++){
        cJSON* entry = cJSON_GetArrayItem(values,i);
        if (!cJSON_IsArray(entry) || cJSON_GetArraySize(entry) != 3){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values[%d] must be [configuration, \"e\", value]\n",
                    df->dfname,i);
            ccexRefuse();
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values[%d] configuration must be an integer in 1..%d\n",
                    df->dfname,i,nconfig);
            ccexRefuse();
        }

        cJSON* spin_item = cJSON_GetArrayItem(entry,1);
        if (!cJSON_IsString(spin_item) || spin_item->valuestring == NULL
            || strcasecmp(spin_item->valuestring,"e") != 0){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values[%d] spin label must be \"e\"\n",
                    df->dfname,i);
            ccexRefuse();
        }

        char what[180];
        snprintf(what,sizeof(what),"Defect.species[%s].detuning.values[%d][2]",df->dfname,i);
        df->detuning[iax_item->valueint] =
            scale_to_mhz * readFiniteNumber(cJSON_GetArrayItem(entry,2),what);
    }
}

static double readDefectGyroScaleToRadkHzPerG(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].gyros.unit must be a string\n",dfname);
        ccexRefuse();
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"radkHz/G") == 0 || strcasecmp(unit,"rad/ms/G") == 0){ return 1.0; }
    if (strcasecmp(unit,"kHz/G") == 0){ return 2.0*M_PI; }
    if (strcasecmp(unit,"MHz/G") == 0){ return 2.0*M_PI*1.0e3; }
    if (strcasecmp(unit,"MHz/T") == 0){ return 2.0*M_PI*1.0e-1; }
    if (strcasecmp(unit,"Hz/T") == 0){ return 2.0*M_PI*1.0e-7; }
    if (strcasecmp(unit,"radMHz/T") == 0 || strcasecmp(unit,"rad/us/T") == 0){ return 1.0e-1; }
    if (strcasecmp(unit,"rad/s/T") == 0){ return 1.0e-7; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].gyros.unit = \"%s\" is not supported\n",dfname,unit);
    ccexRefusalPrint(stderr,"Supported gyro units are radkHz/G, rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T, radMHz/T and rad/s/T.\n");
    ccexRefuse();
}

static double readDefectEqScaleTo1eMinus30m2(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].quadrupole_moments.unit must be a string\n",dfname);
        ccexRefuse();
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"1e-30 m^2") == 0 || strcasecmp(unit,"1e-30m^2") == 0
        || strcasecmp(unit,"10^-30 m^2") == 0 || strcasecmp(unit,"10^-30m^2") == 0
        || strcasecmp(unit,"fm^2") == 0 || strcasecmp(unit,"fm2") == 0
        || strcasecmp(unit,"10 mbarn") == 0){ return 1.0; }
    if (strcasecmp(unit,"m^2") == 0 || strcasecmp(unit,"m2") == 0){ return 1.0e30; }
    if (strcasecmp(unit,"cm^2") == 0 || strcasecmp(unit,"cm2") == 0){ return 1.0e26; }
    if (strcasecmp(unit,"barn") == 0 || strcasecmp(unit,"b") == 0){ return 1.0e2; }
    if (strcasecmp(unit,"mbarn") == 0 || strcasecmp(unit,"millibarn") == 0
        || strcasecmp(unit,"mb") == 0){ return 1.0e-1; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].quadrupole_moments.unit = \"%s\" is not supported\n",dfname,unit);
    ccexRefusalPrint(stderr,"Supported quadrupole-moment units are 1e-30 m^2, fm^2, m^2, cm^2, barn and mbarn.\n");
    ccexRefuse();
}

static double readDefectEfgScaleToHartreePerBohr2(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].efg.unit must be a string\n",dfname);
        ccexRefuse();
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"Hartree/Bohr^2") == 0 || strcasecmp(unit,"Ha/Bohr^2") == 0
        || strcasecmp(unit,"Ha/a0^2") == 0 || strcasecmp(unit,"a.u.") == 0
        || strcasecmp(unit,"au") == 0 || strcasecmp(unit,"atomic_unit") == 0){ return 1.0; }

    // Keep these constants identical to BathSpin_setQuad_fromEFG. Numerically, one
    // eV of electron potential energy per m^2 corresponds to one V/m^2; the electron
    // charge is implicit in the stored eQ convention used by that routine.
    const double hartree_to_ev = 27.211386;
    const double bohr_to_m = 0.5291772e-10;
    const double volt_per_m2_to_au = bohr_to_m*bohr_to_m/hartree_to_ev;

    if (strcasecmp(unit,"V/m^2") == 0 || strcasecmp(unit,"V/m2") == 0){
        return volt_per_m2_to_au;
    }
    if (strcasecmp(unit,"V/angstrom^2") == 0 || strcasecmp(unit,"V/angstrom2") == 0
        || strcasecmp(unit,"V/A^2") == 0 || strcasecmp(unit,"V/A2") == 0){
        return 1.0e20*volt_per_m2_to_au;
    }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].efg.unit = \"%s\" is not supported\n",dfname,unit);
    ccexRefusalPrint(stderr,"Supported EFG units are Hartree/Bohr^2 (a.u.), V/m^2 and V/angstrom^2.\n");
    ccexRefuse();
}

static cJSON* readDefectObjectValues(cJSON* object, const char* dfname,
                                     const char* field_name, int expected_count){

    cJSON* values = cJSON_GetObjectItem(object,"values");
    if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != expected_count){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values must be an array of exactly %d numbers\n",
                dfname,field_name,expected_count);
        ccexRefuse();
    }
    return values;
}

static void readDefectGyros(cJSON* defect, DefectArray* dfa, int idf,
                            int naddspin, bool have_default){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"gyros");

    if (item == NULL || cJSON_IsArray(item)){
        double* values = cJSON_ReadDouble1d(defect,(char*)"gyros",have_default,NULL,naddspin);
        DefectArray_setDefect_idf_gyros(dfa,idf,values);
        freeDouble1d(&values);
        return;
    }
    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].gyros must be an object or the legacy value array\n",df->dfname);
        ccexRefuse();
    }

    cJSON* values = readDefectObjectValues(item,df->dfname,"gyros",naddspin);
    double scale = readDefectGyroScaleToRadkHzPerG(item,df->dfname);
    for (int isp=0; isp<naddspin; isp++){
        char what[180];
        snprintf(what,sizeof(what),"Defect.species[%s].gyros.values[%d]",df->dfname,isp);
        df->gyros[isp] = scale * readFiniteNumber(cJSON_GetArrayItem(values,isp),what);
    }
}

static void readDefectEqs(cJSON* defect, DefectArray* dfa, int idf, int naddspin){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"quadrupole_moments");

    if (item == NULL || cJSON_IsArray(item)){
        double* values = cJSON_ReadDouble1d(defect,(char*)"quadrupole_moments",true,NULL,naddspin);
        if (values != NULL){
            DefectArray_setDefect_idf_eqs(dfa,idf,values);
            freeDouble1d(&values);
        }else{
            double* defaults = allocDouble1d(naddspin);
            DefectArray_setDefect_idf_eqs(dfa,idf,defaults);
            freeDouble1d(&defaults);
        }
        return;
    }
    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].quadrupole_moments must be an object or the legacy value array\n",df->dfname);
        ccexRefuse();
    }

    cJSON* values = readDefectObjectValues(item,df->dfname,"quadrupole_moments",naddspin);
    double scale = readDefectEqScaleTo1eMinus30m2(item,df->dfname);
    for (int isp=0; isp<naddspin; isp++){
        char what[180];
        snprintf(what,sizeof(what),"Defect.species[%s].quadrupole_moments.values[%d]",df->dfname,isp);
        df->eqs[isp] = scale * readFiniteNumber(cJSON_GetArrayItem(values,isp),what);
    }
}

static int findDefectSpinType(Defect* df, const char* spin_name){
    for (int isp=0; isp<df->naddspin; isp++){
        if (strcasecmp(df->types[isp],spin_name) == 0){ return isp; }
    }
    return -1;
}

static void readDefectIndexedTensorWithUnit(cJSON* defect, DefectArray* dfa, int idf,
                                            const char* field_name, MatrixXcd*** target,
                                            int nconfig, int naddspin,
                                            bool coordinate_frame_explicit){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,field_name);

    if (item == NULL || cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharMatrixXcd2d(defect,(char*)field_name,9,target,
                                                df->types,nconfig+1,naddspin);
        return;
    }
    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s must be an object or the legacy indexed array\n",
                df->dfname,field_name);
        ccexRefuse();
    }
    if (!coordinate_frame_explicit){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s] uses object-form %s but has no coordinate_frame\n",
                df->dfname,field_name);
        ccexRefusalPrint(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the tensor basis is explicit.\n");
        ccexRefuse();
    }

    cJSON* values = cJSON_GetObjectItem(item,"values");
    if (!cJSON_IsArray(values)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values must be an indexed tensor array\n",
                df->dfname,field_name);
        ccexRefuse();
    }

    double scale = (strcasecmp(field_name,"hyperfine") == 0)
                 ? readDefectFrequencyScaleToMHz(item,df->dfname,"hyperfine")
                 : readDefectEfgScaleToHartreePerBohr2(item,df->dfname);

    for (int iax=0; iax<=nconfig; iax++){
        for (int isp=0; isp<naddspin; isp++){ (*target)[iax][isp] = MatrixXcd::Zero(3,3); }
    }

    int nentry = cJSON_GetArraySize(values);
    for (int i=0; i<nentry; i++){
        cJSON* entry = cJSON_GetArrayItem(values,i);
        if (!cJSON_IsArray(entry) || cJSON_GetArraySize(entry) != 3){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] must be [configuration, spin, tensor]\n",
                    df->dfname,field_name,i);
            ccexRefuse();
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] configuration must be an integer in 1..%d\n",
                    df->dfname,field_name,i,nconfig);
            ccexRefuse();
        }

        cJSON* spin_item = cJSON_GetArrayItem(entry,1);
        if (!cJSON_IsString(spin_item) || spin_item->valuestring == NULL){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] spin label must be a string\n",
                    df->dfname,field_name,i);
            ccexRefuse();
        }
        int isp = findDefectSpinType(df,spin_item->valuestring);
        if (isp < 0){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] spin type \"%s\" is not in types\n",
                    df->dfname,field_name,i,spin_item->valuestring);
            ccexRefuse();
        }

        char tensor_label[220];
        snprintf(tensor_label,sizeof(tensor_label),"Defect.species[%s].%s.values[%d][2]",
                 df->dfname,field_name,i);
        (*target)[iax_item->valueint][isp] =
            readDefectTensor3x3(cJSON_GetArrayItem(entry,2),tensor_label,scale);
    }
}

// hamiltonian.<section>.jump_operators (adding_MECCE had a top-level jump_operators): the master-equation dissipators,
// [{"bath_name": str, "operator": "-"|"+"|"z"|"x"|"y", "rate": MHz >= 0}], stored as C = sqrt(rate) S_op with the rate in
// rad*kHz. An empty list is kept (rate 0 everywhere) so that the section still names its method.
static void readJumpOperators(Config* cnf, cJSON* section, const char* name){
    for (cJSON* key = section->child; key; key = key->next)
        if (strcmp(key->string,"jump_operators")){
            ccexRefusalPrint(stderr,"Error: hamiltonian.%s has no key \"%s\" (jump_operators)\n",name,key->string);
            ccexRefuse();
        }
    cJSON* arr = cJSON_GetObjectItem(section, "jump_operators");
    if (!arr) return;
    if (cnf->jump_operators) {
        ccexRefusalPrint(stderr, "Error: jump operators are given in both hamiltonian.meCCE and hamiltonian.megCCE\n"); ccexRefuse();
    }
    if (!cJSON_IsArray(arr)) { ccexRefusalPrint(stderr, "Error: hamiltonian.%s.jump_operators must be an array\n", name); ccexRefuse(); }
    const int n = cJSON_GetArraySize(arr);
    JumpOperatorArray* joa = JumpOperatorArray_init();
    if (n > 0) JumpOperatorArray_alloc(joa, n);   // none: njump 0, ops NULL (no zero-length allocation)
    for (int i = 0; i < n; ++i) {
        cJSON* it = cJSON_GetArrayItem(arr, i);
        cJSON* nm = cJSON_GetObjectItem(it, "bath_name");
        cJSON* op = cJSON_GetObjectItem(it, "operator");
        cJSON* rt = cJSON_GetObjectItem(it, "rate");
        if (!cJSON_IsString(nm) || !nm->valuestring[0] || !cJSON_IsString(op) || !rt) {
            ccexRefusalPrint(stderr, "Error: hamiltonian.%s.jump_operators[%d] needs bath_name, operator and rate\n", name, i);
            ccexRefuse();
        }
        if (!isJumpOperatorType(op->valuestring)) {   // before setOp, which keeps 3 characters of op_type
            const bool old = !strcmp(op->valuestring, "+-") || !strcmp(op->valuestring, "-+");
            ccexRefusalPrint(stderr, "Error: hamiltonian.%s.jump_operators[%d].operator must be one of - + z x y%s\n", name, i,
                    old ? " (the old names +- and -+ are now - and +)" : "");
            ccexRefuse();
        }
        const double rate = cJSON_IsNumber(rt) ? rt->valuedouble : -1.0;
        if (!(std::isfinite(rate) && rate >= 0.0)) {
            ccexRefusalPrint(stderr, "Error: hamiltonian.%s.jump_operators[%d].rate must be a finite number >= 0 [MHz]\n", name, i);
            ccexRefuse();
        }
        JumpOperatorArray_setOp(joa, i, nm->valuestring, op->valuestring, MHZ_TO_RADKHZ(rate));
    }
    cnf->jump_operators = joa;
    copyBounded(cnf->jump_operators_section, name, sizeof(cnf->jump_operators_section), "hamiltonian section");
}

/**
 * @brief Read the option from the input file
 * @details Read &General tag  options
 * @note
 *      free the 1d values from cJSON should be free if it is not string
*/
void cJSON_readOptionConfig(Config* cnf, char* fccein){

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }else{
		;//printf("%s",cJSON_Print(root));
	}
    refuseOldInputKeys(root);   // the first reader: an old key name stops the run here (A47)

    ////////////////////////////////////////////////////////////////////////
    // General Options
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printMessage("Read Config Options ...\n");
    }

    if (rank==0){
        printMessage("  - General option-related keys : ");
        printMessage("    [ method, cluster_mode, quantity, propagator, evolution, evolution_update, backend, profile, ");
        printMessage("      cluster_contributions { enabled, max_mb }, ");
        printMessage("      gpu_options { min_dim, packed_max_dim, warp_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb }, ");
        printMessage("      probe_memory_gb, gcce_refinement { enabled, max_sweeps }, ");
        printMessage("      pcce { partition, connections, print_limit, cache { mode, file } }, ");
        printMessage("      bundle { max_dim, merge_cutoff, connections, connection_budget, max_clusters, max_work }, ");
        printMessage("      order, bfield, rbath, rdip, deltat, nstep, nstate, seed ] \n");
    }
    char method[MAX_CHARARRAY_LENGTH], clusterMode[16];
    readMethod(root,method,clusterMode);
    Config_setMethod(cnf,method); // inside: CCE, gCCE, pCCE, dsj, dsjitb, itb (from method and cluster_mode)
    strcpy(cnf->cluster_mode,clusterMode);   // inside: standard, or adaptive for cluster_mode bundle

    char* quantity = cJSON_ReadString(root,"quantity",true,"coherence");
    Config_setQuantity(cnf,quantity); // Current possible options : coherence, noise, dm

    // eigen (default) | expm (gCCE, legacy) | sparse (A2), as a string or as {"method": ..., "min_dim": ..., "tol": ...};
    // min_dim and tol belong to sparse. Config_validatePropagator checks the paths once every source is in.
    cJSON* propagatorItem = cJSON_GetObjectItem(root,"propagator");
    if (cJSON_IsObject(propagatorItem)){
        for (cJSON* key = propagatorItem->child; key; key = key->next)
            if (strcasecmp(key->string,"method") && strcasecmp(key->string,"min_dim") && strcasecmp(key->string,"tol")){
                ccexRefusalPrint(stderr,"Error: propagator has no key \"%s\" (method, min_dim, tol)\n",key->string);
                ccexRefuse();
            }
        cJSON* method = cJSON_GetObjectItem(propagatorItem,"method");
        if (!cJSON_IsString(method) || !method->valuestring){
            ccexRefusalPrint(stderr,"Error: propagator.method must be a string (eigen, expm or sparse)\n");
            ccexRefuse();
        }
        Config_setPropagator(cnf,method->valuestring);
        cJSON* minDim = cJSON_GetObjectItem(propagatorItem,"min_dim");
        cJSON* tol = cJSON_GetObjectItem(propagatorItem,"tol");
        if ((minDim || tol) && strcasecmp(cnf->propagator,"sparse")){
            ccexRefusalPrint(stderr,"Error: propagator.min_dim and propagator.tol belong to method sparse (this has %s)\n",cnf->propagator);
            ccexRefuse();
        }
        if (minDim){
            const double value = readFiniteNumber(minDim,"propagator.min_dim");
            if (value < 1 || value > 2147483647.0 || std::floor(value) != value){
                ccexRefusalPrint(stderr,"Error: propagator.min_dim must be a positive 32-bit integer\n");
                ccexRefuse();
            }
            cnf->sparse_min_dim = (long)value;
        }
        if (tol){
            const double value = readFiniteNumber(tol,"propagator.tol");
            if (!(value > 0.0 && value < 1.0)){
                ccexRefusalPrint(stderr,"Error: propagator.tol must be a number between 0 and 1 (exclusive)\n");
                ccexRefuse();
            }
            cnf->sparse_tol = value;
        }
    }else{
        if (propagatorItem && !cJSON_IsString(propagatorItem)){
            ccexRefusalPrint(stderr,"Error: propagator must be a string (eigen, expm or sparse) or an object {\"method\": ..., "
                           "\"min_dim\": ..., \"tol\": ...}\n");
            ccexRefuse();
        }
        char* propagator = cJSON_ReadString(root,"propagator",true,"eigen");
        Config_setPropagator(cnf,propagator);
    }

    // vector is the DEFAULT, but it only exists where its preconditions hold: gCCE, a pure
    // rho0 (nstate > 0) and propagator=eigen. Do NOT decide that here -- this function
    // runs from inside main.cpp's getopt loop, so -m and -N have not been applied yet.
    // Only record whether the key was written; Config_resolveEvolution settles it once
    // every source is in.
    char* evolution = cJSON_ReadString(root,"evolution",true,"vector");
    Config_setEvolution(cnf,evolution); // gCCE only : vector (default) | matrix (legacy)
    cnf->evolution_isdefault = (cJSON_GetObjectItem(root,"evolution") == NULL);
    // exact | recursive (A16 b); without the key Config_validateEvolutionUpdate picks the default once every source is in (A29)
    cJSON* evolutionUpdate = cJSON_GetObjectItem(root,"evolution_update");
    if (evolutionUpdate){
        if (!cJSON_IsString(evolutionUpdate) || !evolutionUpdate->valuestring ||
            (strcasecmp(evolutionUpdate->valuestring,"exact") &&
             strcasecmp(evolutionUpdate->valuestring,"recursive"))){
            ccexRefusalPrint(stderr,"Error: evolution_update must be exact or recursive\n");
            ccexRefuse();
        }
        strcpy(cnf->evolution_update,evolutionUpdate->valuestring);
    }

    cJSON* backend = cJSON_GetObjectItem(root,"backend");
    if (backend){
        if (!cJSON_IsString(backend) || !backend->valuestring ||
            (strcasecmp(backend->valuestring,"cpu") &&
             strcasecmp(backend->valuestring,"cuda") &&
             strcasecmp(backend->valuestring,"auto"))){
            ccexRefusalPrint(stderr,"Error: backend must be cpu, cuda, or auto\n");
            ccexRefuse();
        }
        strcpy(cnf->backend,backend->valuestring);
    }
    // The GPU settings are one block, "gpu_options" (user, 2026-10-12, improvements.md A26); the flat keys before it are
    // refused with the form that replaces them.
    const char* oldGpuKeys[] = {"gpu_min_dim", "gpu_batch_steps", "gpu_cache_mb", "gpu_batch_clusters", "gpu_cluster_batch_mb"};
    const char* gpuKeys[] = {"min_dim", "evolution_batch_steps", "cache_mb", "evolution_batch_size", "cluster_batch_mb", "packed_max_dim",
                             "warp_max_dim"};
    int* gpuValues[] = {&cnf->gpu_min_dim, &cnf->gpu_batch_steps, &cnf->gpu_cache_mb, &cnf->gpu_batch_clusters, &cnf->gpu_cluster_batch_mb,
                        &cnf->gpu_packed_max_dim, &cnf->gpu_warp_max_dim};
    for (int k=0; k<5; ++k)
        if (cJSON_GetObjectItem(root,oldGpuKeys[k])){
            ccexRefusalPrint(stderr,"Error: %s is now gpu_options.%s (\"gpu_options\": {\"%s\": ...})\n",oldGpuKeys[k],gpuKeys[k],gpuKeys[k]);
            ccexRefuse();
        }
    if (cJSON* block = cJSON_GetObjectItem(root,"gpu_options")){
        if (!cJSON_IsObject(block)){
            ccexRefusalPrint(stderr,"Error: gpu_options must be an object, e.g. \"gpu_options\": {\"min_dim\": 256}\n");
            ccexRefuse();
        }
        for (cJSON* item = block->child; item; item = item->next){
            int k = 0;
            while (k<7 && strcasecmp(item->string,gpuKeys[k])) ++k;
            if (k==7){
                ccexRefusalPrint(stderr,"Error: unknown key gpu_options.%s (min_dim, packed_max_dim, warp_max_dim, evolution_batch_size, "
                               "evolution_batch_steps, cluster_batch_mb, cache_mb)\n",item->string);
                ccexRefuse();
            }
        }
        for (int k=0; k<7; ++k){
            cJSON* item = cJSON_GetObjectItem(block,gpuKeys[k]);   // recorded as read (cJSON_checkInputKeys)
            if (!item) continue;
            char name[64];
            snprintf(name,sizeof name,"gpu_options.%s",gpuKeys[k]);
            const double value = readFiniteNumber(item,name);
            if (k==6){
                // the warp kernel (A40): clusters of a compiled total dimension up to this value; 0 switches it off
                if (value < 0 || value > 32 || std::floor(value) != value){
                    ccexRefusalPrint(stderr,"Error: gpu_options.warp_max_dim must be 0 (off) or an integer up to 32\n");
                    ccexRefuse();
                }
            } else if (k==5){
                // the packed kernel exists for bath dimension 2 and 4 (small_cluster_math.h); 0 switches it off
                if (value!=0 && value!=2 && value!=4){
                    ccexRefusalPrint(stderr,"Error: gpu_options.packed_max_dim must be 0 (off), 2 or 4\n");
                    ccexRefuse();
                }
            } else if (value < (k==2 ? 0 : 1) || value > 2147483647.0 || std::floor(value) != value){
                ccexRefusalPrint(stderr,"Error: %s must be a %s 32-bit integer\n",name,k==2 ? "nonnegative" : "positive");
                ccexRefuse();
            }
            *gpuValues[k] = (int)value;
        }
    }
    // A1: the memory the probes of an adaptive probe step may hold at once on a node, in GB (0: automatic)
    if (cJSON* item = cJSON_GetObjectItem(root,"probe_memory_gb")){
        const double value = readFiniteNumber(item,"probe_memory_gb");
        if (value < 0){ ccexRefusalPrint(stderr,"Error: probe_memory_gb must be nonnegative\n"); ccexRefuse(); }
        cnf->probe_memory_gb = value;
    }
    // gcce_refinement and cluster_contributions are blocks (A47): {"enabled", ...}; a block without "enabled" is on, as an
    // auto_convergence block is (main.cpp, inputAsksAutoConv)
    if (cJSON* block = readBlock(root,"gcce_refinement")){
        cnf->gcce_refinement=readStrictBool(block,"gcce_refinement","enabled",true);
        if(auto* item=cJSON_GetObjectItem(block,"max_sweeps")) {
            const double value=readFiniteNumber(item,"gcce_refinement.max_sweeps");
            if(value<1 || value>INT_MAX || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: gcce_refinement.max_sweeps must be a positive 32-bit integer\n");
                ccexRefuse();
            }
            cnf->gcce_refinement_max_sweeps=int(value);
        }
    }
    cnf->profile = readStrictBool(root,"config","profile",false);
    if (cJSON* block = readBlock(root,"cluster_contributions")){
        cnf->cluster_contributions = readStrictBool(block,"cluster_contributions","enabled",true);
        if(auto* item=cJSON_GetObjectItem(block,"max_mb")) {
            const double value=readFiniteNumber(item,"cluster_contributions.max_mb");
            if(value<1 || value>INT_MAX || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: cluster_contributions.max_mb must be a positive 32-bit integer\n");ccexRefuse();
            }
            cnf->cluster_contribution_max_mb=int(value);
        }
    }
    // checked against rdip and the method in rdip_ladder.h. Every rank reads the input, so a refusal here stops the ranks together
    // (rdipladder::refuse: one message, MPI_Finalize, exit 1); a plain exit on each rank let mpirun kill the slower rank and
    // report 1 or 255 at random, sometimes with one message (A54)
    if(auto* item=cJSON_GetObjectItem(root,"rdip_ladder")) {
        const int count=cJSON_IsArray(item) ? cJSON_GetArraySize(item) : -1;
        if(count<1 || count>64) {
            rdipladder::refuse("must be an array of 1 to 64 numbers");
        }
        for(int k=0;k<count;++k) {
            const double value=readFiniteNumber(cJSON_GetArrayItem(item,k),"rdip_ladder");
            if(value<=0 || (k>0 && float(value)<=cnf->rdip_ladder[k-1])) {
                rdipladder::refuse("must be positive and strictly increasing");
            }
            cnf->rdip_ladder[k]=float(value);
        }
        cnf->nrdip_ladder=count;
    }
    // the bundle block: cluster_mode bundle only (validateAdaptiveConfig refuses it elsewhere)
    cJSON* bundle=readBlock(root,"bundle");
    if(bundle)cnf->adaptive_options_present=true;
    const char* bundleKeys[]={"max_dim","merge_cutoff"};
    for(int k=0;k<2 && bundle;++k) {
        cJSON* item=cJSON_GetObjectItem(bundle,bundleKeys[k]);
        if(!item)continue;
        const std::string name="bundle."+std::string(bundleKeys[k]);
        double value=readFiniteNumber(item,name.c_str());
        if(k<1) {
            if(value<1 || value>2147483647.0 || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());ccexRefuse();
            }
            cnf->bundle_max_dim=int(value);
        } else {
            if(value<0) { ccexRefusalPrint(stderr,"Error: bundle.merge_cutoff must be nonnegative\n");ccexRefuse(); }
            cnf->bundle_merge_cutoff=value;
        }
    }
    const char* expansionKeys[]={"max_clusters","max_work"};
    int* expansionValues[]={&cnf->expansion_max_clusters,&cnf->expansion_max_work};
    for(int k=0;k<2 && bundle;++k)if(auto* item=cJSON_GetObjectItem(bundle,expansionKeys[k])) {
        const std::string name="bundle."+std::string(expansionKeys[k]);
        const double value=readFiniteNumber(item,name.c_str());
        if(value<1 || value>INT_MAX || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());ccexRefuse();
        }
        *expansionValues[k]=int(value);
    }
    cJSON* connections=bundle ? cJSON_GetObjectItem(bundle,"connections") : NULL;
    if(connections) {
        if(!cJSON_IsString(connections) || !connections->valuestring ||
           (strcasecmp(connections->valuestring,"distance") && strcasecmp(connections->valuestring,"auto"))) {
            ccexRefusalPrint(stderr,"Error: bundle.connections must be distance or auto\n");ccexRefuse();
        }
        strcpy(cnf->bundle_connections,connections->valuestring);
    }
    cJSON* connectionBudget=bundle ? cJSON_GetObjectItem(bundle,"connection_budget") : NULL;
    if(connectionBudget) {
        cnf->bundle_connection_budget=readFiniteNumber(connectionBudget,"bundle.connection_budget");
        if(cnf->bundle_connection_budget<0) {
            ccexRefusalPrint(stderr,"Error: bundle.connection_budget must be nonnegative\n");ccexRefuse();
        }
    }
    // adaptive.work_file: the request the refinement driver writes into a step's input; the rest of the adaptive block is
    // the refinement request itself, read in main.cpp (adaptiveArgsFromInput)
    cJSON* adaptiveBlock=readBlock(root,"adaptive");
    cJSON* adaptiveWork=adaptiveBlock ? cJSON_GetObjectItem(adaptiveBlock,"work_file") : NULL;
    if(adaptiveWork) {
        cnf->adaptive_options_present=true;
        if(!cJSON_IsString(adaptiveWork) || !adaptiveWork->valuestring || !*adaptiveWork->valuestring ||
           strlen(adaptiveWork->valuestring)>=MAX_FILEPATH) {
            ccexRefusalPrint(stderr,"Error: adaptive.work_file must be a nonempty bounded path\n");ccexRefuse();
        }
        strcpy(cnf->adaptive_work_file,adaptiveWork->valuestring);
    }
    // the pcce block: cluster_mode pcce only (partition_size, trials, max_iterations, kmeans_pp and print_trials are read
    // with the clusters, cJSON_readOptionCluster)
    cJSON* pcce=readBlock(root,"pcce");
    if(pcce && !methodUsesPcceClusters(cnf->method)) {   // the rest: standard, bundle (adaptive inside), dsj, itb, dsjitb
        ccexRefusalPrint(stderr,"Error: the pcce block belongs to cluster_mode pcce (this input has cluster_mode %s)\n",
                !strcasecmp(cnf->cluster_mode,"adaptive") ? "bundle" : !strcasecmp(cnf->method,"gcce") || !strcasecmp(cnf->method,"cce")
                || !strcasecmp(cnf->method,"mecce") || !strcasecmp(cnf->method,"megcce") ? "standard" : cnf->method);
        ccexRefuse();
    }
    if(auto* item=pcce ? cJSON_GetObjectItem(pcce,"print_limit") : NULL) {
        const double value=readFiniteNumber(item,"pcce.print_limit");
        if(value<0 || value>INT_MAX || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: pcce.print_limit must be a nonnegative 32-bit integer\n");
            ccexRefuse();
        }
        cnf->pcce_detail_limit=int(value);
    }
    cJSON* partition = pcce ? cJSON_GetObjectItem(pcce,"partition") : NULL;
    if (partition){
        if (!cJSON_IsString(partition) || !partition->valuestring ||
            (strcasecmp(partition->valuestring,"deterministic") &&
             strcasecmp(partition->valuestring,"legacy"))){
            ccexRefusalPrint(stderr,"Error: pcce.partition must be deterministic or legacy\n");
            ccexRefuse();
        }
        strcpy(cnf->pcce_partition,partition->valuestring);
    }
    cJSON* pcceConnections = pcce ? cJSON_GetObjectItem(pcce,"connections") : NULL;
    if (pcceConnections){
        if (!cJSON_IsString(pcceConnections) || !pcceConnections->valuestring ||
            (strcasecmp(pcceConnections->valuestring,"centers") &&
             strcasecmp(pcceConnections->valuestring,"spins"))){
            ccexRefusalPrint(stderr,"Error: pcce.connections must be centers or spins\n");
            ccexRefuse();
        }
        strcpy(cnf->pcce_connections,pcceConnections->valuestring);
    }

    cJSON* cache = pcce ? readBlock(pcce,"cache") : NULL;
    cJSON* cacheMode = cache ? cJSON_GetObjectItem(cache,"mode") : NULL;
    if (cacheMode){
        if (!cJSON_IsString(cacheMode) || !cacheMode->valuestring ||
            (strcasecmp(cacheMode->valuestring,"off") &&
             strcasecmp(cacheMode->valuestring,"read") &&
             strcasecmp(cacheMode->valuestring,"readwrite"))){
            ccexRefusalPrint(stderr,"Error: pcce.cache.mode must be off, read, or readwrite\n");
            ccexRefuse();
        }
        strcpy(cnf->pcce_cache_mode,cacheMode->valuestring);
    }
    cJSON* cacheFile = cache ? cJSON_GetObjectItem(cache,"file") : NULL;
    if (cacheFile){
        if (!cJSON_IsString(cacheFile) || !cacheFile->valuestring ||
            !cacheFile->valuestring[0] || strlen(cacheFile->valuestring)>=MAX_FILEPATH){
            ccexRefusalPrint(stderr,"Error: pcce.cache.file must be a nonempty path shorter than %d bytes\n",MAX_FILEPATH);
            ccexRefuse();
        }
        strcpy(cnf->pcce_cache_file,cacheFile->valuestring);
    }

    int order = cJSON_ReadInt(root,"order",false,1);
    if(!strcasecmp(cnf->cluster_mode,"adaptive")) {
        const double value=readFiniteNumber(cJSON_GetObjectItem(root,"order"),"order");
        if(value<1 || value>=INT_MAX || value>cnf->expansion_max_clusters || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: cluster_mode bundle requires a positive integer order (of bundles) up to bundle.max_clusters\n");
            ccexRefuse();
        }
    }
    Config_setOrder(cnf,order);

    // Magnetic field : This value can be set by [bx,by,bz] or bz
    // (the array reader refuses a number, so it is asked only when bfield is an array)
    float* bfield = cJSON_IsArray(cJSON_GetObjectItem(root,"bfield")) ? cJSON_ReadFloat1d(root,"bfield",true,NULL,3) : NULL;

    if (bfield != NULL){
        // set bfield as [bx,by,bz]
        Config_setBfield(cnf,bfield);
        freeFloat1d(&bfield);
    }else{
        // set bfield as [0,0,bz]
        float bfieldtry[3] = {0.0,0.0,0.0};
        float bfield_z = cJSON_ReadFloat(root,"bfield",false,-1);
        bfieldtry[2] = bfield_z;
        Config_setBfield(cnf,bfieldtry);
    }

    ////////////////////////////////////////////////////////////////////////
    // magnetic_field_axis : removed
    ////////////////////////////////////////////////////////////////////////
    // It only ever existed in unreleased development work. Silently ignoring it would
    // leave an old input running with a field direction it no longer describes, so the
    // key is refused outright rather than skipped.
    if (cJSON_GetObjectItem(root,"magnetic_field_axis") != NULL){
        ccexRefusalPrint(stderr,"Error: \"magnetic_field_axis\" was removed.\n");
        ccexRefusalPrint(stderr,"Use \"bfield\": [0,0,Bz] with coordinate_frame_rotation.\n");
        ccexRefuse();
    }
    ////////////////////////////////////////////////////////////////////////
    
    // rbath and rdip are a number (the outer radius) or [inner, outer] (A47: the inner radius was rbathcut, rdipcut)
    cJSON* rbathInner = NULL;
    cJSON* rbathOuter = radiusItem(root,"rbath",&rbathInner);
    if (!rbathOuter){ ccexRefusalPrint(stderr,"Error: rbath is not found in the input file\n"); ccexRefuse(); }
    float rbath = readFloatItem(rbathOuter,rbathInner ? "rbath[1]" : "rbath");
    Config_setRbath(cnf,rbath);

    const bool adaptive=!strcasecmp(cnf->cluster_mode,"adaptive");
    cJSON* rdipInner=NULL;
    cJSON* radius=radiusItem(root,"rdip",&rdipInner);
    const char* rdipName=rdipInner ? "rdip[1]" : "rdip";
    // Simple inputs select the policy by the presence of rdip. Explicit old
    // distance/auto controls retain their existing semantics and validation.
    if(adaptive && !radius && !connections)strcpy(cnf->bundle_connections,"checked");
    const bool autoConnections=adaptive && ( !strcasecmp(cnf->bundle_connections,"auto") ||
                                             !strcasecmp(cnf->bundle_connections,"checked"));
    cnf->adaptive_rdip_is_set=radius!=nullptr;
    float rdip;
    if(adaptive && radius) {
        rdip=float(readFiniteNumber(radius,rdipName));
        if(!std::isfinite(rdip) || rdip<=0) {
            ccexRefusalPrint(stderr,"Error: explicit cluster_mode bundle rdip must be finite and positive; omit it for unlimited auto connections\n");ccexRefuse();
        }
    } else if(radius) rdip=readFloatItem(radius,rdipName);
    else if(autoConnections) rdip=0;
    else { ccexRefusalPrint(stderr,"Error: rdip is not found in the input file\n"); ccexRefuse(); }
    Config_setRdip(cnf,rdip);

    double deltat = cJSON_ReadDouble(root,"deltat",false,-1);
    Config_setDeltat(cnf,deltat);

    int nstep = cJSON_ReadInt(root,"nstep",false,-1);
    Config_setNstep(cnf,nstep);

    float rbathcut = rbathInner ? readFloatItem(rbathInner,"rbath[0]") : 0.0f;
    Config_setRbathcut(cnf,rbathcut);

    float rdipcut = rdipInner ? readFloatItem(rdipInner,"rdip[0]") : 0.0f;
    Config_setRdipcut(cnf,rdipcut);

    int nstate = cJSON_ReadInt(root,"nstate",true,0);
    Config_setNstate(cnf,nstate);

    const int requestedSeed=cJSON_ReadInt(root,"seed",true,-1);
    const int seed=resolveMpiSeed(requestedSeed,MPI_COMM_WORLD);
    if(requestedSeed==-1 && rank==0)
        fprintf(stderr,"[seed] generated on rank 0 and shared with all ranks: %d; set seed=%d to replay this run.\n",seed,seed);
    Config_setSeed(cnf,seed);

    ////////////////////////////////////////////////////////////////////////
    // Qubit and BathFiles
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - File-related keys :");
        printMessage("    [ Qubit.file, gyrofile, bathfile, bath_shift, coordinate_frame_rotation, ");
        printMessage("      defect_axis_reference, Defect.axis_file, bath_state_file, Defect.state_file ] \n");
    }

    // Qubit.file (the qubit of a defect file; the Qubit block reader checks the rest of the block)
    cJSON* qubitBlock = readBlock(root,"Qubit");
    char* qubitfile = qubitBlock ? cJSON_ReadFilePath(qubitBlock,(char*)"file",true,NULL) : NULL;
    if (qubitfile != NULL) {
        Config_allocQubitfile(cnf);
        Config_setQubitfile(cnf,qubitfile);
    }

    // gyrofile
    char* gyrofile = cJSON_ReadFilePath(root,"gyrofile",true,NULL);
    if (gyrofile != NULL) {
        Config_allocGyrofile(cnf);
        Config_setGyrofile(cnf,gyrofile);
    }

    // bathfiles
    int length = 0;
    char** bathfiles = cJSON_ReadFilePath1d(&length, root, "bathfile", true, NULL); // get length
    if (bathfiles != NULL) {
        Config_setNbathfiles(cnf,length);
        Config_allocBathfiles(cnf);
        Config_allocBathadjust(cnf);
        for (int i = 0; i < length; i++){
            Config_setBathfiles_i(cnf,bathfiles[i],i);
        }        
        freeChar2d(&bathfiles,length);
    }
    
    // bath_shift (bathadjust until A47)
    double** bathadjustDefault = allocDouble2d(length,3); // all elements are 0.0
    double** bathadjust = cJSON_ReadDouble2d(root,(char*)"bath_shift",true,bathadjustDefault,length,3);
    
    for (int i = 0; i < length; i++){
        Config_setBathadjust_i(cnf,bathadjust[i],i);
    }

    freeDouble2d(&bathadjustDefault,length);
    freeDouble2d(&bathadjust,length);

    ////////////////////////////////////////////////////////////////////////
    // Coordinate frame rotation
    ////////////////////////////////////////////////////////////////////////
    // Absent section, or "enabled": false, leaves Config at the values set by
    // Config_init -- rotation off, both axes on +z -- so the run is bit-for-bit
    // the legacy one.
    //
    // "bath_coordinate_rotation" was the name while the feature only moved bath
    // POSITIONS. It now also decides the magnetic-field direction, so the name is
    // wrong; it is still accepted, with a warning, so that inputs written against
    // the earlier version keep working.
    const char* ROTKEY = "coordinate_frame_rotation";
    cJSON* RotSection    = cJSON_GetObjectItem(root,ROTKEY);
    cJSON* RotSectionOld = cJSON_GetObjectItem(root,"bath_coordinate_rotation");

    if (RotSection != NULL && RotSectionOld != NULL){
        ccexRefusalPrint(stderr,"Error: both \"coordinate_frame_rotation\" and \"bath_coordinate_rotation\" are present\n");
        ccexRefusalPrint(stderr,"They are the same option under two names. Keep \"coordinate_frame_rotation\".\n");
        ccexRefuse();
    }

    if (RotSection == NULL && RotSectionOld != NULL){
        RotSection = RotSectionOld;
        ROTKEY = "bath_coordinate_rotation";
        if (rank==0){
            fprintf(stderr,"Warning: \"bath_coordinate_rotation\" is deprecated -- rename it to\n"
                           "         \"coordinate_frame_rotation\". The behaviour is unchanged.\n");
        }
    }

    if (RotSection != NULL){

        // A non-object -- "coordinate_frame_rotation": true, or a path string -- has no
        // children, so every sub-key below would read as missing and the run would go on
        // silently WITHOUT the rotation the input clearly asked for.
        if (!cJSON_IsObject(RotSection)){
            ccexRefusalPrint(stderr,"Error: %s must be an object, e.g.\n",ROTKEY);
            ccexRefusalPrint(stderr,"  \"%s\" : { \"enabled\" : true, \"bath_axis\" : [0,0,1], \"qubit_axis\" : [1,1,1] }\n",ROTKEY);
            ccexRefuse();
        }

        bool rot_enabled = readStrictBool(RotSection,ROTKEY,"enabled",false);
        Config_setRot_enabled(cnf,rot_enabled);

        if (rot_enabled){
            double bath_axis[3]  = {0.0,0.0,1.0};
            double qubit_axis[3] = {0.0,0.0,1.0};

            readAxis3(RotSection,ROTKEY,"bath_axis" ,bath_axis ,false);
            readAxis3(RotSection,ROTKEY,"qubit_axis",qubit_axis,false);

            Config_setRot_bath_axis(cnf,bath_axis);
            Config_setRot_qubit_axis(cnf,qubit_axis);

            // Which qubit the rotation is taken about. Optional while nqubit is
            // restricted to 1, where there is nothing to disambiguate; when given it
            // must name an existing qubit exactly (resolved in readBathfiles, where
            // the QubitArray exists).
            // Which basis the external tensor COMPONENTS are written in. No silent
            // default: reading a bath-frame tensor as if it were already rotated (or
            // the reverse) is a double / missing rotation that nothing downstream can
            // detect, so it has to be stated. Required only where it has an effect --
            // enforced below, once hf_readmode and qd_readmode have been read too.
            cJSON* hfFrame = cJSON_GetObjectItem(RotSection,"hf_tensor_frame");
            if (hfFrame != NULL){
                if (!cJSON_IsString(hfFrame) || hfFrame->valuestring == NULL){
                    ccexRefusalPrint(stderr,"Error: %s.hf_tensor_frame must be \"bath\" or \"qubit\"\n",ROTKEY);
                    ccexRefuse();
                }
                Config_setRot_hf_tensor_frame(cnf,hfFrame->valuestring);
            }

            cJSON* qdFrame = cJSON_GetObjectItem(RotSection,"qd_tensor_frame");
            if (qdFrame != NULL){
                if (!cJSON_IsString(qdFrame) || qdFrame->valuestring == NULL){
                    ccexRefusalPrint(stderr,"Error: %s.qd_tensor_frame must be \"bath\" or \"qubit\"\n",ROTKEY);
                    ccexRefuse();
                }
                Config_setRot_qd_tensor_frame(cnf,qdFrame->valuestring);
            }

            // Which frame every Qubit.xyz is written in. One value for the whole
            // QubitArray -- per-qubit frames are not supported, and would not make
            // sense: a single R and r0 move every qubit and every bath spin together.
            cJSON* posFrame = cJSON_GetObjectItem(RotSection,"qubit_position_frame");
            if (posFrame != NULL){
                if (!cJSON_IsString(posFrame) || posFrame->valuestring == NULL){
                    ccexRefusalPrint(stderr,"Error: %s.qubit_position_frame must be the string \"bath\"\n",ROTKEY);
                    ccexRefuse();
                }
                Config_setRot_qubit_position_frame(cnf,posFrame->valuestring);
            }

            cJSON* refItem = cJSON_GetObjectItem(RotSection,"reference_qubit");
            if (refItem != NULL){
                if (!cJSON_IsString(refItem) || refItem->valuestring == NULL){
                    ccexRefusalPrint(stderr,"Error: %s.reference_qubit must be a qubit name (a string)\n",ROTKEY);
                    ccexRefuse();
                }
                Config_setRot_reference_qubit(cnf,refItem->valuestring);
            }
        }
    }

    // Defect.axis_file and Defect.state_file (avaaxfile, exstatefile until A47) sit beside the species list
    cJSON* defectBlock = readBlock(root,"Defect");
    char* avaaxfile = defectBlock ? cJSON_ReadFilePath(defectBlock,(char*)"axis_file",true,NULL) : NULL;
    if (avaaxfile != NULL) {
        Config_allocAvaaxfile(cnf);
        Config_setAvaaxfile(cnf,avaaxfile);
    }

    char* statefile = cJSON_ReadFilePath(root,(char*)"bath_state_file",true,NULL);
    if (statefile != NULL) {
        Config_allocStatefile(cnf);
        Config_setStatefile(cnf,statefile);
    }

    char* exstatefile = defectBlock ? cJSON_ReadFilePath(defectBlock,(char*)"state_file",true,NULL) : NULL;
    if (exstatefile != NULL) {
        Config_allocExstatefile(cnf);
        Config_setExstatefile(cnf,exstatefile);
    }

    ////////////////////////////////////////////////////////////////////////
    // The hamiltonian block (A47): common to both methods, and the terms of one method
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - Hamiltonian-related keys :");
        printMessage("    hamiltonian.common : [ qubit_bath_mean_field, bath_bath_mean_field, ");
        printMessage("      hyperfine { dipole, contact, file, missing_row }, ");
        printMessage("      quadrupole (\"off\" or { structure, file, file_without_defect, cell }) ] ");
        printMessage("    hamiltonian.cCCE : [ knight, hfmedi, delta_A ] ");
        printMessage("    hamiltonian.meCCE / megCCE : [ jump_operators ] \n");
    }
    cJSON* hamiltonian = readBlock(root,"hamiltonian");
    cJSON* common = hamiltonian ? readBlock(hamiltonian,"common") : NULL;
    cJSON* cCCE = hamiltonian ? readBlock(hamiltonian,"cCCE") : NULL;
    // gCCE and the reserved method noise (A49) have no settings yet: an empty section is allowed
    for (const char* name : {"gCCE", "noise"}){
        cJSON* section = hamiltonian ? readBlock(hamiltonian,name) : NULL;
        if (section && section->child){
            ccexRefusalPrint(stderr,"Error: hamiltonian.%s has no settings yet (\"%s\")\n",name,section->child->string);
            ccexRefuse();
        }
    }
    // the master-equation methods: jump_operators only; which method the section belongs to is checked once -m
    // has been applied (main.cpp, validateHamiltonianForMethod)
    for (const char* name : {"meCCE", "megCCE"})
        if (cJSON* section = hamiltonian ? readBlock(hamiltonian,name) : NULL) readJumpOperators(cnf, section, name);
    // a cCCE setting in a gCCE run is refused once -m has been applied (main.cpp, validateHamiltonianForMethod)
    cnf->hamiltonian_ccce_given = cCCE && cCCE->child;

    // the static field of the other bath spins (nstate > 0; with nstate 0 it is not used, as vector_mean_field was not):
    // off, its z part (the default), or the vector J.m_j
    if (cJSON* item = common ? cJSON_GetObjectItem(common,"bath_bath_mean_field") : NULL){
        const char* v = cJSON_IsString(item) && item->valuestring ? item->valuestring : "";
        if (strcasecmp(v,"off") && strcasecmp(v,"z") && strcasecmp(v,"xyz")){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.bath_bath_mean_field must be \"off\", \"z\" or \"xyz\"\n");
            ccexRefuse();
        }
        cnf->vector_mean_field = !strcasecmp(v,"xyz");
        cnf->bath_mean_field_off = !strcasecmp(v,"off");
    }

    ////////////////////////////////////////////////////////////////////////
    // DFT Hyperfine tensor: hamiltonian.common.hyperfine {dipole: point | dft, contact: false | true} is the readmode of
    // the code (0 point/false, 1 point/true, 2 dft/false, 3 dft/true)
    ////////////////////////////////////////////////////////////////////////
    cJSON* hyperfine = common ? readBlock(common,"hyperfine") : NULL;
    int hf_readmode = 0;
    if (hyperfine){
        const char* dipole = cJSON_ReadString(hyperfine,(char*)"dipole",true,(char*)"point");
        if (strcasecmp(dipole,"point") && strcasecmp(dipole,"dft")){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.hyperfine.dipole must be \"point\" or \"dft\"\n");
            ccexRefuse();
        }
        const bool contact = readStrictBool(hyperfine,"hamiltonian.common.hyperfine","contact",false);
        hf_readmode = (strcasecmp(dipole,"dft") ? 0 : 2) + (contact ? 1 : 0);
    }
    Config_setHf_readmode(cnf,hf_readmode);

    if (hf_readmode != 0){
        char* hf_tensorfile = cJSON_ReadFilePath(hyperfine,(char*)"file",false,NULL);
        Config_allocHf_tensorfile(cnf);
        Config_setHf_tensorfile(cnf,hf_tensorfile);

        // a bath spin outside the file's range: stop (0), or take the point-dipole tensor (1)
        const char* missing = cJSON_ReadString(hyperfine,(char*)"missing_row",true,(char*)"stop");
        if (strcasecmp(missing,"stop") && strcasecmp(missing,"point_dipole")){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.hyperfine.missing_row must be \"stop\" or \"point_dipole\"\n");
            ccexRefuse();
        }
        Config_setHf_ignore_oor(cnf,strcasecmp(missing,"stop") ? 1 : 0);
    }

    ////////////////////////////////////////////////////////////////////////
    // DFT Quadrupole tensor: hamiltonian.common.quadrupole is "off" or {structure: plain | strained | bubble, ...}, the
    // readmode of the code 0, or 2, 3, 4
    ////////////////////////////////////////////////////////////////////////
    cJSON* quadrupole = common ? cJSON_GetObjectItem(common,"quadrupole") : NULL;
    int qd_readmode = 0;
    if (quadrupole && !(cJSON_IsString(quadrupole) && !strcasecmp(quadrupole->valuestring,"off"))){
        if (!cJSON_IsObject(quadrupole)){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.quadrupole must be \"off\" or an object {\"structure\": \"plain\"|\"strained\"|"
                           "\"bubble\", \"file\": ..., \"file_without_defect\": ..., \"cell\": [a, b, c]}\n");
            ccexRefuse();
        }
        const char* structure = cJSON_ReadString(quadrupole,(char*)"structure",false,NULL);
        static const char* structures[] = {"plain", "strained", "bubble"};
        for (int k=0; k<3; ++k) if (!strcasecmp(structure,structures[k])) qd_readmode = 2+k;
        if (qd_readmode == 0){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.quadrupole.structure must be \"plain\", \"strained\" or \"bubble\"\n");
            ccexRefuse();
        }
    }
    Config_setQd_readmode(cnf,qd_readmode);

    char* qd_tensorfile = NULL;
    char* qd_tensorfile_woqubit = NULL;
    double* qd_cellpara = NULL;

    switch (qd_readmode){
        case 2:
            qd_tensorfile = cJSON_ReadFilePath(quadrupole,(char*)"file",false,NULL);
            Config_allocQd_tensorfile(cnf);
            Config_setQd_tensorfile(cnf,qd_tensorfile);

            break;
        case 3: case 4:
            qd_tensorfile = cJSON_ReadFilePath(quadrupole,(char*)"file",false,NULL);
            Config_allocQd_tensorfile(cnf);
            Config_setQd_tensorfile(cnf,qd_tensorfile);

            // Both keys are REQUIRED here, and were meant to be. Passing them through as
            // optional handed Config_setQd_tensorfile_woqubit and Config_setQd_cellpara a
            // NULL, which reach strcpy and memcpy -- so an input that omitted either one
            // segfaulted instead of being told what was missing. "false" as the _default
            // argument was a string literal, which is a non-null pointer and therefore
            // true, so qd_cellpara asked for the opposite of what it meant.
            qd_tensorfile_woqubit = cJSON_ReadFilePath(quadrupole,(char*)"file_without_defect",false,NULL);
            Config_allocQd_tensorfile_woqubit(cnf);
            Config_setQd_tensorfile_woqubit(cnf,qd_tensorfile_woqubit);

            double* qd_cellpara = cJSON_ReadDouble1d(quadrupole,(char*)"cell", false, NULL,3);
            Config_setQd_cellpara(cnf,qd_cellpara);
            freeDouble1d(&qd_cellpara);

            break;

    }

    ////////////////////////////////////////////////////////////////////////
    // hamiltonian.cCCE: knight, hfmedi, and delta_A (false: the hyperfine tensors of a cluster's spins are replaced by
    // their mean, so their Knight fields no longer differ; A47)
    ////////////////////////////////////////////////////////////////////////
    Config_setHfmedi(cnf, cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","hfmedi",false) : false);
    Config_setKnight(cnf, cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","knight",false) : false);
    cnf->delta_A = cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","delta_A",true) : true;

    ////////////////////////////////////////////////////////////////////////
    // Defect axis reference
    ////////////////////////////////////////////////////////////////////////
    // Read here, after the rotation section, because the only value it can take names
    // an axis that only that section defines.
    cJSON* DefAxisItem = cJSON_GetObjectItem(root,"defect_axis_reference");

    if (DefAxisItem != NULL){

        if (!cJSON_IsString(DefAxisItem) || DefAxisItem->valuestring == NULL){
            ccexRefusalPrint(stderr,"Error: defect_axis_reference must be a string\n");
            ccexRefusalPrint(stderr,"The only accepted value is \"qubit_axis\".\n");
            ccexRefuse();
        }
        if (strcasecmp(DefAxisItem->valuestring,"qubit_axis") != 0){
            ccexRefusalPrint(stderr,"Error: defect_axis_reference = \"%s\" is not available\n",
                    DefAxisItem->valuestring);
            ccexRefusalPrint(stderr,"The only accepted value is \"qubit_axis\".\n");
            ccexRefuse();
        }
        if (!Config_getRot_enabled(cnf)){
            ccexRefusalPrint(stderr,"Error: defect_axis_reference = \"qubit_axis\" needs a qubit axis,\n");
            ccexRefusalPrint(stderr,"and coordinate_frame_rotation is what defines one in the source frame.\n");
            ccexRefusalPrint(stderr,"Enable coordinate_frame_rotation, or drop defect_axis_reference.\n");
            ccexRefuse();
        }
        Config_setDefect_axis_reference(cnf,DEFECT_AXIS_REFERENCE_QUBIT_AXIS);
    }

    ////////////////////////////////////////////////////////////////////////
    // Cross-checks that need both the rotation section and the readmodes
    ////////////////////////////////////////////////////////////////////////
    if (Config_getRot_enabled(cnf)){

        // readmode 1 discards the file's anisotropic tensor and rebuilds it as a
        // point-dipole tensor from the source geometry, keeping only the isotropic
        // Fermi-contact term, which has no frame. So hf_tensor_frame decides nothing
        // there -- not required, and said out loud if it is given anyway.
        if (Config_getHf_readmode(cnf) == 1 && Config_getRot_hf_tensor_frame(cnf)[0] != '\0' && rank==0){
            printMessage("  Note: hf_tensor_frame is not applicable to hyperfine dipole point with contact true.");
            printMessage("        Only the isotropic Fermi-contact term is taken from the file, and the");
            printMessage("        point-dipole tensor is transformed from the bath frame regardless.\n");
        }

        if ((Config_getHf_readmode(cnf) == 2 || Config_getHf_readmode(cnf) == 3)
            && Config_getRot_hf_tensor_frame(cnf)[0] == '\0'){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.hyperfine dipole \"dft\" with the coordinate frame rotation on, but\n");
            ccexRefusalPrint(stderr,"coordinate_frame_rotation.hf_tensor_frame is not set.\n");
            ccexRefusalPrint(stderr,"State the basis the tensor COMPONENTS in hamiltonian.common.hyperfine.file are written in:\n");
            ccexRefusalPrint(stderr,"  \"bath\"  : the original bath Cartesian basis -- they get transformed\n");
            ccexRefusalPrint(stderr,"  \"qubit\" : already the qubit-aligned basis -- they are left alone\n");
            ccexRefusalPrint(stderr,"(The POSITION columns of that file are always read in the original frame.)\n");
            ccexRefuse();
        }
        if (Config_getQd_readmode(cnf) != 0 && Config_getRot_qd_tensor_frame(cnf)[0] == '\0'){
            ccexRefusalPrint(stderr,"Error: hamiltonian.common.quadrupole (structure %s) with the coordinate frame rotation on, but\n",
                    Config_getQd_readmode(cnf) == 2 ? "plain" : Config_getQd_readmode(cnf) == 3 ? "strained" : "bubble");
            ccexRefusalPrint(stderr,"coordinate_frame_rotation.qd_tensor_frame is not set.\n");
            ccexRefusalPrint(stderr,"State the basis the tensor COMPONENTS in hamiltonian.common.quadrupole.file are written in:\n");
            ccexRefusalPrint(stderr,"  \"bath\"  : the original bath Cartesian basis -- they get transformed\n");
            ccexRefusalPrint(stderr,"  \"qubit\" : already the qubit-aligned basis -- they are left alone\n");
            ccexRefuse();
        }
    }

    cJSON_Delete(root);    
    freeChar1d(&data);
}

/**
 * @brief Read the option from the input file
 * @details Read &Qubit tag  options
*/
void cJSON_readOptionQubitArray(QubitArray* qa, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Qubit Options ...\n");
        printMessage("  - Read values of main-key 'Qubit'");
        printMessage("    sub-key of 'Qubit' : [ qubit, intmap, alpha_state, beta_state, psi0 ] or [ file, spin, zfs, alphams, betams ] \n");
        printMessage("  - Read values of sub-key 'qubit'");
        printMessage("    sub-sub-key : [ name, spin, gyro, xyz, detuning, alphams, betams ] \n");
        printMessage("  - Read values of sub-key 'intmap'");
        printMessage("    sub-sub-key : [ between, tensor, unit, tensor_frame ] \n");
        printMessage("  - Read hamiltonian.common.qubit_bath_mean_field");
    }

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }

    int nqubitDefault = 1;

    char nameDefault[MAX_CHARARRAY_LENGTH] = "q";
    float spinDefault = 1.0;
    double gyroDefault = GAMMA_ELECTRON;
    double detuningDefault = 0.0;
    float alphaMsDefault = 1.0;
    float betaMsDefault = 0.0;
    
    MatrixXcd intmapDefault = MatrixXcd::Zero(3,3);
    MatrixXcd psiDefault;
    // the Overhauser field of the bath on the qubits (Qubit.overhaus until A47); both forms of the Qubit block read it
    cJSON* hamiltonian = readBlock(root,"hamiltonian");
    cJSON* common = hamiltonian ? readBlock(hamiltonian,"common") : NULL;
    const bool overhaus = common ? readStrictBool(common,"hamiltonian.common","qubit_bath_mean_field",false) : false;
    QubitArray_setOverhaus(qa,overhaus);

    // Qubit.file and Qubit.qubit are two descriptions of the qubits (qubitfile and the Qubit block until A47); one was
    // ignored without a word when both were given (review 20261013-readme R01: a "two-qubit" tutorial ran one qubit)
    cJSON* QubitSection = readBlock(root,"Qubit");
    if (QubitSection && cJSON_GetObjectItem(QubitSection,"file") != NULL){
        for (const char* key : {"qubit", "intmap", "alpha_state", "beta_state", "psi0"})
            if (cJSON_GetObjectItem(QubitSection,key) != NULL){
                if (rank==0){
                    ccexRefusalPrint(stderr,"Error: the Qubit block has both \"file\" and \"%s\"; keep one form (\"file\" for one "
                                   "qubit from a defect file, \"qubit\" for one or more qubits written out)\n",key);
                }
                ccexRefuse();
            }
    }

    // read Qubit.file
    char* _qubitfile = QubitSection ? cJSON_ReadFilePath(QubitSection,(char*)"file",true,NULL) : NULL;
    
    if (_qubitfile != NULL) {
        // when the qubitfile exist : doesn't read "Qubit section"
        // It is possible to read only one qubit

        ////////////////////////////////////////////////////////////////////////
        // set nqubit
        ////////////////////////////////////////////////////////////////////////
        // set one qubit
        QubitArray_setNqubit(qa,nqubitDefault);
        QubitArray_allocQubit(qa);
        
        ////////////////////////////////////////////////////////////////////////
        // set qubit properties by nqubit
        ////////////////////////////////////////////////////////////////////////
        // set qubit properties
        snprintf(nameDefault,MAX_CHARARRAY_LENGTH,"q%d",0);
        QubitArray_setQubit_i_name(qa,nameDefault,0);
        QubitArray_setQubit_i_gyro(qa,gyroDefault,0);
        QubitArray_setQubit_i_detuning(qa,detuningDefault,0);
        float spin = cJSON_ReadFloat(QubitSection,(char*)"spin",true,spinDefault);
        QubitArray_setQubit_i_spin(qa,spin,0);

        // set the qubit state (alpha, beta)
        float alphams = cJSON_ReadFloat(QubitSection,(char*)"alphams",true,alphaMsDefault);
        float betams = cJSON_ReadFloat(QubitSection,(char*)"betams",true,betaMsDefault);
        QubitArray_setQubit_i_alpha_fromMs(qa,alphams,0);
        QubitArray_setQubit_i_beta_fromMs(qa,betams,0);

        ////////////////////////////////////////////////////////////////////////
        //  Intmap
        ////////////////////////////////////////////////////////////////////////
        // set the options of qubit-related Hamiltonian tensor
        //
        // alloc Interaction map
        QubitArray_allocIntmap(qa);

        // read interaction map
        MatrixXcd tensor = cJSON_ReadTensor(QubitSection,(char*)"zfs",true,intmapDefault);
        tensor = KHZ_TO_RADKHZ(tensor);

        // Set ZFS in qubit intmap
        QubitArray_setIntmap_i_j(qa,tensor,0,0); // qubit index:0,0
        // Qubit.zfs (qzfs until A47) has no tensor_frame sub-key, and this path is single-qubit only, so it
        // stays where the legacy code left it: read as written, never rotated.
        QubitArray_setIntmap_i_j_frame(qa,INTMAP_EXPLICIT_UNSPECIFIED,0,0);

        // mediatedTerm IO
        cJSON_Delete(root);
        freeChar1d(&data);
        return;
    }else{
        // without Qubit.file : the qubits are written out in Qubit.qubit
        // It is possible to read multiple qubits
        int iqubit = 0;

        ////////////////////////////////////////////////////////////////////////
        // set nqubit: the length of Qubit.qubit (Qubit.nqubit until A47)
        ////////////////////////////////////////////////////////////////////////
        cJSON* qubitArray = QubitSection ? cJSON_GetObjectItem(QubitSection,"qubit") : NULL;
        if (!cJSON_IsArray(qubitArray) || cJSON_GetArraySize(qubitArray) < 1){
            ccexRefusalPrint(stderr, "Error: Qubit.qubit must be an array of one entry per qubit (or Qubit.file for one qubit from a "
                            "defect file)\n");
            ccexRefuse();
        }
        int nqubit = cJSON_GetArraySize(qubitArray);
        QubitArray_setNqubit(qa,nqubit);
        QubitArray_allocQubit(qa);

        ////////////////////////////////////////////////////////////////////////
        // set qubit properties by nqubit
        ////////////////////////////////////////////////////////////////////////
        cJSON* qubit;

        cJSON_ArrayForEach(qubit, qubitArray){
            if (qubit == NULL) {
                ccexRefusalPrint(stderr, "Error: The length of \"qubit\" is different from nqubit(%d) \n", nqubit);
                ccexRefuse();
            }
            // read qubit properties from the input file
            snprintf(nameDefault,MAX_CHARARRAY_LENGTH,"q%d",iqubit);
            char* name = cJSON_ReadString(qubit,"name",true,nameDefault);
            float spin = cJSON_ReadFloat(qubit,"spin",true,spinDefault);
            double gyro = cJSON_ReadDouble(qubit,"gyro",true,gyroDefault);
            double* xyz = cJSON_ReadDouble1d(qubit,"xyz",false,NULL,3);
            // detuning: a number, or the one-element array [x] the Defect form suggests
            // (13 tutorial and check_numerics inputs write [0.0]); it used to be ignored for
            // the default, so [x] with x != 0 was silently 0 before 2026-10-02 (D18)
            cJSON* detuningItem = cJSON_GetObjectItem(qubit,"detuning");
            if (cJSON_IsArray(detuningItem)){
                if (cJSON_GetArraySize(detuningItem) != 1 || !cJSON_IsNumber(cJSON_GetArrayItem(detuningItem,0))){
                    ccexRefusalPrint(stderr,"Error: Qubit.qubit[%d].detuning must be a number (or [number])\n",iqubit);
                    ccexRefuse();
                }
            }
            double detuning = cJSON_IsArray(detuningItem) ? cJSON_GetArrayItem(detuningItem,0)->valuedouble
                                                          : cJSON_ReadDouble(qubit,"detuning",true,detuningDefault); //kHz
            float alphams = cJSON_ReadFloat(qubit,"alphams",true,alphaMsDefault); 
            float betams = cJSON_ReadFloat(qubit,"betams",true,betaMsDefault); 
            // written, not defaulted, per qubit and side: alphaidx/betaidx would override them, and the run checks that what
            // was written agrees (D43; review b467fe R1)
            if ((int)qa->_alphaGiven.size() <= iqubit){ qa->_alphaGiven.resize(iqubit+1,0); qa->_betaGiven.resize(iqubit+1,0); }
            qa->_alphaGiven[iqubit] = cJSON_GetObjectItem(qubit,"alphams") != NULL;
            qa->_betaGiven[iqubit]  = cJSON_GetObjectItem(qubit,"betams") != NULL;

            // set qubit properties
            QubitArray_setQubit_i_name(qa,name,iqubit);
            QubitArray_setQubit_i_spin(qa,spin,iqubit);
            QubitArray_setQubit_i_xyz(qa,xyz,iqubit);
            QubitArray_setQubit_i_gyro(qa,gyro,iqubit);
            QubitArray_setQubit_i_detuning(qa,KHZ_TO_RADKHZ(detuning),iqubit); // radkHz
            QubitArray_setQubit_i_alpha_fromMs(qa,alphams,iqubit);
            QubitArray_setQubit_i_beta_fromMs(qa,betams,iqubit);
            iqubit++;
            freeDouble1d(&xyz);
        }
        
        ////////////////////////////////////////////////////////////////////////
        //  Intmap
        ////////////////////////////////////////////////////////////////////////
        // alloc Interaction map
        QubitArray_allocIntmap(qa);
        QubitArray_setIntmap_dipAuto(qa);

        // Set the initial interaction map
        // as dipolar interaction tensor in radkHz
        QubitArray_setIntmap_dipAuto(qa);

        // read interaction map
        cJSON* intmapArray = cJSON_GetObjectItem(QubitSection,"intmap");
        cJSON* intmap;
        cJSON_ArrayForEach(intmap, intmapArray){

            int qubit1_idx; char* qubit1_name;
            int qubit2_idx; char* qubit2_name;
            cJSON* between = cJSON_GetObjectItem(intmap,"between");
            
            // read which qubits are interacting
            if (cJSON_IsArray(between)) {
                qubit1_name = cJSON_GetStringValue(cJSON_GetArrayItem(between,0));
                qubit2_name = cJSON_GetStringValue(cJSON_GetArrayItem(between,1));
                qubit1_idx = QubitArray_getQubitIdx_fromName(qa,qubit1_name);
                qubit2_idx = QubitArray_getQubitIdx_fromName(qa,qubit2_name);
                // check if the qubit name exists
                if (qubit1_idx == -1 || qubit2_idx == -1 ) {
                    ccexRefusalPrint(stderr, "Error: The qubit name is not found in the input file\n");
                    ccexRefuse();
                }
                if (qubit1_idx > qubit2_idx){
                    ccexRefusalPrint(stderr, "Error: First qubit index is greater than second qubit index\n");
                    ccexRefusalPrint(stderr, "       qubit1_name = %s, qubit2_name = %s\n", qubit1_name, qubit2_name);
                    ccexRefusalPrint(stderr, "       qubit1_idx = %d, qubit2_idx = %d\n", qubit1_idx, qubit2_idx);
                    ccexRefusalPrint(stderr, "       qubit1_idx should be less than qubit2_idx\n");
                    ccexRefuse();
                }
            }else{
                ccexRefusalPrint(stderr, "Error: \"between\" is not found in the input file\n");
                ccexRefuse();
            }

            // read tensor properties from the input file
            // Read tensor properties from the input file. The absent-unit path is
            // intentionally kHz for backward compatibility; explicit units are
            // normalized to kHz here and then to CCEX's internal rad*kHz below.
            MatrixXcd tensor = cJSON_ReadTensor(intmap,"tensor",true,intmapDefault);
            tensor *= readIntmapFrequencyScaleToKHz(intmap,qubit1_name,qubit2_name);
            tensor = KHZ_TO_RADKHZ(tensor);

            // Which basis the components are written in. Recorded, never guessed: a ZFS
            // tensor and a pair tensor look the same, and so does one that has already
            // been rotated. Whether the omission is fatal depends on nqubit and on
            // whether the rotation is on, which validateCoordinateFrameRotationInputs
            // decides once the whole input is in.
            IntmapProvenance frame = INTMAP_EXPLICIT_UNSPECIFIED;
            cJSON* tframe = cJSON_GetObjectItem(intmap,"tensor_frame");
            if (tframe != NULL){
                if (!cJSON_IsString(tframe) || tframe->valuestring == NULL){
                    ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].tensor_frame must be \"bath\" or \"qubit\"\n",
                            qubit1_name,qubit2_name);
                    ccexRefuse();
                }
                if (strcasecmp(tframe->valuestring,"bath") == 0){
                    frame = INTMAP_EXPLICIT_BATH;
                }else if (strcasecmp(tframe->valuestring,"qubit") == 0){
                    frame = INTMAP_EXPLICIT_QUBIT;
                }else{
                    ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].tensor_frame = \"%s\" is not available\n",
                            qubit1_name,qubit2_name,tframe->valuestring);
                    ccexRefusalPrint(stderr,"  \"bath\"  : components are in the original bath Cartesian basis\n");
                    ccexRefusalPrint(stderr,"  \"qubit\" : components are already in the common computational basis\n");
                    ccexRefusalPrint(stderr,"            (reference-qubit-aligned; NOT an individual local frame)\n");
                    ccexRefuse();
                }
            }

            // set Interaction map
            QubitArray_setIntmap_i_j(qa,tensor,qubit1_idx,qubit2_idx);
            QubitArray_setIntmap_i_j_frame(qa,frame,qubit1_idx,qubit2_idx);
        }

        ////////////////////////////////////////////////////////////////////////
        //  psia, psib, psi0
        ////////////////////////////////////////////////////////////////////////
        // Qubit.alpha_state and beta_state (alphaidx/betaidx and psia/psib until A47): {"index": n}, the n-th eigenstate
        // of the qubit Hamiltonian by descending energy, or {"vector": [...]}, the state itself. A side not given is the
        // product of the qubits' alphams (betams); each side is decided on its own (QubitArray_setPsiaPsib), with no
        // order of precedence between the forms. Not both in one state.
        int dim = QubitArray_dim(qa);
        const char* sides[] = {"alpha_state", "beta_state"};
        for (int side=0; side<2; ++side){
            cJSON* state = readBlock(QubitSection,sides[side]);
            MatrixXcd vector = psiDefault;   // psiDefault is the 0 dimension matrix: not given
            if (state){
                cJSON* indexItem = cJSON_GetObjectItem(state,"index");
                cJSON* vectorItem = cJSON_GetObjectItem(state,"vector");
                if ((indexItem != NULL) == (vectorItem != NULL)){
                    ccexRefusalPrint(stderr,"Error: Qubit.%s must have one of \"index\" and \"vector\" (%s)\n",sides[side],
                            indexItem ? "it has both" : "it has neither");
                    ccexRefuse();
                }
                if (indexItem){
                    const int index = cJSON_ReadInt(state,(char*)"index",false,0);
                    if (index < 0 || index >= dim){
                        ccexRefusalPrint(stderr,"Error: Qubit.%s.index must be 0 to %d (the qubit eigenstates by descending energy)\n",
                                sides[side],dim-1);
                        ccexRefuse();
                    }
                    int*& slot = side == 0 ? qa->_alphaidx : qa->_betaidx;
                    slot = (int*)allocArray1d(1,sizeof(int));
                    *slot = index;
                }else{
                    double* values = cJSON_ReadDouble1d(state,(char*)"vector",false,NULL,dim);
                    vector = Double1dToMatrixXcd(values,dim);
                    freeDouble1d(&values);
                }
            }
            if (side == 0) QubitArray_setPsia(qa,vector);
            else           QubitArray_setPsib(qa,vector);
        }
        double* psi0mat = cJSON_ReadDouble1d(QubitSection,"psi0",true,NULL,dim);

        if (psi0mat != NULL) {
            QubitArray_setPsi0(qa,Double1dToMatrixXcd(psi0mat,dim));
            freeDouble1d(&psi0mat);
        }else{
            QubitArray_setPsi0(qa,psiDefault);
        }
    }
    cJSON_Delete(root);
    freeChar1d(&data);
}

/**
 * @brief Read the option from the input file
 * @details Read &Cluster tag  options

*/
void cJSON_readOptionCluster(Cluster* clus, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Cluster Options ...");
        printMessage("  [ order, method, cluster_mode, add_subclusters, max_clusters_per_order ] \n");
    }
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }

    int order = cJSON_ReadInt(root,"order",false,1);
    Cluster_setOrder(clus,order);

    char method[MAX_CHARARRAY_LENGTH], clusterMode[16];
    readMethod(root,method,clusterMode); // read twice (cJSON_readOptionConfig)
    Cluster_setMethod(clus,method);

    // add_subclusters and max_clusters_per_order (addsubclus and nk until A47)
    bool addsubclus = true;
    if (cJSON* item = cJSON_GetObjectItem(root,"add_subclusters")){
        if (!cJSON_IsBool(item)){
            ccexRefusalPrint(stderr,"Error: add_subclusters must be true or false\n");
            ccexRefuse();
        }
        addsubclus = cJSON_IsTrue(item);
    }
    Cluster_setAddsubclus(clus,addsubclus);

    Cluster_allocNk(clus);
    int* nk = cJSON_ReadInt1d(root,(char*)"max_clusters_per_order",true,NULL,order+1);
    if (nk!=NULL){
        Cluster_setNk(clus,nk);
    }else{;} 
    freeInt1d(&nk);
    // all element would be "0" in which we consider all pairs within rdip


    if (methodUsesPcceClusters(Cluster_getMethod(clus))){

        if (rank==0){
            printMessage("  pcce : [ partition_size, trials, max_iterations, kmeans_pp, print_trials ] \n");
        }

        // Do not dereference a null default argument even when partition_size is supplied.
        // Zero or fractional sizes also make the partition divisions invalid.
        cJSON* pcce = readBlock(root,"pcce");
        const char* keys[] = {"partition_size", "trials", "max_iterations"};
        int values[] = {0, 30000, 30000};
        for (int k=0; k<3; ++k){
            cJSON* item = pcce ? cJSON_GetObjectItem(pcce,keys[k]) : NULL;
            if (!item){
                if (k==0){
                    ccexRefusalPrint(stderr,"Error: cluster_mode pcce requires pcce.partition_size\n");
                    ccexRefuse();
                }
                continue;
            }
            const std::string name = "pcce." + std::string(keys[k]);
            double value = readFiniteNumber(item,name.c_str());
            if (value < 1 || value > 2147483647.0 || std::floor(value) != value){
                ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());
                ccexRefuse();
            }
            values[k] = (int)value;
        }
        int sK = values[0], max_trial = values[1], max_iter = values[2];
        bool kmeans_pp   = pcce ? readStrictBool(pcce,"pcce","kmeans_pp",true) : true;
        bool iter_detail = pcce ? readStrictBool(pcce,"pcce","print_trials",false) : false;

        Cluster_setSk(clus,sK);
        Cluster_setMax_trial(clus,max_trial);
        Cluster_setMax_iter(clus,max_iter);
        Cluster_setKmeans_pp(clus,kmeans_pp);
        Cluster_setIter_detail(clus,iter_detail);
    }

    cJSON_Delete(root);
    freeChar1d(&data);
}

////////////////////////// 
// Pulse HS version  //
////////////////////////// 
// One fraction of a Manual sequence (review 5f7834 R2): a JSON number, or a string "a/b" of two integers or a
// number, all of it; finite and not negative. Until 2026-10-02 a number fraction read a null string (an exception)
// and "1/2x" read as 1/2.
static double readSequenceFraction(cJSON* item, int row) {
    double value = 0;
    if (cJSON_IsNumber(item)) value = item->valuedouble;
    else if (cJSON_IsString(item) && item->valuestring) {
        const char* text = item->valuestring; char* end = nullptr;
        const char* slash = strchr(text, '/');
        if (slash) {
            errno = 0; const long num = strtol(text, &end, 10);
            const bool numOk = end == slash && slash != text;
            const long den = strtol(slash+1, &end, 10);
            if (!numOk || end == slash+1 || *end != '\0' || errno == ERANGE || den == 0) {
                ccexRefusalPrint(stderr, "Error: sequence[%d] fraction \"%s\" must be a/b with integers and b != 0\n", row, text);
                ccexRefuse();
            }
            value = double(num)/double(den);
        } else {
            errno = 0; value = strtod(text, &end);
            if (end == text || *end != '\0' || errno == ERANGE) {
                ccexRefusalPrint(stderr, "Error: sequence[%d] fraction \"%s\" is not a number\n", row, text);
                ccexRefuse();
            }
        }
    } else {
        ccexRefusalPrint(stderr, "Error: sequence[%d] fraction must be a number or a string such as \"1/2\"\n", row);
        ccexRefuse();
    }
    if (!std::isfinite(value) || value < 0) {
        ccexRefusalPrint(stderr, "Error: sequence[%d] fraction must be finite and not negative\n", row);
        ccexRefuse();
    }
    return value;
}

char parse_axis(cJSON* axis_json, int row_index) {
//char parse_axis(cJSON* axis_json) {
    if (!cJSON_IsString(axis_json)) {
        ccexRefusalPrint(stderr, "[Error] Axis at row %d is not a string.\n", row_index);
        ccexRefuse();
    }

    const char* axis_str = axis_json->valuestring;
    if (!axis_str || strlen(axis_str) == 0) {
        ccexRefusalPrint(stderr, "[Error] Axis at row %d is empty.\n", row_index);
        ccexRefuse();
    }

    // Check if the string is one of the allowed axes: X, Y, Z, I
    if (strlen(axis_str) > 1) {
        ccexRefusalPrint(stderr, "[Error] Axis at row %d is too long. Only single characters are allowed.\n", row_index);
        ccexRefuse();
    }

    char c = toupper(axis_str[0]);

    if (c == 'X' || c == 'Y' || c == 'Z' || c == 'I') {
        return c;
    }

    ccexRefusalPrint(stderr, "[Error] Invalid axis '%s'. Only 'X', 'Y', 'Z', 'I' are allowed.\n", axis_str);
    ccexRefuse();
}

std::string read_json_file(const std::string& fccein) {
    std::ifstream file(fccein, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open ccein file: " + fccein);
    file.seekg(0, std::ios::end);
    std::string contents(file.tellg(), '\0');
    file.seekg(0, std::ios::beg);
    file.read(&contents[0], contents.size());
    return contents;
}

void cJSON_readOptionPulse(Pulse* pulse, char* fccein) {

    if (rank == 0) {
        printf("\n");
        printMessage("Read Pulse Options ...");
        printMessage("  pulse : [ count, name, sequence, per_qubit, duration, detuning_factor ] \n");
    }

    char*  data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank == 0) {
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }

    // the pulse block (npulse, pulsename, sequence, pulseiter, pulse_time and detuning_factor until A47)
    cJSON* block = readBlock(root,"pulse");
    if (block == NULL){
        ccexRefusalPrint(stderr,"Error: pulse is not found in the input file (\"pulse\": {\"count\": 1} for a Hahn echo, 0 for a "
                       "Ramsey, n > 1 for CPMG)\n");
        ccexRefuse();
    }

    bool pulseiter = readStrictBool(block,"pulse","per_qubit",false);
    Pulse_setPulseiter(pulse, pulseiter);

    bool made = false;

    // Read "pulse number" (required)
    int npulse = cJSON_ReadInt(block,(char*)"count",false, -1);
    Pulse_setNpulse(pulse,npulse);

    // ================= //
    // New Parameters !! //
    // Input pulse duration time (unit : ns)
    double pulse_time     = cJSON_ReadDouble(block, (char*)"duration", true, 0.0);
    Pulse_setPulseTime(pulse,pulse_time);
    // Input pulse detuning offset (default : 0)
    double detuning_factor = cJSON_ReadDouble(block, (char*)"detuning_factor", true, 1.0);
    Pulse_setPulseDetuningFactor(pulse,detuning_factor);
    // ================= //

    char* pulsename = cJSON_ReadString(block,(char*)"name",true,(char*)"None");
    Pulse_setPulsename(pulse,pulsename);

    Pulse_allocSequence(pulse);
    Pulse_allocAxes(pulse);
    Pulse_allocAngles(pulse);
    Pulse_allocSequenceIndices(pulse);

    cJSON* sequence_array = cJSON_GetObjectItem(block, "sequence");
    if (sequence_array == NULL){
        allocateDefaultSequence(pulse);
        made = true;
    }
    if (strcasecmp(pulse->pulsename, "Manual") == 0) {
        if (!(sequence_array && cJSON_IsArray(sequence_array))) {
            ccexRefusalPrint(stderr, "[Error] pulse.sequence is missing or not a valid JSON array. [case pulse name = manual]\n");
            ccexRefuse();
        }

        int count = cJSON_GetArraySize(sequence_array);
        if (count != (pulse->npulse) + 1) {
            throw std::runtime_error("Mismatch between pulse.sequence length and pulse.count + 1");
        }

        double total_frac = 0.0;
        for (int i = 0; i < (pulse->npulse)+1; ++i) {
            cJSON* item = cJSON_GetArrayItem(sequence_array, i);
            // every row is [fraction, axis, angle]; a malformed row used to be skipped and leave zeros (R2)
            if (!cJSON_IsArray(item) || cJSON_GetArraySize(item) != 3) {
                ccexRefusalPrint(stderr, "Error: sequence[%d] must be [fraction, axis, angle]\n", i);
                ccexRefuse();
            }
            cJSON* frac  = cJSON_GetArrayItem(item, 0);
            cJSON* axis  = cJSON_GetArrayItem(item, 1);
            cJSON* angle = cJSON_GetArrayItem(item, 2);

            pulse->sequence[i][0] = total_frac;
            double fraction       = readSequenceFraction(frac, i);
            total_frac            = total_frac + fraction;
            pulse->sequence[i][1] = total_frac ;
            pulse->sequence[i][2] = fraction;                               // fraction
            pulse->pulse_axes[i]    = parse_axis(axis, i);                  // axis
            // degrees as given: valueint cut 90.5 to 90 (R2)
            char label[64]; snprintf(label, sizeof(label), "sequence[%d] angle", i);
            pulse->pulse_angles[i]  = readFiniteNumber(angle, label);       // angle
        }
        assign_sequence_indices(pulse);
        // the fractions make up the whole time: a sum short of 1 used to pass (R2)
        if (std::fabs(total_frac - 1.0) > 1e-9){
            ccexRefusalPrint(stderr, "Error: the fractions of the Manual sequence sum to %.12g, not 1\n", total_frac);
            ccexRefuse();
        }
    }

    // DEER bath pulses, the block pulse.bath (adding_MECCE src/option.cpp:677-730 had them as top-level keys)
    if (cJSON* bath = readBlock(block,"bath")) {
        for (cJSON* key = bath->child; key; key = key->next)
            if (strcmp(key->string,"count") && strcmp(key->string,"sequence") && strcmp(key->string,"defect") &&
                strcmp(key->string,"energy_shift") && strcmp(key->string,"spin") && strcmp(key->string,"alphams") &&
                strcmp(key->string,"betams")){
                ccexRefusalPrint(stderr,"Error: pulse.bath has no key \"%s\" (count, sequence, defect, energy_shift, spin, alphams, betams)\n",
                        key->string);
                ccexRefuse();
            }
        cJSON* jseq = cJSON_GetObjectItem(bath, "sequence");     // every key looked up, so the input check knows them
        cJSON* jdef = cJSON_GetObjectItem(bath, "defect");
        cJSON* jsh  = cJSON_GetObjectItem(bath, "energy_shift");
        cJSON* jbs  = cJSON_GetObjectItem(bath, "spin");
        cJSON* jba  = cJSON_GetObjectItem(bath, "alphams");
        cJSON* jbb  = cJSON_GetObjectItem(bath, "betams");
        pulse->bnpulse = cJSON_ReadInt(bath, (char*)"count", true, 0);    // type and integer checks as pulse.count
        if (pulse->bnpulse < 0) { ccexRefusalPrint(stderr, "Error: pulse.bath.count must be >= 0\n"); ccexRefuse(); }
        if (jbs) pulse->bspin    = readFiniteNumber(jbs, "pulse.bath.spin");
        if (jba) pulse->balphams = readFiniteNumber(jba, "pulse.bath.alphams");
        if (jbb) pulse->bbetams  = readFiniteNumber(jbb, "pulse.bath.betams");
        if (pulse->bnpulse > 0) {
            if (!(jdef && cJSON_IsString(jdef) && jdef->valuestring && jdef->valuestring[0])) {
                ccexRefusalPrint(stderr, "Error: pulse.bath.defect is required when pulse.bath.count > 0\n"); ccexRefuse();
            }
            copyBounded(pulse->bpulse_defect, jdef->valuestring, sizeof(pulse->bpulse_defect), "pulse.bath.defect");
            std::vector<double> tones;
            if (jsh && cJSON_IsNumber(jsh)) tones.push_back(readFiniteNumber(jsh, "pulse.bath.energy_shift"));
            else if (jsh && cJSON_IsArray(jsh))
                for (int k = 0; k < cJSON_GetArraySize(jsh); ++k)
                    tones.push_back(readFiniteNumber(cJSON_GetArrayItem(jsh, k), "pulse.bath.energy_shift"));
            if (tones.empty()) {
                ccexRefusalPrint(stderr, "Error: pulse.bath.energy_shift needs at least one finite number [MHz]\n"); ccexRefuse();
            }
            pulse->n_bpulse_shift = (int)tones.size();
            pulse->bpulse_energy_shifts = (double*)calloc(tones.size(), sizeof(double));
            for (size_t k = 0; k < tones.size(); ++k) pulse->bpulse_energy_shifts[k] = MHZ_TO_RADKHZ(tones[k]);
            if (!(jseq && cJSON_IsArray(jseq)) || cJSON_GetArraySize(jseq) != pulse->bnpulse + 1) {
                ccexRefusalPrint(stderr, "Error: pulse.bath.sequence must have pulse.bath.count+1 rows of [fraction, axis, angle]\n");
                ccexRefuse();
            }
            const int nrow = pulse->bnpulse + 1;
            pulse->bsequence_fracs = (double*)calloc(nrow, sizeof(double));
            pulse->bpulse_axes     = (char*)calloc(nrow + 1, sizeof(char));
            pulse->bpulse_angles   = (double*)calloc(nrow, sizeof(double));
            double total = 0.0;
            for (int i = 0; i < nrow; ++i) {
                cJSON* row = cJSON_GetArrayItem(jseq, i);
                if (!cJSON_IsArray(row) || cJSON_GetArraySize(row) != 3) {
                    ccexRefusalPrint(stderr, "Error: pulse.bath.sequence[%d] must be [fraction, axis, angle]\n", i); ccexRefuse();
                }
                pulse->bsequence_fracs[i] = readSequenceFraction(cJSON_GetArrayItem(row, 0), i);   // its messages say sequence[i]
                pulse->bpulse_axes[i]     = parse_axis(cJSON_GetArrayItem(row, 1), i);
                char label[64]; snprintf(label, sizeof(label), "pulse.bath.sequence[%d] angle", i);
                pulse->bpulse_angles[i]   = readFiniteNumber(cJSON_GetArrayItem(row, 2), label);
                total += pulse->bsequence_fracs[i];
            }
            if (std::fabs(total - 1.0) > 1e-9) {
                ccexRefusalPrint(stderr, "Error: the fractions of pulse.bath.sequence sum to %.12g, not 1\n", total); ccexRefuse();
            }
        }   // what the merge needs of the two sequences is checked where it runs (Pulse_resolveBathSequence)
    }

    //////////////////////////
    cJSON_Delete(root);
    freeChar1d(&data);
}

// savemode as a string (the mode) or as {"mode": ..., "nodiv": true|false} (2026-10-13). nodiv (default true) writes
// _noDiv beside _wiDiv: the group analyzer's error correction, the standing checks and the paper scripts read it; false
// writes the coherence (_wiDiv) alone (improvements.md C, user decision (b), 2026-10-08). With readNodiv false only the
// mode is taken (a refinement step input, whose nodiv the run keeps from the user's input).
void cJSON_readSavemode(cJSON* item, Output* op, bool readNodiv){
    if (cJSON_IsString(item) && item->valuestring){
        Output_setSavemode(op,item->valuestring);
        if (readNodiv) op->save_nodiv = true;
        return;
    }
    if (!cJSON_IsObject(item)){
        ccexRefusalPrint(stderr,"Error: savemode must be a string (coherence, coherence_tagged, cluster_coherence, cluster_density_matrix) "
                       "or an object {\"mode\": ..., \"nodiv\": true|false}\n");
        ccexRefuse();
    }
    for (cJSON* key = item->child; key; key = key->next)
        if (strcasecmp(key->string,"mode") && strcasecmp(key->string,"nodiv")){
            ccexRefusalPrint(stderr,"Error: savemode has no key \"%s\" (mode, nodiv)\n",key->string);
            ccexRefuse();
        }
    cJSON* mode = cJSON_GetObjectItem(item,"mode");
    if (mode && (!cJSON_IsString(mode) || !mode->valuestring)){
        ccexRefusalPrint(stderr,"Error: savemode.mode must be a string (coherence, coherence_tagged, cluster_coherence, cluster_density_matrix)\n");
        ccexRefuse();
    }
    Output_setSavemode(op,mode ? mode->valuestring : (char*)"coherence");
    cJSON* nodiv = cJSON_GetObjectItem(item,"nodiv");
    if (nodiv && !cJSON_IsBool(nodiv)){
        ccexRefusalPrint(stderr,"Error: savemode.nodiv must be a boolean (true or false)\n");
        ccexRefuse();
    }
    if (readNodiv) op->save_nodiv = nodiv ? (bool)nodiv->valueint : true;
}

void cJSON_readOptionOutput(Output* op, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Output Options ...");
        printMessage("  [ savemode { mode, nodiv }, outfile ] ");
    }
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }

    // save_nodiv is savemode.nodiv since 2026-10-13 (refused with the new form, as the flat gpu_* keys)
    if (cJSON_GetObjectItem(root,"save_nodiv")){
        ccexRefusalPrint(stderr,"Error: save_nodiv is now savemode.nodiv (\"savemode\": {\"mode\": \"coherence\", \"nodiv\": false})\n");
        ccexRefuse();
    }
    Output_setSavemode(op,(char*)"coherence");
    op->save_nodiv = true;
    if (cJSON* item = cJSON_GetObjectItem(root,"savemode")) cJSON_readSavemode(item,op,true);

    char* outfileDefault = NULL;
    char* outfile = cJSON_ReadString(root,"outfile",true,outfileDefault);
	if (outfile != NULL){
	    Output_allocOutfile(op);
    	Output_setOutfile(op,outfile);
	}

    cJSON_Delete(root);
    freeChar1d(&data);
}


/**
 * @brief Read the option from the input file
 * @details Read &Defect tag  options
 * @note Unit of the input file :
 * - rxyzs      : [Angstrom]
 * - hypf       : [MHz]
 * - efg        : [Hartree/Bohr^2]
 * - zfs        : [MHz]
 * - gyros      : [radkHz/G]
 * - eqs        : [10^-30 m^2]
 * - detuning   : [MHz]
*/
void cJSON_readOptionDefectArray(DefectArray* dfa, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Defect Options ...\n");
        printMessage("  - Read values of main-key 'Defect' { state_file, axis_file, species }");
        printMessage("    Defect.species[] : [ name, n_nuclear_spins, n_configurations, coordinate_frame, ( type : single value) ");
        printMessage("                axis, types, spins, gyros, quadrupole_moments,    ( type : array ) ");
        printMessage("                nuclear_positions, hyperfine, efg, zfs, detuning ]  ( type : [ axis, spname, array ] ) \n");
    }

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            ccexRefusalPrint(stdout,"Error before: %s\n", cJSON_GetErrorPtr());
        }
        ccexRefuse();
        freeChar1d(&data);
    }

    int ndefect = 0;

    // Defect is {"state_file", "axis_file", "species": [...]} (a list until A47; refuseOldInputKeys)
    cJSON* defectBlock = readBlock(root,"Defect");
    cJSON* defectArray = defectBlock ? cJSON_GetObjectItem(defectBlock,"species") : NULL;
    if (defectArray && !cJSON_IsArray(defectArray)){
        ccexRefusalPrint(stderr,"Error: Defect.species must be an array of the defect species\n");
        ccexRefuse();
    }
    cJSON* defect;

    // Read the number of defect information
    cJSON_ArrayForEach(defect, defectArray){ndefect++;}

    // // Allocate the defect array
    DefectArray_setNdefect(dfa,ndefect);
    DefectArray_allocDefect(dfa);

    // Default values
    bool apprxDefault = true; // apprx is no longer an input (A47): it changed nothing
    int naddspinDefault = 0;
    int navaaxDefault = 0;

    // Set the each defect information
    int idf = 0;
    int length = 0;
    cJSON_ArrayForEach(defect,defectArray){

        ////////////////////////////////////////////////////////////////////////
        // General information && Allocation
        ////////////////////////////////////////////////////////////////////////
        int naddspin = cJSON_ReadInt(defect,(char*)"n_nuclear_spins",true,naddspinDefault);
        int navaax = cJSON_ReadInt(defect,(char*)"n_configurations",true,navaaxDefault);
        
        DefectArray_setDefect_idf_naddspin(dfa,idf,naddspin);
        DefectArray_setDefect_idf_navaax(dfa,idf,navaax+1); // 0th : main spin
        DefectArray_allocDefect_idf(dfa,idf,navaax+1,naddspin);

        char* dfname = cJSON_ReadString(defect,(char*)"name",false,NULL);
        DefectArray_setDefect_idf_dfname(dfa,idf,dfname);

        DefectArray_setDefect_idf_apprx(dfa,idf,apprxDefault);

        bool coordinate_frame_explicit = false;
        DefectCoordinateFrame coordinate_frame =
            readDefectCoordinateFrame(defect,dfname,&coordinate_frame_explicit);
        DefectArray_setDefect_idf_coordinate_frame(dfa,idf,coordinate_frame);

        cJSON* axis_item = cJSON_GetObjectItem(defect,"axis");
        if (axis_item != NULL){
            if (!coordinate_frame_explicit){
                ccexRefusalPrint(stderr,"Error: Defect.species[%s].axis is present but coordinate_frame is missing\n",dfname);
                ccexRefusalPrint(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the axis basis is explicit.\n");
                ccexRefuse();
            }
            double axis_raw[3] = {0.0,0.0,0.0};
            double axis_normalized[3];
            char section_label[160];
            char axis_label[160];
            snprintf(section_label,sizeof(section_label),"Defect.species[%s]",dfname);
            snprintf(axis_label,sizeof(axis_label),"Defect.species[%s].axis",dfname);
            readAxis3(defect,section_label,"axis",axis_raw,true);
            normalizeAxisOrDie(axis_raw,axis_normalized,axis_label);
            DefectArray_setDefect_idf_axis(dfa,idf,axis_normalized);
        }

        ////////////////////////////////////////////////////////////////////////
        // Spin information
        ////////////////////////////////////////////////////////////////////////
        bool HaveDefault = true;
        if (naddspin >= 1){HaveDefault=false;}

        // If there is nuclear spin on defect site, 
        // then the following values should be included in your cce.json:
        // types, spins, gyros (So the "HaveDefault" value would become false)
        char** types = cJSON_ReadString1d(defect,"types",HaveDefault,NULL, naddspin);
        DefectArray_setDefect_idf_types(dfa,idf,types);
        freeChar2d(&types,naddspin);
        
        float* spins = cJSON_ReadFloat1d(defect,"spins",HaveDefault,NULL,naddspin);    
        DefectArray_setDefect_idf_spins(dfa,idf,spins);
        freeFloat1d(&spins);


        readDefectGyros(defect,dfa,idf,naddspin,HaveDefault);
        readDefectEqs(defect,dfa,idf,naddspin);

        ////////////////////////////////////////////////////////////////////////
        // Relative position && tensor information for each axis/additional spins
        ////////////////////////////////////////////////////////////////////////

        // Read nuclear_positions (rxyzs until A47) other wise set 0.0
        cJSON_ReadDefectInfo_IntCharDoubleArray(defect,(char*)"nuclear_positions",3,&(dfa->defect[idf]->rxyzs),dfa->defect[idf]->types,navaax+1,naddspin);

        // Unit-aware objects are converted to the legacy internal units while parsing.
        // Their tensor basis is explicit for the same reason as object-form zfs.
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"hyperfine",&(dfa->defect[idf]->hypf),
                                        navaax,naddspin,coordinate_frame_explicit);
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"efg",&(dfa->defect[idf]->efg),
                                        navaax,naddspin,coordinate_frame_explicit);

        // Object forms provide one physical tensor shared by every indexed
        // configuration. The legacy list remains per-index and is unchanged.
        readDefectZfs(defect,dfa,idf,navaax,coordinate_frame_explicit);

        // The object form adds an explicit frequency unit; the legacy array stays MHz.
        readDefectDetuning(defect,dfa,idf,navaax);

        // ////////////////////////////////////////////////////////////////////////
        idf++;
    }
    cJSON_Delete(root);
    freeChar1d(&data);
}

////////////////////////////////////////////////////////////////
char* cJSON_ReadFccein(char* fccein){

    FILE* inputFile = fopen(fccein, "r");
    
    if (inputFile == NULL) {
        ccexRefusalPerror("Error opening file");
        ccexRefuse();
    }

    // '#' and '!' start a comment to the end of the line, but only outside a JSON string: a key such as "#note" or a
    // path such as "run#1" used to be cut at the '#' and the input refused or misread (2026-10-08). The text is read
    // whole, so a line is no longer limited to 1023 characters either.
    std::string text;
    char chunk[4096];
    size_t got;
    while ((got = fread(chunk, 1, sizeof(chunk), inputFile)) > 0) text.append(chunk, got);
    fclose(inputFile);

    std::string kept;
    kept.reserve(text.size());
    bool inString = false, escaped = false;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (inString) {
            kept += c;
            if (escaped) escaped = false;
            else if (c == '\\') escaped = true;
            else if (c == '"') inString = false;
        } else if (c == '#' || c == '!') {
            while (i < text.size() && text[i] != '\n') ++i;   // drop the comment, keep the newline
            if (i < text.size()) kept += '\n';
        } else {
            if (c == '"') inString = true;
            kept += c;
        }
    }

    char* data = allocChar1d(kept.size() + 1);
    memcpy(data, kept.c_str(), kept.size() + 1);
    return data;
}

// A key that is present has to have the reader's type. The defaulting readers used to
// take a present key of another type for a missing one and run with the default, so
// "knight": "true" was false and "gyro": "-17608" the default gyro, with only "use
// default value" in the log (D18, 2026-10-02). Absent keys keep their defaults.
static void refuseWrongType(cJSON* item, const char* key, bool right, const char* what){
    if (item != NULL && !right){
        ccexRefusalPrint(stderr, "Error: %s must be %s in the input file\n", key, what);
        ccexRefuse();
    }
}

// An integer option given as 2.5 or 1e10 was truncated or saturated by cJSON's valueint.
static int integerValue(cJSON* item, const char* key){
    const double v = item->valuedouble;
    if (!(v >= INT_MIN && v <= INT_MAX) || std::floor(v) != v){
        ccexRefusalPrint(stderr, "Error: %s must be an integer in the input file\n", key);
        ccexRefuse();
    }
    return (int)v;
}

char* cJSON_ReadFilePath(cJSON* root, char* key, bool _default, char* default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsString(item), "a string (a file path)");

    if (cJSON_IsString(item) && (item->valuestring != NULL)) {
        if (access(item->valuestring, R_OK) == 0) {
            return item->valuestring;
        }else{
            // The path is returned anyway, so an unreadable file does not stop the run --
            // whatever the caller meant to fill from it just keeps its default. Silence
            // here is how a qubit ends up at the cell origin and the run still "works".
            if (rank==0){
                fprintf(stderr,"Warning: %s : cannot read \"%s\" -- continuing with defaults\n",
                        name, item->valuestring);
            }
            return item->valuestring;
        }
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

// An array element that is not a string has no valuestring at all, and strcpy
// from NULL was a segfault (`"bathfile": [123]`, found 2026-09-19).
static const char* stringElement(cJSON* element, const char* key, int i){
    if (!cJSON_IsString(element) || element->valuestring == NULL){
        ccexRefusalPrint(stderr, "Error: %s[%d] is not a string in the input file\n",key,i);
        ccexRefuse();
    }
    return element->valuestring;
}

char** cJSON_ReadFilePath1d(int* length, cJSON* root, char* key, bool _default, char** default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsArray(item), "an array of file paths");

    if (cJSON_IsArray(item)) {
        char** array = NULL;
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){

            ////access 
            //if (access(itemElement->valuestring, R_OK) != 0) {
            //    fprintf(stderr, "Error: %s is not found in the input file\n",key);
            //    fprintf(stderr, "Current path: %s\n",itemElement->valuestring);
            //    exit(EXIT_FAILURE);
            //}

            if (i==0){
                array = allocChar2d(1,MAX_FILEPATH);
            }else{
                reallocChar2d(&array,i,i+1,MAX_FILEPATH);
            }
            copyBounded(array[i],stringElement(itemElement,name,i),MAX_FILEPATH,name);
            i++;
        }
        *length = i;
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

char* cJSON_ReadString(cJSON* root, char* key, bool _default, char* default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsString(item), "a string");
    if (cJSON_IsString(item) && (item->valuestring != NULL)) {
        return item->valuestring;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

char** cJSON_ReadString1d(cJSON* root, char* key, bool _default, char** default_value, int size){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages

    cJSON* item = cJSON_GetObjectItem(root, key);

    refuseWrongType(item, name, cJSON_IsArray(item), "an array of strings");
    if (cJSON_IsArray(item)) {
        char** array = allocChar2d(size,MAX_CHARARRAY_LENGTH);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if (i>=size){
                ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                ccexRefuse();
            }else{
                copyBounded(array[i],stringElement(itemElement,name,i),MAX_CHARARRAY_LENGTH,name);
            }
            i++;
        }
        if (i < size){   // a short array was filled with empty strings and used (D42)
            ccexRefusalPrint(stderr, "Error: %s has %d entries in the input file, %d are needed\n",name,i,size);
            ccexRefuse();
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

int cJSON_ReadInt(cJSON* root, char* key, bool _default, int default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages

    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsNumber(item), "an integer");
    if (cJSON_IsNumber(item)) {
        return integerValue(item, name);
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

int* cJSON_ReadInt1d(cJSON* root, char* key, bool _default, int* default_value, int size){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    refuseWrongType(item, name, cJSON_IsArray(item), "an array of integers");
    if (cJSON_IsArray(item)) {
        int* array = allocInt1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                ccexRefuse();
            }else{
                if (!cJSON_IsNumber(itemElement)){
                    ccexRefusalPrint(stderr, "Error: %s[%d] must be an integer in the input file\n",name,i);
                    ccexRefuse();
                }
                array[i] = integerValue(itemElement, name);
            }
            i++;
        }
        if (i != size) {
            ccexRefusalPrint(stderr, "Error: %s is too short(%d) in the input file (expected size = %d)\n",name,i,size);
            ccexRefuse();
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

// The number readers below take a finite value only, and for a float one that fits (review fa7f39 R3): cJSON reads
// 1e999 as inf, and "bfield": 1e40 became inf when stored as a float.
static double finiteNumber(cJSON* item, const char* key, int i, int j, bool asFloat){
    const double v = item->valuedouble;
    if (!std::isfinite(v) || (asFloat && std::fabs(v) > std::numeric_limits<float>::max())){
        char where[64] = "";
        if (i >= 0) snprintf(where, sizeof(where), j >= 0 ? "[%d][%d]" : "[%d]", i, j);
        ccexRefusalPrint(stderr, "Error: %s%s must be a finite number%s in the input file\n", key, where,
                asFloat ? " within the float range" : "");
        ccexRefuse();
    }
    return v;
}

double cJSON_ReadDouble(cJSON* root, char* key, bool _default, double default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsNumber(item), "a number");
    if (cJSON_IsNumber(item)) {
        return finiteNumber(item, name, -1, -1, false);
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;    
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

double* cJSON_ReadDouble1d(cJSON* root, char* key, bool _default, double* default_value, int size){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    refuseWrongType(item, name, cJSON_IsArray(item), "an array of numbers");
    if (cJSON_IsArray(item)) {
        double* array = allocDouble1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                ccexRefuse();
            }
            // A non-number element used to be read anyway: cJSON leaves valuedouble at 0
            // for a string or a bool, so "bfield": [0,0,"500"] became [0,0,0] silently.
            if (!cJSON_IsNumber(itemElement)){
                ccexRefusalPrint(stderr, "Error: %s[%d] must be a number in the input file\n",name,i);
                ccexRefuse();
            }
            array[i] = finiteNumber(itemElement, name, i, -1, false);
            i++;
        }
        if (i != size) {
            ccexRefusalPrint(stderr, "Error: %s is too short(%d) in the input file (expected size = %d)\n",name,i,size);
            ccexRefuse();
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}


double** cJSON_ReadDouble2d(cJSON* root, char* key, bool _default, double** default_value, int row, int col){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages

    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsArray(item), "an array of arrays of numbers");
    double** array = allocDouble2d(row,col);

    if (cJSON_IsArray(item)) {
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= row){
                ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                ccexRefuse();
            }else{
                if (!cJSON_IsArray(itemElement)){
                    ccexRefusalPrint(stderr, "Error: %s[%d] must be an array of numbers in the input file\n",name,i);
                    ccexRefuse();
                }
                int j = 0;
                cJSON* itemElement_j;
                cJSON_ArrayForEach(itemElement_j, itemElement){
                    if(j >= col){
                        ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                        ccexRefuse();
                    }else{
                        if (!cJSON_IsNumber(itemElement_j)){
                            ccexRefusalPrint(stderr, "Error: %s[%d][%d] must be a number in the input file\n",name,i,j);
                            ccexRefuse();
                        }
                        array[i][j] = finiteNumber(itemElement_j, name, i, j, false);
                    }
                    j++;
                }
                if (j < col){   // a short row was filled with 0 and used (D42)
                    ccexRefusalPrint(stderr, "Error: %s[%d] has %d numbers in the input file, %d are needed\n",name,i,j,col);
                    ccexRefuse();
                }
            }
            i++;
        }
        if (i < row){
            ccexRefusalPrint(stderr, "Error: %s has %d rows in the input file, %d are needed\n",name,i,row);
            ccexRefuse();
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            copyDouble2d(array,(const double**)default_value,row,col);
            return array;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

float cJSON_ReadFloat(cJSON* root, char* key, bool _default, float default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    refuseWrongType(item, name, cJSON_IsNumber(item), "a number");
    if (cJSON_IsNumber(item)) {
        return finiteNumber(item, name, -1, -1, true);
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;    
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

float* cJSON_ReadFloat1d(cJSON* root, char* key, bool _default, float* default_value, int size){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    refuseWrongType(item, name, cJSON_IsArray(item), "an array of numbers");
    if (cJSON_IsArray(item)) {
        float* array = allocFloat1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                ccexRefusalPrint(stderr, "Error: %s is too long in the input file\n",name);
                ccexRefuse();
            }
            // A non-number element used to be read anyway: cJSON leaves valuedouble at 0
            // for a string or a bool, so "bfield": [0,0,"500"] became [0,0,0] silently.
            if (!cJSON_IsNumber(itemElement)){
                ccexRefusalPrint(stderr, "Error: %s[%d] must be a number in the input file\n",name,i);
                ccexRefuse();
            }
            array[i] = finiteNumber(itemElement, name, i, -1, true);
            i++;
        }
        if (i != size) {
            ccexRefusalPrint(stderr, "Error: %s is too short(%d) in the input file (expected size = %d)\n",name,i,size);
            ccexRefuse();
        }
        return array;
    }else{
        if(_default){
            if (rank==0 && strcasecmp(name,"bfield")!=0){printf("           %s : use default value .. \n",name);}
            return default_value;
            
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

bool cJSON_ReadBool(cJSON* root, char* key, bool _default, bool default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages
    cJSON* item = cJSON_GetObjectItem(root, key);
    // 0 and 1 are taken as false and true: 29 order-2 inputs write "addsubclus": 0, which
    // was ignored for the default true until 2026-10-02 and changes nothing at order 2
    const bool cBool = cJSON_IsNumber(item) && (item->valuedouble == 0.0 || item->valuedouble == 1.0);
    refuseWrongType(item, name, cJSON_IsBool(item) || cBool, "true or false (or 0 or 1)");
    if (cBool) {
        return item->valuedouble == 1.0;
    }
    if (cJSON_IsBool(item)) {
        return item->valueint;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",name);}
            return default_value;    
        }else{
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
            ccexRefuse();
        }
    }
}

MatrixXcd cJSON_ReadTensor(cJSON* root, char* key, bool _default, MatrixXcd default_value){
    const std::string shownName = shownKey(root,key); const char* name = shownName.c_str();   // the path, for the messages

    cJSON* tensor = cJSON_GetObjectItem(root,key);
    refuseWrongType(tensor, name, cJSON_IsArray(tensor), "a 3x3 array of numbers");
    if (!cJSON_IsArray(tensor)){
        if(_default) return default_value;
        ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",name);
        ccexRefuse();
    }
    // Exactly 3 rows of 3 finite numbers (D28, 2026-10-02): a shorter array used to be padded with
    // zeros ("tensor": [[1]] ran as diag(1,0,0)), as the object form of Defect never allowed.
    if (cJSON_GetArraySize(tensor) != 3){
        ccexRefusalPrint(stderr, "Error: %s must be a 3x3 array of numbers (it has %d rows)\n",name,cJSON_GetArraySize(tensor));
        ccexRefuse();
    }
    MatrixXcd mat = MatrixXcd::Zero(3,3);
    for (int i = 0; i < 3; ++i){
        cJSON* tensor_i = cJSON_GetArrayItem(tensor, i);
        if (!cJSON_IsArray(tensor_i) || cJSON_GetArraySize(tensor_i) != 3){
            ccexRefusalPrint(stderr, "Error: %s row %d must contain exactly 3 numbers\n",name,i);
            ccexRefuse();
        }
        for (int j = 0; j < 3; ++j){
            char label[160]; snprintf(label,sizeof(label),"%s[%d][%d]",name,i,j);
            mat(i,j) = readFiniteNumber(cJSON_GetArrayItem(tensor_i, j),label);
        }
    }
    return mat;
}

// D28 (2026-10-02): the [index, name, values] entries of the indexed Defect form are checked as the
// object form checks them. Before, a string index or name read as 0 or crashed, a fractional index was
// truncated, a string value read as 0, and a field that was not an array was skipped.
static void refuseNonArrayField(cJSON* item, const char* key){
    if (item && !cJSON_IsArray(item)){
        ccexRefusalPrint(stderr,"Error: Defect %s must be an array of [index, name, values] entries\n",key);
        ccexRefuse();
    }
}
static int indexedEntryIndex(cJSON* entry, const char* key, int i, int navaax){
    cJSON* item = cJSON_GetArrayItem(entry,0);
    if (!cJSON_IsNumber(item) || !std::isfinite(item->valuedouble) || item->valuedouble != std::floor(item->valuedouble)
        || item->valuedouble < 0 || item->valuedouble >= navaax){
        ccexRefusalPrint(stderr,"Error: %s[%d] index is out of range: it must be an integer in 0..%d\n",key,i,navaax-1);
        ccexRefuse();
    }
    return (int)item->valuedouble;
}
static const char* indexedEntryName(cJSON* entry, const char* key, int i){
    cJSON* item = cJSON_GetArrayItem(entry,1);
    if (!cJSON_IsString(item)){
        ccexRefusalPrint(stderr,"Error: %s[%d] spin name must be a string\n",key,i);
        ccexRefuse();
    }
    return item->valuestring;
}
static double indexedEntryValue(cJSON* item, const char* key, int i, int j){
    char label[160]; snprintf(label,sizeof(label),"%s[%d] value %d",key,i,j);
    return readFiniteNumber(item,label);
}

void cJSON_ReadDefectInfo_IntCharDoubleArray(cJSON* root, char* key, int valuecount, double**** array, char** types, int navaax, int naddspin){

    // itemArray2d : (int, char, doubleArray) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);
    refuseNonArrayField(itemArray2d, key);

    if (itemArray2d == NULL && (navaax == 0 && naddspin == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        for (int isp=0; isp<naddspin; isp++){
            for (int j=0; j<valuecount; j++){
                (*array)[iax][isp][j] = 0.0;
            }
        }
    }

    // set the array from the input file
    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, doubleArray
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            ccexRefusalPrint(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            ccexRefusalPrint(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            ccexRefuse();
        }

        // Find axis index.
        //
        // NOTE the parameter called navaax here is the caller's navaax + 1: index 0 is
        // reserved for the main spin, so DefectArray_allocDefect_idf allocates navaax+1
        // entries and the valid range is 0 .. navaax.
        //
        // The bound was `iax > navaax` before, which accepted one index past the end of
        // every array this reads into -- rxyzs, hypf, efg, zfs and detuning alike. An
        // input addressing configuration navaax+1 therefore ran, and corrupted memory
        // while it did. It is refused now, which is the one way this branch can stop a
        // v1.1.0 input that used to complete.
        int iax = indexedEntryIndex(itemArray1d, key, i, navaax);
        if (iax < 0 || iax >= navaax) {
            ccexRefusalPrint(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            ccexRefusalPrint(stderr, "Error: %s[%d] is out of range\n", key, i);
            ccexRefuse();
        }
        
        // Find spin index
        const char* spname = indexedEntryName(itemArray1d, key, i);
        int isp = findIndexChar(types,0,naddspin-1,(char*)spname);
        if (isp == -1) {
            ccexRefusalPrint(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",spname);
            ccexRefuse();
        }
        
        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            ccexRefusalPrint(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            ccexRefusalPrint(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            ccexRefuse();
        }
    
        for (int j = 0; j < valuecount; j++) {
            (*array)[iax][isp][j] = indexedEntryValue(cJSON_GetArrayItem(values, j), key, i, j);
        }
    }
}

void cJSON_ReadDefectInfo_IntCharMatrixXcd2d(cJSON* root, char* key, int valuecount, MatrixXcd*** array, char** types, int navaax, int naddspin){

    if (valuecount != 9){
        ccexRefusalPrint(stderr, "Error: valuecount should be 9\n");
        ccexRefuse();
    }

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);
    refuseNonArrayField(itemArray2d, key);

    if (itemArray2d == NULL && (navaax == 0 && naddspin == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        for (int isp=0; isp<naddspin; isp++){
            (*array)[iax][isp] = MatrixXcd::Zero(3,3);
        }
    }



    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            ccexRefusalPrint(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            ccexRefuse();
        }

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = indexedEntryIndex(itemArray1d, key, i, navaax);
        if (iax < 0 || iax >= navaax) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            ccexRefusalPrint(stderr, "Error: %s[%d] is out of range\n", key, i);
            ccexRefuse();
        }

        // Find spin index
        const char* spname = indexedEntryName(itemArray1d, key, i);
        int isp = findIndexChar(types,0,naddspin-1,(char*)spname);
        if (isp == -1) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            ccexRefusalPrint(stderr, "Error: %s is not found in the input file\n",spname);
            ccexRefuse();
        }

        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            ccexRefusalPrint(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            ccexRefuse();
        }

        for (int j = 0; j < valuecount; j++) {
            int row = j / 3;
            int col = j % 3;
            double value = indexedEntryValue(cJSON_GetArrayItem(values, j), key, i, j);
            (*array)[iax][isp](row,col) = doublec(value,0.0);
        }
    }
}

void cJSON_ReadDefectInfo_IntCharMatrixXcd1d(cJSON* root, char* key, int valuecount, MatrixXcd** array, int navaax){

    if (valuecount != 9){
        ccexRefusalPrint(stderr, "Error: valuecount should be 9\n");
        ccexRefuse();
    }

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);
    refuseNonArrayField(itemArray2d, key);

    if (itemArray2d == NULL && (navaax == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        (*array)[iax] = MatrixXcd::Zero(3,3);
    }

    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            ccexRefusalPrint(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            ccexRefuse();
        }

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = indexedEntryIndex(itemArray1d, key, i, navaax);
        if (iax < 0 || iax >= navaax) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            ccexRefusalPrint(stderr, "Error: %s[%d] is out of range\n", key, i);
            ccexRefuse();
        }

        // Find spin index
        const char* spname = indexedEntryName(itemArray1d, key, i);
        if (strcasecmp(spname,"e") != 0) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            ccexRefusalPrint(stderr, "Error: key : %s, %s is not found in the input file\n",key,spname);
            ccexRefuse();
        }

        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            ccexRefusalPrint(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            ccexRefuse();
        }

        for (int j = 0; j < valuecount; j++) {
            int row = j / 3;
            int col = j % 3;
            double value = indexedEntryValue(cJSON_GetArrayItem(values, j), key, i, j);
            (*array)[iax](row,col) = doublec(value,0.0);
        }
    }
}



void cJSON_ReadDefectInfo_IntCharDouble(cJSON* root, char* key, double** array, int navaax){

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);
    refuseNonArrayField(itemArray2d, key);

    if (itemArray2d == NULL && (navaax == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        (*array)[iax] = 0.0;
    }

    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            ccexRefusalPrint(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            ccexRefuse();
        }

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = indexedEntryIndex(itemArray1d, key, i, navaax);
        if (iax < 0 || iax >= navaax) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            ccexRefusalPrint(stderr, "Error: %s[%d] is out of range\n", key, i);
            ccexRefuse();
        }

        // Find spin index
        const char* spname = indexedEntryName(itemArray1d, key, i);
        if (strcasecmp(spname,"e") != 0) {
            ccexRefusalPrint(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            ccexRefusalPrint(stderr, "Error: key : %s, %s is not found in the input file\n",key,spname);
            ccexRefuse();
        }

        // Find double array values
        (*array)[iax] = indexedEntryValue(cJSON_GetArrayItem(itemArray1d, 2), key, i, 0);
    }
}

// ---------------------------------------------------------------------------------------------
// D24 (2026-10-02): a key no reader looked up is either a typo or not used by this run. A typo of a
// key that was looked up but is absent from the same object took that key's default without a word
// ("nstates" for nstate, "alphams" at the top of the qubitfile form for qalphams), so it stops the
// run; any other key nobody read is reported. Call after every cJSON_readOption* above has run. The
// "Adaptive" block (and its old name "Refine") is checked by its own reader and is skipped here.
static size_t keyDistance(const std::string& a, const std::string& b){
    std::vector<size_t> row(b.size()+1);
    for (size_t j = 0; j <= b.size(); ++j) row[j] = j;
    for (size_t i = 1; i <= a.size(); ++i){
        size_t diagonal = row[0]; row[0] = i;
        for (size_t j = 1; j <= b.size(); ++j){
            const size_t up = row[j];
            row[j] = std::min({row[j]+1, row[j-1]+1, diagonal+(a[i-1] == b[j-1] ? 0 : 1)});
            diagonal = up;
        }
    }
    return row[b.size()];
}

static bool objectHasKey(const cJSON* object, const std::string& key){
    for (const cJSON* c = object->child; c; c = c->next)
        if (c->string && lowerKey(c->string) == key) return true;
    return false;
}

static void collectUnreadKeys(const cJSON* object, const std::string& shown, const std::string& block,
                              std::vector<std::string>& typos, std::vector<std::string>& unread){
    // the keys looked up in this block: a key is read when block.key was looked up, and a typo is
    // judged only against the keys of the same block
    std::vector<std::string> blockKeys;
    for (const std::string& k : lookedUpKeys()){
        const size_t dot = k.rfind('.');
        const std::string parent = k.substr(0, dot == std::string::npos ? 0 : dot);
        const std::string name = dot == std::string::npos ? k : k.substr(dot+1);
        if (parent == block && (dot != std::string::npos || block.empty())) blockKeys.push_back(name);
    }
    for (const cJSON* c = object->child; c; c = c->next){
        if (!c->string) continue;
        const std::string key = lowerKey(c->string);
        if (block.empty() && (key == "adaptive" || key == "refine")) continue;
        const std::string where = shown.empty() ? c->string : shown + "." + c->string;
        if (!lookedUpKeys().count(joinKey(block, key))){
            std::string best; size_t bestDistance = 3;
            if (key.size() >= 4)
                for (const std::string& known : blockKeys){
                    const size_t d = keyDistance(key, known);
                    if (d < bestDistance && !objectHasKey(object, known)){ bestDistance = d; best = known; }
                }
            if (!best.empty()) typos.push_back("\"" + where + "\" is not read; did you mean \"" + best
                                               + "\"? (it is absent, so its default was used)");
            else unread.push_back(where);
            continue;
        }
        const std::string inner = joinKey(block, key);
        if (cJSON_IsObject(c)) collectUnreadKeys(c, where, inner, typos, unread);
        else if (cJSON_IsArray(c)){
            int i = 0;
            for (const cJSON* e = c->child; e; e = e->next, ++i)
                if (cJSON_IsObject(e)) collectUnreadKeys(e, where + "[" + std::to_string(i) + "]", inner + "[]", typos, unread);
        }
    }
}

void cJSON_checkInputKeys(char* fccein){
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);
    if (root == NULL){ freeChar1d(&data); return; }   // the readers above have reported a parse failure
    std::vector<std::string> typos, unread;
    collectUnreadKeys(root, "", "", typos, unread);
    cJSON_Delete(root); freeChar1d(&data);
    // on stdout, with the readers' own report: mpirun forwards stdout and stderr in chunks, so a
    // warning on stderr could land in the middle of a report line that a check reads back
    if (rank == 0){
        for (const std::string& k : unread)
            printf("Warning: input key \"%s\" is not read by this run and is ignored\n", k.c_str());
        fflush(stdout);
    }
    if (!typos.empty()){
        if (rank == 0)
            for (const std::string& t : typos) ccexRefusalPrint(stderr, "Error: input key %s\n", t.c_str());
        ccexRefuse();
    }
}
