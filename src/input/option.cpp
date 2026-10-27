/**
 * @file option.cpp
 * @brief The input file reader: parsing the JSON input, checking its keys and refusing removed ones.
 */
#include "../../include/input/mpi_seed.h"
#include "../../include/input/option.h"
#include "../../include/input/json.h"
#include "../../include/input/reader.h"
#include "../../include/common/memory.h"
#include "../../include/calculation/superoperator.h"
#include "../../include/auto_convergence/rdip_ladder.h"
#include "../../include/calculation/bath_drive.h"

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
static void applyInputUnits(cJSON* root);
static cJSON* trackedParse(const char* text){
    cJSON* root = cJSON_Parse(text);
    if (root) applyInputUnits(root);    // A67: {"value", "unit"} to the plain number in the key's own unit, before any reader
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
static int integerValue(cJSON* item, const char* key);
static void refuseUnreadable(const std::string& key, const std::string& file);
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

////////////////////////////////////////////////////////////////////////////////////////////////
// A67 (2026-10-15, user: a unit per key, not a units block): a key below takes a plain number (or list) in its own unit,
// as before, or {"value": number or list, "unit": "..."}. The object is replaced by the plain value in the key's own unit
// right after the parse, so every reader and the unread-key check see what they always saw. The conversion is one exact
// operation: a multiplication by an exact power of ten going to a smaller unit, a division by one going to a larger
// (1/1000 has no exact binary form; x/1000 is the double nearest the decimal x/1000, so "5" us is 0.005 ms bit for bit).
// The key's own unit, written out, leaves the value untouched.
// the keys the input wrote with a unit object (shown path), for the note on a plain number in a unit that is easy to
// mistake (A94, R4: Qubit.qubit[].detuning is kHz, Defect.species[].detuning MHz)
static std::set<std::string>& unitWritten(){ static std::set<std::string> keys; return keys; }
namespace {
struct UnitFactor { const char* unit; double times; double over; };
struct UnitKind { const char* what; std::vector<UnitFactor> units; };
const UnitKind kTimeMs   {"ms",  {{"s",1e3,1}, {"ms",1,1}, {"us",1,1e3}, {"ns",1,1e6}}};
const UnitKind kTimeNs   {"ns",  {{"s",1e9,1}, {"ms",1e6,1}, {"us",1e3,1}, {"ns",1,1}}};
const UnitKind kFieldG   {"G",   {{"G",1,1}, {"mT",10,1}, {"T",1e4,1}}};
const UnitKind kLengthA  {"A",   {{"A",1,1}, {"nm",10,1}}};
const UnitKind kFreqMHz  {"MHz", {{"Hz",1,1e6}, {"kHz",1,1e3}, {"MHz",1,1}, {"GHz",1e3,1}}};
const UnitKind kFreqKHz  {"kHz", {{"Hz",1,1e3}, {"kHz",1,1}, {"MHz",1e3,1}, {"GHz",1e6,1}}};
const UnitKind kRateRadMs{"rad/ms", {{"rad/s",1,1e3}, {"rad/ms",1,1}, {"rad/us",1e3,1}}};   // A94 (G20): angular frequency

std::string allowedUnits(const UnitKind& kind){
    std::string out;
    for (const UnitFactor& f : kind.units) out += (out.empty() ? "" : ", ") + std::string(f.unit);
    return out;
}
// every leaf of a number or a (nested) list of numbers
bool numbersOnly(const cJSON* v){
    if (cJSON_IsNumber(v)) return true;
    if (!cJSON_IsArray(v) || !v->child) return false;
    for (const cJSON* c = v->child; c; c = c->next) if (!numbersOnly(c)) return false;
    return true;
}
void scaleNumbers(cJSON* v, const UnitFactor& f){
    if (cJSON_IsNumber(v)){
        if (f.times != 1) cJSON_SetNumberHelper(v, v->valuedouble * f.times);
        else if (f.over != 1) cJSON_SetNumberHelper(v, v->valuedouble / f.over);
        return;
    }
    for (cJSON* c = v->child; c; c = c->next) scaleNumbers(c, f);
}
// parent.key (looked up without regard to case, as the readers do): an object there becomes its value in the key's unit
void applyUnit(cJSON* parent, const char* key, const std::string& shown, const UnitKind& kind){
    cJSON* item = untrackedItem(parent, key);
    if (!cJSON_IsObject(item)) return;                  // a plain value, or none: the reader's as before
    cJSON* value = NULL; cJSON* unit = NULL;
    int n = 0;
    for (cJSON* c = item->child; c; c = c->next, ++n){
        if (!strcasecmp(c->string,"value")) value = c;
        else if (!strcasecmp(c->string,"unit")) unit = c;
    }
    const std::string allowed = allowedUnits(kind);
    if (n != 2 || !value || !unit || !cJSON_IsString(unit) || !numbersOnly(value)){
        ccexRefusalPrint(stderr,"Error: %s must be a number (in %s) or {\"value\": number or list of numbers, \"unit\": \"...\"} "
                         "with the unit one of %s\n", shown.c_str(), kind.what, allowed.c_str());
        ccexRefuse();
    }
    const UnitFactor* factor = NULL;
    for (const UnitFactor& f : kind.units) if (!strcmp(f.unit, unit->valuestring)) factor = &f;
    if (!factor){
        ccexRefusalPrint(stderr,"Error: %s.unit \"%s\" is not known, the unit is one of %s (a plain number is in %s)\n",
                         shown.c_str(), unit->valuestring, allowed.c_str(), kind.what);
        ccexRefuse();
    }
    unitWritten().insert(shown);
    cJSON* plain = cJSON_Duplicate(value, 1);
    scaleNumbers(plain, *factor);
    const std::string name = item->string;              // the key as written: the replacement keeps it
    cJSON_ReplaceItemInObjectCaseSensitive(parent, name.c_str(), plain);
}
}   // namespace

// The keys with a unit (A67). Not here: what already takes {"values", "unit"} (the Defect block, Qubit.intmap), the
// numbers inside the bath and state files, the gyromagnetic ratios; auto_convergence.time_window.value is the runner's
// (example/tools/AutoConvergence).
static void applyInputUnits(cJSON* root){
    if (!cJSON_IsObject(root)) return;
    applyUnit(root, "deltat", "deltat", kTimeMs);
    applyUnit(root, "bfield", "bfield", kFieldG);
    applyUnit(root, "rbath", "rbath", kLengthA);
    applyUnit(root, "rdip", "rdip", kLengthA);
    cJSON* conf = untrackedItem(root, "configuration");
    if (cJSON_IsObject(conf)) applyUnit(conf, "bath_shift", "configuration.bath_shift", kLengthA);
    cJSON* qubitBlock = untrackedItem(root, "Qubit");
    cJSON* qubits = cJSON_IsObject(qubitBlock) ? untrackedItem(qubitBlock, "qubit") : NULL;
    if (cJSON_IsArray(qubits)){
        int i = 0;
        for (cJSON* q = qubits->child; q; q = q->next, ++i){
            if (!cJSON_IsObject(q)) continue;
            const std::string at = "Qubit.qubit[" + std::to_string(i) + "].";
            applyUnit(q, "xyz", at + "xyz", kLengthA);
            applyUnit(q, "detuning", at + "detuning", kFreqKHz);
        }
    }
    cJSON* pulse = untrackedItem(root, "pulse");
    if (cJSON_IsObject(pulse)){
        applyUnit(pulse, "duration", "pulse.duration", kTimeNs);
        cJSON* bath = untrackedItem(pulse, "bath");
        if (cJSON_IsObject(bath)){
            applyUnit(bath, "energy_shift", "pulse.bath.energy_shift", kFreqMHz);
            applyUnit(bath, "bandwidth", "pulse.bath.bandwidth", kFreqMHz);
        }
    }
    // A94 (G20): cluster.bundle.merge_cutoff (rad/ms) and hamiltonian.common.quadrupole.cell (A) take a unit too
    cJSON* bundle = untrackedItem(untrackedItem(root, "cluster"), "bundle");
    if (cJSON_IsObject(bundle)) applyUnit(bundle, "merge_cutoff", "cluster.bundle.merge_cutoff", kRateRadMs);
    cJSON* quadrupole = untrackedItem(untrackedItem(untrackedItem(root, "hamiltonian"), "common"), "quadrupole");
    if (cJSON_IsObject(quadrupole)) applyUnit(quadrupole, "cell", "hamiltonian.common.quadrupole.cell", kLengthA);
    cJSON* hamiltonian = untrackedItem(root, "hamiltonian");
    if (cJSON_IsObject(hamiltonian))
        for (cJSON* section = hamiltonian->child; section; section = section->next){
            cJSON* jumps = cJSON_IsObject(section) ? untrackedItem(section, "jump_operators") : NULL;
            if (!cJSON_IsArray(jumps)) continue;
            int i = 0;
            for (cJSON* j = jumps->child; j; j = j->next, ++i)
                if (cJSON_IsObject(j))
                    applyUnit(j, "rate", "hamiltonian." + std::string(section->string) + ".jump_operators[" + std::to_string(i) + "].rate",
                              kFreqMHz);
        }
}

static void refuseOldInputKeys(cJSON* root){
    struct Renamed { const char* block; const char* key; const char* now; };
    static const Renamed renamed[] = {
        {"", "qubitfile", "is now an entry of Qubit.qubit: \"Qubit\": {\"qubit\": [{\"file\": ...}]}"},
        {"", "qspin", "is now the \"spin\" of the Qubit.qubit entry: [{\"file\": ..., \"spin\": ...}]"},
        {"", "qzfs", "is now a Qubit.intmap self entry: \"intmap\": [{\"between\": [\"q0\", \"q0\"], "
                     "\"tensor\": ..., \"unit\": \"kHz\"}]"},
        {"", "qalphams", "is now \"Qubit.alpha_state\": {\"ms\": [m]}"},
        {"", "qbetams", "is now \"Qubit.beta_state\": {\"ms\": [m]}"},
        {"", "rbathcut", "is now the inner radius of \"rbath\": [inner, outer]"},
        {"", "rdipcut", "is now the inner radius of \"rdip\": [inner, outer]"},
        {"", "bathadjust", "is now \"configuration.bath_shift\""},
        {"", "statefile", "is now \"configuration.bath_state_file\": [one name per bath file]"},
        // A64: the bath configuration's files are one block, "configuration" {path, k_range, gyrofile, bathfile, bath_shift,
        // bath_state_file}
        {"", "gyrofile", "is now \"configuration.gyrofile\" (\"configuration\": {\"gyrofile\": ...})"},
        {"", "bathfile", "is now \"configuration.bathfile\" (\"configuration\": {\"bathfile\": [...]})"},
        {"", "bath_shift", "is now \"configuration.bath_shift\""},
        {"", "bath_state_file", "is now \"configuration.bath_state_file\": [one name per bath file]"},
        {"", "nstate", "is now \"configuration.s_range\": [1, N] (the sampled states 1 .. N), leave it out for the ensemble"},
        // A94 (2026-10-16): the sampled states are a range, as the configurations are
        {"configuration", "nstate", "is now \"configuration.s_range\": [first, last] (nstate N: [1, N]; nstate 0, the "
                                    "ensemble: leave s_range out)"},
        {"", "seed", "is now \"configuration.seed\""},
        {"", "exstatefile", "is now \"configuration.defect_state_file\": [one name per bath file]"},
        {"", "avaaxfile", "is now \"configuration.defect_axis_file\": [one file per bath file]"},
        // A64: the Defect block holds the species (the model) alone; the files of a configuration are the configuration's
        {"Defect", "axis_file", "is now \"configuration.defect_axis_file\": [one file per bath file]"},
        {"Defect", "state_file", "is now \"configuration.defect_state_file\": [one name per bath file]"},
        {"", "addsubclus", "is now \"cluster.selection.add_subclusters\" (true or false)"},
        {"", "nk", "is now \"cluster.selection.max_per_order\""},
        {"", "vector_mean_field", "is now \"hamiltonian.common.bath_bath_mean_field\" (\"xyz\" for true, \"z\" for false)"},
        {"", "hf_readmode", "is now \"hamiltonian.common.hyperfine\" {\"dipole\": \"point\"|\"dft\", \"contact\": true|false} "
                            "(0 point/false, 1 point/true, 2 dft/false, 3 dft/true)"},
        {"", "hf_tensorfile", "is now \"hamiltonian.common.hyperfine.file\""},
        {"", "hf_ignore_oor", "is now \"hamiltonian.common.hyperfine.missing_row\" (0 \"stop\", 1 \"point_dipole\")"},
        {"", "qd_readmode", "is now \"hamiltonian.common.quadrupole\" (0 \"off\", 2, 3, 4 {\"structure\": \"plain\", "
                            "\"strained\", \"bubble\"})"},
        {"", "qd_tensorfile", "is now \"hamiltonian.common.quadrupole.file\""},
        {"", "qd_tensorfile_woqubit", "is now \"hamiltonian.common.quadrupole.file_without_defect\""},
        {"", "qd_cellpara", "is now \"hamiltonian.common.quadrupole.cell\""},
        {"", "knight", "is now \"hamiltonian.cCCE.knight\""},
        {"", "hfmedi", "is now \"hamiltonian.cCCE.hfmedi\""},
        {"", "DefectTotSpin", "was removed: the qubit spin is used (also for the S of a v3 tensor file)"},
        {"", "CorrTotSpin", "was removed: it served old-format VASP tensor files only, remake such a file with "
                            "example/tools/SpinTensorReader (format v3)"},
        {"", "sK", "is now \"cluster.pcce.partition_size\""},
        {"", "max_trial", "is now \"cluster.pcce.trials\""},
        {"", "max_iter", "is now \"cluster.pcce.max_iterations\""},
        {"", "kmeans_pp", "is now \"cluster.pcce.kmeans_pp\""},
        {"", "iter_detail", "is now \"cluster.pcce.print_trials\" (the number of trials printed, 0 none)"},
        {"", "pcce_detail_limit", "is now \"cluster.pcce.print_trials\" (the number of trials printed)"},
        {"", "pcce_partition", "is now \"cluster.pcce.partition\""},
        {"", "pcce_connections", "is now \"cluster.pcce.connections\""},
        {"", "pcce_cache_mode", "is now \"cluster.pcce.cache.mode\""},
        {"", "pcce_cache_file", "is now \"cluster.pcce.cache.file\""},
        {"", "bundle_max_dim", "is now \"cluster.bundle.max_dim\""},
        {"", "bundle_merge_cutoff", "is now \"cluster.bundle.merge_cutoff\""},
        {"", "bundle_connections", "is now \"cluster.bundle.connections\""},
        {"", "bundle_connection_budget", "is now \"cluster.bundle.connection_budget\""},
        {"", "expansion_max_clusters", "is now \"cluster.bundle.max_clusters\""},
        {"", "expansion_max_work", "is now \"cluster.bundle.max_work\""},
        {"", "adaptive_work_file", "is now \"cluster.adaptive.work_file\""},
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
        {"", "gcce_refinement_max_sweeps", "is now \"execution.gcce_refinement.max_sweeps\""},
        {"", "cluster_contribution_max_mb", "is now \"output.cluster_contributions.max_mb\""},
        // A65 (2026-10-15): how it is computed in "execution", the clusters in "cluster", what is written in "output",
        // the frame rotation in "frame"
        {"", "backend", "is now \"execution.backend\" (\"execution\": {\"backend\": ...})"},
        {"", "gpu_options", "is now \"execution.gpu_options\""},
        {"", "propagator", "is now \"execution.propagator\""},
        {"", "evolution", "is now \"execution.evolution\""},
        {"", "evolution_update", "is now \"execution.evolution_update\""},
        {"", "gcce_refinement", "is now \"execution.gcce_refinement\""},
        {"", "column_evolution", "is now \"execution.column_evolution\""},
        {"", "profile", "is now \"execution.profile\""},
        {"", "cluster_mode", "is now \"cluster.mode\" (\"cluster\": {\"mode\": ...})"},
        {"", "add_subclusters", "is now \"cluster.selection.add_subclusters\""},
        {"", "max_clusters_per_order", "is now \"cluster.selection.max_per_order\""},
        {"", "pcce", "is now \"cluster.pcce\""},
        {"", "bundle", "is now \"cluster.bundle\""},
        {"", "adaptive", "is now \"cluster.adaptive\""},
        {"", "probe_memory_gb", "is now \"cluster.adaptive.probe_memory_gb\""},
        {"", "outfile", "is now \"output.file\" (\"output\": {\"file\": ...})"},
        {"", "savemode", "is now \"output.savemode\""},
        {"", "cluster_contributions", "is now \"output.cluster_contributions\""},
        {"", "cluster_only", "is now \"output.cluster_only\""},
        {"", "coordinate_frame_rotation", "is now \"frame.rotation\" (on when the block is given, \"enabled\" defaults to true)"},
        {"", "bath_coordinate_rotation", "is now \"frame.rotation\" (on when the block is given, \"enabled\" defaults to true)"},
        {"", "defect_axis_reference", "was removed: it only named the qubit axis of frame.rotation in the log, remove it"},
        {"cluster", "max_clusters_per_order", "is now \"cluster.selection.max_per_order\""},
        {"cluster", "add_subclusters", "is now \"cluster.selection.add_subclusters\""},
        {"cluster", "cluster_mode", "is now \"cluster.mode\""},
        {"Qubit", "nqubit", "was removed: the number of qubits is the length of Qubit.qubit"},
        {"Qubit", "overhaus", "is now \"hamiltonian.common.qubit_bath_mean_field\""},
        {"Qubit", "alphaidx", "is now \"Qubit.alpha_state\": {\"index\": n}"},
        {"Qubit", "betaidx", "is now \"Qubit.beta_state\": {\"index\": n}"},
        {"Qubit", "psia", "is now \"Qubit.alpha_state\": {\"vector\": [...]}"},
        {"Qubit", "psib", "is now \"Qubit.beta_state\": {\"vector\": [...]}"},
        {"Qubit", "zfs", "is now a Qubit.intmap self entry: \"intmap\": [{\"between\": [\"q0\", \"q0\"], "
                         "\"tensor\": ..., \"unit\": \"kHz\"}] (\"q0\" is the qubit's name, its default), with \"frame\" there "
                         "when frame.rotation is on (\"qubit\": as Qubit.zfs was read)"},                                   // A93
        {"Qubit", "frame", "was removed: the frame of the ZFS is the \"frame\" of its Qubit.intmap self entry"},
        {"Qubit", "alphams", "is now \"Qubit.alpha_state\": {\"ms\": [m]}"},        // A93
        {"Qubit", "betams", "is now \"Qubit.beta_state\": {\"ms\": [m]}"},
        // A94 (G12): the Qubit.file form is an entry of Qubit.qubit, the position read from the file in place of "xyz"
        {"Qubit", "file", "is now an entry of Qubit.qubit: \"qubit\": [{\"file\": ..., \"spin\": ...}] (\"file\" in place of "
                          "\"xyz\", with the entry's name, spin, gyro, detuning)"},
        {"Qubit", "spin", "is now the \"spin\" of the Qubit.qubit entry: \"qubit\": [{\"file\": ..., \"spin\": ...}]"},
        {"Qubit", "name", "is now the \"name\" of the Qubit.qubit entry: \"qubit\": [{\"file\": ..., \"name\": ...}]"},
    };
    static const Renamed species[] = {
        {"", "dfname", "is now \"name\""},
        {"", "apprx", "was removed (it changed nothing)"},
        {"", "naddspin", "is now \"n_onsite_spins\""},
        {"", "navaax", "is now \"n_axis_states\""},
        {"", "n_nuclear_spins", "is now \"n_onsite_spins\""},      // A65
        {"", "n_configurations", "is now \"n_axis_states\""},
        {"", "eqs", "is now \"quadrupole_moments\""},
        {"", "rxyzs", "is now \"nuclear_positions\""},
        {"", "hypf", "is now \"hyperfine\""},
        {"", "coordinate_frame", "is now \"frame\" (\"bath\" is now \"file\")"},     // A90
        {"", "axis", "is now \"zfs.axis\": \"zfs\": {\"D\": ..., \"axis\": [...], \"unit\": ...} (only the D form of zfs used "
                     "it)"},                                                                // A94 (G5)
    };
    // A90 (2026-10-16): frame.rotation keeps qubit_axis and enabled; the frame of a set of directional numbers is the key
    // "frame" beside them ("file" | "qubit")
    static const Renamed rotation[] = {
        {"", "bath_axis", "was removed: it only turned the computational x and y axes about qubit_axis (x is now "
                          "qubit_axis x [0,0,1], the old default)"},
        {"", "reference_qubit", "was removed: it only chose the rotation origin, a translation no result depends on (the "
                                "origin is the first qubit)"},
        {"", "qubit_position_frame", "was removed: its one value was \"bath\" (Qubit.qubit[].xyz is always in the bath "
                                     "file's frame)"},
        {"", "hf_tensor_frame", "is now \"hamiltonian.common.hyperfine.frame\" (\"bath\" is now \"file\")"},
        {"", "qd_tensor_frame", "is now \"hamiltonian.common.quadrupole.frame\" (\"bath\" is now \"file\")"},
    };
    std::vector<std::string> found;
    for (const Renamed& r : renamed){
        cJSON* block = *r.block ? untrackedItem(root,r.block) : root;
        if (untrackedItem(block,r.key))
            found.push_back(std::string("input key \"") + (*r.block ? std::string(r.block) + "." : "") + r.key + "\" " + r.now);
    }
    if (cJSON* rot = untrackedItem(untrackedItem(root,"frame"),"rotation"))
        for (const Renamed& r : rotation)
            if (untrackedItem(rot,r.key))
                found.push_back(std::string("input key \"frame.rotation.") + r.key + "\" " + r.now);
    // A93 (2026-10-16): the qubit states are Qubit.alpha_state/beta_state only
    if (cJSON* qubits = untrackedItem(untrackedItem(root,"Qubit"),"qubit")){
        int i = 0;
        for (cJSON* entry = cJSON_IsArray(qubits) ? qubits->child : NULL; entry; entry = entry->next, ++i)
            for (const char* side : {"alpha", "beta"})
                if (untrackedItem(entry,(std::string(side) + "ms").c_str()))
                    found.push_back("input key \"Qubit.qubit[" + std::to_string(i) + "]." + side + "ms\" is now \"Qubit." + side
                                    + "_state\": {\"ms\": [m_1, m_2, ...]}, one m_s per qubit in the order of Qubit.qubit");
    }
    static const Renamed intmapEntry[] = {
        {"", "tensor_frame", "is now \"frame\" (\"bath\" is now \"file\")"},
    };
    if (cJSON* intmap = untrackedItem(untrackedItem(root,"Qubit"),"intmap")){
        int i = 0;
        for (cJSON* entry = cJSON_IsArray(intmap) ? intmap->child : NULL; entry; entry = entry->next, ++i)
            for (const Renamed& r : intmapEntry)
                if (untrackedItem(entry,r.key))
                    found.push_back("input key \"Qubit.intmap[" + std::to_string(i) + "]." + r.key + "\" " + r.now);
    }
    if (cJSON* method = untrackedItem(root,"method")){
        const char* v = cJSON_IsString(method) ? method->valuestring : "";
        if (!strcasecmp(v,"cce")) found.push_back(std::string("method \"") + v + "\" is now \"cCCE\"");
        for (const char* shape : {"pcce", "dsj", "itb", "dsjitb"})
            if (!strcasecmp(v,shape))
                found.push_back(std::string("method \"") + v + "\" is now \"method\": \"cCCE\" with \"cluster\": {\"mode\": \""
                                + shape + "\"}");
    }
    if (cJSON* mode = untrackedItem(untrackedItem(root,"cluster"),"mode"))
        if (cJSON_IsString(mode) && !strcasecmp(mode->valuestring,"adaptive"))
            found.push_back("cluster.mode \"adaptive\" is now \"bundle\" (the bundle partition, the cluster.adaptive block "
                            "asks for the refinement)");
    if (cJSON* item = untrackedItem(untrackedItem(root,"execution"),"gcce_refinement"))
        if (cJSON_IsBool(item))
            found.push_back("\"execution.gcce_refinement\": true|false is now \"execution.gcce_refinement\": {\"enabled\": true|false}");
    if (cJSON* item = untrackedItem(untrackedItem(root,"output"),"cluster_contributions"))
        if (cJSON_IsBool(item))
            found.push_back("\"output.cluster_contributions\": true|false is now \"output.cluster_contributions\": {\"enabled\": true|false}");
    if (cJSON* defect = untrackedItem(root,"Defect")){
        cJSON* list = defect;
        if (cJSON_IsArray(defect))
            found.push_back("\"Defect\": [...] is now \"Defect\": {\"species\": [...]} (with \"defect_state_file\" and \"defect_axis_file\" in configuration)");
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

////////////////////////////////////////////////////////////////////////////////////////////////
// A64 (2026-10-15): the bath configuration is one block, "configuration": {"path", "k_range", "gyrofile", "bathfile",
// "bath_shift", "bath_state_file", "defect_axis_file", "defect_state_file", "nstate", "seed"}. path is put in front of the
// block's relative file names (an absolute name is used as it is). bath_state_file, defect_axis_file and defect_state_file
// (Defect.axis_file and Defect.state_file until 2026-10-15) are lists of one name per bath file, read
// in the order of the bath files. With "k_range": [first, last] every k from first to last is a configuration of its own,
// run one after another by main.cpp, and "{k}" in the entries of configuration.bathfile, configuration.bath_state_file,
// configuration.defect_axis_file, configuration.defect_state_file and in pcce.cache.file is k ("{{k}}" is a literal {k}); without k_range nothing is
// substituted. cJSON_readConfigurationRange checks k_range and every file of every k before anything is calculated;
// cJSON_readOptionConfig reads the configuration k set by cJSON_setConfigurationIndex.
static bool configurationRanged = false;
static long configurationK = 0;
void cJSON_setConfigurationIndex(bool ranged, long k){ configurationRanged = ranged; configurationK = k; }

// name with {k} as k and {{k}} as {k} (with k_range); *found: a {k} was there
static std::string configurationName(const std::string& name, bool ranged, long k, bool* found = nullptr){
    if (found) *found = false;
    if (!ranged) return name;
    std::string out;
    for (size_t i = 0; i < name.size();){
        if (!name.compare(i,5,"{{k}}")){ out += "{k}"; i += 5; }
        else if (!name.compare(i,3,"{k}")){ out += std::to_string(k); i += 3; if (found) *found = true; }
        else out += name[i++];
    }
    return out;
}
// configuration.path: "" for the working directory (absent or ".")
static std::string configurationFolder(cJSON* block){
    cJSON* item = block ? cJSON_GetObjectItem(block,"path") : NULL;
    if (!item) return "";
    if (!cJSON_IsString(item) || !item->valuestring || !item->valuestring[0] || strlen(item->valuestring) >= MAX_FILEPATH/2){
        ccexRefusalPrint(stderr,"Error: configuration.path must be a nonempty folder name (relative to the working directory, "
                       "or absolute)\n");
        ccexRefuse();
    }
    return strcmp(item->valuestring,".") ? item->valuestring : "";
}
static std::string inFolder(const std::string& folder, const std::string& name){
    if (folder.empty() || name.empty() || name[0] == '/') return name;
    return folder.back() == '/' ? folder + name : folder + "/" + name;
}
// a file name of the configuration block as it is read: {k} substituted (with k_range), in folder
static std::string configurationPath(const std::string& folder, const std::string& name, bool ranged, long k){
    return inFolder(folder,configurationName(name,ranged,k));
}
// A list of file names, one per bath file (configuration.bath_state_file, defect_axis_file, defect_state_file): {k}
// substituted, in folder; empty when the key is absent. "" is a bath file without such a file (kept as ""). A string, or
// a list of another length, is refused with the form.
static std::vector<std::string> perBathFileNames(cJSON* block, const char* key, int nbathfiles, const std::string& folder,
                                                 bool ranged, long k){
    std::vector<std::string> out;
    if (!block) return out;
    const std::string shown = shownKey(block,key);
    cJSON* item = cJSON_GetObjectItem(block,key);
    if (!item){
        if (rank==0 && verbosity){printf("        %s : use default value ..\n",shown.c_str());}
        return out;
    }
    if (!cJSON_IsArray(item) || cJSON_GetArraySize(item) != nbathfiles){
        ccexRefusalPrint(stderr,"Error: %s must be a list of one name per bath file (%d, as configuration.bathfile), e.g. "
                       "[\"name\"] for one bath file\n",shown.c_str(),nbathfiles);
        ccexRefuse();
    }
    int i = 0;
    for (cJSON* e = item->child; e; e = e->next, ++i){
        if (!cJSON_IsString(e) || !e->valuestring || strlen(e->valuestring) >= MAX_FILEPATH/2){
            ccexRefusalPrint(stderr,"Error: %s[%d] is not a string (a file name)\n",shown.c_str(),i);
            ccexRefuse();
        }
        out.push_back(configurationPath(folder,e->valuestring,ranged,k));
    }
    return out;
}
// A file name of the configuration block: {k} substituted, in folder; false when absent. The checks and messages of
// cJSON_ReadFilePath, on the name that is read.
static bool configurationFile(cJSON* block, const char* key, const std::string& folder, bool ranged, long k, std::string* out){
    if (!block) return false;
    const std::string shownName = shownKey(block,key);
    cJSON* item = cJSON_GetObjectItem(block,key);
    refuseWrongType(item, shownName.c_str(), cJSON_IsString(item), "a string (a file path)");
    if (!item){
        if (rank==0 && verbosity){printf("        %s : use default value ..\n",shownName.c_str());}
        return false;
    }
    *out = configurationPath(folder,item->valuestring,ranged,k);
    if (!out->empty()) refuseUnreadable(shownName,*out);   // A94: refused, not "continuing with defaults"; "" is no file
    return true;
}
// configuration.k_range or s_range: [first, last], two integers >= lowest, first <= last; false when the key is absent
static bool readIndexRange(cJSON* block, const char* key, long lowest, long* first, long* last){
    cJSON* range = block ? cJSON_GetObjectItem(block,key) : NULL;
    if (!range) return false;
    long value[2] = {0, 0};
    bool ok = cJSON_IsArray(range) && cJSON_GetArraySize(range) == 2;
    for (int i = 0; ok && i < 2; ++i){
        cJSON* item = cJSON_GetArrayItem(range,i);
        ok = cJSON_IsNumber(item) && std::isfinite(item->valuedouble) && item->valuedouble >= lowest &&
             item->valuedouble <= INT_MAX && std::floor(item->valuedouble) == item->valuedouble;
        if (ok) value[i] = (long)item->valuedouble;
    }
    if (!ok || value[0] > value[1]){
        if (lowest == 0)
            ccexRefusalPrint(stderr,"Error: configuration.%s must be [first, last], two nonnegative integers with first <= last "
                           "(both included)\n",key);
        else
            ccexRefusalPrint(stderr,"Error: configuration.%s must be [first, last], two integers of %ld or more with first <= "
                           "last (both included)%s\n",key,lowest,cJSON_IsNumber(range) ? ", a list also for one ([n, n])" : "");
        ccexRefuse();
    }
    *first = value[0];
    *last = value[1];
    return true;
}
// configuration.k_range: [first, last], two nonnegative integers, first <= last; false without k_range
static bool readConfigurationRange(cJSON* block, long* first, long* last){
    return readIndexRange(block,"k_range",0,first,last);
}
// configuration.s_range (A94; nstate before): [first, last] of the sampled bath states, 1 <= first <= last; false (the
// ensemble) without s_range
static bool readStateRange(cJSON* block, long* first, long* last){
    return readIndexRange(block,"s_range",1,first,last);
}
// the input has Defect species (their axes and on-site states are read from the configuration's files)
static bool inputHasDefectSpecies(cJSON* root){
    cJSON* species = untrackedItem(untrackedItem(root,"Defect"),"species");
    return cJSON_IsArray(species) && cJSON_GetArraySize(species) > 0;
}
// A94: a file the input names and the run reads has to open; there is no fallback any more (a missing state or axis file
// was drawn from the seed, a missing tensor file was a warning here and an error later). "" is no file (A64): passed.
static void refuseUnreadable(const std::string& key, const std::string& file){
    if (file.empty() || access(file.c_str(),R_OK) == 0) return;
    ccexRefusalPrint(stderr,"Error: %s: cannot read \"%s\" (%s)\n",key.c_str(),file.c_str(),strerror(errno));
    ccexRefuse();
}
// the pcce partition cache (A64): its file takes {k} like the configuration's files; NULL when there is no enabled cache
static cJSON* enabledPcceCacheFile(cJSON* pcce){
    cJSON* cache = cJSON_IsObject(pcce) ? untrackedItem(pcce,"cache") : NULL;
    cJSON* mode = cJSON_IsObject(cache) ? untrackedItem(cache,"mode") : NULL;
    if (!cJSON_IsString(mode) || !strcasecmp(mode->valuestring,"off")) return NULL;
    cJSON* file = untrackedItem(cache,"file");
    return cJSON_IsString(file) ? file : NULL;
}

bool cJSON_readConfigurationRange(char* fccein, long* first, long* last){
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);
    freeChar1d(&data);
    if (!root) return false;                    // the readers report it
    cJSON* block = untrackedItem(root,"configuration");
    if (!cJSON_IsObject(block) || !readConfigurationRange(block,first,last)){ cJSON_Delete(root); return false; }
    const std::string folder = configurationFolder(block);
    // the names {k} can be in, with the folder each is read from and whether it is a prefix of state files
    struct Named { std::string key, name, folder; bool prefix; };
    std::vector<Named> names;
    const auto add = [&](cJSON* list, const std::string& key, const std::string& from, bool prefix){
        int i = 0;
        for (cJSON* e = cJSON_IsArray(list) ? list->child : NULL; e; e = e->next, ++i)
            if (cJSON_IsString(e) && e->valuestring[0])   // "": that bath file has none
                names.push_back({key + "[" + std::to_string(i) + "]", e->valuestring, from, prefix});
    };
    add(untrackedItem(block,"bathfile"),"configuration.bathfile",folder,false);
    add(untrackedItem(block,"bath_state_file"),"configuration.bath_state_file",folder,true);
    add(untrackedItem(block,"defect_axis_file"),"configuration.defect_axis_file",folder,false);
    add(untrackedItem(block,"defect_state_file"),"configuration.defect_state_file",folder,true);
    bool any = false;
    for (const Named& n : names){ bool found = false; configurationName(n.name,true,0,&found); any = any || found; }
    if (!any){
        ccexRefusalPrint(stderr,"Error: configuration.k_range is given, but no file name has {k} (configuration.bathfile, "
                       "configuration.bath_state_file, defect_axis_file, defect_state_file), without {k} every k would be the "
                       "same configuration\n");
        ccexRefuse();
    }
    // an enabled pcce partition cache: one file per configuration, or the second would read the first one's partition
    if (cJSON* cacheFile = enabledPcceCacheFile(untrackedItem(untrackedItem(root,"cluster"),"pcce"))){   // cluster.pcce (A65)
        bool found = false;
        configurationName(cacheFile->valuestring,true,0,&found);
        if (!found){
            ccexRefusalPrint(stderr,"Error: configuration.k_range with the pcce partition cache: cluster.pcce.cache.file \"%s\" needs {k} "
                           "(one cache file per configuration, e.g. \"partition_{k}.cache\")\n",cacheFile->valuestring);
            ccexRefuse();
        }
    }
    // every file of every k, before anything is calculated; a state file is its name with the state s of configuration.
    // s_range ({s}, or appended, A94). The state files are read with s_range only, the Defect files with Defect species.
    long sFirst = 0, sLast = 0;
    const bool sampled = readStateRange(block,&sFirst,&sLast);
    const bool defects = inputHasDefectSpecies(root);
    for (long k = *first; k <= *last; ++k)
        for (const Named& n : names){
            bool found = false;
            const std::string name = inFolder(n.folder,configurationName(n.name,true,k,&found));
            if (!found) continue;               // the same file for every k: the ordinary reader checks it
            if ((n.prefix && !sampled) || (n.key.rfind("configuration.defect",0) == 0 && !defects)) continue;   // not read
            for (long s = n.prefix ? sFirst : 0; s <= (n.prefix ? sLast : 0); ++s){
                const std::string file = n.prefix ? stateFileName(name,s) : name;
                if (access(file.c_str(),R_OK) != 0){
                    ccexRefusalPrint(stderr,"Error: configuration k = %ld: cannot read %s (%s \"%s\")\n",k,file.c_str(),
                                     n.key.c_str(),n.name.c_str());
                    ccexRefuse();
                }
            }
        }
    cJSON_Delete(root);
    return true;
}

// The calculation method and the cluster shape (A47): "method" is cCCE, gCCE, meCCE or megCCE (the master-equation
// methods, adding_MECCE), "cluster_mode" standard, pcce, bundle, dsj, itb or dsjitb. Inside, the code keeps the names it
// compares against - CCE, gCCE, pCCE, dsj, itb, dsjitb, meCCE, megCCE, pgCCE, pmeCCE, pmegCCE, and cluster_mode adaptive for
// bundle - so nothing past the readers changed (and the cluster_only inventory writes what it wrote).
// cluster_mode pcce runs with cCCE, gCCE, meCCE and megCCE: inside pCCE, pgCCE, pmeCCE and pmegCCE (methodUsesPcceClusters
// and methodKernel, general.h); the adding_MECCE branch's method names pgCCE, pmeCCE and pmegCCE are refused with the new
// form. gCCE also runs with bundle; dsj, itb and dsjitb run with cCCE only (gCCE with them is not implemented yet;
// improvements.md A47). noise is a reserved method name (A49): accepted by the reader, refused as not implemented yet.
// A method name: cCCE, gCCE, meCCE or megCCE (true for gCCE); reserved, old and unknown names are refused.
static bool checkedMethodName(const char* method, const char* where){
    const bool gcce = !strcasecmp(method,"gcce");
    const bool mecce = !strcasecmp(method,"mecce"), megcce = !strcasecmp(method,"megcce");
    if (!strcasecmp(method,"noise")){
        ccexRefusalPrint(stderr,"Error: %s \"%s\" is not implemented yet (a reserved name, the methods that run are cCCE, gCCE, meCCE "
                       "and megCCE)\n",where,method);
        ccexRefuse();
    }
    if (!gcce && !mecce && !megcce && strcasecmp(method,"ccce")){
        // CCE, pCCE, dsj, itb and dsjitb never come here: refuseOldInputKeys refuses them with the new form, before any reader
        static const char* branch[][2] = {{"pgcce","gCCE"}, {"pmecce","meCCE"}, {"pmegcce","megCCE"}};   // adding_MECCE
        const char* now = nullptr;
        for (const auto& b : branch) if (!strcasecmp(method,b[0])) now = b[1];
        if (now)
            ccexRefusalPrint(stderr,"Error: %s \"%s\" is now method \"%s\" with \"cluster.mode\": \"pcce\" in the input\n",where,method,now);
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
    if (k == 6){   // "adaptive" (bundle now) never comes here: refuseOldInputKeys refuses it first
        ccexRefusalPrint(stderr,"Error: cluster.mode must be standard, pcce, bundle, dsj, itb or dsjitb, not \"%s\"\n",shape);
        ccexRefuse();
    }
    if ((mecce || megcce) && k != 0 && k != 1){
        ccexRefusalPrint(stderr,"Error: method %s with cluster.mode \"%s\" is not implemented (cluster.mode standard or pcce)\n",
                mecce ? "meCCE" : "megCCE", shapes[k]);
        ccexRefuse();
    }
    if (gcce && k != 0 && k != 1 && k != 2){
        ccexRefusalPrint(stderr,"Error: method gCCE with cluster.mode \"%s\" is not implemented yet, cluster.mode %s runs with method cCCE\n",
                shapes[k],shapes[k]);
        ccexRefuse();
    }
    static const char* internal[] = {"CCE", "pCCE", "CCE", "dsj", "itb", "dsjitb"};
    if (k == 1 && (gcce || mecce || megcce)) strcpy(methodOut,gcce ? "pgCCE" : mecce ? "pmeCCE" : "pmegCCE");
    else strcpy(methodOut,gcce ? "gCCE" : mecce ? "meCCE" : megcce ? "megCCE" : internal[k]);
    strcpy(modeOut,k == 2 ? "adaptive" : "standard");
}

static cJSON* readBlockOrEmpty(cJSON* parent, const char* key);
static void readMethod(cJSON* root, char* methodOut, char* modeOut){
    cJSON* mode = cJSON_GetObjectItem(readBlockOrEmpty(root,"cluster"),"mode");   // cluster.mode (A65; cluster_mode before)
    if (mode && (!cJSON_IsString(mode) || !mode->valuestring)){
        ccexRefusalPrint(stderr,"Error: cluster.mode must be a string (standard, pcce, bundle, dsj, itb, dsjitb)\n");
        ccexRefuse();
    }
    const char* method = cJSON_ReadString(root,(char*)"method",true,(char*)"cCCE");
    internalMethod(method,mode ? mode->valuestring : "standard","method",methodOut,modeOut);
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
// A65: a block whose keys all have defaults (execution, cluster, output, frame): an absent one is read as an empty object,
// added to parent (freed with it) under its path, so its readers look keys up as in a written one
static cJSON* readBlockOrEmpty(cJSON* parent, const char* key){
    if (cJSON* item = readBlock(parent,key)) return item;
    cJSON* empty = cJSON_CreateObject();
    cJSON_AddItemToObject(parent,key,empty);
    const auto path = objectPaths().find(parent);
    const auto shown = shownPaths().find(parent);
    recordPaths(empty, joinKey(path == objectPaths().end() ? "" : path->second, lowerKey(key)),
                joinKey(shown == shownPaths().end() ? "" : shown->second, key));
    return empty;
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

////////////////////////////////////////////////////////////////////////////////////////////////
// A90 (2026-10-16): every place that holds directional numbers says which frame they are written in, with the one key
// "frame": "file" (the bath file's frame, before the rotation: turned with the positions) or "qubit" (the computational
// frame, qubit_axis = z: used as written). With frame.rotation on it is required (no default); with the rotation off it
// is not used, and a written one is warned about once per key.
enum FrameValue { FRAME_ABSENT, FRAME_FILE, FRAME_QUBIT };
static bool frameRotationOn = false;   // frame.rotation, set by cJSON_readOptionConfig before the other readers run
static const char* const kFrameValues =
    "\"file\" (written in the bath file's frame, before the rotation: turned with the positions) or \"qubit\" (written in "
    "the computational frame, qubit_axis = z: used as written)";

// parent.frame, checked: absent, "file" or "qubit" (case ignored); "bath" (the value until A90) is refused with "file"
static FrameValue readFrameKey(cJSON* parent, const std::string& shown){
    cJSON* item = cJSON_GetObjectItem(parent,"frame");
    if (item == NULL) return FRAME_ABSENT;
    const char* v = (cJSON_IsString(item) && item->valuestring) ? item->valuestring : NULL;
    if (v && !strcasecmp(v,"file")) return FRAME_FILE;
    if (v && !strcasecmp(v,"qubit")) return FRAME_QUBIT;
    if (v && !strcasecmp(v,"bath"))
        ccexRefusalPrint(stderr,"Error: %s.frame \"bath\" is now \"file\" (written in the bath file's frame, turned with the "
                         "positions)\n%s\n",shown.c_str(),kConvertHint);
    else
        ccexRefusalPrint(stderr,"Error: %s.frame must be %s\n",shown.c_str(),kFrameValues);
    ccexRefuse();
}

// A frame given while the rotation is off: once per key (generic, e.g. "Defect.species[].frame")
static void warnFrameUnused(FrameValue f, const char* generic){
    static std::set<std::string> warned;
    if (f == FRAME_ABSENT || frameRotationOn || !warned.insert(generic).second) return;
    warnNotUsed(generic,"frame.rotation is off, so the numbers are used as written");   // A94 (G3): one wording
}

// parent.frame of numbers that need it: required with the rotation on, warned about with it off
static FrameValue frameOf(cJSON* parent, const std::string& shown, const char* generic, const char* what){
    const FrameValue f = readFrameKey(parent,shown);
    if (!frameRotationOn){ warnFrameUnused(f,generic); return f; }
    if (f == FRAME_ABSENT){
        ccexRefusalPrint(stderr,"Error: frame.rotation is on, so %s.frame must say the frame of %s: %s\n",
                         shown.c_str(),what,kFrameValues);
        ccexRefuse();
    }
    return f;
}

// A94 (R3): a unit is spelled as written in the table, case and all, as the A67 units are (mHz is not MHz, mb not Mb).
// Until A94 the Defect and Qubit.intmap units were compared without regard to case; a spelling that differs in case only
// is refused with the table's.
static bool unitIs(const char* given, const char* name, const std::string& where){
    if (!strcmp(given,name)) return true;
    if (strcasecmp(given,name)) return false;
    ccexRefusalPrint(stderr,"Error: %s.unit \"%s\" is written \"%s\" (units are case-sensitive, as mHz and MHz differ)\n",
                     where.c_str(),given,name);
    ccexRefuse();
}

static double readDefectFrequencyScaleToMHz(cJSON* object, const char* dfname,
                                            const char* field_name){
    const std::string where = std::string("Defect.species[") + dfname + "]." + field_name;

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; } // documented default
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.unit must be a string\n",dfname,field_name);
        ccexRefuse();
    }

    if (unitIs(item->valuestring,"Hz",where)){ return 1.0e-6; }
    if (unitIs(item->valuestring,"kHz",where)){ return 1.0e-3; }
    if (unitIs(item->valuestring,"MHz",where)){ return 1.0; }
    if (unitIs(item->valuestring,"GHz",where)){ return 1.0e3; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.unit = \"%s\" is not supported\n",
            dfname,field_name,item->valuestring);
    ccexRefusalPrint(stderr,"Supported frequency units are Hz, kHz, MHz and GHz.\n");
    ccexRefuse();
}

static double readIntmapFrequencyScaleToKHz(cJSON* intmap, const char* qubit1_name,
                                            const char* qubit2_name){
    const std::string where = std::string("Qubit.intmap[\"") + qubit1_name + "\",\"" + qubit2_name + "\"]";

    cJSON* item = cJSON_GetObjectItem(intmap,"unit");
    if (item == NULL){ return 1.0; } // legacy Qubit.intmap tensors are in kHz
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].unit must be a string\n",
                qubit1_name,qubit2_name);
        ccexRefuse();
    }

    if (unitIs(item->valuestring,"Hz",where)){ return 1.0e-3; }
    if (unitIs(item->valuestring,"kHz",where)){ return 1.0; }
    if (unitIs(item->valuestring,"MHz",where)){ return 1.0e3; }
    if (unitIs(item->valuestring,"GHz",where)){ return 1.0e6; }

    ccexRefusalPrint(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].unit = \"%s\" is not supported\n",
            qubit1_name,qubit2_name,item->valuestring);
    ccexRefusalPrint(stderr,"Supported frequency units are Hz, kHz, MHz and GHz.\n");
    ccexRefuse();
}

// A94 (R3): a Defect value object takes "values" (its key until A94) or "value" (the A67 keys' word), one of the two
static cJSON* defectValues(cJSON* object, const char* dfname, const char* field){
    cJSON* values = cJSON_GetObjectItem(object,"values");
    cJSON* value = cJSON_GetObjectItem(object,"value");
    if (values && value){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s has both \"value\" and \"values\" (one key, either spelling)\n",
                         dfname,field);
        ccexRefuse();
    }
    return values ? values : value;
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

static void readDefectZfs(cJSON* defect, DefectArray* dfa, int idf, int nconfig){

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
    if (nconfig < 1){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s] has an object-form zfs but n_axis_states = %d\n",df->dfname,nconfig);
        ccexRefusalPrint(stderr,"At least one axis state (n_axis_states) is required.\n");
        ccexRefuse();
    }

    cJSON* D_item = cJSON_GetObjectItem(item,"D");
    cJSON* tensor_item = cJSON_GetObjectItem(item,"tensor");
    cJSON* axis_item = cJSON_GetObjectItem(item,"axis");
    if ((D_item != NULL) == (tensor_item != NULL)){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs must contain exactly one of D or tensor\n",df->dfname);
        ccexRefuse();
    }
    // A94 (G5): E took 0 only (an axis fixes z, not the x and y an E needs), so it was removed; the axis belongs to the D
    // form (Defect.species[].axis until A94, refused by refuseOldInputKeys)
    if (untrackedItem(item,"E")){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs.E was removed (only 0 was taken): leave it out of the D form, or "
                         "write zfs.tensor for an E other than 0\n%s\n",df->dfname,kConvertHint);
        ccexRefuse();
    }
    if (tensor_item != NULL && axis_item != NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs.axis belongs to the D form (zfs.tensor has its own axes)\n",df->dfname);
        ccexRefuse();
    }
    if (D_item != NULL){
        char section_label[160];
        snprintf(section_label,sizeof(section_label),"Defect.species[%s].zfs",df->dfname);
        double axis_raw[3] = {0.0,0.0,0.0};
        double axis_normalized[3];
        if (axis_item == NULL){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].zfs uses the D form but zfs.axis is missing (the symmetry axis "
                             "D is along)\n",df->dfname);
            ccexRefuse();
        }
        readAxis3(item,section_label,"axis",axis_raw,true);
        normalizeAxisOrDie(axis_raw,axis_normalized,(std::string(section_label) + ".axis").c_str());
        DefectArray_setDefect_idf_axis(dfa,idf,axis_normalized);
    }

    double scale_to_mhz = readDefectFrequencyScaleToMHz(item,df->dfname,"zfs");
    MatrixXcd zfs = MatrixXcd::Zero(3,3);
    DefectZfsInputMode mode;

    if (D_item != NULL){
        char what[160];
        snprintf(what,sizeof(what),"Defect.species[%s].zfs.D",df->dfname);
        double D_mhz = scale_to_mhz * readFiniteNumber(D_item,what);

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

    cJSON* values = defectValues(item,df->dfname,"detuning");
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
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values[%d] must be [axis state, \"e\", value]\n",
                    df->dfname,i);
            ccexRefuse();
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].detuning.values[%d] axis state must be an integer in 1..%d\n",
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

// the gyromagnetic ratio units of Defect.species[].gyros and (A94, G20) Qubit.qubit[].gyro, to rad kHz/G
static double gyroScaleToRadkHzPerG(const char* unit, const std::string& where){
    if (unitIs(unit,"radkHz/G",where) || unitIs(unit,"rad/ms/G",where)){ return 1.0; }
    if (unitIs(unit,"kHz/G",where)){ return 2.0*M_PI; }
    if (unitIs(unit,"MHz/G",where)){ return 2.0*M_PI*1.0e3; }
    if (unitIs(unit,"MHz/T",where)){ return 2.0*M_PI*1.0e-1; }
    if (unitIs(unit,"Hz/T",where)){ return 2.0*M_PI*1.0e-7; }
    if (unitIs(unit,"radMHz/T",where) || unitIs(unit,"rad/us/T",where)){ return 1.0e-1; }
    if (unitIs(unit,"rad/s/T",where)){ return 1.0e-7; }

    ccexRefusalPrint(stderr,"Error: %s.unit = \"%s\" is not supported\n",where.c_str(),unit);
    ccexRefusalPrint(stderr,"Supported gyro units are radkHz/G, rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T, radMHz/T and rad/s/T.\n");
    ccexRefuse();
}

static double readDefectGyroScaleToRadkHzPerG(cJSON* object, const char* dfname){
    const std::string where = std::string("Defect.species[") + dfname + "].gyros";

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].gyros.unit must be a string\n",dfname);
        ccexRefuse();
    }
    return gyroScaleToRadkHzPerG(item->valuestring,where);
}

// Qubit.qubit[].gyro (A94, G20): a number in rad kHz/G as before, or {"value": number, "unit": ...} with the units of
// Defect.species[].gyros (the unit rad kHz/G leaves the number as written)
static double readQubitGyro(cJSON* qubit, int iqubit, double gyroDefault){
    cJSON* item = cJSON_GetObjectItem(qubit,"gyro");
    if (!cJSON_IsObject(item)) return cJSON_ReadDouble(qubit,(char*)"gyro",true,gyroDefault);
    const std::string where = "Qubit.qubit[" + std::to_string(iqubit) + "].gyro";
    cJSON* value = NULL; cJSON* unit = NULL;
    int n = 0;
    for (cJSON* c = item->child; c; c = c->next, ++n){
        if (!strcasecmp(c->string,"value")) value = c;
        else if (!strcasecmp(c->string,"unit")) unit = c;
    }
    if (n != 2 || !cJSON_IsNumber(value) || !cJSON_IsString(unit)){
        ccexRefusalPrint(stderr,"Error: %s must be a number (rad kHz/G) or {\"value\": number, \"unit\": \"...\"}\n",where.c_str());
        ccexRefuse();
    }
    cJSON_GetObjectItem(item,"value"); cJSON_GetObjectItem(item,"unit");   // read (cJSON_checkInputKeys)
    const double scale = gyroScaleToRadkHzPerG(unit->valuestring,where);
    const double v = readFiniteNumber(value,(where + ".value").c_str());
    return scale == 1.0 ? v : scale * v;
}

static double readDefectEqScaleTo1eMinus30m2(cJSON* object, const char* dfname){
    const std::string where = std::string("Defect.species[") + dfname + "].quadrupole_moments";

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].quadrupole_moments.unit must be a string\n",dfname);
        ccexRefuse();
    }

    const char* unit = item->valuestring;
    if (unitIs(unit,"1e-30 m^2",where) || unitIs(unit,"1e-30m^2",where)
        || unitIs(unit,"10^-30 m^2",where) || unitIs(unit,"10^-30m^2",where)
        || unitIs(unit,"fm^2",where) || unitIs(unit,"fm2",where)
        || unitIs(unit,"10 mbarn",where)){ return 1.0; }
    if (unitIs(unit,"m^2",where) || unitIs(unit,"m2",where)){ return 1.0e30; }
    if (unitIs(unit,"cm^2",where) || unitIs(unit,"cm2",where)){ return 1.0e26; }
    if (unitIs(unit,"barn",where) || unitIs(unit,"b",where)){ return 1.0e2; }
    if (unitIs(unit,"mbarn",where) || unitIs(unit,"millibarn",where)
        || unitIs(unit,"mb",where)){ return 1.0e-1; }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].quadrupole_moments.unit = \"%s\" is not supported\n",dfname,unit);
    ccexRefusalPrint(stderr,"Supported quadrupole-moment units are 1e-30 m^2, fm^2, m^2, cm^2, barn and mbarn.\n");
    ccexRefuse();
}

static double readDefectEfgScaleToHartreePerBohr2(cJSON* object, const char* dfname){
    const std::string where = std::string("Defect.species[") + dfname + "].efg";

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        ccexRefusalPrint(stderr,"Error: Defect.species[%s].efg.unit must be a string\n",dfname);
        ccexRefuse();
    }

    const char* unit = item->valuestring;
    if (unitIs(unit,"Hartree/Bohr^2",where) || unitIs(unit,"Ha/Bohr^2",where)
        || unitIs(unit,"Ha/a0^2",where) || unitIs(unit,"a.u.",where)
        || unitIs(unit,"au",where) || unitIs(unit,"atomic_unit",where)){ return 1.0; }

    // Keep these constants identical to BathSpin_setQuad_fromEFG. Numerically, one
    // eV of electron potential energy per m^2 corresponds to one V/m^2; the electron
    // charge is implicit in the stored eQ convention used by that routine.
    const double hartree_to_ev = 27.211386;
    const double bohr_to_m = 0.5291772e-10;
    const double volt_per_m2_to_au = bohr_to_m*bohr_to_m/hartree_to_ev;

    if (unitIs(unit,"V/m^2",where) || unitIs(unit,"V/m2",where)){
        return volt_per_m2_to_au;
    }
    if (unitIs(unit,"V/angstrom^2",where) || unitIs(unit,"V/angstrom2",where)
        || unitIs(unit,"V/A^2",where) || unitIs(unit,"V/A2",where)){
        return 1.0e20*volt_per_m2_to_au;
    }

    ccexRefusalPrint(stderr,"Error: Defect.species[%s].efg.unit = \"%s\" is not supported\n",dfname,unit);
    ccexRefusalPrint(stderr,"Supported EFG units are Hartree/Bohr^2 (a.u.), V/m^2 and V/angstrom^2.\n");
    ccexRefuse();
}

static cJSON* readDefectObjectValues(cJSON* object, const char* dfname,
                                     const char* field_name, int expected_count){

    cJSON* values = defectValues(object,dfname,field_name);
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
                                            int nconfig, int naddspin){

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

    cJSON* values = defectValues(item,df->dfname,field_name);
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
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] must be [axis state, spin, tensor]\n",
                    df->dfname,field_name,i);
            ccexRefuse();
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            ccexRefusalPrint(stderr,"Error: Defect.species[%s].%s.values[%d] axis state must be an integer in 1..%d\n",
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
    // A65: how it is computed, the clusters, what is written, the frame rotation (absent: read as empty, every default)
    cJSON* execution = readBlockOrEmpty(root,"execution");
    cJSON* clusterBlock = readBlockOrEmpty(root,"cluster");
    cJSON* output = readBlockOrEmpty(root,"output");
    cJSON* frame = readBlockOrEmpty(root,"frame");
    // the configuration block (A64): {k} of the configuration main.cpp runs (cJSON_setConfigurationIndex)
    cJSON* configuration = readBlock(root,"configuration");
    long rangeFirst = 0, rangeLast = 0;
    const bool ranged = readConfigurationRange(configuration,&rangeFirst,&rangeLast);
    const long configurationIndex = configurationRanged ? configurationK : rangeFirst;
    const std::string folder = configurationFolder(configuration);

    ////////////////////////////////////////////////////////////////////////
    // General Options
    ////////////////////////////////////////////////////////////////////////
    if (rank==0 && verbosity){
        printf("    >> Read Config Options\n\n");
    }

    if (rank==0 && verbosity){
        printf("      - General option-related keys :\n");
        printf("        [ method, quantity, order, bfield, rbath, rdip, deltat, nstep ]\n");
        printf("        execution : [ backend, propagator, evolution, evolution_update, profile,\n");
        printf("          gpu_options { min_dim, packed_max_dim, warp_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb },\n");
        printf("          gcce_refinement { enabled, max_sweeps }, column_evolution { enabled, min_dim, columns } ]\n");
        printf("        cluster : [ mode, pcce { partition, connections, print_trials, cache { mode, file } },\n");
        printf("          bundle { max_dim, merge_cutoff, connections, connection_budget, max_clusters, max_work }, adaptive { probe_memory_gb } ]\n");
        printf("        output : [ cluster_only, cluster_contributions { enabled, max_mb } ]\n\n");
    }
    char method[MAX_CHARARRAY_LENGTH], clusterMode[16];
    readMethod(root,method,clusterMode);
    Config_setMethod(cnf,method); // inside: CCE, gCCE, pCCE, dsj, dsjitb, itb (from method and cluster_mode)
    strcpy(cnf->cluster_mode,clusterMode);   // inside: standard, or adaptive for cluster_mode bundle

    char* quantity = cJSON_ReadString(root,"quantity",true,"coherence");
    Config_setQuantity(cnf,quantity); // Current possible options : coherence, noise, dm

    // eigen (default) | expm (gCCE, legacy) | sparse (A2), as a string or as {"method": ..., "min_dim": ..., "tol": ...};
    // min_dim and tol belong to sparse. Config_validatePropagator checks the paths once every source is in.
    cJSON* propagatorItem = cJSON_GetObjectItem(execution,"propagator");
    if (cJSON_IsObject(propagatorItem)){
        for (cJSON* key = propagatorItem->child; key; key = key->next)
            if (strcasecmp(key->string,"method") && strcasecmp(key->string,"min_dim") && strcasecmp(key->string,"tol")){
                ccexRefusalPrint(stderr,"Error: execution.propagator has no key \"%s\" (method, min_dim, tol)\n",key->string);
                ccexRefuse();
            }
        cJSON* method = cJSON_GetObjectItem(propagatorItem,"method");
        if (!cJSON_IsString(method) || !method->valuestring){
            ccexRefusalPrint(stderr,"Error: execution.propagator.method must be a string (eigen, expm or sparse)\n");
            ccexRefuse();
        }
        Config_setPropagator(cnf,method->valuestring);
        cJSON* minDim = cJSON_GetObjectItem(propagatorItem,"min_dim");
        cJSON* tol = cJSON_GetObjectItem(propagatorItem,"tol");
        if ((minDim || tol) && strcasecmp(cnf->propagator,"sparse")){
            ccexRefusalPrint(stderr,"Error: execution.propagator.min_dim and execution.propagator.tol belong to method sparse (this has %s)\n",cnf->propagator);
            ccexRefuse();
        }
        if (minDim){
            const double value = readFiniteNumber(minDim,"execution.propagator.min_dim");
            if (value < 1 || value > 2147483647.0 || std::floor(value) != value){
                ccexRefusalPrint(stderr,"Error: execution.propagator.min_dim must be a positive 32-bit integer\n");
                ccexRefuse();
            }
            cnf->sparse_min_dim = (long)value;
        }
        if (tol){
            const double value = readFiniteNumber(tol,"execution.propagator.tol");
            if (!(value > 0.0 && value < 1.0)){
                ccexRefusalPrint(stderr,"Error: execution.propagator.tol must be a number between 0 and 1 (exclusive)\n");
                ccexRefuse();
            }
            cnf->sparse_tol = value;
        }
    }else{
        if (propagatorItem && !cJSON_IsString(propagatorItem)){
            ccexRefusalPrint(stderr,"Error: execution.propagator must be a string (eigen, expm or sparse) or an object {\"method\": ..., "
                           "\"min_dim\": ..., \"tol\": ...}\n");
            ccexRefuse();
        }
        char* propagator = cJSON_ReadString(execution,"propagator",true,"eigen");
        Config_setPropagator(cnf,propagator);
    }

    // vector is the DEFAULT, but it only exists where its preconditions hold: gCCE, a pure
    // rho0 (nstate > 0) and propagator=eigen. Do NOT decide that here: only record whether
    // the key was written; Config_resolveEvolution settles it once every reader has run.
    char* evolution = cJSON_ReadString(execution,"evolution",true,"vector");
    Config_setEvolution(cnf,evolution); // gCCE only : vector (default) | matrix (legacy)
    cnf->evolution_isdefault = (cJSON_GetObjectItem(execution,"evolution") == NULL);
    // exact | recursive (A16 b); without the key Config_validateEvolutionUpdate picks the default once every source is in (A29)
    cJSON* evolutionUpdate = cJSON_GetObjectItem(execution,"evolution_update");
    if (evolutionUpdate){
        if (!cJSON_IsString(evolutionUpdate) || !evolutionUpdate->valuestring ||
            (strcasecmp(evolutionUpdate->valuestring,"exact") &&
             strcasecmp(evolutionUpdate->valuestring,"recursive"))){
            ccexRefusalPrint(stderr,"Error: execution.evolution_update must be exact or recursive\n");
            ccexRefuse();
        }
        strcpy(cnf->evolution_update,evolutionUpdate->valuestring);
    }

    cJSON* backend = cJSON_GetObjectItem(execution,"backend");
    if (backend){
        if (!cJSON_IsString(backend) || !backend->valuestring ||
            (strcasecmp(backend->valuestring,"cpu") &&
             strcasecmp(backend->valuestring,"cuda") &&
             strcasecmp(backend->valuestring,"auto"))){
            ccexRefusalPrint(stderr,"Error: execution.backend must be cpu, cuda, or auto\n");
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
        if (untrackedItem(root,oldGpuKeys[k]) || untrackedItem(execution,oldGpuKeys[k])){
            ccexRefusalPrint(stderr,"Error: %s is now execution.gpu_options.%s (\"execution\": {\"gpu_options\": {\"%s\": ...}})\n",oldGpuKeys[k],gpuKeys[k],gpuKeys[k]);
            ccexRefuse();
        }
    if (cJSON* block = cJSON_GetObjectItem(execution,"gpu_options")){
#ifndef CCEX_USE_CUDA
        // A94 (G3): main.out is built without CUDA (main_gpu.out reads them); the values are still checked
        warnNotUsed("execution.gpu_options","this main.out is built without CUDA, every cluster runs on the CPU (main_gpu.out uses them)");
#endif
        if (!cJSON_IsObject(block)){
            ccexRefusalPrint(stderr,"Error: execution.gpu_options must be an object, e.g. \"gpu_options\": {\"min_dim\": 256}\n");
            ccexRefuse();
        }
        for (cJSON* item = block->child; item; item = item->next){
            int k = 0;
            while (k<7 && strcasecmp(item->string,gpuKeys[k])) ++k;
            if (k==7){
                ccexRefusalPrint(stderr,"Error: unknown key execution.gpu_options.%s (min_dim, packed_max_dim, warp_max_dim, evolution_batch_size, "
                               "evolution_batch_steps, cluster_batch_mb, cache_mb)\n",item->string);
                ccexRefuse();
            }
        }
        for (int k=0; k<7; ++k){
            cJSON* item = cJSON_GetObjectItem(block,gpuKeys[k]);   // recorded as read (cJSON_checkInputKeys)
            if (!item) continue;
            char name[64];
            snprintf(name,sizeof name,"execution.gpu_options.%s",gpuKeys[k]);
            const double value = readFiniteNumber(item,name);
            if (k==6){
                // the warp kernel (A40): clusters of a compiled total dimension up to this value; 0 switches it off
                if (value < 0 || value > 32 || std::floor(value) != value){
                    ccexRefusalPrint(stderr,"Error: execution.gpu_options.warp_max_dim must be 0 (off) or an integer up to 32\n");
                    ccexRefuse();
                }
            } else if (k==5){
                // the packed kernel exists for bath dimension 2 and 4 (small_cluster_math.h); 0 switches it off
                if (value!=0 && value!=2 && value!=4){
                    ccexRefusalPrint(stderr,"Error: execution.gpu_options.packed_max_dim must be 0 (off), 2 or 4\n");
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
    if (cJSON* item = cJSON_GetObjectItem(readBlockOrEmpty(clusterBlock,"adaptive"),"probe_memory_gb")){
        const double value = readFiniteNumber(item,"cluster.adaptive.probe_memory_gb");
        if (value < 0){ ccexRefusalPrint(stderr,"Error: cluster.adaptive.probe_memory_gb must be nonnegative\n"); ccexRefuse(); }
        cnf->probe_memory_gb = value;
    }
    // gcce_refinement and cluster_contributions are blocks (A47): {"enabled", ...}; a block without "enabled" is on, as an
    // auto_convergence block is (main.cpp, inputAsksAutoConv)
    if (cJSON* block = readBlock(execution,"gcce_refinement")){
        cnf->gcce_refinement=readStrictBool(block,"execution.gcce_refinement","enabled",true);
        if(auto* item=cJSON_GetObjectItem(block,"max_sweeps")) {
            const double value=readFiniteNumber(item,"execution.gcce_refinement.max_sweeps");
            if(value<1 || value>INT_MAX || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: execution.gcce_refinement.max_sweeps must be a positive 32-bit integer\n");
                ccexRefuse();
            }
            cnf->gcce_refinement_max_sweeps=int(value);
        }
    }
    // A62: the one-state cCCE time points as matrix columns for large clusters; a block without "enabled" is on, as above
    if (cJSON* block = readBlock(execution,"column_evolution")){
        cnf->column_evolution=readStrictBool(block,"execution.column_evolution","enabled",true);
        const char* keys[]={"min_dim","columns"};
        int* values[]={&cnf->column_evolution_min_dim,&cnf->column_evolution_columns};
        for(int k=0;k<2;++k){
            cJSON* item=cJSON_GetObjectItem(block,keys[k]);
            if(!item)continue;
            const std::string name="execution.column_evolution."+std::string(keys[k]);
            const double value=readFiniteNumber(item,name.c_str());
            if(value<1 || value>INT_MAX || std::floor(value)!=value){
                ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());ccexRefuse();
            }
            *values[k]=int(value);
        }
    }
    cnf->profile = readStrictBool(execution,"execution","profile",false);
    // A63: the plan of the old -P option: <outfile>.clusters.json, then stop before the dynamics (checked in main.cpp)
    if (cJSON* item = cJSON_GetObjectItem(output,"cluster_only")){
        if (!cJSON_IsBool(item)){
            ccexRefusalPrint(stderr,"Error: output.cluster_only must be a boolean (true or false)\n");
            ccexRefuse();
        }
        cnf->cluster_only = cJSON_IsTrue(item);
    }
    if (cJSON* block = readBlock(output,"cluster_contributions")){
        cnf->cluster_contributions = readStrictBool(block,"output.cluster_contributions","enabled",true);
        if(auto* item=cJSON_GetObjectItem(block,"max_mb")) {
            const double value=readFiniteNumber(item,"output.cluster_contributions.max_mb");
            if(value<1 || value>INT_MAX || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: output.cluster_contributions.max_mb must be a positive 32-bit integer\n");ccexRefuse();
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
    cJSON* bundle=readBlock(clusterBlock,"bundle");
    if(bundle)cnf->adaptive_options_present=true;
    const char* bundleKeys[]={"max_dim","merge_cutoff"};
    for(int k=0;k<2 && bundle;++k) {
        cJSON* item=cJSON_GetObjectItem(bundle,bundleKeys[k]);
        if(!item)continue;
        const std::string name="cluster.bundle."+std::string(bundleKeys[k]);
        double value=readFiniteNumber(item,name.c_str());
        if(k<1) {
            if(value<1 || value>2147483647.0 || std::floor(value)!=value) {
                ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());ccexRefuse();
            }
            cnf->bundle_max_dim=int(value);
        } else {
            if(value<0) { ccexRefusalPrint(stderr,"Error: cluster.bundle.merge_cutoff must be nonnegative\n");ccexRefuse(); }
            cnf->bundle_merge_cutoff=value;
        }
    }
    const char* expansionKeys[]={"max_clusters","max_work"};
    int* expansionValues[]={&cnf->expansion_max_clusters,&cnf->expansion_max_work};
    for(int k=0;k<2 && bundle;++k)if(auto* item=cJSON_GetObjectItem(bundle,expansionKeys[k])) {
        const std::string name="cluster.bundle."+std::string(expansionKeys[k]);
        const double value=readFiniteNumber(item,name.c_str());
        if(value<1 || value>INT_MAX || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());ccexRefuse();
        }
        *expansionValues[k]=int(value);
    }
    cJSON* connections=bundle ? cJSON_GetObjectItem(bundle,"connections") : NULL;
    if(connections) {
        // A94 (G8): checked, the mode a bundle input without rdip and connections runs, can be written too
        if(!cJSON_IsString(connections) || !connections->valuestring ||
           (strcasecmp(connections->valuestring,"distance") && strcasecmp(connections->valuestring,"auto") &&
            strcasecmp(connections->valuestring,"checked"))) {
            ccexRefusalPrint(stderr,"Error: cluster.bundle.connections must be distance, auto or checked\n");ccexRefuse();
        }
        strcpy(cnf->bundle_connections,connections->valuestring);
    }
    cJSON* connectionBudget=bundle ? cJSON_GetObjectItem(bundle,"connection_budget") : NULL;
    if(connectionBudget) {
        cnf->bundle_connection_budget=readFiniteNumber(connectionBudget,"cluster.bundle.connection_budget");
        if(cnf->bundle_connection_budget<0) {
            ccexRefusalPrint(stderr,"Error: cluster.bundle.connection_budget must be nonnegative\n");ccexRefuse();
        }
    }
    // adaptive.work_file: the request the refinement driver writes into a step's input; the rest of the adaptive block is
    // the refinement request itself, read in main.cpp (adaptiveArgsFromInput)
    cJSON* adaptiveBlock=readBlock(clusterBlock,"adaptive");
    cJSON* adaptiveWork=adaptiveBlock ? cJSON_GetObjectItem(adaptiveBlock,"work_file") : NULL;
    if(adaptiveWork) {
        cnf->adaptive_options_present=true;
        if(!cJSON_IsString(adaptiveWork) || !adaptiveWork->valuestring || !*adaptiveWork->valuestring ||
           strlen(adaptiveWork->valuestring)>=MAX_FILEPATH) {
            ccexRefusalPrint(stderr,"Error: cluster.adaptive.work_file must be a nonempty bounded path\n");ccexRefuse();
        }
        strcpy(cnf->adaptive_work_file,adaptiveWork->valuestring);
    }
    // the pcce block: cluster_mode pcce only (partition_size, trials, max_iterations, kmeans_pp and print_trials are read
    // with the clusters, cJSON_readOptionCluster)
    cJSON* pcce=readBlock(clusterBlock,"pcce");
    if(pcce && !methodUsesPcceClusters(cnf->method)) {   // the rest: standard, bundle (adaptive inside), dsj, itb, dsjitb
        ccexRefusalPrint(stderr,"Error: the cluster.pcce block belongs to cluster.mode pcce (this input has cluster.mode %s)\n",
                !strcasecmp(cnf->cluster_mode,"adaptive") ? "bundle" : !strcasecmp(cnf->method,"gcce") || !strcasecmp(cnf->method,"cce")
                || !strcasecmp(cnf->method,"mecce") || !strcasecmp(cnf->method,"megcce") ? "standard" : cnf->method);
        ccexRefuse();
    }
    // A94 (G7): print_trials is the number of rank 0's trials printed (0, the default: none). It was a switch beside
    // print_limit, the most printed, two keys for one thing (print_trials false and print_limit 0 both printed none).
    if(untrackedItem(pcce,"print_limit")) {
        ccexRefusalPrint(stderr,"Error: cluster.pcce.print_limit is now cluster.pcce.print_trials: the number of trials printed "
                       "(0: none)\n%s\n",kConvertHint);
        ccexRefuse();
    }
    if(auto* item=pcce ? cJSON_GetObjectItem(pcce,"print_trials") : NULL) {
        if(cJSON_IsBool(item)) {
            ccexRefusalPrint(stderr,"Error: cluster.pcce.print_trials is a number, the trials printed (0: none; true with the "
                           "old default print_limit is 32)\n%s\n",kConvertHint);
            ccexRefuse();
        }
        const double value=readFiniteNumber(item,"cluster.pcce.print_trials");
        if(value<0 || value>INT_MAX || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: cluster.pcce.print_trials must be a nonnegative 32-bit integer (0: none)\n");
            ccexRefuse();
        }
        cnf->pcce_detail_limit=int(value);
    }
    cJSON* partition = pcce ? cJSON_GetObjectItem(pcce,"partition") : NULL;
    if (partition){
        if (!cJSON_IsString(partition) || !partition->valuestring ||
            (strcasecmp(partition->valuestring,"deterministic") &&
             strcasecmp(partition->valuestring,"legacy"))){
            ccexRefusalPrint(stderr,"Error: cluster.pcce.partition must be deterministic or legacy\n");
            ccexRefuse();
        }
        strcpy(cnf->pcce_partition,partition->valuestring);
    }
    cJSON* pcceConnections = pcce ? cJSON_GetObjectItem(pcce,"connections") : NULL;
    if (pcceConnections){
        if (!cJSON_IsString(pcceConnections) || !pcceConnections->valuestring ||
            (strcasecmp(pcceConnections->valuestring,"centers") &&
             strcasecmp(pcceConnections->valuestring,"spins"))){
            ccexRefusalPrint(stderr,"Error: cluster.pcce.connections must be centers or spins\n");
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
            ccexRefusalPrint(stderr,"Error: cluster.pcce.cache.mode must be off, read, or readwrite\n");
            ccexRefuse();
        }
        strcpy(cnf->pcce_cache_mode,cacheMode->valuestring);
    }
    cJSON* cacheFile = cache ? cJSON_GetObjectItem(cache,"file") : NULL;
    if (cacheFile){
        if (!cJSON_IsString(cacheFile) || !cacheFile->valuestring ||
            !cacheFile->valuestring[0] || strlen(cacheFile->valuestring)>=MAX_FILEPATH){
            ccexRefusalPrint(stderr,"Error: cluster.pcce.cache.file must be a nonempty path shorter than %d bytes\n",MAX_FILEPATH);
            ccexRefuse();
        }
        // {k} as in the configuration's files (A64: one cache file per configuration of configuration.k_range)
        const std::string name = configurationName(cacheFile->valuestring,ranged,configurationIndex);
        if (name.size() >= MAX_FILEPATH){
            ccexRefusalPrint(stderr,"Error: cluster.pcce.cache.file must be a nonempty path shorter than %d bytes\n",MAX_FILEPATH);
            ccexRefuse();
        }
        strcpy(cnf->pcce_cache_file,name.c_str());
        // A94: read mode takes the stored partition only, so its file has to be there (readwrite makes a missing one)
        if (!strcasecmp(cnf->pcce_cache_mode,"read")) refuseUnreadable("cluster.pcce.cache.file",name);
    }

    int order = cJSON_ReadInt(root,"order",false,1);
    if(!strcasecmp(cnf->cluster_mode,"adaptive")) {
        const double value=readFiniteNumber(cJSON_GetObjectItem(root,"order"),"order");
        if(value<1 || value>=INT_MAX || value>cnf->expansion_max_clusters || std::floor(value)!=value) {
            ccexRefusalPrint(stderr,"Error: cluster.mode bundle requires a positive integer order (of bundles) up to cluster.bundle.max_clusters\n");
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
        ccexRefusalPrint(stderr,"Use \"bfield\": [0,0,Bz] with frame.rotation.\n");
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
    if(adaptive && radius && !strcasecmp(cnf->bundle_connections,"checked")) {
        ccexRefusalPrint(stderr,"Error: cluster.bundle.connections checked takes no rdip (the run checks the links itself), "
                       "remove rdip or write connections distance or auto\n");
        ccexRefuse();
    }
    const bool autoConnections=adaptive && ( !strcasecmp(cnf->bundle_connections,"auto") ||
                                             !strcasecmp(cnf->bundle_connections,"checked"));
    cnf->adaptive_rdip_is_set=radius!=nullptr;
    float rdip;
    if(adaptive && radius) {
        rdip=float(readFiniteNumber(radius,rdipName));
        if(!std::isfinite(rdip) || rdip<=0) {
            ccexRefusalPrint(stderr,"Error: explicit cluster.mode bundle rdip must be finite and positive, omit it for unlimited auto connections\n");ccexRefuse();
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

    // the sampled states and the seed: the configuration's (A64). A94: configuration.s_range [first, last] (nstate N until
    // A94 was [1, N]); absent: the ensemble. Inside, nstate is last and state_first first: the states 1 .. first-1 are
    // drawn (simulator.cpp) and not calculated, so state s is the same state as in a run of [1, last].
    long stateFirst = 0, stateLast = 0;
    const bool sampled = readStateRange(configuration,&stateFirst,&stateLast);
    int nstate = sampled ? int(stateLast) : 0;
    Config_setNstate(cnf,nstate);
    cnf->state_first = sampled ? int(stateFirst) : 1;

    // configuration.seed (A94): a number, the seed of the first configuration, configuration k of k_range taking seed +
    // (k - first) as a generated seed does (A73); or a list of one seed per configuration (k_range last - first + 1 of
    // them, one without k_range). Until A94 one number was the seed of every configuration (the converter writes such an
    // input as a list of the same number). A73: a generated seed is generated once, for the first configuration, and
    // configuration k takes generated + (k - first): a different, reproducible sequence of states for every k.
    cJSON* seedItem=configuration ? cJSON_GetObjectItem(configuration,"seed") : NULL;
    int requestedSeed=-1;
    bool seedListed=false;
    if(cJSON_IsArray(seedItem)) {
        const long count=ranged ? rangeLast-rangeFirst+1 : 1;
        if(cJSON_GetArraySize(seedItem)!=count) {
            ccexRefusalPrint(stderr,"Error: configuration.seed as a list has one seed per configuration, %ld here (%s), it has %d\n",
                             count,ranged ? "configuration.k_range last - first + 1" : "one configuration without k_range",
                             cJSON_GetArraySize(seedItem));
            ccexRefuse();
        }
        cJSON* entry=cJSON_GetArrayItem(seedItem,int(configurationIndex-rangeFirst));
        if(!cJSON_IsNumber(entry)) {
            ccexRefusalPrint(stderr,"Error: configuration.seed[%ld] must be an integer\n",configurationIndex-rangeFirst);
            ccexRefuse();
        }
        requestedSeed=integerValue(entry,"configuration.seed[]");
        seedListed=true;
    } else if(seedItem) {
        refuseWrongType(seedItem,"configuration.seed",cJSON_IsNumber(seedItem),"an integer or a list of one integer per configuration");
        requestedSeed=integerValue(seedItem,"configuration.seed");
    }
    static int generatedSeed=-1;
    const bool reuse=requestedSeed==-1 && generatedSeed!=-1 && !seedListed;
    const int base=reuse ? generatedSeed : resolveMpiSeed(requestedSeed,MPI_COMM_WORLD);
    if(requestedSeed==-1 && !seedListed) generatedSeed=base;
    const int seed=!seedListed && ranged ? int(unsigned(base)+unsigned(configurationIndex-rangeFirst)) : base;
    // A76: the seed is reported only for a run that draws from it: sampled bath states (nstate > 0), the axes and onsite-spin
    // states of Defect species (drawn when their files are absent), the pCCE partition (its trials and its cache key). An
    // ensemble run of a bath without Defect species draws nothing.
    cJSON* defectBlock=cJSON_GetObjectItem(root,"Defect");
    cJSON* defectSpecies=cJSON_IsObject(defectBlock) ? cJSON_GetObjectItem(defectBlock,"species") : NULL;
    cnf->seed_used=nstate>0 || methodUsesPcceClusters(method) ||
                   (cJSON_IsArray(defectSpecies) && cJSON_GetArraySize(defectSpecies)>0);
    if(requestedSeed==-1 && !ranged && rank==0 && cnf->seed_used)
        fprintf(stderr,"[seed] generated on rank 0 and shared with all ranks: %d, set seed=%d to replay this run.\n",seed,seed);
    if(requestedSeed==-1 && ranged && rank==0 && cnf->seed_used)
        fprintf(stderr,"[seed] configuration k=%ld: seed %d (generated %d + k - %ld), set seed=%d to replay this configuration alone.\n",
                configurationIndex,seed,base,rangeFirst,seed);
    Config_setSeed(cnf,seed);

    ////////////////////////////////////////////////////////////////////////
    // Qubit and BathFiles
    ////////////////////////////////////////////////////////////////////////
    if (rank==0 && verbosity){
        printf("\n");
        printf("      - File-related keys :\n");
        printf("        [ configuration { path, k_range, gyrofile, bathfile, bath_shift, bath_state_file,\n");
        printf("          defect_axis_file, defect_state_file, s_range, seed }, frame.rotation ]\n\n");
    }

    // the configuration block (A64): gyrofile, bathfile, bath_shift and bath_state_file, the names in configuration.path
    std::string gyrofile;
    if (configurationFile(configuration,"gyrofile",folder,ranged,configurationIndex,&gyrofile)) {
        Config_allocGyrofile(cnf);
        Config_setGyrofile(cnf,const_cast<char*>(gyrofile.c_str()));
    }

    // bathfiles
    int length = 0;
    char** bathfiles = configuration ? cJSON_ReadFilePath1d(&length, configuration, "bathfile", true, NULL) : NULL; // get length
    if (bathfiles != NULL) {
        Config_setNbathfiles(cnf,length);
        Config_allocBathfiles(cnf);
        Config_allocBathadjust(cnf);
        for (int i = 0; i < length; i++){
            const std::string name = configurationPath(folder,bathfiles[i],ranged,configurationIndex);
            refuseUnreadable(shownKey(configuration,"bathfile") + "[" + std::to_string(i) + "]",name);   // A94
            Config_setBathfiles_i(cnf,const_cast<char*>(name.c_str()),i);
        }
        freeChar2d(&bathfiles,length);
    }

    // bath_shift (bathadjust until A47)
    double** bathadjustDefault = allocDouble2d(length,3); // all elements are 0.0
    double** bathadjust = configuration ? cJSON_ReadDouble2d(configuration,(char*)"bath_shift",true,bathadjustDefault,length,3)
                                        : allocDouble2d(length,3);

    for (int i = 0; i < length; i++){
        Config_setBathadjust_i(cnf,bathadjust[i],i);
    }

    freeDouble2d(&bathadjustDefault,length);
    freeDouble2d(&bathadjust,length);

    ////////////////////////////////////////////////////////////////////////
    // Coordinate frame rotation
    ////////////////////////////////////////////////////////////////////////
    // frame.rotation (A65; coordinate_frame_rotation before, refused): an absent block leaves Config at the values set by
    // Config_init -- rotation off, both axes on +z -- so the run is bit-for-bit the legacy one. A given block is on unless
    // it says "enabled": false.
    const char* ROTKEY = "frame.rotation";
    cJSON* RotSection    = cJSON_GetObjectItem(frame,"rotation");

    if (RotSection != NULL){

        // A non-object -- "coordinate_frame_rotation": true, or a path string -- has no
        // children, so every sub-key below would read as missing and the run would go on
        // silently WITHOUT the rotation the input clearly asked for.
        if (!cJSON_IsObject(RotSection)){
            ccexRefusalPrint(stderr,"Error: %s must be an object, e.g.\n",ROTKEY);
            ccexRefusalPrint(stderr,"  \"%s\" : { \"qubit_axis\" : [1,1,1] }\n",ROTKEY);
            ccexRefuse();
        }

        bool rot_enabled = readStrictBool(RotSection,ROTKEY,"enabled",true);
        Config_setRot_enabled(cnf,rot_enabled);

        if (rot_enabled){
            // A90: qubit_axis only (bath_axis, reference_qubit and qubit_position_frame were removed; the frames of the
            // tensors are the "frame" keys beside them: refuseOldInputKeys)
            double qubit_axis[3] = {0.0,0.0,1.0};
            readAxis3(RotSection,ROTKEY,"qubit_axis",qubit_axis,false);
            Config_setRot_qubit_axis(cnf,qubit_axis);
        }
    }
    frameRotationOn = Config_getRot_enabled(cnf);   // A90: the frame keys of this reader and the Qubit and Defect readers

    // configuration.defect_axis_file and defect_state_file (avaaxfile, exstatefile until A47, then Defect.axis_file and
    // Defect.state_file until A64)
    // one name per bath file each (A64), in the order of configuration.bathfile
    const int nbath = Config_getNbathfiles(cnf);
    const struct { const char* key; void (*alloc)(Config*); void (*set)(Config*,char*,int); } perBath[] = {
        {"defect_axis_file",  Config_allocAvaaxfile,   Config_setAvaaxfile_i},
        {"bath_state_file",   Config_allocStatefile,   Config_setStatefile_i},
        {"defect_state_file", Config_allocExstatefile, Config_setExstatefile_i}};
    const bool defects = inputHasDefectSpecies(root);
    for (const auto& f : perBath) {
        const std::vector<std::string> names = perBathFileNames(configuration,f.key,nbath,folder,ranged,configurationIndex);
        if (names.empty()) continue;
        // A94: the files the run reads have to open (the states of s_range, the axes with Defect species); a key the
        // run does not read is reported (G3) and left unread
        const bool stateFile = strcmp(f.key,"defect_axis_file") != 0;
        const bool defectFile = strcmp(f.key,"bath_state_file") != 0;
        const std::string shown = shownKey(configuration,f.key);
        if (stateFile && !sampled){
            warnNotUsed(shown.c_str(),"the ensemble (no configuration.s_range) reads no state files");
            continue;
        }
        if (defectFile && !defects){
            warnNotUsed(shown.c_str(),"the input has no Defect species");
            continue;
        }
        for (int i = 0; i < nbath; i++){
            if (names[i].empty()) continue;           // "": this bath file's states or axes are drawn from the seed
            const std::string key = shown + "[" + std::to_string(i) + "]";
            if (!stateFile) refuseUnreadable(key,names[i]);
            else for (long st = stateFirst; st <= stateLast; ++st) refuseUnreadable(key,stateFileName(names[i],st));
        }
        f.alloc(cnf);
        for (int i = 0; i < nbath; i++) f.set(cnf,const_cast<char*>(names[i].c_str()),i);
    }

    ////////////////////////////////////////////////////////////////////////
    // The hamiltonian block (A47): common to both methods, and the terms of one method
    ////////////////////////////////////////////////////////////////////////
    if (rank==0 && verbosity){
        printf("\n");
        printf("      - Hamiltonian-related keys :\n");
        printf("        hamiltonian.common : [ qubit_bath_mean_field, bath_bath_mean_field,\n");
        printf("          hyperfine { dipole, contact, file, missing_row },\n");
        printf("          quadrupole (\"off\" or { structure, file, file_without_defect, cell }) ]\n");
        printf("        hamiltonian.cCCE : [ knight, hfmedi, delta_A ]\n");
        printf("        hamiltonian.meCCE / megCCE : [ jump_operators ]\n\n");
    }
    cJSON* hamiltonian = readBlock(root,"hamiltonian");
    // A94 (G1) named the sections by role, corrections and dissipation; reverted 2026-10-17, the sections stay by method
    if (hamiltonian){
        static const char* const renamed[][2] = {{"corrections", "\"hamiltonian.cCCE\" (knight, hfmedi, delta_A; meCCE reads it too)"},
                                                 {"dissipation", "\"hamiltonian.meCCE\" or \"hamiltonian.megCCE\" (jump_operators, "
                                                                 "the section of the run's method)"}};
        for (const auto& o : renamed)
            if (untrackedItem(hamiltonian,o[0])){
                ccexRefusalPrint(stderr,"Error: hamiltonian.%s is now %s\n%s\n",o[0],o[1],kConvertHint);
                ccexRefuse();
            }
    }
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
    // the master-equation methods: jump_operators only; which method the section belongs to is checked once every
    // reader has run (main.cpp, validateHamiltonianForMethod)
    for (const char* name : {"meCCE", "megCCE"})
        if (cJSON* section = hamiltonian ? readBlock(hamiltonian,name) : NULL) readJumpOperators(cnf, section, name);
    // a cCCE setting in a gCCE run is refused once every reader has run (main.cpp, validateHamiltonianForMethod)
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
        if (nstate == 0)   // A94 (G2, G3)
            warnNotUsed("hamiltonian.common.bath_bath_mean_field","the ensemble (no configuration.s_range) has no sampled bath "
                        "state to make a static field from");
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

    // hamiltonian.common.hyperfine.frame (frame.rotation.hf_tensor_frame until A90): the frame of the tensor COMPONENTS
    // in hyperfine.file (its position columns are always read in the file frame). Needed with "dipole": "dft" only; with
    // "point" the tensors are point dipoles from the positions, turned with them, and the contact term has no direction.
    if (hyperfine){
        const FrameValue f = (frameRotationOn && hf_readmode >= 2)
            ? frameOf(hyperfine,"hamiltonian.common.hyperfine","hamiltonian.common.hyperfine.frame",
                      "the tensor components of hamiltonian.common.hyperfine.file")
            : readFrameKey(hyperfine,"hamiltonian.common.hyperfine");
        warnFrameUnused(f,"hamiltonian.common.hyperfine.frame");
        if (f != FRAME_ABSENT) Config_setHf_frame(cnf,f == FRAME_FILE ? "file" : "qubit");
        if (f != FRAME_ABSENT && frameRotationOn && hf_readmode < 2 && rank==0){
            printMessage("  Note: hamiltonian.common.hyperfine.frame is not used with dipole \"point\": the tensors are point");
            printMessage("        dipoles from the positions, turned with them, and the contact term has no direction.\n");
        }
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
    // hamiltonian.common.quadrupole.frame (frame.rotation.qd_tensor_frame until A90): the frame of the tensor components
    // in quadrupole.file (and file_without_defect)
    if (qd_readmode != 0){
        const FrameValue f = frameOf(quadrupole,"hamiltonian.common.quadrupole","hamiltonian.common.quadrupole.frame",
                                     "the tensor components of hamiltonian.common.quadrupole.file");
        if (f != FRAME_ABSENT) Config_setQd_frame(cnf,f == FRAME_FILE ? "file" : "qubit");
    }

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
    // hamiltonian.cCCE: knight, hfmedi, and delta_A (the second-order delta_A(m_s).I of Maze 2008 eq. A3, A83)
    ////////////////////////////////////////////////////////////////////////
    Config_setHfmedi(cnf, cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","hfmedi",false) : false);
    Config_setKnight(cnf, cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","knight",false) : false);
    cnf->delta_A = cCCE ? readStrictBool(cCCE,"hamiltonian.cCCE","delta_A",false) : false;

    cJSON_Delete(root);    
    freeChar1d(&data);
}

// Qubit.intmap: the qubit-qubit tensors, a self entry (the same name twice) is the zero-field splitting. Both forms of
// the Qubit block read it (A93, 2026-10-16: Qubit.zfs of the Qubit.file form, fixed kHz, and its Qubit.frame are now
// a self entry of the Qubit.file qubit, named by Qubit.name).
static void readQubitIntmap(cJSON* QubitSection, QubitArray* qa){
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
        MatrixXcd tensor = cJSON_ReadTensor(intmap,"tensor",true,MatrixXcd::Zero(3,3));
        // A94 (R4): without "unit" the tensor is kHz (a Defect zfs is MHz); the log says so for a tensor that is not zero
        if (!cJSON_GetObjectItem(intmap,"unit") && !tensor.isZero(0.0) && rank==0)
            printf("      Note: Qubit.intmap[\"%s\",\"%s\"].tensor has no unit and is read in kHz (a Defect zfs is MHz), "
                   "\"unit\" names it\n",qubit1_name,qubit2_name);
        tensor *= readIntmapFrequencyScaleToKHz(intmap,qubit1_name,qubit2_name);
        tensor = KHZ_TO_RADKHZ(tensor);

        // Which frame the components are written in (frame, tensor_frame until A90). Recorded, never guessed: a ZFS
        // tensor and a pair tensor look the same, and so does one that has already been turned. Required with the
        // rotation on (A90; until then a single qubit's entry without it was read as written).
        const std::string where = std::string("Qubit.intmap[\"") + qubit1_name + "\",\"" + qubit2_name + "\"]";
        const FrameValue f = frameOf(intmap,where,"Qubit.intmap[].frame","its tensor");
        const IntmapProvenance frame = f == FRAME_FILE ? INTMAP_EXPLICIT_FILE : f == FRAME_QUBIT ? INTMAP_EXPLICIT_QUBIT
                                                                                  : INTMAP_EXPLICIT_UNSPECIFIED;

        // set Interaction map
        QubitArray_setIntmap_i_j(qa,tensor,qubit1_idx,qubit2_idx);
        QubitArray_setIntmap_i_j_frame(qa,frame,qubit1_idx,qubit2_idx);
    }
}

// A93 (2026-10-16): the qubit states are Qubit.alpha_state and Qubit.beta_state only, each exactly one of
//   {"ms": [m_1, m_2, ...]}  one m_s per qubit in the order of Qubit.qubit ([m] for the qubit of Qubit.file): the product
//                            of the |m_s> states (Qubit.qubit[].alphams/betams, Qubit.alphams/betams with Qubit.file
//                            until A93)
//   {"index": n}             the n-th eigenstate of the qubit Hamiltonian by descending energy (not with Qubit.file)
//   {"vector": [...]}        the state itself (not with Qubit.file)
// Absent, a side is the product of m_s +1 (alpha) or 0 (beta) on every qubit, as before. The m_s of each qubit (used by
// cCCE, the mean fields, knight/hfmedi/delta_A, DEER, ...) are the "ms" list, or the defaults +1 and 0 when the side is
// an index or a vector (as before, when alphams/betams were not written beside one).
// A94 (G13): a state vector, dim components, each a number (real) or [re, im]
static MatrixXcd readStateVector(cJSON* parent, const char* key, const std::string& shown, int dim){
    cJSON* item = cJSON_GetObjectItem(parent,key);
    if (!cJSON_IsArray(item) || cJSON_GetArraySize(item) != dim){
        ccexRefusalPrint(stderr,"Error: %s must be a list of %d components (the qubit dimension), each a number or [re, im]\n",
                         shown.c_str(),dim);
        ccexRefuse();
    }
    MatrixXcd vector = MatrixXcd::Zero(dim,1);
    for (int i = 0; i < dim; ++i){
        cJSON* c = cJSON_GetArrayItem(item,i);
        const std::string label = shown + "[" + std::to_string(i) + "]";
        if (cJSON_IsArray(c)){
            if (cJSON_GetArraySize(c) != 2){
                ccexRefusalPrint(stderr,"Error: %s must be a number or [re, im]\n",label.c_str());
                ccexRefuse();
            }
            vector(i,0) = doublec(readFiniteNumber(cJSON_GetArrayItem(c,0),(label + "[0]").c_str()),
                                  readFiniteNumber(cJSON_GetArrayItem(c,1),(label + "[1]").c_str()));
        }else{
            vector(i,0) = doublec(readFiniteNumber(c,label.c_str()),0.0);
        }
    }
    return vector;
}

static void readQubitStates(cJSON* QubitSection, QubitArray* qa){
    const int nqubit = QubitArray_getNqubit(qa);
    const int dim = QubitArray_dim(qa);
    const char* sides[] = {"alpha_state", "beta_state"};
    const float msDefault[] = {1.0f, 0.0f};
    for (int side=0; side<2; ++side){
        std::vector<float> ms(nqubit, msDefault[side]);
        MatrixXcd vector;                   // the 0 dimension matrix: not given
        cJSON* state = readBlock(QubitSection,sides[side]);
        if (state){
            cJSON* msItem = cJSON_GetObjectItem(state,"ms");
            cJSON* indexItem = cJSON_GetObjectItem(state,"index");
            cJSON* vectorItem = cJSON_GetObjectItem(state,"vector");
            const int given = (msItem != NULL) + (indexItem != NULL) + (vectorItem != NULL);
            if (given != 1){
                ccexRefusalPrint(stderr,"Error: Qubit.%s must have exactly one of \"ms\", \"index\" and \"vector\" (%s)\n",
                                 sides[side], given == 0 ? "it has none" : "it has more than one");
                ccexRefuse();
            }
            if (msItem){
                const bool listOk = cJSON_IsArray(msItem) && cJSON_GetArraySize(msItem) == nqubit;
                if (!listOk){
                    ccexRefusalPrint(stderr,"Error: Qubit.%s.ms must be a list of one m_s per qubit, %d numbers in the "
                                     "order of Qubit.qubit%s\n", sides[side], nqubit,
                                     cJSON_IsArray(msItem) ? (std::string(" (it has ") + std::to_string(cJSON_GetArraySize(msItem))
                                                              + ")").c_str() : " (a list also for one qubit)");
                    ccexRefuse();
                }
                for (int i=0; i<nqubit; ++i){
                    const std::string label = std::string("Qubit.") + sides[side] + ".ms[" + std::to_string(i) + "]";
                    ms[i] = (float)readFiniteNumber(cJSON_GetArrayItem(msItem,i),label.c_str());
                }
            }else if (indexItem){
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
                vector = readStateVector(state,"vector",std::string("Qubit.") + sides[side] + ".vector",dim);   // G13: [re, im] too
            }
        }
        // with {"index"} or {"vector"} nothing names a level per qubit: a placeholder level is stored, the default, or,
        // where the default is no level of the qubit's spin (a half-integer spin), +1/2 (alpha) and -1/2 (beta).
        // pulse.per_qubit, which turns each qubit between these levels, refuses {"index"} and {"vector"} (A94).
        if (state && !cJSON_GetObjectItem(state,"ms"))
            for (int i=0; i<nqubit; ++i){
                const double S = QubitArray_getQubit_i_spin(qa,i);
                const double m = ms[i];
                const bool level = std::fabs(m) <= S + 1e-9 && std::fabs((S - m) - std::round(S - m)) < 1e-9;
                if (!level) ms[i] = side == 0 ? 0.5f : -0.5f;
            }
        for (int i=0; i<nqubit; ++i){
            if (side == 0) QubitArray_setQubit_i_alpha_fromMs(qa,ms[i],i);
            else           QubitArray_setQubit_i_beta_fromMs(qa,ms[i],i);
        }
        if (side == 0) QubitArray_setPsia(qa,vector);
        else           QubitArray_setPsib(qa,vector);
    }
}

/**
 * @brief Read the option from the input file
 * @details Read &Qubit tag  options
*/
void cJSON_readOptionQubitArray(QubitArray* qa, char* fccein){

    if (rank==0 && verbosity){
        printf("\n");
        printf("    >> Read Qubit Options\n\n");
        printf("      - Read values of main-key 'Qubit'\n");
        printf("        sub-key of 'Qubit' : [ qubit, intmap, alpha_state, beta_state, psi0 ]\n\n");
        printf("      - Read values of sub-key 'qubit'\n");
        printf("        sub-sub-key : [ name, spin, gyro, xyz or file, detuning ]\n\n");
        printf("      - Read values of sub-key 'intmap'\n");
        printf("        sub-sub-key : [ between, tensor, unit, frame ]\n\n");
        printf("      - Read hamiltonian.common.qubit_bath_mean_field\n");
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

    char nameDefault[MAX_CHARARRAY_LENGTH] = "q";
    float spinDefault = 1.0;
    double gyroDefault = GAMMA_ELECTRON;
    double detuningDefault = 0.0;

    MatrixXcd psiDefault;
    // the Overhauser field of the bath on the qubits (Qubit.overhaus until A47); both forms of the Qubit block read it
    cJSON* hamiltonian = readBlock(root,"hamiltonian");
    cJSON* common = hamiltonian ? readBlock(hamiltonian,"common") : NULL;
    const bool overhaus = common ? readStrictBool(common,"hamiltonian.common","qubit_bath_mean_field",false) : false;
    QubitArray_setOverhaus(qa,overhaus);

    // the qubits are the entries of Qubit.qubit (A94, G12: the Qubit.file form, one qubit from a defect file, is an entry
    // {"file": ...} in place of "xyz"; its old keys are refused by refuseOldInputKeys)
    cJSON* QubitSection = readBlock(root,"Qubit");
    {
        int iqubit = 0;

        ////////////////////////////////////////////////////////////////////////
        // set nqubit: the length of Qubit.qubit (Qubit.nqubit until A47)
        ////////////////////////////////////////////////////////////////////////
        cJSON* qubitArray = QubitSection ? cJSON_GetObjectItem(QubitSection,"qubit") : NULL;
        if (!cJSON_IsArray(qubitArray) || cJSON_GetArraySize(qubitArray) < 1){
            ccexRefusalPrint(stderr, "Error: Qubit.qubit must be an array of one entry per qubit, each with \"xyz\" or \"file\" "
                            "(a defect file whose first line is the position)\n");
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
                ccexRefusalPrint(stderr, "Error: The length of \"qubit\" is different from nqubit(%d)\n", nqubit);
                ccexRefuse();
            }
            // read qubit properties from the input file
            snprintf(nameDefault,MAX_CHARARRAY_LENGTH,"q%d",iqubit);
            char* name = cJSON_ReadString(qubit,"name",true,nameDefault);
            float spin = cJSON_ReadFloat(qubit,"spin",true,spinDefault);
            double gyro = readQubitGyro(qubit,iqubit,gyroDefault);
            // the position: "xyz", or the first line of "file" (A94, G12)
            double* xyz = NULL;
            if (cJSON_GetObjectItem(qubit,"file")){
                if (cJSON_GetObjectItem(qubit,"xyz")){
                    ccexRefusalPrint(stderr,"Error: Qubit.qubit[%d] has both \"file\" and \"xyz\", the position is one of them\n",iqubit);
                    ccexRefuse();
                }
                const char* path = cJSON_ReadString(qubit,(char*)"file",false,NULL);   // opened below, refused with its index
                xyz = allocDouble1d(3);
                readQubitPositionFile(path,xyz,iqubit);
            }else{
                xyz = cJSON_ReadDouble1d(qubit,"xyz",false,NULL,3);
            }
            // detuning: a number, or the one-element array [x] the Defect form suggests
            // (13 tutorial and check_numerics inputs write [0.0]); it used to be ignored for
            // the default, so [x] with x != 0 was silently 0 before 2026-10-02 (D18)
            // A94 (G21): a number only ([x] was a second form of the same value)
            cJSON* detuningItem = cJSON_GetObjectItem(qubit,"detuning");
            if (cJSON_IsArray(detuningItem)){
                ccexRefusalPrint(stderr,"Error: Qubit.qubit[%d].detuning is a number (kHz) or {\"value\": ..., \"unit\": ...}, "
                                 "not a list: write [x] as x\n%s\n",iqubit,kConvertHint);
                ccexRefuse();
            }
            double detuning = cJSON_ReadDouble(qubit,"detuning",true,detuningDefault); //kHz

            // A94 (R4): a plain number is kHz here and MHz in Defect.species[].detuning; the log says which was read
            if (detuning != 0.0 && rank==0 && !unitWritten().count("Qubit.qubit[" + std::to_string(iqubit) + "].detuning"))
                printf("      Note: Qubit.qubit[%d].detuning %g is read in kHz (= %g MHz; a Defect detuning is MHz), "
                       "{\"value\": ..., \"unit\": ...} names the unit\n",iqubit,detuning,detuning*1e-3);

            // set qubit properties
            QubitArray_setQubit_i_name(qa,name,iqubit);
            QubitArray_setQubit_i_spin(qa,spin,iqubit);
            QubitArray_setQubit_i_xyz(qa,xyz,iqubit);
            QubitArray_setQubit_i_gyro(qa,gyro,iqubit);
            QubitArray_setQubit_i_detuning(qa,KHZ_TO_RADKHZ(detuning),iqubit); // radkHz
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
        readQubitIntmap(QubitSection,qa);

        ////////////////////////////////////////////////////////////////////////
        //  psia, psib, psi0
        ////////////////////////////////////////////////////////////////////////
        // Qubit.alpha_state and beta_state (alphaidx/betaidx and psia/psib until A47; the qubits' alphams/betams until A93)
        readQubitStates(QubitSection,qa);
        int dim = QubitArray_dim(qa);
        if (cJSON_GetObjectItem(QubitSection,"psi0")){
            QubitArray_setPsi0(qa,readStateVector(QubitSection,"psi0","Qubit.psi0",dim));
            // A94 (G3): cCCE and meCCE take the qubit's two levels only
            char method[MAX_CHARARRAY_LENGTH], mode[16];
            readMethod(root,method,mode);
            const char* kernel = methodKernel(method);
            if (!strcmp(kernel,"cce") || !strcmp(kernel,"mecce"))
                warnNotUsed("Qubit.psi0","cCCE and meCCE use alpha_state and beta_state only (the qubit is projected out)");
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

    if (rank==0 && verbosity){
        printf("\n");
        printf("    >> Read Cluster Options\n");
        printf("      [ order, method, cluster.mode, cluster.selection { add_subclusters, max_per_order } ]\n\n");
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

    // cluster.selection {add_subclusters, max_per_order} (A65; addsubclus and nk until A47, then top-level add_subclusters
    // and max_clusters_per_order)
    cJSON* clusterBlock = readBlockOrEmpty(root,"cluster");
    cJSON* selection = readBlockOrEmpty(clusterBlock,"selection");
    bool addsubclus = true;
    if (cJSON* item = cJSON_GetObjectItem(selection,"add_subclusters")){
        if (!cJSON_IsBool(item)){
            ccexRefusalPrint(stderr,"Error: cluster.selection.add_subclusters must be true or false\n");
            ccexRefuse();
        }
        addsubclus = cJSON_IsTrue(item);
    }
    Cluster_setAddsubclus(clus,addsubclus);

    Cluster_allocNk(clus);
    int* nk = cJSON_ReadInt1d(selection,(char*)"max_per_order",true,NULL,order+1);
    if (nk!=NULL){
        Cluster_setNk(clus,nk);
        // A89: only entries 3 and up cut anything (makeHashClusterOn, standard and pcce). Entries 0-2 are never compared
        // (1 and 2 are overwritten by the counts), and dsj, itb, dsjitb never read the list. A nonzero entry that is not
        // applied is reported once.
        const char* m = Cluster_getMethod(clus);
        const bool noneApplied = !strcasecmp(m,"dsj") || !strcasecmp(m,"itb") || !strcasecmp(m,"dsjitb");
        std::string unused;
        for (int i=0; i<=order; ++i)
            if (nk[i] != 0 && (i < 3 || noneApplied)) unused += (unused.empty() ? "" : ", ") + std::to_string(i);
        if (!unused.empty() && strcmp(clusterMode,"adaptive") != 0){   // bundle refuses a nonzero entry
            const bool one = unused.find(',') == std::string::npos;
            const std::string why = std::string(one ? "entry " : "entries ") + unused + (one ? " is" : " are") + " ignored, "
                                  + (noneApplied ? "this cluster.mode keeps every cluster" : "only orders 3 and up are cut");
            warnNotUsed("cluster.selection.max_per_order",why.c_str());   // A94 (G3): one wording
        }
    }
    freeInt1d(&nk);
    // all element would be "0" in which we consider all pairs within rdip


    if (methodUsesPcceClusters(Cluster_getMethod(clus))){

        if (rank==0 && verbosity){
            printf("      cluster.pcce : [ partition_size, trials, max_iterations, kmeans_pp, print_trials ]\n\n");
        }

        // Do not dereference a null default argument even when partition_size is supplied.
        // Zero or fractional sizes also make the partition divisions invalid.
        cJSON* pcce = readBlock(clusterBlock,"pcce");
        const char* keys[] = {"partition_size", "trials", "max_iterations"};
        int values[] = {0, 30000, 30000};
        for (int k=0; k<3; ++k){
            cJSON* item = pcce ? cJSON_GetObjectItem(pcce,keys[k]) : NULL;
            if (!item){
                if (k==0){
                    ccexRefusalPrint(stderr,"Error: cluster.mode pcce requires cluster.pcce.partition_size\n");
                    ccexRefuse();
                }
                continue;
            }
            const std::string name = "cluster.pcce." + std::string(keys[k]);
            double value = readFiniteNumber(item,name.c_str());
            if (value < 1 || value > 2147483647.0 || std::floor(value) != value){
                ccexRefusalPrint(stderr,"Error: %s must be a positive 32-bit integer\n",name.c_str());
                ccexRefuse();
            }
            values[k] = (int)value;
        }
        int sK = values[0], max_trial = values[1], max_iter = values[2];
        bool kmeans_pp   = pcce ? readStrictBool(pcce,"cluster.pcce","kmeans_pp",true) : true;
        // print_trials N > 0 (A94, read with the limit in cJSON_readOptionConfig): the first N trials of rank 0
        cJSON* printTrials = pcce ? cJSON_GetObjectItem(pcce,"print_trials") : NULL;
        bool iter_detail = cJSON_IsNumber(printTrials) && printTrials->valuedouble > 0;

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
static double readSequenceFraction(cJSON* item, int row, const char* key = "sequence") {
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
                ccexRefusalPrint(stderr, "Error: %s[%d] fraction \"%s\" must be a/b with integers and b != 0\n", key, row, text);
                ccexRefuse();
            }
            value = double(num)/double(den);
        } else {
            errno = 0; value = strtod(text, &end);
            if (end == text || *end != '\0' || errno == ERANGE) {
                ccexRefusalPrint(stderr, "Error: %s[%d] fraction \"%s\" is not a number\n", key, row, text);
                ccexRefuse();
            }
        }
    } else {
        ccexRefusalPrint(stderr, "Error: %s[%d] fraction must be a number or a string such as \"1/2\"\n", key, row);
        ccexRefuse();
    }
    if (!std::isfinite(value) || value < 0) {
        ccexRefusalPrint(stderr, "Error: %s[%d] fraction must be finite and not negative\n", key, row);
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

    if (rank == 0 && verbosity) {
        printf("\n");
        printf("    >> Read Pulse Options\n");
        printf("      pulse : [ count, name, sequence, per_qubit, duration, detuning_factor ]\n\n");
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
    // A94 (user, 2026-10-16): pulse.per_qubit turns each qubit between its own two levels, the m_s of Qubit.alpha_state and
    // beta_state {"ms": [...]}. With {"index"} or {"vector"} no level per qubit is given (A93 took a default for it).
    if (pulseiter) {
        cJSON* qubitBlock = untrackedItem(root,"Qubit");
        for (const char* side : {"alpha_state", "beta_state"}) {
            cJSON* state = untrackedItem(qubitBlock,side);
            if (cJSON_IsObject(state) && !untrackedItem(state,"ms")) {
                ccexRefusalPrint(stderr,"Error: pulse.per_qubit true needs Qubit.alpha_state and Qubit.beta_state as {\"ms\": "
                                 "[m_1, m_2, ...]} (the two levels each qubit is turned between), Qubit.%s is {\"%s\"}\n",
                                 side, untrackedItem(state,"index") ? "index" : "vector");
                ccexRefuse();
            }
        }
    }

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

    // A94 (R1): pulse.name is one of five names, and a name that fixes the number of pulses has to agree with pulse.count
    // (an unknown name was kept, and the count decided without a word: {"count": 1, "name": "Hahn"} ran a Hahn echo,
    // {"count": 0, "name": "CPMG"} a Ramsey). Checked before the arrays (count + 1 rows) are allocated: WAHUHA always
    // writes five rows (A86).
    if (cJSON_GetObjectItem(block,"name")) {
        static const char* const names[] = {"Ramsey", "HahnEcho", "CPMG", "WAHUHA", "Manual"};
        const char* known = nullptr;
        for (const char* n : names) if (!strcasecmp(pulse->pulsename, n)) known = n;
        if (!known) {
            ccexRefusalPrint(stderr, "Error: pulse.name must be Ramsey, HahnEcho, CPMG, WAHUHA or Manual (not \"%s\"), or left "
                             "out: pulse.count alone gives Ramsey (0), HahnEcho (1) or CPMG (more)\n", pulse->pulsename);
            ccexRefuse();
        }
        const char* needs = !strcmp(known,"Ramsey") && npulse != 0   ? "0"
                          : !strcmp(known,"HahnEcho") && npulse != 1 ? "1"
                          : !strcmp(known,"CPMG") && npulse <= 1     ? "2 or more"
                          : !strcmp(known,"WAHUHA") && npulse != 4   ? "4 (the four-pulse sequence)" : nullptr;
        if (needs) {
            ccexRefusalPrint(stderr, "Error: pulse.name %s needs pulse.count %s (this input has %d)\n", known, needs, npulse);
            ccexRefuse();
        }
    }
    // a sequence is read only for Manual, with another name it left the arrays zero (A86)
    cJSON* sequence_array = cJSON_GetObjectItem(block, "sequence");
    if (sequence_array != NULL && strcasecmp(pulse->pulsename, "Manual") != 0) {
        ccexRefusalPrint(stderr, "Error: pulse.sequence is read only with pulse.name Manual (pulse.name %s), write "
                         "\"name\": \"Manual\" or leave the sequence out\n", pulse->pulsename);
        ccexRefuse();
    }

    Pulse_allocSequence(pulse);
    Pulse_allocAxes(pulse);
    Pulse_allocAngles(pulse);
    Pulse_allocSequenceIndices(pulse);

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
            ccexRefusalPrint(stderr, "Error: pulse.sequence has %d rows, pulse.count + 1 = %d (one row per pulse and the last "
                             "free interval)\n", count, (pulse->npulse) + 1);
            ccexRefuse();
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
        // A94 (G6): spin was the S of the driven spins written again (refused when it differed), it is taken from them;
        // alphams and betams are one key, levels
        if (untrackedItem(bath,"spin")) {
            ccexRefusalPrint(stderr,"Error: pulse.bath.spin was removed: the spin of the driven bath spins (pulse.bath.defect) is "
                             "used, remove it\n%s\n",kConvertHint);
            ccexRefuse();
        }
        if (untrackedItem(bath,"alphams") || untrackedItem(bath,"betams")) {
            ccexRefusalPrint(stderr,"Error: pulse.bath.alphams and betams are now pulse.bath.levels: [alphams, betams] (default "
                             "[0.5, -0.5])\n%s\n",kConvertHint);
            ccexRefuse();
        }
        for (cJSON* key = bath->child; key; key = key->next)
            if (strcmp(key->string,"count") && strcmp(key->string,"sequence") && strcmp(key->string,"defect") &&
                strcmp(key->string,"energy_shift") && strcmp(key->string,"bandwidth") && strcmp(key->string,"levels")){
                ccexRefusalPrint(stderr,"Error: pulse.bath has no key \"%s\" (count, sequence, defect, energy_shift, bandwidth, "
                        "levels)\n", key->string);
                ccexRefuse();
            }
        cJSON* jseq = cJSON_GetObjectItem(bath, "sequence");     // every key looked up, so the input check knows them
        cJSON* jdef = cJSON_GetObjectItem(bath, "defect");
        cJSON* jsh  = cJSON_GetObjectItem(bath, "energy_shift");
        cJSON* jlv  = cJSON_GetObjectItem(bath, "levels");
        pulse->bnpulse = cJSON_ReadInt(bath, (char*)"count", true, 0);    // type and integer checks as pulse.count
        if (pulse->bnpulse < 0) { ccexRefusalPrint(stderr, "Error: pulse.bath.count must be >= 0\n"); ccexRefuse(); }
        if (jlv) {
            if (!cJSON_IsArray(jlv) || cJSON_GetArraySize(jlv) != 2) {
                ccexRefusalPrint(stderr, "Error: pulse.bath.levels must be [m_a, m_b], the two levels the bath pulse turns "
                                 "between\n");
                ccexRefuse();
            }
            pulse->balphams = readFiniteNumber(cJSON_GetArrayItem(jlv,0), "pulse.bath.levels[0]");
            pulse->bbetams  = readFiniteNumber(cJSON_GetArrayItem(jlv,1), "pulse.bath.levels[1]");
        }
        // A67: the half width of the window around a tone (MHz); absent: kBathPulseBandwidthRadkHz itself, not a round trip
        if (cJSON* jbw = cJSON_GetObjectItem(bath, "bandwidth")) {
            const double mhz = readFiniteNumber(jbw, "pulse.bath.bandwidth");
            if (!(mhz > 0.0)) {
                ccexRefusalPrint(stderr, "Error: pulse.bath.bandwidth must be a positive number [MHz] (%g)\n", mhz); ccexRefuse();
            }
            pulse->bpulse_bandwidth = MHZ_TO_RADKHZ(mhz);
        }
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
                pulse->bsequence_fracs[i] = readSequenceFraction(cJSON_GetArrayItem(row, 0), i, "pulse.bath.sequence");
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

    // Continuous bath drive, the block pulse.bath_drive (port of pyCCEX bdrive; calculation/bath_drive.h): tones of
    // lines [MHz], rabi [MHz], offset [MHz], phase [deg]; shifts, qubit_transverse, line_tolerance [MHz]. What it needs of
    // the rest of the run (method, bath, sequence) is checked where it runs (BathDrive_check)
    if (cJSON* drv = readBlock(block,"bath_drive")) {
        for (cJSON* key = drv->child; key; key = key->next)
            if (strcmp(key->string,"tones") && strcmp(key->string,"shifts") && strcmp(key->string,"qubit_transverse") &&
                strcmp(key->string,"line_tolerance")){
                ccexRefusalPrint(stderr,"Error: pulse.bath_drive has no key \"%s\" (tones, shifts, qubit_transverse, "
                                 "line_tolerance)\n", key->string);
                ccexRefuse();
            }
        BathDrive* dr = new BathDrive();
        cJSON* jt = cJSON_GetObjectItem(drv, "tones");
        if (!cJSON_IsArray(jt)){     // [] is allowed: every spin in the frame of its own line, no field (the undriven reference)
            ccexRefusalPrint(stderr,"Error: pulse.bath_drive.tones must be a list of {\"lines\", \"rabi\", "
                             "\"offset\", \"phase\"}\n");
            ccexRefuse();
        }
        for (int it = 0; it < cJSON_GetArraySize(jt); ++it){
            cJSON* t = cJSON_GetArrayItem(jt, it);
            if (!cJSON_IsObject(t)){
                ccexRefusalPrint(stderr,"Error: pulse.bath_drive.tones[%d] must be an object\n", it); ccexRefuse();
            }
            for (cJSON* key = t->child; key; key = key->next)
                if (strcmp(key->string,"lines") && strcmp(key->string,"rabi") && strcmp(key->string,"offset") &&
                    strcmp(key->string,"phase")){
                    ccexRefusalPrint(stderr,"Error: pulse.bath_drive.tones[%d] has no key \"%s\" (lines, rabi, offset, "
                                     "phase)\n", it, key->string);
                    ccexRefuse();
                }
            BathDriveTone tone;
            char what[96];
            cJSON* jl = cJSON_GetObjectItem(t, "lines");
            if (jl){
                snprintf(what, sizeof what, "pulse.bath_drive.tones[%d].lines", it);
                if (cJSON_IsNumber(jl)) tone.lines.push_back(MHZ_TO_RADKHZ(readFiniteNumber(jl, what)));
                else if (cJSON_IsArray(jl)){
                    for (int k = 0; k < cJSON_GetArraySize(jl); ++k)
                        tone.lines.push_back(MHZ_TO_RADKHZ(readFiniteNumber(cJSON_GetArrayItem(jl, k), what)));
                }else{
                    ccexRefusalPrint(stderr,"Error: %s must be a number or a list of numbers [MHz]\n", what); ccexRefuse();
                }
            }
            if (cJSON* v = cJSON_GetObjectItem(t, "rabi")){
                snprintf(what, sizeof what, "pulse.bath_drive.tones[%d].rabi", it);
                tone.rabi = MHZ_TO_RADKHZ(readFiniteNumber(v, what));
                if (tone.rabi < 0.0){
                    ccexRefusalPrint(stderr,"Error: %s must be >= 0 (a sign goes into phase)\n", what); ccexRefuse();
                }
            }
            if (cJSON* v = cJSON_GetObjectItem(t, "offset")){
                snprintf(what, sizeof what, "pulse.bath_drive.tones[%d].offset", it);
                tone.offset = MHZ_TO_RADKHZ(readFiniteNumber(v, what));
            }
            if (cJSON* v = cJSON_GetObjectItem(t, "phase")){
                snprintf(what, sizeof what, "pulse.bath_drive.tones[%d].phase", it);
                tone.phase = readFiniteNumber(v, what) * M_PI / 180.0;
            }
            dr->tones.push_back(tone);
        }
        dr->shifts = readStrictBool(drv, "pulse.bath_drive", "shifts", false);
        dr->qubit_transverse = readStrictBool(drv, "pulse.bath_drive", "qubit_transverse", true);
        dr->line_tolerance = MHZ_TO_RADKHZ(1e-3);
        if (cJSON* v = cJSON_GetObjectItem(drv, "line_tolerance"))
            dr->line_tolerance = MHZ_TO_RADKHZ(readFiniteNumber(v, "pulse.bath_drive.line_tolerance"));
        delete pulse->bdrive;
        pulse->bdrive = dr;
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
        ccexRefusalPrint(stderr,"Error: output.savemode must be a string (coherence, coherence_tagged, cluster_coherence, cluster_density_matrix) "
                       "or an object {\"mode\": ..., \"nodiv\": true|false}\n");
        ccexRefuse();
    }
    for (cJSON* key = item->child; key; key = key->next)
        if (strcasecmp(key->string,"mode") && strcasecmp(key->string,"nodiv")){
            ccexRefusalPrint(stderr,"Error: output.savemode has no key \"%s\" (mode, nodiv)\n",key->string);
            ccexRefuse();
        }
    cJSON* mode = cJSON_GetObjectItem(item,"mode");
    if (mode && (!cJSON_IsString(mode) || !mode->valuestring)){
        ccexRefusalPrint(stderr,"Error: output.savemode.mode must be a string (coherence, coherence_tagged, cluster_coherence, cluster_density_matrix)\n");
        ccexRefuse();
    }
    Output_setSavemode(op,mode ? mode->valuestring : (char*)"coherence");
    cJSON* nodiv = cJSON_GetObjectItem(item,"nodiv");
    if (nodiv && !cJSON_IsBool(nodiv)){
        ccexRefusalPrint(stderr,"Error: output.savemode.nodiv must be a boolean (true or false)\n");
        ccexRefuse();
    }
    if (readNodiv) op->save_nodiv = nodiv ? (bool)nodiv->valueint : true;
}

void cJSON_readOptionOutput(Output* op, char* fccein){

    if (rank==0 && verbosity){
        printf("\n");
        printf("    >> Read Output Options\n");
        printf("      output : [ path, file, savemode { mode, nodiv } ]\n");
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
        ccexRefusalPrint(stderr,"Error: save_nodiv is now output.savemode.nodiv (\"output\": {\"savemode\": {\"mode\": \"coherence\", \"nodiv\": false}})\n");
        ccexRefuse();
    }
    Output_setSavemode(op,(char*)"coherence");
    op->save_nodiv = true;
    cJSON* output = readBlockOrEmpty(root,"output");   // A65: savemode and outfile were top-level keys
    if (cJSON* item = cJSON_GetObjectItem(output,"savemode")) cJSON_readSavemode(item,op,true);

    // output.path (A65): the folder of the outputs, below -o DIR (main.cpp, resolveOutputName)
    op->outpath[0] = '\0';
    if (cJSON* item = cJSON_GetObjectItem(output,"path")){
        if (!cJSON_IsString(item) || !item->valuestring || !item->valuestring[0] || strlen(item->valuestring) >= sizeof op->outpath){
            ccexRefusalPrint(stderr,"Error: output.path must be a nonempty folder name (relative to the working directory, or "
                           "absolute)\n");
            ccexRefuse();
        }
        strcpy(op->outpath,item->valuestring);
    }

    char* outfileDefault = NULL;
    char* outfile = cJSON_ReadString(output,"file",true,outfileDefault);
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

    if (rank==0 && verbosity){
        printf("\n");
        printf("    >> Read Defect Options\n\n");
        printf("      - Read values of main-key 'Defect' { species }\n");
        printf("        Defect.species[] : [ name, n_onsite_spins, n_axis_states, frame, ( type : single value)\n");
        printf("                    axis, types, spins, gyros, quadrupole_moments,    ( type : array )\n");
        printf("                    nuclear_positions, hyperfine, efg, zfs, detuning ]  ( type : [ axis, spname, array ] )\n\n");
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

    // Defect is {"species": [...]} (a list until A47; its axis and state files are the configuration's since A64)
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
        int naddspin = cJSON_ReadInt(defect,(char*)"n_onsite_spins",true,naddspinDefault);   // n_nuclear_spins until A65
        int navaax = cJSON_ReadInt(defect,(char*)"n_axis_states",true,navaaxDefault);        // n_configurations until A65
        
        DefectArray_setDefect_idf_naddspin(dfa,idf,naddspin);
        DefectArray_setDefect_idf_navaax(dfa,idf,navaax+1); // 0th : main spin
        DefectArray_allocDefect_idf(dfa,idf,navaax+1,naddspin);

        char* dfname = cJSON_ReadString(defect,(char*)"name",false,NULL);
        DefectArray_setDefect_idf_dfname(dfa,idf,dfname);

        DefectArray_setDefect_idf_apprx(dfa,idf,apprxDefault);

        // Defect.species[].frame (coordinate_frame until A90): the frame of nuclear_positions, hyperfine, efg and zfs (with
        // its axis); required with the rotation on (until A90 an absent one was "qubit"). A94 (G4): only for a species
        // that writes one of them, as Qubit.intmap and hyperfine ask only where numbers have a direction (a P1 of
        // detunings alone has none)
        bool directional = false;
        for (const char* key : {"nuclear_positions", "hyperfine", "efg", "zfs"})
            directional = directional || untrackedItem(defect,key) != NULL;
        const FrameValue frame = directional
            ? frameOf(defect,std::string("Defect.species[") + dfname + "]","Defect.species[].frame",
                      "its nuclear_positions, tensors and zfs axis")
            : readFrameKey(defect,std::string("Defect.species[") + dfname + "]");
        if (!directional && frame != FRAME_ABSENT)
            warnNotUsed("Defect.species[].frame","a species without nuclear_positions, hyperfine, efg or zfs has no direction to "
                        "turn");
        DefectArray_setDefect_idf_coordinate_frame(dfa,idf,frame == FRAME_FILE ? DEFECT_COORDINATE_FRAME_FILE
                                                                                : DEFECT_COORDINATE_FRAME_QUBIT);

        // the symmetry axis is zfs.axis of the D form (A94, G5; Defect.species[].axis before), read with the zfs

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
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"hyperfine",&(dfa->defect[idf]->hypf),
                                        navaax,naddspin);
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"efg",&(dfa->defect[idf]->efg),
                                        navaax,naddspin);

        // Object forms provide one physical tensor shared by every indexed
        // configuration. The legacy list remains per-index and is unchanged.
        readDefectZfs(defect,dfa,idf,navaax);

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

    // '#', '!' and '//' start a comment to the end of the line, but only outside a JSON string: a key such as "#note" or
    // a path such as "run#1" or "http://x" used to be cut at the '#' and the input refused or misread (2026-10-08). The
    // text is read whole, so a line is no longer limited to 1023 characters either. '//' (2026-10-16) lets an input be
    // a .jsonc file that vim and VS Code colour as JSON with comments.
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
        } else if (c == '#' || c == '!' || (c == '/' && i + 1 < text.size() && text[i + 1] == '/')) {
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
        // A94: a file the input names has to open (it was a warning, "continuing with defaults", and the readers refused
        // it later or took a default)
        refuseUnreadable(shownName,item->valuestring);
        return item->valuestring;
    }else{
        if(_default){
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity && strcasecmp(name,"bfield")!=0){printf("        %s : use default value ..\n",name);}
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
            if (rank==0 && verbosity){printf("        %s : use default value ..\n",name);}
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
        if ((block == "cluster" && key == "adaptive") || (block.empty() && key == "refine")) continue;
        const std::string where = shown.empty() ? c->string : shown + "." + c->string;
        if (!lookedUpKeys().count(joinKey(block, key))){
            std::string best; size_t bestDistance = 3;
            if (key.size() >= 4)
                for (const std::string& known : blockKeys){
                    const size_t d = keyDistance(key, known);
                    if (d < bestDistance && !objectHasKey(object, known)){ bestDistance = d; best = known; }
                }
            if (!best.empty()) typos.push_back("\"" + where + "\" is not read, did you mean \"" + best
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
