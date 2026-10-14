#!/usr/bin/env bash
#
# Refused combinations, on the PHYSICAL fixtures.
# ---------------------------------------------------------------------------
# test/code_verification/check_rotation/run.sh [10] already covers the refusals on its own toy systems.
# What is checked here is that the same guards fire on the real inputs -- a DFT
# tensor file, two qubits with explicit intmap tensors -- and that they fire
# BEFORE any tensor file is opened, which is the part that would regress quietly
# if the check ever drifted later into the run.
#
# A guard that passes here is one that cannot be reached by accident: the run
# stops, says why, and has not read anything it would have had to interpret in
# the wrong frame.
#
# Usage:  bash test/code_verification/check_refusals/run.sh
# Needs:  bin/main.out; ../check_input_forms/legacy.sh is run first (it writes the legacy cceins).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
# the legacy-form inputs the refused variants are edited from
bash "$HERE/../check_input_forms/legacy.sh" > /dev/null || { echo "legacy inputs could not be written"; exit 1; }

# reject <sys> <name> <expected message fragment> <python edit> [read]
# "read": the refusal is about a tensor file, so it comes after a reader has started
reject(){
    local sys="$1" name="$2" want="$3" edit="$4" when="${5:-}"
    local w="$WORK/$sys"
    [ -f "$w/legacy.json" ] || { echo "  FAIL  $name : no legacy.json"; FAIL=$((FAIL+1)); return; }
    python3 - "$w" "$name" <<PY
import json, sys
w, name = sys.argv[1:3]
d = json.load(open(f"{w}/legacy.json"))
$edit
d["outfile"] = f"./{name}"
json.dump(d, open(f"{w}/{name}.json","w"), indent=4)
PY
    local extra=""
    [ "$sys" = "S3" ] && extra="-s $FIX/single_data/state_DiaP1_1ppm_ -S $FIX/single_data/stateEx_DiaP1_1ppm_ -N 1"
    [ "$sys" = "S5" ] && extra="-s $EX/tutorial/03_Diamond_multiNV_natab/NV2_13C2/state -N 1"
    ( cd "$w" && run_ccex "$BIN_DEV" "$name.json" "proc_$name" $extra )
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return
    fi
    if ! grep -qF "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    # the guard must run before either tensor reader announces itself
    if [ "$when" != read ] && grep -qE "Read the (Hyperfine|Quadrupole) interaction from DFT inputfile" "$w/proc_$name"; then
        echo "  FAIL  $name : refused, but a DFT tensor file was read first"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}

R111='{"enabled": True, "bath_axis": [0.0,0.0,1.0], "qubit_axis": [1.0,1.0,1.0]}'

# A valid input with a command-line option that is not a number (review 78d583 R5): -N and -B were read by
# atoi/atof, so "-N nope" ran an ensemble (nstate 0) and "-N 1.5" one state.
reject_cli(){
    local sys="$1" name="$2" want="$3"; shift 3
    local w="$WORK/$sys"
    ( cd "$w" && run_ccex "$BIN_DEV" "${CLI_INPUT:-legacy.json}" "proc_$name" "$@" )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
echo "=== S1a : command-line numbers"
reject_cli S1a rj_cli_nstate_word "-N must be an integer" -N nope
reject_cli S1a rj_cli_nstate_fraction "-N must be an integer" -N 1.5
reject_cli S1a rj_cli_bfield_word "-B must be a finite number" -B nope
reject_cli S1a rj_cli_bfield_trailing "-B must be a finite number" -B 500G
# -P writes its plan under the output name: with neither outfile nor -o the inventory writer was handed a null name and
# stopped with "basic_string: construction from null" (A55)
python3 - "$WORK/S1a" <<'PY'
import json, sys
d = json.load(open(f"{sys.argv[1]}/legacy.json")); d.pop("outfile", None)
json.dump(d, open(f"{sys.argv[1]}/no_outfile.json", "w"), indent=4)
PY
CLI_INPUT=no_outfile.json reject_cli S1a rj_cli_plan_no_outfile "-P needs an output name" -P
# The command line wins over the input file whatever the order (review b9eb09 R1): -N before -f took nstate from the file.
order_check(){
    local w="$WORK/S1a" a b
    a="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -N 3 -f legacy.json -P -o ./order_a < /dev/null 2>&1 | grep -E '^ +nstate +:')"
    b="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f legacy.json -N 3 -P -o ./order_b < /dev/null 2>&1 | grep -E '^ +nstate +:')"
    if [ -n "$a" ] && [ "$a" = "$b" ] && grep -Eq ':[[:space:]]+3([[:space:]]|$)' <<< "$a"; then echo "  PASS  cli_order_nstate"; PASS=$((PASS+1))
    else echo "  FAIL  cli_order_nstate : '-N 3 -f' gave [$a], '-f -N 3' gave [$b]"; FAIL=$((FAIL+1)); fi
}
order_check
# a Qubit.file (qubitfile until A47) the input names but that is not there ran the qubit at (0,0,0) (review 4dd5fd R1)
reject S1a rj_qubitfile_missing "cannot open the qubit file" "d['Qubit']['file'] = d['Qubit']['file'] + '.missing'"
# Qubit.file and Qubit.qubit together (qubitfile and a Qubit block until A47): the block was ignored without a word, so
# a "two-qubit" tutorial ran one qubit
# (review 20261013-readme R01)
reject S1a rj_qubitfile_and_qubit_block "the Qubit block has both \"file\" and \"qubit\"; keep one form" \
  "d['Qubit']['qubit'] = [{'name': 'q1', 'spin': 1.0, 'gyro': -17608.59705, 'xyz': [0.0, 0.0, 0.0]}]"
# quantity noise was accepted and computed the coherence (review 20261013-doxygen D22), in the file and on the command line
reject S1a rj_quantity_noise "quantity \"noise\" is not supported" "d['quantity'] = 'noise'"
reject_cli S1a rj_cli_quantity_noise "quantity \"noise\" is not supported" -q noise
# --help is -h (review 20261013-doxygen D21): the same text, exit 0
help_alias(){
    local w="$WORK/S1a" a b ra rb
    # from "Usage:" on: the banner above it carries the start time
    a="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -h < /dev/null 2>&1)"; ra=$?; a="$(sed -n '/^Usage:/,$p' <<< "$a")"
    b="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" --help < /dev/null 2>&1)"; rb=$?; b="$(sed -n '/^Usage:/,$p' <<< "$b")"
    if [ $ra -eq 0 ] && [ $rb -eq 0 ] && [ "$a" = "$b" ] && grep -q "^Usage:" <<< "$a"; then echo "  PASS  help_alias (--help = -h, exit 0)"; PASS=$((PASS+1))
    else echo "  FAIL  help_alias : -h exit $ra, --help exit $rb, same text: $([ "$a" = "$b" ] && echo yes || echo no)"; FAIL=$((FAIL+1)); fi
}
help_alias
# nstep 0 left the result arrays empty and the output read result[0] (review 7eed38 R1)
reject S1a rj_nstep_zero "nstep (0 < 1)" "d['nstep'] = 0"
# A Manual sequence (review 5f7834 R2): a sum of fractions short of 1, a fraction that is not a/b or a number all
# through, a row that is not [fraction, axis, angle] used to pass (or skip the row); a number fraction used to throw
# and an angle of 90.5 became 90 - sequence_accepted checks those two now go through as given.
MANUAL="d['method'] = 'gCCE'; d['pulse'] = {'count': 2, 'name': 'Manual'}"   # gCCE: cCCE refuses an I pulse on its own
reject S1a rj_seq_sum_short "sum to 0.75" "$MANUAL; d['pulse']['sequence'] = [['1/4','I',0],['1/4','X',180],['1/4','I',0]]"
reject S1a rj_seq_fraction_text "must be a/b" "$MANUAL; d['pulse']['sequence'] = [['1/2x','I',0],['1/4','X',180],['1/4','I',0]]"
reject S1a rj_seq_row_short "must be [fraction, axis, angle]" "$MANUAL; d['pulse']['sequence'] = [['1/2','I'],['1/4','X',180],['1/4','I',0]]"
sequence_accepted(){
    local w="$WORK/S1a"
    python3 - "$w" <<'PY'
import json, sys
w = sys.argv[1]; d = json.load(open(f"{w}/legacy.json"))
d.update(method='gCCE', pulse={'count': 2, 'name': 'Manual', 'sequence': [[0.25,'I',0],['1/2','X',90.5],[0.25,'I',0]]}, outfile='./seq_ok')
json.dump(d, open(f"{w}/seq_ok.json", "w"))
PY
    if (cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f seq_ok.json -P -o ./seq_ok < /dev/null > proc_seq_ok 2>&1) &&
       grep -Eq '0\.250 +0\.750 +0\.500 +X +90\.50' "$w/proc_seq_ok"; then
        echo "  PASS  sequence_accepted (a number fraction, an angle of 90.5)"; PASS=$((PASS+1))
    else echo "  FAIL  sequence_accepted : see $w/proc_seq_ok"; FAIL=$((FAIL+1)); fi
}
sequence_accepted
# a finite double outside the float it is stored in became inf (review fa7f39 R3)
reject S1a rj_bfield_float_overflow "within the float range" "d['bfield'] = 1e40"

echo "=== S4 : DFT tensor file + rotation"
reject S4 rj_hf_frame_missing "hamiltonian.common.hyperfine dipole \"dft\" with the coordinate frame rotation on" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'"
reject S4 rj_qd_frame_missing "hamiltonian.common.quadrupole (structure plain) with the coordinate frame rotation on" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'"
reject S4 rj_tilted_bfield "bfield" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'
d['bfield'] = [100.0, 0.0, 30000.0]"
reject S4 rj_pos_frame_qubit "qubit_position_frame" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'qubit'"

echo
echo "=== S5 : two qubits + rotation"
reject S5 rj_refq_missing "reference_qubit is required" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'
for e in d['Qubit']['intmap']: e['tensor_frame'] = 'bath'"
reject S5 rj_intmap_frame_missing "intmap tensor" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['reference_qubit'] = 'q1'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'"
# D43: a state index (Qubit.alpha_state/beta_state {"index"}, alphaidx/betaidx until A47) counts the qubit eigenstates by
# descending energy and overrides alphams/betams; written together and naming other states, the run stops (S5 writes
# alphams 1 / betams 0 for both NVs; index 2 is not |0,0>, the lowest level)
reject S5 rj_alphaidx_vs_alphams "Keep one of the two" \
  "d['Qubit']['alpha_state'] = {'index': 1}; d['Qubit']['beta_state'] = {'index': 2}" read
# only what is written constrains (review b467fe R1): q1's betams 0 alone against beta index 0, the highest level |+1,+1>
reject S5 rj_alphaidx_one_side "Keep one of the two" \
  "d['Qubit']['alpha_state'] = {'index': 1}; d['Qubit']['beta_state'] = {'index': 0}
for q in d['Qubit']['qubit']: q.pop('alphams', None); q.pop('betams', None)
d['Qubit']['qubit'][0]['betams'] = 0" read
# and the same kind of input that agrees runs: q1's alphams 1 alone with alpha index 0 (|+1,+1>); beta index 1 is not written
# against, so it is free. A short run (nstep 2), since the check sits where the run picks psia/psib.
idx_ms_accepted(){
    local w="$WORK/S5"
    python3 - "$w" <<'PY'
import json, sys
w = sys.argv[1]; d = json.load(open(f"{w}/legacy.json"))
for q in d['Qubit']['qubit']: q.pop('alphams', None); q.pop('betams', None)
d['Qubit']['qubit'][0]['alphams'] = 1
d['Qubit'].update(alpha_state={'index': 0}, beta_state={'index': 1}); d.update(nstep=2, outfile='./idx_ms_ok')
json.dump(d, open(f"{w}/idx_ms_ok.json", "w"))
PY
    if (cd "$w" && run_ccex "$BIN_DEV" idx_ms_ok.json proc_idx_ms_ok -s "$EX/tutorial/03_Diamond_multiNV_natab/NV2_13C2/state" -N 1) &&
       grep -Eq 'psia is mostly +: +m_s \|\+1, \+1>' "$w/proc_idx_ms_ok"; then
        echo "  PASS  idx_ms_accepted (an alpha index with a partly written alphams that agrees)"; PASS=$((PASS+1))
    else echo "  FAIL  idx_ms_accepted : see $w/proc_idx_ms_ok"; FAIL=$((FAIL+1)); fi
}
idx_ms_accepted
# D28: Qubit.intmap[].tensor shorter than 3x3 was padded with zeros
reject S5 rj_intmap_tensor_short "3x3" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0]]"
reject S5 rj_intmap_tensor_string "finite JSON number" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0,0,0],[0,'a',0],[0,0,1.0]]"
# review 20261008-autoconv-d42 #5: a qubit array shorter than nqubit left the last qubit allocated but never set; the
# number of qubits is the length of Qubit.qubit since A47 (Qubit.nqubit removed), so what is left to refuse is no qubit
reject S5 rj_qubit_array_short "Qubit.qubit must be an array of one entry per qubit" \
  "d['Qubit']['qubit'] = []"
reject S5 rj_refq_unknown "reference_qubit" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['reference_qubit'] = 'q9'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'
for e in d['Qubit']['intmap']: e['tensor_frame'] = 'bath'"

echo
echo "=== S3 : removed key, and the Defect object forms"
reject S3 rj_magnetic_field_axis "magnetic_field_axis" \
  "d['magnetic_field_axis'] = [1.0,1.0,1.0]"
reject S3 rj_defaxis_without_rot "defect_axis_reference" \
  "d['defect_axis_reference'] = 'qubit_axis'"
reject S3 rj_zfs_de_no_axis "axis is missing" \
  "df = d['Defect']['species'][0]
df['coordinate_frame'] = 'qubit'
df['zfs'] = {'D': 2870.0, 'unit': 'MHz'}"
reject S3 rj_zfs_E_nonzero "zfs.E is non-zero" \
  "df = d['Defect']['species'][0]
df['coordinate_frame'] = 'qubit'
df['axis'] = [1.0,1.0,1.0]
df['zfs'] = {'D': 2870.0, 'E': 5.0, 'unit': 'MHz'}"
# axis is validated before zfs, so these are two different guards, not one.
reject S3 rj_axis_no_frame "axis is present but coordinate_frame is missing" \
  "df = d['Defect']['species'][0]
df['axis'] = [1.0,1.0,1.0]"
reject S3 rj_zfs_object_no_frame "object-form zfs but has no coordinate_frame" \
  "df = d['Defect']['species'][0]
df['zfs'] = {'D': 2870.0, 'unit': 'MHz'}"
# The one behaviour change that can stop a v1.1.0 input which used to complete. The
# arrays hold n_configurations+1 entries (index 0 is the main spin), so n_configurations+1 is one past the
# end; it used to be accepted and written to. S3's own fixture uses 1..n_configurations, and
# axis 1 shows that side is untouched.
reject S3 rj_iax_out_of_range "out of range" \
  "df = d['Defect']['species'][0]
df['hyperfine'] = df['hyperfine'] + [[df['n_configurations'] + 1, '14N', [0.0]*9]]"

# D28 (2026-10-02): the entries of the indexed form are checked as the object form checks them. A string
# or fractional index, a non-string spin name, a non-number value and a field that is not an array were
# read as 0, truncated, crashed or skipped.
reject S3 rj_iax_string "out of range" \
  "df = d['Defect']['species'][0]
df['hyperfine'] = df['hyperfine'] + [['1', '14N', [0.0]*9]]"
reject S3 rj_iax_fraction "out of range" \
  "df = d['Defect']['species'][0]
df['hyperfine'] = df['hyperfine'] + [[1.5, '14N', [0.0]*9]]"
reject S3 rj_spin_not_string "spin name must be a string" \
  "df = d['Defect']['species'][0]
df['hyperfine'] = df['hyperfine'] + [[1, 7, [0.0]*9]]"
reject S3 rj_value_string "must be a finite JSON number" \
  "df = d['Defect']['species'][0]
df['hyperfine'] = df['hyperfine'] + [[1, '14N', ['x'] + [0.0]*8]]"
reject S3 rj_rxyzs_not_array "must be an array of [index, name, values] entries" \
  "df = d['Defect']['species'][0]
df['nuclear_positions'] = 5"

# D42 (2026-10-08): an array shorter than the input needs was filled with empty strings or 0 and used. The bath_shift cases
# ran to the end before (zero offsets); the types case was refused later anyway (the hypf entries name 14N, which the empty
# type then lacked), so it checks the earlier, clearer refusal.
reject S3 rj_types_short "types has 0 entries" \
  "d['Defect']['species'][0]['types'] = []"
reject S3 rj_bathadjust_row_short "bath_shift[0] has 2 numbers" \
  "d['bath_shift'] = [[0.0, 0.0]]"
reject S3 rj_bathadjust_rows_short "bath_shift has 0 rows" \
  "d['bath_shift'] = []"
# n_configurations 4 means axes 1-4: an axis file with 5 was read one past the tables (2026-10-09); refused while the file is read
reject S3 rj_avaax_past_navaax "The principal axis is out of range" \
  "lines = open(d['Defect']['axis_file']).read().split('\n')
open(f'{w}/avaax_axis5','w').write('\n'.join([lines[0]] + ['5' if l.strip() else l for l in lines[1:]]))
d['Defect']['axis_file'] = f'{w}/avaax_axis5'" read

reject S3 rj_bad_unit "unit" \
  "df = d['Defect']['species'][0]
df['coordinate_frame'] = 'qubit'
df['quadrupole_moments'] = {'values': [2.044], 'unit': 'furlongs'}"
reject S3 rj_rot_alias_both "both" \
  "d['coordinate_frame_rotation'] = $R111
d['bath_coordinate_rotation'] = $R111"

echo
echo "=== S4 : a tensor file without a condition the reader needs"
# READ_Tensor_const / READ_Tensor_etc said so and then ended with exit 0 until 2026-10-01 (D12)
TENSOR_EDIT='
import re
key, cond = ARGS
block = d["hamiltonian"]["common"][key]
src = block["file"]
dst = f"{w}/{name}_" + src.rsplit("/", 1)[1]
open(dst, "w").write("".join(l for l in open(src) if cond not in l))
block["file"] = dst'
reject S4 rj_hf_no_gfactor "Tensor Condition (g-factor___) can not read" \
  "ARGS = ('hyperfine', 'g-factor___')
$TENSOR_EDIT" read
reject S4 rj_qd_no_etc "Tensor Condition (etc____) can not read" \
  "ARGS = ('quadrupole', 'etc____')
$TENSOR_EDIT" read

echo
echo "=== a present key of the wrong type"
# The defaulting readers took a present key of another type for a missing one and ran
# with the default until 2026-10-02 (D18): "knight": "true" was false, a quoted gyro the
# default gyro, nstep 10.5 was 10. Absent keys keep their defaults; a qubit's detuning may be [x]. (The switches that
# took 0 or 1 - "addsubclus": 0 in 29 order-2 inputs - are true or false only since A47.)
reject S1a rj_type_knight "hamiltonian.cCCE.knight must be a boolean (true or false)" "d['hamiltonian'] = {'cCCE': {'knight': 'true'}}"
reject S1a rj_type_hfmedi "hamiltonian.cCCE.hfmedi must be a boolean (true or false)" "d['hamiltonian'] = {'cCCE': {'hfmedi': 'true'}}"
reject S1a rj_type_bool_number "add_subclusters must be true or false" "d['add_subclusters'] = 2"
reject S1a rj_type_bool_zero "add_subclusters must be true or false" "d['add_subclusters'] = 0"
reject S5 rj_type_qubit_gyro "gyro must be a number" "d['Qubit']['qubit'][0]['gyro'] = str(d['Qubit']['qubit'][0]['gyro'])"
reject S1a rj_type_nstep_fraction "nstep must be an integer" "d['nstep'] = d['nstep'] + 0.5"
reject S1a rj_type_bfield_string "bfield must be a number" "d['bfield'] = str(d['bfield'])"
reject S1a rj_type_nk_element "max_clusters_per_order[1] must be an integer" "d['max_clusters_per_order'] = [0, '1', 0]"
reject S1a rj_type_bathadjust_element "bath_shift[0][2] must be a number" "d['bath_shift'] = [[0.0, 0.0, '1.0']]"
reject S1a rj_type_method "method must be a string" "d['method'] = 1"

echo
echo "=== a key nobody reads, close to one that took its default"
# Until 2026-10-02 (D24) such a key was ignored without a word: "nstates" left nstate at 0, and
# the top-level "alphams" of the qubitfile form (which read qalphams) left the NV branch at +1; the same inside a block
# (Qubit.alphms for Qubit.alphams, A47).
reject S1a rj_typo_nstates "did you mean \"nstate\"" "d['nstates'] = d.pop('nstate')"
reject S1a rj_typo_qubitfile_alphams "did you mean \"alphams\"" "d['Qubit']['alphms'] = d['Qubit'].pop('alphams')"

echo
echo "=== a result that cannot be written in full"
# Output_save and Output_save_info checked only that the files opened (and
# Output_save_info only the first of the two): with a full disk, a quota or a
# file-size limit the run left a truncated result and ended with exit 0 until
# 2026-10-01 (D10). On a small case whose result is 18 KB (47 KB in info mode):
#   full     the result files are links to /dev/full, so every write fails with
#            ENOSPC, as on a full disk. Until 2026-10-02 a 4 KB file-size limit
#            (ulimit -f) did this, which Open MPI's mpirun did not survive: its own
#            shared-memory files met the limit and MPI_Init failed (idm, job 16853)
#   nodiv    a directory where the _noDiv file has to go, so only the second file
#            cannot be opened
CASE="$REPO/test/code_verification/check_numerics/cases/natab_cce2_hahn"
write_fail(){
    local name="$1" savemode="$2" how="$3" want="$4"
    local outw="$WORK/write_failure/$name"
    rm -rf "${outw:?}"; mkdir -p "$outw"
    python3 - "$CASE/ccein.json" "$outw/ccein.json" "$savemode" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); d["savemode"] = sys.argv[3]
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
    [ "$how" = nodiv ] && mkdir "$outw/out_noDiv"
    [ "$how" = full ] && ln -s /dev/full "$outw/out_noDiv" && ln -s /dev/full "$outw/out_wiDiv"
    # relative paths in the case's input are read from the case's folder
    ( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$outw/ccein.json" -o "$outw/out" 2>&1 < /dev/null ) > "$outw/proc"
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "  FAIL  $name : ended with exit 0"; FAIL=$((FAIL+1))
    elif ! grep -qF "$want $outw/out_" "$outw/proc"; then
        echo "  FAIL  $name : exit=$rc but no message naming the file"; FAIL=$((FAIL+1))
    else
        echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
    fi
}
write_fail write_full_normal  coherence        full  "cannot write file"
write_fail write_full_info    coherence_tagged full  "cannot write file"
write_fail open_nodiv_normal  coherence        nodiv "cannot open file"
write_fail open_nodiv_info    coherence_tagged nodiv "cannot open file"

echo
echo "=== the factors and except work modes of method D's regroup"
# adaptive-regroup-stage3.md 4절, -stage3b.md 1절: a request CCEX cannot answer for these modes is refused at the
# request, before any state is calculated. The driver checks the same before it asks (an in-process refusal would end
# the job), so these requests are written by hand: a P1 input of the regroup regression (six sites, five states).
P1IN="$REPO/test/cluster_verification/check_adaptive/regression/inputs/regroup-p1-0.json"
reject_work(){
    local name="$1" want="$2" edit="$3" request="$4"; shift 4
    local w="$WORK/adaptive_work"; mkdir -p "$w"
    python3 - "$P1IN" "$w" "$name" "$request" <<PY
import json, sys
src, w, name, request = sys.argv[1:5]
d = json.load(open(src))
$edit
r = json.loads(request); r["response"] = f"{w}/{name}.response"
json.dump(r, open(f"{w}/{name}.request.json", "w"))
d["adaptive"] = {"work_file": f"{w}/{name}.request.json"}
d["outfile"] = f"{w}/{name}"
json.dump(d, open(f"{w}/{name}.json", "w"), indent=1)
PY
    # relative paths in the input are the repository's
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" "$@" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    if ls "$w/$name.response".*state-* > /dev/null 2>&1; then
        echo "  FAIL  $name : refused, but after a state was calculated"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
B='"bundles": [[0, 3], [1, 5], [2, 4]]'
F="{\"mode\": \"factors\", $B, \"state\": 1, \"step\": 3}"
E="{\"mode\": \"except\", $B, \"except\": [[0, 1, 3, 5]]}"
reject_work rj_factors_state_zero "factors state is not a state of this run" "" "${F/\"state\": 1/\"state\": 0}"
reject_work rj_factors_state_past "factors state is not a state of this run" "" "${F/\"state\": 1/\"state\": 6}"
reject_work rj_factors_step "factors step is outside 0..nstep-1" "" "${F/\"step\": 3/\"step\": 30}"
reject_work rj_factors_step_text "factors work requires integer state and step" "" "${F/\"step\": 3/\"step\": \"3\"}"
reject_work rj_factors_gcce "factors/except work requires method cce" "d['method'] = 'gCCE'" "$F"
reject_work rj_factors_order1 "factors/except work requires bundle order 2" "d['order'] = 1" "$F"
reject_work rj_factors_probes "probes require probe mode" "" "${F%\}}, \"probes\": [[0]]}"
reject_work rj_factors_contributions "cannot be combined with cluster_contributions" "d['cluster_contributions'] = {'enabled': True}" "$F"
reject_work rj_factors_plan_only "probe, factors and except work require dynamics" "" "$F" -P
reject_work rj_state_without_factors "state and step require factors mode" "" "{\"mode\": \"run\", $B, \"state\": 1, \"step\": 3}"
reject_work rj_except_gcce "factors/except work requires method cce" "d['method'] = 'gCCE'" "$E"
reject_work rj_except_order1 "factors/except work requires bundle order 2" "d['order'] = 1" "$E"
reject_work rj_except_probes "probes require probe mode" "" "${E%\}}, \"probes\": [[0]]}"
reject_work rj_except_empty "except mode requires a nonempty except list" "" "{\"mode\": \"except\", $B, \"except\": []}"
reject_work rj_except_twice "except lists a term twice" "" "{\"mode\": \"except\", $B, \"except\": [[0, 3], [0, 3]]}"
reject_work rj_except_not_a_term "an except entry is not a term of this partition's expansion" "" "{\"mode\": \"except\", $B, \"except\": [[0, 1]]}"
reject_work rj_except_plan_only "probe, factors and except work require dynamics" "" "$E" -P
reject_work rj_except_without_mode "except requires except mode" "" "{\"mode\": \"run\", $B, \"except\": [[0, 3]]}"

echo
echo "=== evolution_update"
# exact | recursive (2026-10-11): recursive (A16 b, A18) exists on the CPU for cCCE (single state and ensemble) and for
# gCCE with propagator=eigen (evolution vector and matrix); on the GPU the run would compute as exact and gCCE expm has no
# phases to carry over, so a given recursive there is refused (without the key those run exact, A29). On check_numerics
# cases, refused at the start. (The ensemble cCCE and gCCE matrix were refused until A18; they run now,
# ../check_input_forms/evolution_update.sh.) Since A36 (2026-10-12) backend auto in this CPU-only binary computes every
# cluster on the CPU and a given recursive runs (check_input_forms E5); backend cuda stays refused here, and recursive
# with backend auto or cuda in bin/main_gpu.out is refused in test/gpu_verification/check_gpu/check_options.py.
reject_case(){    # reject_case <case> <name> <expected message fragment> <python edit> [command-line options]
    local case="$1" name="$2" want="$3" edit="$4"; shift 4
    local w="$WORK/evolution_update_refusals"; mkdir -p "$w"
    python3 - "$REPO/test/code_verification/check_numerics/cases/$case/ccein.json" "$w/$name.json" <<PY
import json, sys
d = json.load(open(sys.argv[1]))
$edit
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
    # relative paths in the case's input are read from the case's folder
    ( cd "$REPO/test/code_verification/check_numerics/cases/$case" &&
      "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" -o "$w/$name" "$@" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
reject_case p1_cce2_single rj_update_value "evolution_update must be exact or recursive" "d['evolution_update'] = 'fast'"
reject_case p1_cce2_single rj_update_type "evolution_update must be exact or recursive" "d['evolution_update'] = 1"
reject_case p1_manual_hahn rj_update_gcce_expm "this run has method=gCCE with propagator=expm" \
  "d['evolution_update'] = 'recursive'; d['evolution'] = 'matrix'; d['propagator'] = 'expm'"
reject_case p1_cce2_single rj_update_backend "this run has backend=cuda" \
  "d['evolution_update'] = 'recursive'; d['backend'] = 'cuda'"

echo "=== gpu_options"
# The GPU settings are one block since 2026-10-12 (A26, docs/user/InputFile_Backend.md#tag_gpu_options): the flat keys
# before it are refused with the new form named, and so are an unknown key of the block, a block that is not an object,
# and a packed_max_dim other than 0, 2 or 4. packed_max_dim needs method=CCE and the GPU binary (this check runs
# bin/main.out), both refused before the bath is read. The accepted forms are in ../check_input_forms/gpu_options.sh.
reject_case p1_cce2_single rj_gpu_min_dim "gpu_min_dim is now gpu_options.min_dim" "d['gpu_min_dim'] = 256"
reject_case p1_cce2_single rj_gpu_batch_steps "gpu_batch_steps is now gpu_options.evolution_batch_steps" "d['gpu_batch_steps'] = 64"
reject_case p1_cce2_single rj_gpu_batch_clusters "gpu_batch_clusters is now gpu_options.evolution_batch_size" "d['gpu_batch_clusters'] = 8"
reject_case p1_cce2_single rj_gpu_cluster_batch_mb "gpu_cluster_batch_mb is now gpu_options.cluster_batch_mb" "d['gpu_cluster_batch_mb'] = 64"
reject_case p1_cce2_single rj_gpu_cache_mb "gpu_cache_mb is now gpu_options.cache_mb" "d['gpu_cache_mb'] = 64"
reject_case p1_cce2_single rj_gpu_options_type "gpu_options must be an object" "d['gpu_options'] = 256"
reject_case p1_cce2_single rj_gpu_options_key "unknown key gpu_options.packed_dim" "d['gpu_options'] = {'packed_dim': 4}"
reject_case p1_cce2_single rj_packed_3 "gpu_options.packed_max_dim must be 0 (off), 2 or 4" "d['gpu_options'] = {'packed_max_dim': 3}"
reject_case p1_cce2_single rj_packed_8 "gpu_options.packed_max_dim must be 0 (off), 2 or 4" "d['gpu_options'] = {'packed_max_dim': 8}"
reject_case p1_cce2_single rj_packed_string "finite JSON number" "d['gpu_options'] = {'packed_max_dim': '4'}"
reject_case p1_cce2_single rj_min_dim_zero "gpu_options.min_dim must be a positive 32-bit integer" "d['gpu_options'] = {'min_dim': 0}"
reject_case p1_cce2_single rj_packed_cpu_binary "gpu_options.packed_max_dim needs the GPU binary" "d['gpu_options'] = {'packed_max_dim': 4}"
reject_case p1_manual_hahn rj_packed_gcce "gpu_options.packed_max_dim works with method=CCE" \
  "d['gpu_options'] = {'packed_max_dim': 4}; d['method'] = 'gCCE'"

echo "=== savemode"
# savemode is a string or {"mode": ..., "nodiv": ...} since 2026-10-13 (docs/user/InputFile_Output.md#tag_savemode): the
# names before it and the top-level save_nodiv are refused with the new form named, and so are an unknown key of the block,
# a block value of the wrong type and an unknown mode. The accepted block is in ../check_input_forms/nodiv.sh.
reject_case p1_cce2_single rj_savemode_normal 'savemode "normal" is now "coherence"' "d['savemode'] = 'normal'"
reject_case p1_cce2_single rj_savemode_info 'savemode "info" is now "coherence_tagged"' "d['savemode'] = 'info'"
reject_case p1_cce2_single rj_savemode_all 'savemode "all" is now "cluster_coherence"' "d['savemode'] = 'all'"
reject_case p1_cce2_single rj_savemode_allfull 'savemode "allfull" is now "cluster_density_matrix"' "d['savemode'] = 'allfull'"
reject_case p1_cce2_single rj_savemode_old_in_block 'savemode "all" is now "cluster_coherence"' "d['savemode'] = {'mode': 'all'}"
reject_case p1_cce2_single rj_save_nodiv "save_nodiv is now savemode.nodiv" "d['save_nodiv'] = False"
reject_case p1_cce2_single rj_savemode_key 'savemode has no key "nodvi" (mode, nodiv)' "d['savemode'] = {'mode': 'coherence', 'nodvi': False}"
reject_case p1_cce2_single rj_savemode_type "savemode must be a string" "d['savemode'] = 1"
reject_case p1_cce2_single rj_savemode_mode_type "savemode.mode must be a string" "d['savemode'] = {'mode': 1}"
reject_case p1_cce2_single rj_savemode_nodiv_type "savemode.nodiv must be a boolean" "d['savemode'] = {'mode': 'coherence', 'nodiv': 0}"
reject_case p1_cce2_single rj_savemode_unknown 'savemode "averaged" is not one of' "d['savemode'] = 'averaged'"

echo
echo "=== propagator"
# propagator (A2, 2026-10-12): expm is gCCE only; sparse runs for one sampled state of the cCCE and gCCE evolution=vector
# with one qubit, backend=cpu, evolution_update=exact, without knight and hfmedi - everything else is refused at the
# start (the paths that run, and the malformed values, are in ../check_input_forms/propagator.sh).
SP="d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_cce2_single rj_prop_expm_cce "propagator=expm is implemented for method=gCCE only" "d['propagator'] = 'expm'"
reject_case p1_cce2_ensemble rj_prop_sparse_ensemble "this run has nstate=0 (the ensemble)" "$SP"
reject_case p1_manual_hahn rj_prop_sparse_matrix "this run has method=gCCE with evolution=matrix" "$SP; d['evolution'] = 'matrix'"
reject_case p1_cce2_single rj_prop_sparse_cuda "this run has backend=cuda" "$SP; d['backend'] = 'cuda'"
reject_case p1_manual_hahn rj_prop_sparse_auto "this run has backend=auto" "$SP; d['backend'] = 'auto'"
reject_case p1_cce2_single rj_prop_sparse_recursive "this run has evolution_update=recursive" "$SP; d['evolution_update'] = 'recursive'"
reject_case p1_manual_hahn rj_prop_sparse_recursive_gcce "this run has evolution_update=recursive" "$SP; d['evolution_update'] = 'recursive'"
reject_case natab_knight_shift rj_prop_sparse_knight "this run has the knight term" "$SP; d['nstate'] = 1"
reject_case hbn_hfmedi rj_prop_sparse_hfmedi "this run has the hfmedi term" "$SP; d['nstate'] = 1"
# a second NV 20 A away added to the one-qubit Qubit block (the check_numerics cases hold one qubit each; the old
# natab_gcce2_two_qubits held one too, through its qubitfile, and is natab_gcce2_bathadjust since 2026-10-13)
TWO="q = d['Qubit']; q2 = dict(q['qubit'][0], name='qubit2', xyz=[x + 20.0 for x in q['qubit'][0]['xyz']]); q['qubit'].append(q2); q['intmap'].append(dict(q['intmap'][0], between=['qubit2', 'qubit2']))"
reject_case p1_gcce2_zfs rj_prop_sparse_two_qubits "this run has nqubit=2" "$SP; $TWO"

echo
echo "=== hyperfine tensor file format v3"
# (2026-10-13) A file with the provenance line "# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz),
# source = ..., S = ..." is read x 1, so it must not contradict itself or the input: line 1 has to be the tag of its
# source, the line has to be line 2 and well formed, and S has to be the qubit spin. Each file is the h-BN v3 file of
# check_numerics with one edit; what such files are read as is ../check_input_forms/tensor_format.sh.
TF="$WORK/tensor_format_refusals"; mkdir -p "$TF"
V3="$REPO/test/code_verification/check_numerics/inputs/hbn/hBN_VB_flat_ATensor_HSE_CellPBE_v3"
mk_tensor(){ sed -e "$2" "$V3" > "$TF/$1"; }    # mk_tensor <name> <sed script>
mk_tensor qe_tag       '1s/^VASP$/QUANTUMESPRESSO/'
mk_tensor qe_source    '2s/source = VASP OUTCAR/source = QE hyperfine.out/'
mk_tensor no_tag       '1d'
mk_tensor spin_half    '2s/S = 1$/S = 0.5/'
mk_tensor spin_fraction '2s/S = 1$/S = 0.7/'
mk_tensor no_spin      '2s/, S = 1$//'
mk_tensor version4     '2s/CCEX-tensor v3/CCEX-tensor v4/'
mk_tensor bad_source   '2s/VASP OUTCAR/VASP vasprun.xml/'
mk_tensor twice        '2p'
HB=hbn_cce2_hf3_qd0
reject_case $HB rj_tensor_tag_vs_source_qe "line 1 is 'QUANTUMESPRESSO' but the provenance line says source = VASP OUTCAR" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/qe_tag'"
reject_case $HB rj_tensor_source_vs_tag_vasp "line 1 is 'VASP' but the provenance line says source = QE hyperfine.out" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/qe_source'"
reject_case $HB rj_tensor_no_tag "it must be line 2" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/no_tag'"
reject_case $HB rj_tensor_spin_vs_defect "spin S = 0.5 defect (provenance line) but the qubit spin is 1" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/spin_half'"
reject_case $HB rj_tensor_defect_vs_spin "spin S = 1 defect (provenance line) but the qubit spin is 1.5" "d['hamiltonian']['common']['hyperfine']['file'] = '$V3'; d['Qubit'].update(spin=1.5, alphams=1.5, betams=0.5)"
reject_case $HB rj_tensor_spin_fraction "S must be a positive half-integer" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/spin_fraction'"
reject_case $HB rj_tensor_no_spin "has no \", S = <S>\"" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/no_spin'"
reject_case $HB rj_tensor_version4 "is not a tensor format this version reads" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/version4'"
reject_case $HB rj_tensor_unknown_source "names source 'VASP vasprun.xml'" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/bad_source'"
reject_case $HB rj_tensor_two_lines "more than one '# CCEX-tensor' line" "d['hamiltonian']['common']['hyperfine']['file'] = '$TF/twice'"

echo
echo "=== input keys renamed (A47)"
# The input keys were renamed and grouped on 2026-10-13 ((internal notes), docs/user/InputFile.md). An
# old name is refused with its new one, as savemode's old values are, and the message names the converter; each case puts
# one old key (or old value, or old form) back into a converted check_numerics input.
OLDK='input key "'
reject_case natab_cce2_hahn rj_old_npulse "${OLDK}npulse\" is now \"pulse.count\"" "d['npulse'] = 1"
reject_case natab_cce2_hahn rj_old_hint "example/tools/InputKeys/convert_input_keys.py" "d['npulse'] = 1"
reject_case natab_cce2_hahn rj_old_qubitfile "${OLDK}qubitfile\" is now \"Qubit.file\"" "d['qubitfile'] = d['Qubit'].pop('file')"
reject_case natab_cce2_hahn rj_old_rbathcut "${OLDK}rbathcut\" is now the inner radius of \"rbath\": [inner, outer]" "d['rbathcut'] = 1.0"
reject_case natab_cce2_hahn rj_old_addsubclus "${OLDK}addsubclus\" is now \"add_subclusters\"" "d['addsubclus'] = 0"
reject_case natab_cce2_hahn rj_old_nk "${OLDK}nk\" is now \"max_clusters_per_order\"" "d['nk'] = [0, 0, 0]"
reject_case natab_cce2_hahn rj_old_statefile "${OLDK}statefile\" is now \"bath_state_file\"" "d['statefile'] = 'state'"
reject_case natab_cce2_hahn rj_old_vector_mean_field "${OLDK}vector_mean_field\" is now \"hamiltonian.common.bath_bath_mean_field\"" "d['vector_mean_field'] = True"
reject_case natab_cce2_hahn rj_old_knight "${OLDK}knight\" is now \"hamiltonian.cCCE.knight\"" "d['knight'] = True"
reject_case hbn_cce2_hf3_qd0 rj_old_hf_readmode "${OLDK}hf_readmode\" is now \"hamiltonian.common.hyperfine\"" "d['hf_readmode'] = 3"
reject_case hbn_cce2_hf3_qd0 rj_old_DefectTotSpin "${OLDK}DefectTotSpin\" was removed: the qubit spin is used" "d['DefectTotSpin'] = 1.0"
reject_case hbn_cce2_hf3_qd0 rj_old_CorrTotSpin "${OLDK}CorrTotSpin\" was removed" "d['CorrTotSpin'] = 0.0"
reject_case p1_pcce2_sk2 rj_old_sK "${OLDK}sK\" is now \"pcce.partition_size\"" "d['sK'] = d['pcce'].pop('partition_size')"
reject_case p1_pcce2_sk2 rj_old_method_pcce "method \"pCCE\" is now \"method\": \"cCCE\" with \"cluster_mode\": \"pcce\"" "d['method'] = 'pCCE'; d.pop('cluster_mode')"
reject_case natab_cce2_hahn rj_old_method_cce "method \"CCE\" is now \"cCCE\"" "d['method'] = 'CCE'"
reject_case natab_cce2_hahn rj_old_cluster_mode "cluster_mode \"adaptive\" is now \"bundle\"" "d['cluster_mode'] = 'adaptive'"
reject_case natab_cce2_hahn rj_old_bundle_max_dim "${OLDK}bundle_max_dim\" is now \"bundle.max_dim\"" "d['bundle_max_dim'] = 8"
reject_case natab_cce2_hahn rj_old_gcce_refinement "\"gcce_refinement\": true|false is now \"gcce_refinement\": {\"enabled\": true|false}" "d['gcce_refinement'] = False"
reject_case hbn_gcce2_oh rj_old_nqubit "${OLDK}Qubit.nqubit\" was removed" "d['Qubit']['nqubit'] = 1"
reject_case hbn_gcce2_oh rj_old_psia "${OLDK}Qubit.psia\" is now \"Qubit.alpha_state\": {\"vector\": [...]}" "d['Qubit']['psia'] = d['Qubit'].pop('alpha_state')['vector']"
reject_case hbn_gcce2_oh rj_old_overhaus "${OLDK}Qubit.overhaus\" is now \"hamiltonian.common.qubit_bath_mean_field\"" "d['Qubit']['overhaus'] = True"
reject_case p1_cce2_single rj_old_defect_list "\"Defect\": [...] is now \"Defect\": {\"species\": [...]}" "d['Defect'] = d['Defect']['species']"
reject_case p1_cce2_single rj_old_dfname "${OLDK}Defect.species[0].dfname\" is now \"name\"" "d['Defect']['species'][0]['dfname'] = d['Defect']['species'][0].pop('name')"
reject_case p1_cce2_single rj_old_exstatefile "${OLDK}exstatefile\" is now \"Defect.state_file\"" "d['exstatefile'] = 'x'"
reject_case natab_cce2_hahn rj_old_cli_method "-m \"CCE\" is now \"cCCE\"" "" -m CCE
reject_case natab_cce2_hahn rj_old_cli_method_pcce "-m \"pCCE\" is now method \"cCCE\" with \"cluster_mode\": \"pcce\" in the input" "" -m pCCE
# cluster_mode pcce runs with gCCE, meCCE and megCCE too (internal pgCCE, pmeCCE, pmegCCE); the partition cache stays with
# cCCE + pcce, and the branch's method names (adding_MECCE) are refused with the new form (from the input, and from -m)
reject_case p1_pcce2_sk2 rj_pg_cache "enabled cache requires method cCCE with cluster_mode pcce" \
  "d['method']='gCCE'; d['pcce']['cache']={'mode':'readwrite','file':'x'}"
reject_case p1_pcce2_sk2 rj_old_pgcce "method \"pgCCE\" is now method \"gCCE\" with \"cluster_mode\": \"pcce\" in the input" "d['method']='pgCCE'"
reject_case natab_cce2_hahn rj_old_pmecce_cli "-m \"pmeCCE\" is now method \"meCCE\" with \"cluster_mode\": \"pcce\" in the input" "" -m pmeCCE
# gCCE with the cluster shapes that run with cCCE only: not implemented yet (from the input)
reject_case natab_dsj2 rj_gcce_dsj "method gCCE with cluster_mode \"dsj\" is not implemented yet" "d['method'] = 'gCCE'"
reject_case natab_dsj2 rj_gcce_itb "method gCCE with cluster_mode \"itb\" is not implemented yet" "d['method'] = 'gCCE'; d['cluster_mode'] = 'itb'"
reject_case natab_dsj2 rj_gcce_dsjitb "method gCCE with cluster_mode \"dsjitb\" is not implemented yet" "d['method'] = 'gCCE'; d['cluster_mode'] = 'dsjitb'"
# a hamiltonian setting of the other method: hamiltonian.cCCE in a gCCE run (written so, or after -m), hamiltonian.gCCE
reject_case p1_manual_hahn rj_hamiltonian_ccce_in_gcce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE; this run is gCCE" \
  "d.setdefault('hamiltonian', {})['cCCE'] = {'knight': True}"
reject_case natab_knight_shift rj_hamiltonian_ccce_cli_gcce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE; this run is gCCE" "" -m gCCE
reject_case natab_cce2_hahn rj_hamiltonian_gcce_key "hamiltonian.gCCE has no settings yet" "d['hamiltonian'] = {'gCCE': {'knight': True}}"
# the reserved method noise (A49): refused as not implemented yet, from the input and from -m
reject_case natab_cce2_hahn rj_method_noise "method \"noise\" is not implemented yet" "d['method'] = 'noise'"
reject_case natab_cce2_hahn rj_method_noise_cli "\"noise\" is not implemented yet" "" -m noise
# the adaptive block belongs to cluster_mode pcce or bundle; the pcce block to pcce
reject_case natab_cce2_hahn rj_adaptive_standard "the adaptive block belongs to cluster_mode pcce or bundle (this input has cluster_mode standard)" \
  "d['adaptive'] = {'output': 'x', 'atol': 1e-8, 'rtol': 1e-6}"
reject_case natab_cce2_hahn rj_pcce_block_standard "the pcce block belongs to cluster_mode pcce" "d['pcce'] = {'partition_size': 2}"
# the two new settings (A47; what they compute is in ../check_analytic/): delta_A false is for a bath of nuclear spins,
# and bath_bath_mean_field takes off, z or xyz
reject_case p1_cce2_single rj_delta_A_defect "hamiltonian.cCCE.delta_A false is for a bath of nuclear spins" \
  "d['hamiltonian'] = {'cCCE': {'delta_A': False}}"
reject_case natab_cce2_hahn rj_delta_A_overhaus "hamiltonian.cCCE.delta_A false is for a bath of nuclear spins" \
  "d['hamiltonian'] = {'cCCE': {'delta_A': False}, 'common': {'qubit_bath_mean_field': True}}"
reject_case natab_cce2_hahn rj_mean_field_value "hamiltonian.common.bath_bath_mean_field must be \"off\", \"z\" or \"xyz\"" \
  "d['hamiltonian'] = {'common': {'bath_bath_mean_field': 'vector'}}"
# one qubit state, one form
reject_case hbn_gcce2_oh rj_state_index_and_vector "Qubit.alpha_state must have one of \"index\" and \"vector\" (it has both)" \
  "d['Qubit']['alpha_state']['index'] = 0"
reject_case hbn_gcce2_oh rj_state_neither "Qubit.beta_state must have one of \"index\" and \"vector\" (it has neither)" \
  "d['Qubit']['beta_state'] = {}"

echo
echo "=== DEER bath pulses"
# pulse.bath (DEER, 2026-10-13): pulses on the bath spins of one defect species. Only gCCE and megCCE have the bath in
# their matrices (cluster_mode standard or pcce), and only their CPU paths (gCCE with propagator eigen or expm, megCCE)
# expand a pulse beyond the qubit; the block's form is checked as it is read, and the merged qubit+bath sequence before
# the run (each pulse exactly once). The refusal names the method as the input writes it (meCCE + pcce is pmeCCE inside).
DEER="d['pulse']['bath'] = {'count': 1, 'defect': 'P1', 'energy_shift': [84.946496], 'sequence': [['0.5','X',180.0],['0.5','I',0.0]]}"
reject_case p1_cce2_single rj_deer_cce "bath pulses (pulse.bath.count > 0) are implemented for methods gCCE and megCCE" "$DEER"
reject_case p1_cce2_single rj_deer_mecce "implemented for methods gCCE and megCCE (cluster_mode standard or pcce); this run is method meCCE" \
  "$DEER; d['method']='meCCE'"
reject_case p1_pcce2_sk2 rj_deer_pmecce "implemented for methods gCCE and megCCE (cluster_mode standard or pcce); this run is method meCCE" \
  "$DEER; d['method']='meCCE'"
reject_case p1_pcce2_sk2 rj_deer_pcce_ccce "implemented for methods gCCE and megCCE (cluster_mode standard or pcce); this run is method cCCE" \
  "$DEER"
reject_case p1_manual_hahn rj_deer_cuda "bath pulses (pulse.bath.count > 0) run on the CPU only" "$DEER; d['backend']='cuda'"
reject_case p1_manual_hahn rj_deer_detuning "bath pulses need pulse.detuning_factor = 1" "$DEER; d['pulse']['detuning_factor']=0.5"
reject_case p1_manual_hahn rj_deer_pulseiter "bath pulses cannot be combined with pulse.per_qubit = true" "$DEER; d['pulse']['per_qubit']=True"
reject_case p1_manual_hahn rj_deer_sparse "bath pulses need propagator eigen or expm" "$DEER; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_manual_hahn rj_deer_meg_sparse "bath pulses with method megCCE need propagator eigen (propagator sparse" \
  "$DEER; d['method']='megCCE'; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_manual_hahn rj_deer_warp "bath pulses cannot use gpu_options.warp_max_dim > 0" "$DEER; d['gpu_options'] = {'warp_max_dim': 27}"
reject_case p1_manual_hahn rj_deer_bundle "bath pulses need cluster_mode standard or pcce (this run has cluster_mode bundle)" \
  "$DEER; d['cluster_mode']='bundle'"
reject_case p1_manual_hahn rj_deer_negative "pulse.bath.count must be >= 0" "d['pulse']['bath'] = {'count': -1}"
reject_case p1_manual_hahn rj_deer_rows "pulse.bath.sequence must have pulse.bath.count+1 rows" "$DEER; d['pulse']['bath']['sequence']=[['1','I',0.0]]"
reject_case p1_manual_hahn rj_deer_sum "the fractions of pulse.bath.sequence sum to" "$DEER; d['pulse']['bath']['sequence']=[['0.5','X',180.0],['0.4','I',0.0]]"
reject_case p1_manual_hahn rj_deer_defect "pulse.bath.defect is required when pulse.bath.count > 0" "$DEER; del d['pulse']['bath']['defect']"
reject_case p1_manual_hahn rj_deer_tones "pulse.bath.energy_shift needs at least one finite number" "$DEER; d['pulse']['bath']['energy_shift']=[]"
reject_case p1_manual_hahn rj_deer_unknown "pulse.bath has no key \"tone\"" "$DEER; d['pulse']['bath']['tone']=1"
# the pulse of the last row would come after the last interval: never applied, so a non-identity one is refused
reject_case p1_manual_hahn rj_deer_last_row "the last row of pulse.bath.sequence must be identity" "$DEER; d['pulse']['bath']['sequence']=[['0.5','X',180.0],['0.5','X',180.0]]"
# two pulses of one channel within 1e-9 would share one edge of the merged sequence and one of them be lost; the chained
# case: two bath pulses 1.8e-9 apart, each within 1e-9 of the qubit pi at 0.5
reject_case p1_manual_hahn rj_deer_two_qubit_pulses "pulse.sequence has two qubit pulses at the same time; bath pulses need them apart" \
  "$DEER; d['pulse'].update({'count': 2, 'name': 'Manual', 'sequence': [['0.5','X',90.0],['0','Y',90.0],['0.5','I',0.0]]})"
reject_case p1_manual_hahn rj_deer_two_bath_pulses "pulse.bath.sequence has two bath pulses at the same time" \
  "$DEER; d['pulse']['bath'].update({'count': 2, 'sequence': [['0.5','X',90.0],['0','Y',90.0],['0.5','I',0.0]]})"
reject_case p1_manual_hahn rj_deer_two_bath_pulses_chained "pulse.bath.sequence has two bath pulses at the same time" \
  "$DEER; d['pulse']['bath'].update({'count': 2, 'sequence': [['0.4999999991','X',90.0],['0.0000000018','Y',90.0],['0.4999999991','I',0.0]]})"
# the branch's flat keys (adding_MECCE): refused with the new name (otherwise a warning and a run without the pump)
reject_case natab_cce2_hahn rj_old_bnpulse "${OLDK}bnpulse\" is now \"pulse.bath.count\"" "d['bnpulse'] = 1"
reject_case natab_cce2_hahn rj_old_bsequence "${OLDK}bsequence\" is now \"pulse.bath.sequence\"" "d['bsequence'] = [['1','I',0.0]]"
reject_case natab_cce2_hahn rj_old_bpulse_defect "${OLDK}bpulse_defect\" is now \"pulse.bath.defect\"" "d['bpulse_defect'] = 'P1'"
reject_case natab_cce2_hahn rj_old_bpulse_energy_shift "${OLDK}bpulse_energy_shift\" is now \"pulse.bath.energy_shift\"" "d['bpulse_energy_shift'] = [84.946496]"
reject_case natab_cce2_hahn rj_old_bspin "${OLDK}bspin\" is now \"pulse.bath.spin\"" "d['bspin'] = 0.5"
reject_case natab_cce2_hahn rj_old_balphams "${OLDK}balphams\" is now \"pulse.bath.alphams\"" "d['balphams'] = 0.5"
reject_case natab_cce2_hahn rj_old_bbetams "${OLDK}bbetams\" is now \"pulse.bath.betams\"" "d['bbetams'] = -0.5"

echo
echo "=== master equation (meCCE, megCCE)"
# meCCE and megCCE (adding_MECCE, 2026-10-13): Lindblad dissipation by the jump operators of hamiltonian.<method>, on the
# CPU with a matrix exponential per step. The section belongs to its method; the branch's top-level key is refused with
# the new place; what the ME kernels do not take is refused before the bath is read.
MEJ(){ echo "d.setdefault('hamiltonian', {})['$1'] = {'jump_operators': $2}"; }   # MEJ <section> <python list>
JUMP="$(MEJ meCCE "[{'bath_name':'e','operator':'-','rate':1e-4}]")"
reject_case p1_manual_hahn rj_me_section_gcce "hamiltonian.meCCE (jump_operators) belongs to method meCCE; this run is gCCE" "$JUMP"
reject_case p1_cce2_single rj_me_section_other "hamiltonian.megCCE (jump_operators) belongs to method megCCE; this run is meCCE" \
  "d['method']='meCCE'; $(MEJ megCCE "[{'bath_name':'e','operator':'-','rate':1e-4}]")"
# the run's method as the input writes it: the default cCCE, not the internal name CCE (an old name the reader refuses)
reject_case p1_cce2_single rj_me_section_ccce "hamiltonian.meCCE (jump_operators) belongs to method meCCE; this run is cCCE" \
  "d.pop('method'); $JUMP"
reject_case p1_cce2_single rj_me_section_key "hamiltonian.meCCE has no key \"knight\" (jump_operators)" \
  "d['method']='meCCE'; d['hamiltonian'] = {'meCCE': {'knight': True}}"
reject_case p1_cce2_single rj_old_jump_operators "${OLDK}jump_operators\" is now \"hamiltonian.meCCE.jump_operators\"" \
  "d['jump_operators'] = [{'bath_name':'e','operator':'-','rate':1e-4}]"
reject_case p1_cce2_single rj_me_op "jump_operators[0].operator must be one of - + z x y" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'bath_name':'e','operator':'q','rate':1e-4}]")"
reject_case p1_cce2_single rj_me_old_notation "the old names +- and -+ are now - and +" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'bath_name':'e','operator':'+-','rate':1e-4}]")"
reject_case p1_cce2_single rj_me_rate "jump_operators[0].rate must be a finite number >= 0" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'bath_name':'e','operator':'z','rate':-1}]")"
reject_case p1_cce2_single rj_me_name "jump_operators[0] needs bath_name, operator and rate" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'operator':'z','rate':1e-4}]")"
reject_case natab_manual_hahn rj_me_manual "method meCCE takes its pulses from pulse.count" \
  "d['method']='meCCE'; d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',180],['1/2','I',0]]}"
reject_case p1_cce2_single rj_me_recursive "this run has method=meCCE. Use evolution_update=exact" \
  "d['method']='meCCE'; d['evolution_update']='recursive'"
reject_case p1_cce2_single rj_me_qubit_jump "method meCCE does not take qubit jump operators" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'bath_name':'qubit','operator':'z','rate':1e-4}]")"
# meCCE with cluster_mode pcce (inside pmeCCE) projects the qubit out as well; the message names the method as written
reject_case p1_pcce2_sk2 rj_pme_qubit_jump "method meCCE does not take qubit jump operators" \
  "d['method']='meCCE'; $(MEJ meCCE "[{'bath_name':'qubit','operator':'z','rate':1e-4}]")"
# the method as the input writes it in every refusal a gCCE / meCCE / megCCE + pcce run (inside pgCCE, pmeCCE, pmegCCE) can
# meet: the internal names must not reach the messages
reject_case p1_pcce2_sk2 rj_pme_recursive "this run has method=meCCE" "d['method']='meCCE'; d['evolution_update']='recursive'"
reject_case p1_pcce2_sk2 rj_pme_vector "evolution=vector is implemented for method=gCCE only (this run has method=meCCE)" \
  "d['method']='meCCE'; d['evolution']='vector'"
reject_case p1_pcce2_sk2 rj_pmeg_expm "propagator=expm is implemented for method=gCCE only (this run has method=megCCE)" \
  "d['method']='megCCE'; d['propagator']='expm'"
reject_case p1_pcce2_sk2 rj_pme_pulse "method meCCE calculates every pulse as an ideal pi flip" \
  "d['method']='meCCE'; d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',90],['1/2','I',0]]}"
reject_case p1_pcce2_sk2 rj_pme_rdip_ladder "rdip_ladder needs method cCCE (cluster_mode standard or pcce) or gCCE (cluster_mode standard); this run is method meCCE" \
  "d['method']='meCCE'; d['rdip_ladder']=[600]"
reject_case p1_pcce2_sk2 rj_pg_packed "gpu_options.packed_max_dim works with method=CCE (this run: gCCE, cluster_mode pcce)" \
  "d['method']='gCCE'; d['gpu_options']={'packed_max_dim': 4}"
reject_case p1_pcce2_sk2 rj_pg_warp "gpu_options.warp_max_dim works with method=CCE or gCCE (this run: gCCE, cluster_mode pcce)" \
  "d['method']='gCCE'; d['gpu_options']={'warp_max_dim': 27}"
for name in rj_pme_qubit_jump rj_pme_recursive rj_pme_vector rj_pmeg_expm rj_pme_pulse rj_pme_rdip_ladder rj_pg_packed rj_pg_warp \
            rj_deer_pmecce; do
    if grep -i "error" "$WORK/evolution_update_refusals/proc_$name" | grep -qiE "p(g|me|meg)cce"; then
        echo "  FAIL  $name : an internal method name in the message"; FAIL=$((FAIL+1))
    else
        echo "  PASS  $name names no internal method"; PASS=$((PASS+1))
    fi
done
reject_case p1_cce2_single rj_me_sparse "method meCCE needs propagator eigen" \
  "d['method']='meCCE'; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_cce2_single rj_me_shape "method meCCE with cluster_mode \"dsj\" is not implemented" \
  "d['method']='meCCE'; d['cluster_mode']='dsj'"
reject_case p1_manual_hahn rj_meg_pulseiter "method megCCE needs pulse.per_qubit = false" "d['method']='megCCE'; d['pulse']['per_qubit']=True"
reject_case p1_manual_hahn rj_meg_detuning "method megCCE needs pulse.detuning_factor = 1" "d['method']='megCCE'; d['pulse']['detuning_factor']=0.5"
reject_case p1_manual_hahn rj_meg_density_matrix "method megCCE does not write savemode cluster_density_matrix" \
  "d['method']='megCCE'; d['savemode']='cluster_density_matrix'"
reject_case p1_manual_hahn rj_meg_ccce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE; this run is megCCE" \
  "d['method']='megCCE'; d.setdefault('hamiltonian', {})['cCCE'] = {'knight': True}"
# quantity dm is not implemented for the ME methods (in the input, and from -q after -m)
reject_case p1_cce2_single rj_me_dm "method meCCE computes quantity coherence only" "d['method']='meCCE'; d['quantity']='dm'"
reject_case p1_manual_hahn rj_meg_dm_cli "method megCCE computes quantity coherence only" "" -m megCCE -q dm
# the pcce block in an ME run names the run's cluster_mode (standard), not the method
reject_case natab_cce2_hahn rj_me_pcce_block "the pcce block belongs to cluster_mode pcce (this input has cluster_mode standard)" \
  "d['method']='meCCE'; d['pcce'] = {'partition_size': 2}"
# rdip_ladder assembles cCCE, gCCE and pCCE products only; refused for ME before the bath is read (the ladder's own check
# comes after the clusters)
reject_case p1_cce2_single rj_me_rdip_ladder "rdip_ladder needs method cCCE (cluster_mode standard or pcce) or gCCE (cluster_mode standard); this run is method meCCE" \
  "d['method']='meCCE'; d['rdip_ladder']=[600]"
# guarded already by v2's own whitelists: checked here so a later change cannot open them for ME silently
reject_case p1_manual_hahn rj_meg_cuda "backend=cuda supports method cCCE (cluster_mode standard or pcce), or method gCCE with propagator=eigen" \
  "d['method']='megCCE'; d['backend']='cuda'"
reject_case p1_manual_hahn rj_meg_vector "evolution=vector is implemented for method=gCCE only" "d['method']='megCCE'; d['evolution']='vector'"
reject_case p1_manual_hahn rj_meg_expm "propagator=expm is implemented for method=gCCE only (this run has method=megCCE)" \
  "d['method']='megCCE'; d['propagator']='expm'"
reject_case p1_cce2_single rj_me_plan "-P supports cCCE, gCCE and cluster_mode pcce plans" "d['method']='meCCE'" -P
# the qubit count (two NV, S5): meCCE takes one qubit, and a qubit jump operator needs one qubit
reject S5 rj_me_two_qubits "method meCCE supports one qubit" "d['method']='meCCE'"
reject S5 rj_me_qubit_jump_two "a qubit jump operator needs one qubit" \
  "d['method']='megCCE'; d['hamiltonian'] = {'megCCE': {'jump_operators': [{'bath_name':'qubit','operator':'z','rate':1e-4}]}}"
# a jump whose bath_name is neither the qubit nor a bath spin of the run (bath file and Defect sub-spins) would run at
# rate 0 without a word: refused once the bath is read, with the names the bath has. meCCE on check_master_equation's
# S=1/2 electron-spin system (spins named e; the first entry matches, the second does not), megCCE on the P1 bath
# (nstep and rbath small, so that a build without the check runs these in seconds)
MEI="$REPO/test/code_verification/check_master_equation/inputs"
reject_case p1_cce2_single rj_me_jump_name 'hamiltonian.meCCE.jump_operators[1].bath_name "P1x" names no bath spin of this run (bath spins: e;' \
  "d = {'method': 'meCCE', 'quantity': 'coherence', 'Qubit': {'file': '$MEI/qubit_espin', 'spin': 0.5, 'alphams': 0.5, 'betams': -0.5}, 'gyrofile': '$MEI/gyro_espin', 'bathfile': ['$MEI/bath_espin'], 'order': 2, 'bfield': 300, 'rbath': 500, 'rdip': 300, 'deltat': 0.04, 'nstep': 5, 'nstate': 0, 'pulse': {'count': 1}}; $(MEJ meCCE "[{'bath_name':'e','operator':'-','rate':1e-4},{'bath_name':'P1x','operator':'+','rate':1e-4}]")"
reject_case p1_manual_hahn rj_meg_jump_name 'hamiltonian.megCCE.jump_operators[0].bath_name "P1x" names no bath spin of this run (bath spins: P1;' \
  "d['method']='megCCE'; d['nstep']=5; d['rbath']=200; $(MEJ megCCE "[{'bath_name':'P1x','operator':'z','rate':1e-4}]")"
# a bath_name that is the qubit's name and a bath spin's: the megCCE kernel takes it as a qubit jump and drops the bath
# spins, so it is refused once the bath is read (the qubit of p1_gcce2_zfs named P1, as its bath spins)
reject_case p1_gcce2_zfs rj_meg_jump_both_names 'hamiltonian.megCCE.jump_operators[0].bath_name "P1" names both the qubit and bath spins of this run' \
  "d['method']='megCCE'; d['nstep']=5; d['Qubit']['qubit'][0]['name']='P1'; d['Qubit']['intmap'][0]['between']=['P1','P1']; $(MEJ megCCE "[{'bath_name':'P1','operator':'z','rate':1e-4}]")"
# the jump sections as they are read: one method's section only, and jump_operators a list
reject_case p1_cce2_single rj_me_jump_both_sections "jump operators are given in both hamiltonian.meCCE and hamiltonian.megCCE" \
  "d['method']='meCCE'; d['hamiltonian'] = {'meCCE': {'jump_operators': [{'bath_name':'P1','operator':'-','rate':1e-4}]}, 'megCCE': {'jump_operators': [{'bath_name':'P1','operator':'-','rate':1e-4}]}}"
reject_case p1_cce2_single rj_me_jump_not_array "hamiltonian.meCCE.jump_operators must be an array" \
  "d['method']='meCCE'; d['hamiltonian'] = {'meCCE': {'jump_operators': {'bath_name':'P1','operator':'-','rate':1e-4}}}"
# a rotation other than an ideal pi flip names the method that applies it: megCCE for meCCE, gCCE for cCCE
reject_case p1_cce2_single rj_me_pulse_rotation "is X 90. Use method megCCE for other rotations." \
  "d['method']='meCCE'; d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',90],['1/2','I',0]]}"
reject_case p1_cce2_single rj_ccce_pulse_rotation "is X 90. Use method gCCE for other rotations." \
  "d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',90],['1/2','I',0]]}"

echo
echo "=== -f and -m as getopt reads them, keys without regard to case (review 20261013-180049 I01-I03)"
# I01: the input file was found by a search for an argument equal to "-f", so "-fFILE" (and "-vfFILE") skipped the adaptive
# block: a bundle input with one ran as a plain bundle run (curve off by up to 2.84), and adaptive with cluster_mode
# standard was not refused. The file is found by getopt now; the refusal in both attached forms, then the merge-grow
# regression input with an adaptive block run spaced and attached: both refine, to the same final curve.
reject_attached(){    # reject_attached <option word before the file name> <name>: rj_adaptive_standard's input, the file attached
    local w="$WORK/evolution_update_refusals" rc
    ( cd "$REPO/test/code_verification/check_numerics/cases/natab_cce2_hahn" &&
      "$MPIRUN" -n 1 "$BIN_DEV" "$1$w/rj_adaptive_standard.json" -o "$w/$2" > "$w/proc_$2" 2>&1 < /dev/null )
    rc=$?
    if [ $rc -ne 0 ] && grep -qF "the adaptive block belongs to cluster_mode pcce or bundle (this input has cluster_mode standard)" "$w/proc_$2"
    then echo "  PASS  $2 (exit=$rc)"; PASS=$((PASS+1))
    else echo "  FAIL  $2 : exit=$rc, see $w/proc_$2"; FAIL=$((FAIL+1)); fi
}
reject_attached -f rj_adaptive_standard_attached
reject_attached -vf rj_adaptive_standard_combined
adaptive_attached(){
    local w="$WORK/adaptive_attached" f ok=1; rm -rf "${w:?}"; mkdir -p "$w"
    python3 - "$REPO/test/cluster_verification/check_adaptive/regression/inputs/merge-grow-0.json" "$w" <<'PY'
import json, sys
src, w = sys.argv[1:3]
d = json.load(open(src))
for form in ("spaced", "attached", "upper"):
    # upper: the adaptive keys in another case (I03), read as the lower-case ones
    keys = ("Policy", "Output", "ATOL", "Rtol") if form == "upper" else ("policy", "output", "atol", "rtol")
    d["adaptive"] = dict(zip(keys, ("merge", f"{w}/{form}_adaptive", 1e-8, 2e-6)))
    d["outfile"] = f"{w}/{form}"
    json.dump(d, open(f"{w}/{form}.json", "w"), indent=1)
PY
    # relative paths in the input are the repository's
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/spaced.json" > "$w/proc_spaced" 2>&1 < /dev/null ) || ok=0
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" "-f$w/attached.json" > "$w/proc_attached" 2>&1 < /dev/null ) || ok=0
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/upper.json" > "$w/proc_upper" 2>&1 < /dev/null ) || ok=0
    for f in spaced attached upper; do
        grep -q '"status": "validated_full_bath_model"' "$w/${f}_adaptive/summary.json" 2>/dev/null || ok=0
    done
    for f in attached upper; do
        cmp -s "$w/spaced_adaptive/final-bundles.json" "$w/${f}_adaptive/final-bundles.json" || ok=0
        cmp -s "$w/spaced_adaptive/0009-run/coherence_wiDiv" "$w/${f}_adaptive/0009-run/coherence_wiDiv" || ok=0
    done
    if [ $ok -eq 1 ]; then echo "  PASS  adaptive_attached (-f FILE, -fFILE and Policy/Output/ATOL/Rtol refine alike)"; PASS=$((PASS+1))
    else echo "  FAIL  adaptive_attached : see $w"; FAIL=$((FAIL+1)); fi
}
adaptive_attached
# accept_same <case> <name> <python edit> <python edit of the reference> [command-line options]: the edited input runs,
# and its raw result is the reference's, byte for byte
accept_same(){
    local case="$1" name="$2" edit="$3" ref="$4"; shift 4
    local w="$WORK/accepted_forms" c="$REPO/test/code_verification/check_numerics/cases/$case" f ok=1; mkdir -p "$w"
    python3 - "$c/ccein.json" "$w/$name" <<PY
import json, sys
d = json.load(open(sys.argv[1])); $edit
json.dump(d, open(sys.argv[2] + ".json", "w"), indent=1)
d = json.load(open(sys.argv[1])); $ref
json.dump(d, open(sys.argv[2] + "_ref.json", "w"), indent=1)
PY
    rm -f "$w/$name"_*Div "$w/$name"_state* "$w/$name"_ref_*
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/${name}_ref.json" -o "$w/${name}_ref" > "$w/proc_${name}_ref" 2>&1 < /dev/null ) || ok=0
    # options before a word -f go before the file, the rest after it
    local pre=() post=("$@") i
    for ((i=0; i<$#; i++)); do [ "${post[i]}" = -f ] && { pre=("${post[@]:0:i}"); post=("${post[@]:i+1}"); break; }; done
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" "${pre[@]}" -f "$w/$name.json" -o "$w/$name" "${post[@]}" > "$w/proc_$name" 2>&1 < /dev/null ) || ok=0
    local n=0
    for f in "$w/${name}_ref"_*; do
        n=$((n+1)); cmp -s "$f" "$w/$name${f#$w/${name}_ref}" || ok=0
    done
    [ $n -gt 0 ] || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  $name"; PASS=$((PASS+1))
    else echo "  FAIL  $name : see $w/proc_$name"; FAIL=$((FAIL+1)); fi
}
# I02: -m replaces the file's method before the combination with cluster_mode is judged: a file with gCCE + pcce runs
# under -m cCCE as the cCCE + pcce file does (before or after -f); the other way round runs as gCCE + pcce (pgCCE, the
# master-equation port's cluster_mode pcce for gCCE)
accept_same p1_pcce2_sk2 ac_gcce_pcce_cli_ccce "d['method'] = 'gCCE'" "" -m cCCE
accept_same p1_pcce2_sk2 ac_gcce_pcce_cli_ccce_first "d['method'] = 'gCCE'" "" -m cCCE -f
# the file's method is still read as a method name under -m: a reserved one stays refused (noise; meCCE runs since the
# master-equation port)
reject_case p1_pcce2_sk2 rj_file_method_under_cli "method \"noise\" is not implemented yet" "d['method'] = 'noise'" -m cCCE
# I03: the keys of the savemode and propagator blocks in another case
accept_same natab_cce2_hahn ac_savemode_keys_case "d['savemode'] = {'Mode': 'coherence', 'NoDiv': True}" \
  "d['savemode'] = {'mode': 'coherence', 'nodiv': True}"
accept_same natab_cce2_hahn ac_propagator_keys_case "d['propagator'] = {'Method': 'eigen'}" "d['propagator'] = {'method': 'eigen'}"

echo
echo "==========================================="
echo " refused combinations : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
