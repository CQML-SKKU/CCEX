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

# refused_with <name> <exit status> <expected message fragment> <output file>: the run was refused with that message;
# otherwise the failure is counted and shown, and the status is 1
refused_with(){
    local name="$1" rc="$2" want="$3" proc="$4"
    if [ "$rc" -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return 1; fi
    if ! grep -qF -- "$want" "$proc"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$proc" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return 1
    fi
}

# reject <sys> <name> <expected message fragment> <python edit> [read]
# "read": the refusal is about a tensor file, so it comes after a reader has started
reject(){
    local sys="$1" name="$2" want="$3" edit="$4" when="${5:-}"
    local w="$WORK/$sys"
    [ -f "$w/legacy.json" ] || { echo "  FAIL  $name : no legacy.json"; FAIL=$((FAIL+1)); return; }
    python3 - "$w" "$name" <<PY
import json, sys; from input_layout import flat, nest
w, name = sys.argv[1:3]
d = flat(json.load(open(f"{w}/legacy.json")))
$edit
d["outfile"] = f"./{name}"
json.dump(nest(d), open(f"{w}/{name}.json","w"), indent=4)
PY
    # S3 and S5 passed their state files and -N 1 on the command line as well until A63; legacy.json holds the same keys
    ( cd "$w" && run_ccex "$BIN_DEV" "$name.json" "proc_$name" )
    local rc=$?
    refused_with "$name" $rc "$want" "$w/proc_$name" || return
    # the guard must run before either tensor reader announces itself
    if [ "$when" != read ] && grep -qE "Read the (Hyperfine|Quadrupole) interaction from DFT inputfile" "$w/proc_$name"; then
        echo "  FAIL  $name : refused, but a DFT tensor file was read first"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}

R111='{"enabled": True, "qubit_axis": [1.0,1.0,1.0]}'
# A90: the frames of the tensor files are the "frame" keys of their hamiltonian blocks
HF_FILE="d['hamiltonian']['common']['hyperfine']['frame'] = 'file'"
QD_FILE="d['hamiltonian']['common']['quadrupole']['frame'] = 'file'"

# A valid input with a command-line option that was removed (A63): refused, with the input key that replaces it, before
# the input is read. (Until A63 this section checked -N and -B values that are not numbers, review 78d583 R5; the keys
# nstate and bfield have their own checks.)
reject_cli(){
    local sys="$1" name="$2" want="$3"; shift 3
    local w="$WORK/$sys"
    ( cd "$w" && run_ccex "$BIN_DEV" "${CLI_INPUT:-legacy.json}" "proc_$name" "$@" )
    local rc=$?
    refused_with "$name" $rc "$want" "$w/proc_$name" || return
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
echo "=== S1a : removed command-line options (A63)"
reject_cli S1a rj_cli_removed_bathfile '-I was removed: give the bath files in the input, "configuration": {"bathfile"' -I bath
reject_cli S1a rj_cli_removed_nstate '-N was removed: give the bath states in the input, "configuration": {"s_range": [1, N]}' -N 3
reject_cli S1a rj_cli_removed_method '-m was removed: give the method in the input, "method"' -m gCCE
reject_cli S1a rj_cli_removed_quantity '-q was removed: give the quantity in the input, "quantity"' -q dm
reject_cli S1a rj_cli_removed_bfield '-B was removed: give the field in the input, "bfield"' -B 500
reject_cli S1a rj_cli_removed_statefile '-s was removed: give the bath state files in the input, "configuration": {"bath_state_file"' -s state
reject_cli S1a rj_cli_removed_defect_state "-S was removed: give the defect state files in the input, \"configuration\": {\"defect_state_file\"" -S state
reject_cli S1a rj_cli_removed_defect_axis "-a was removed: give the defect axis files in the input, \"configuration\": {\"defect_axis_file\"" -a axis
reject_cli S1a rj_cli_removed_plan '-P was removed: write "output.cluster_only": true in the input' -P
reject_cli S1a rj_cli_removed_plan_combined '-P was removed' -vP
reject_cli S1a rj_cli_removed_auto_conv '--auto_conv was removed: an "auto_convergence" block in the input' --auto_conv
reject_cli S1a rj_cli_removed_adaptive '--adaptive was removed: write the refinement options as a "cluster.adaptive" block' --adaptive
reject_cli S1a rj_cli_removed_auto_status '--auto_status was removed' --auto_status
reject_cli S1a rj_cli_unknown "unknown option -x" -x
# -o DIR (A63) is an output folder: an absolute outfile beside it names a second place; a file where the folder should be
python3 - "$WORK/S1a" <<'PY'
import json, sys; from input_layout import flat, nest
d = flat(json.load(open(f"{sys.argv[1]}/legacy.json"))); d["outfile"] = f"{sys.argv[1]}/abs"
json.dump(nest(d), open(f"{sys.argv[1]}/abs_outfile.json", "w"), indent=4)
PY
CLI_INPUT=abs_outfile.json reject_cli S1a rj_cli_outfile_absolute "is an absolute path and -o" -o "$WORK/S1a/out_abs"
: > "$WORK/S1a/not_a_folder"
reject_cli S1a rj_cli_output_not_folder "cannot create the output folder" -o "$WORK/S1a/not_a_folder"
# an output folder that cannot be written (A73): refused at the start, not when the first output is written
mkdir -p "$WORK/S1a/read_only"; chmod 555 "$WORK/S1a/read_only"
if [ -w "$WORK/S1a/read_only" ]; then echo "  SKIP  rj_output_not_writable (the folder stays writable: root?)"
else
    reject_cli S1a rj_output_not_writable "cannot write the outputs" -o "$WORK/S1a/read_only"
    if grep -q "Clusterize" "$WORK/S1a/proc_rj_output_not_writable"; then
        echo "  FAIL  rj_output_not_writable : refused after the calculation started"; FAIL=$((FAIL+1)); fi
fi
chmod 755 "$WORK/S1a/read_only"
# -e (A63) writes its example and calculates nothing: not with -f, and never over an existing file
reject_cli S1a rj_cli_example_with_input "-e writes an example input and calculates nothing, leave -f out" -e
example_exists(){
    local w="$WORK/S1a/example_exists" rc; rm -rf "${w:?}"; mkdir -p "$w"; : > "$w/ccex_example.jsonc"
    ( cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -e > proc 2>&1 < /dev/null ); rc=$?
    if [ $rc -ne 0 ] && grep -qF "ccex_example.jsonc exists, -e does not overwrite it" "$w/proc" && [ ! -s "$w/ccex_example.jsonc" ]
    then echo "  PASS  rj_cli_example_exists (exit=$rc)"; PASS=$((PASS+1))
    else echo "  FAIL  rj_cli_example_exists : exit=$rc, see $w/proc"; FAIL=$((FAIL+1)); fi
}
example_exists
# a qubit position file (Qubit.file until A94, qubitfile until A47) the input names but that is not there ran the qubit
# at (0,0,0) (review 4dd5fd R1)
reject S1a rj_qubitfile_missing 'Qubit.qubit[0].file: cannot read' "d['Qubit']['qubit'][0]['file'] += '.missing'"
# A94: a Qubit.qubit entry has its position from "file" or "xyz", not both (Qubit.file beside Qubit.qubit until then:
# the block was ignored without a word, so a "two-qubit" tutorial ran one qubit, review 20261013-readme R01)
reject S1a rj_qubitfile_and_qubit_block 'Qubit.qubit[0] has both "file" and "xyz", the position is one of them' \
  "d['Qubit']['qubit'][0]['xyz'] = [0.0, 0.0, 0.0]"
# A90: with the rotation on a single qubit's Qubit.intmap entry states its frame; the ZFS of the Qubit.file form is such
# an entry since A93 (Qubit.zfs and its Qubit.frame until then)
reject S1a rj_qubit_zfs_frame_missing 'frame.rotation is on, so Qubit.intmap["q0","q0"].frame must say' \
  "d['coordinate_frame_rotation'] = $R111
d['Qubit']['intmap'] = [{'between': ['q0', 'q0'], 'tensor': [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]], 'unit': 'kHz'}]"
reject S1a rj_one_qubit_intmap_frame_missing 'frame.rotation is on, so Qubit.intmap["q0","q0"].frame must say' \
  "d['coordinate_frame_rotation'] = $R111
d['Qubit'] = {'qubit': [{'name': 'q0', 'spin': 1.0, 'gyro': -17608.59705, 'xyz': [0.0, 0.0, 0.0]}],
              'intmap': [{'between': ['q0', 'q0'], 'tensor': [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]}]}"
# quantity noise was accepted and computed the coherence (review 20261013-doxygen D22)
reject S1a rj_quantity_noise "quantity \"noise\" is not supported" "d['quantity'] = 'noise'"
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
import json, sys; from input_layout import flat, nest
w = sys.argv[1]; d = flat(json.load(open(f"{w}/legacy.json")))
d.update(method='gCCE', pulse={'count': 2, 'name': 'Manual', 'sequence': [[0.25,'I',0],['1/2','X',90.5],[0.25,'I',0]]}, outfile='./seq_ok',
         cluster_only=True)
json.dump(nest(d), open(f"{w}/seq_ok.json", "w"))
PY
    if (cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f seq_ok.json < /dev/null > proc_seq_ok 2>&1) &&
       grep -Eq '0\.250 +0\.750 +0\.500 +X +90\.50' "$w/proc_seq_ok"; then
        echo "  PASS  sequence_accepted (a number fraction, an angle of 90.5)"; PASS=$((PASS+1))
    else echo "  FAIL  sequence_accepted : see $w/proc_seq_ok"; FAIL=$((FAIL+1)); fi
}
sequence_accepted
# A86: WAHUHA writes five rows whatever the count while the arrays have count + 1 rows (counts 0-3 overflowed, 5 gave
# another curve without a word); a sequence with a name other than Manual left the arrays zero (gCCE aborted with
# "Invalid axis character"); a Manual row count that does not match ended in an uncaught throw
for c in 1 3 5; do
    reject S1a rj_wahuha_count$c "pulse.name WAHUHA needs pulse.count 4 (the four-pulse sequence) (this input has $c)" \
      "d['method'] = 'gCCE'; d['pulse'] = {'count': $c, 'name': 'WAHUHA'}"
done
WAHUHA_ROWS="[['1/6','X',90],['1/6','Y',270],['2/6','Y',90],['1/6','X',270],['1/6','I',0]]"
reject S1a rj_seq_wahuha "pulse.sequence is read only with pulse.name Manual (pulse.name WAHUHA)" \
  "d['method'] = 'gCCE'; d['pulse'] = {'count': 4, 'name': 'WAHUHA', 'sequence': $WAHUHA_ROWS}"
reject S1a rj_seq_no_name "pulse.sequence is read only with pulse.name Manual (pulse.name None)" \
  "d['pulse'] = {'count': 1, 'sequence': [['1/2','X',180],['1/2','I',0]]}"
reject S1a rj_seq_rows "pulse.sequence has 2 rows, pulse.count + 1 = 3" "$MANUAL; d['pulse']['sequence'] = [['1/2','X',180],['1/2','I',0]]"
# WAHUHA with count 4 still runs to the end (gCCE order 1, a small bath) and logs its five rows
wahuha4_runs(){
    local w="$WORK/S1a"
    python3 - "$w" <<'PY'
import json, sys; from input_layout import flat, nest
w = sys.argv[1]; d = flat(json.load(open(f"{w}/legacy.json")))
d.update(method='gCCE', order=1, rbath=15, nstep=10, pulse={'count': 4, 'name': 'WAHUHA'}, outfile='./wahuha4')
json.dump(nest(d), open(f"{w}/wahuha4.json", "w"))
PY
    rm -f "$w"/wahuha4_wiDiv
    if (cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f wahuha4.json < /dev/null > proc_wahuha4 2>&1) && [ -s "$w/wahuha4_wiDiv" ] &&
       grep -Eq ' Y +270' "$w/proc_wahuha4" && grep -Eq ' X +270' "$w/proc_wahuha4"; then
        echo "  PASS  wahuha4_runs (count 4, exit 0)"; PASS=$((PASS+1))
    else echo "  FAIL  wahuha4_runs : see $w/proc_wahuha4"; FAIL=$((FAIL+1)); fi
}
wahuha4_runs
# a finite double outside the float it is stored in became inf (review fa7f39 R3)
reject S1a rj_bfield_float_overflow "within the float range" "d['bfield'] = 1e40"

echo "=== S4 : DFT tensor file + rotation"
reject S4 rj_hf_frame_missing "frame.rotation is on, so hamiltonian.common.hyperfine.frame must say the frame of" \
  "d['coordinate_frame_rotation'] = $R111
$QD_FILE"
reject S4 rj_qd_frame_missing "frame.rotation is on, so hamiltonian.common.quadrupole.frame must say the frame of" \
  "d['coordinate_frame_rotation'] = $R111
$HF_FILE"
reject S4 rj_tilted_bfield "bfield" \
  "d['coordinate_frame_rotation'] = $R111
$HF_FILE
$QD_FILE
d['bfield'] = [100.0, 0.0, 30000.0]"
# A90: the keys removed from frame.rotation and the old places and value of the tensor-file frames
reject S4 rj_a90_bath_axis '"frame.rotation.bath_axis" was removed' \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['bath_axis'] = [0.0, 0.0, 1.0]
$HF_FILE
$QD_FILE"
reject S4 rj_a90_position_frame '"frame.rotation.qubit_position_frame" was removed' \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'
$HF_FILE
$QD_FILE"
reject S4 rj_a90_hf_tensor_frame '"frame.rotation.hf_tensor_frame" is now "hamiltonian.common.hyperfine.frame"' \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'
$QD_FILE"
reject S4 rj_a90_qd_tensor_frame '"frame.rotation.qd_tensor_frame" is now "hamiltonian.common.quadrupole.frame"' \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'qubit'
$HF_FILE"
reject S4 rj_a90_hf_frame_bath 'hamiltonian.common.hyperfine.frame "bath" is now "file"' \
  "d['coordinate_frame_rotation'] = $R111
d['hamiltonian']['common']['hyperfine']['frame'] = 'bath'
$QD_FILE"
reject S4 rj_a90_qd_frame_unknown 'hamiltonian.common.quadrupole.frame must be "file"' \
  "d['coordinate_frame_rotation'] = $R111
$HF_FILE
d['hamiltonian']['common']['quadrupole']['frame'] = 'crystal'"

echo
echo "=== S5 : two qubits + rotation"
reject S5 rj_a90_reference_qubit '"frame.rotation.reference_qubit" was removed' \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['reference_qubit'] = 'q1'
for e in d['Qubit']['intmap']: e['frame'] = 'file'"
reject S5 rj_intmap_frame_missing "frame.rotation is on, so Qubit.intmap[" \
  "d['coordinate_frame_rotation'] = $R111"
reject S5 rj_a90_intmap_tensor_frame '"Qubit.intmap[0].tensor_frame" is now "frame"' \
  "d['coordinate_frame_rotation'] = $R111
for e in d['Qubit']['intmap']: e['tensor_frame'] = 'bath'"
reject S5 rj_a90_intmap_frame_bath '.frame "bath" is now "file"' \
  "d['coordinate_frame_rotation'] = $R111
for e in d['Qubit']['intmap']: e['frame'] = 'bath'"
# A93 (2026-10-16): the qubit states are Qubit.alpha_state/beta_state only, each one of {"ms": [one m_s per qubit]},
# {"index": n} and {"vector": [...]}; the per-qubit alphams/betams (and Qubit.alphams/betams of the Qubit.file form) are
# refused with the new form. (Until A93 an index beside a written alphams/betams that named another state stopped the
# run, D43; the two can no longer be written together.)
reject S5 rj_a93_old_qubit_alphams 'input key "Qubit.qubit[1].alphams" is now "Qubit.alpha_state": {"ms": [m_1, m_2, ...]}' \
  "d['Qubit']['qubit'][1]['alphams'] = 1"
reject S5 rj_a93_old_qubit_betams 'input key "Qubit.qubit[0].betams" is now "Qubit.beta_state": {"ms": [m_1, m_2, ...]}' \
  "d['Qubit']['qubit'][0]['betams'] = 0"
reject S1a rj_a93_old_file_alphams 'input key "Qubit.alphams" is now "Qubit.alpha_state": {"ms": [m]}' \
  "d['Qubit']['alphams'] = d['Qubit'].pop('alpha_state')['ms'][0]"
reject S1a rj_a93_old_file_betams 'input key "Qubit.betams" is now "Qubit.beta_state": {"ms": [m]}' \
  "d['Qubit']['betams'] = d['Qubit'].pop('beta_state')['ms'][0]"
reject S5 rj_a93_ms_and_index 'Qubit.alpha_state must have exactly one of "ms", "index" and "vector" (it has more than one)' \
  "d['Qubit']['alpha_state'] = {'ms': [1, 1], 'index': 0}"
reject S5 rj_a93_ms_too_short 'Qubit.alpha_state.ms must be a list of one m_s per qubit, 2 numbers in the order of Qubit.qubit (it has 1)' \
  "d['Qubit']['alpha_state'] = {'ms': [1]}"
reject S5 rj_a93_ms_too_long 'Qubit.beta_state.ms must be a list of one m_s per qubit, 2 numbers in the order of Qubit.qubit (it has 3)' \
  "d['Qubit']['beta_state'] = {'ms': [0, 0, 0]}"
reject S5 rj_a93_ms_not_number 'Qubit.alpha_state.ms[1] must be a finite JSON number' \
  "d['Qubit']['alpha_state'] = {'ms': [1, '1']}"
reject S1a rj_a93_file_ms_number 'Qubit.alpha_state.ms must be a list of one m_s per qubit, 1 numbers in the order of Qubit.qubit (a list also for one qubit)' \
  "d['Qubit']['alpha_state'] = {'ms': 1}"
# A94: the qubit of a position file is a Qubit.qubit entry, so {"index"} and {"vector"} are taken (Qubit.file refused them);
# pulse.per_qubit turns each qubit between the m_s of {"ms"} and refuses them (user 2026-10-16: A93 took a guessed level)
reject S1a rj_per_qubit_index 'pulse.per_qubit true needs Qubit.alpha_state and Qubit.beta_state as {"ms": [m_1, m_2, ...]}' \
  "d['Qubit']['beta_state'] = {'index': 1}; d['pulse']['per_qubit'] = True"
reject S1a rj_per_qubit_vector 'Qubit.alpha_state is {"vector"}' \
  "d['Qubit']['alpha_state'] = {'vector': [1, 0, 0]}; d['pulse']['per_qubit'] = True"
# the ZFS of the Qubit.file form is a Qubit.intmap self entry (A93; Qubit.zfs in kHz and Qubit.frame until then), naming
# the qubit by Qubit.name (default "q0")
reject S1a rj_a93_old_file_zfs 'input key "Qubit.zfs" is now a Qubit.intmap self entry: "intmap"' \
  "d['Qubit']['zfs'] = [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]"
reject S1a rj_a93_old_file_frame 'input key "Qubit.frame" was removed: the frame of the ZFS is the "frame" of its Qubit.intmap self entry' \
  "d['Qubit']['frame'] = 'qubit'"
reject S1a rj_a93_file_intmap_name "The qubit name is not found in the input file" \
  "d['Qubit']['qubit'][0]['name'] = 'NV'
d['Qubit']['intmap'] = [{'between': ['q0', 'q0'], 'tensor': [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]}]"
# the phase of an index eigenvector (A93): its largest component real and positive, the first of equal ones. One spin-1/2
# qubit in a field along x (the eigenvectors (1, +-1)/sqrt2, a tie) and along y ((1, +-i)/sqrt2), and one NV in a tilted
# field (mixed states with one largest component); psia and psib as the log prints them (two decimals). A short run.
phase_rule(){
    local name="$1" edit="$2" w="$WORK/S5"
    python3 - "$w" "$name" <<PY
import json, sys; from input_layout import flat, nest
w, name = sys.argv[1:3]; d = flat(json.load(open(f"{w}/legacy.json")))
$edit
d.update(nstep=2, outfile=f'./{name}')
json.dump(nest(d), open(f"{w}/{name}.json", "w"))
PY
    if ! (cd "$w" && run_ccex "$BIN_DEV" "$name.json" "proc_$name"); then
        echo "  FAIL  $name : the run failed, see $w/proc_$name"; FAIL=$((FAIL+1)); return; fi
    if python3 - "$w/proc_$name" <<'PY'
import re, sys
text = open(sys.argv[1]).read()
def vector(name):         # psia : [(re,im), (re,im), ...] as the log prints it, the last one printed
    line = [l for l in text.splitlines() if re.match(rf'\s*{name}\s*:', l)][-1]
    nums = [float(x) for x in re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', line.split(':', 1)[1])]
    return [complex(nums[k], nums[k+1]) for k in range(0, len(nums), 2)]
ok = True
for name in ('psia', 'psib'):
    v = vector(name)
    big = max(abs(c) for c in v)
    first = next(c for c in v if abs(c) >= big*(1 - 1e-6))
    good = first.real > 0 and abs(first.imag) <= 1e-12
    print(f"        {name} = {v}: first largest component {first} {'real > 0' if good else 'NOT real > 0'}")
    ok &= good
sys.exit(0 if ok else 1)
PY
    then echo "  PASS  $name (the largest component of psia and psib is real and positive)"; PASS=$((PASS+1))
    else echo "  FAIL  $name : see $w/proc_$name"; FAIL=$((FAIL+1)); fi
}
SPIN_HALF="d['Qubit'] = {'qubit': [{'name': 'e', 'spin': 0.5, 'xyz': [0.0, 0.0, 0.0]}], 'alpha_state': {'index': 0}, 'beta_state': {'index': 1}}"
phase_rule a93_phase_tie_x "$SPIN_HALF; d['bfield'] = [100.0, 0.0, 0.0]"
phase_rule a93_phase_tie_y "$SPIN_HALF; d['bfield'] = [0.0, 100.0, 0.0]"
phase_rule a93_phase_mixed "q = d['Qubit']; q.update(qubit=q['qubit'][:1], intmap=q['intmap'][:1], alpha_state={'index': 0}, beta_state={'index': 2}); d['bfield'] = [300.0, 200.0, 500.0]"
# D28: Qubit.intmap[].tensor shorter than 3x3 was padded with zeros
reject S5 rj_intmap_tensor_short "3x3" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0]]"
reject S5 rj_intmap_tensor_string "finite JSON number" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0,0,0],[0,'a',0],[0,0,1.0]]"
# review 20261008-autoconv-d42 #5: a qubit array shorter than nqubit left the last qubit allocated but never set; the
# number of qubits is the length of Qubit.qubit since A47 (Qubit.nqubit removed), so what is left to refuse is no qubit
reject S5 rj_qubit_array_short "Qubit.qubit must be an array of one entry per qubit" \
  "d['Qubit']['qubit'] = []"

echo
echo "=== S3 : removed key, and the Defect object forms"
reject S3 rj_magnetic_field_axis "magnetic_field_axis" \
  "d['magnetic_field_axis'] = [1.0,1.0,1.0]"
reject S3 rj_defaxis_without_rot "defect_axis_reference" \
  "d['defect_axis_reference'] = 'qubit_axis'"
reject S3 rj_zfs_de_no_axis "axis is missing" \
  "df = d['Defect']['species'][0]
df['frame'] = 'qubit'
df['zfs'] = {'D': 2870.0, 'unit': 'MHz'}"
reject S3 rj_zfs_E_nonzero "zfs.E was removed (only 0 was taken)" \
  "df = d['Defect']['species'][0]
df['frame'] = 'qubit'
df['zfs'] = {'D': 2870.0, 'E': 5.0, 'axis': [1.0,1.0,1.0], 'unit': 'MHz'}"
# A90: with the rotation on every species states its frame (an absent one was "qubit" until then, and an axis or an
# object-form zfs asked for coordinate_frame whether the rotation was on or not); the old name is refused with the new
reject S3 rj_defect_frame_missing "frame.rotation is on, so Defect.species[" \
  "d['coordinate_frame_rotation'] = $R111
df = d['Defect']['species'][0]
df['zfs'] = {'D': 2870.0, 'axis': [1.0,1.0,1.0], 'unit': 'MHz'}"
reject S3 rj_a90_coordinate_frame '"Defect.species[0].coordinate_frame" is now "frame"' \
  "df = d['Defect']['species'][0]
df['coordinate_frame'] = 'bath'"
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
  "d['configuration']['bath_shift'] = [[0.0, 0.0]]"
reject S3 rj_bathadjust_rows_short "bath_shift has 0 rows" \
  "d['configuration']['bath_shift'] = []"
# n_configurations 4 means axes 1-4: an axis file with 5 was read one past the tables (2026-10-09); refused while the file is read
reject S3 rj_avaax_past_navaax "The principal axis is out of range" \
  "lines = open(d['configuration']['defect_axis_file'][0]).read().split('\n')
open(f'{w}/avaax_axis5','w').write('\n'.join([lines[0]] + ['5' if l.strip() else l for l in lines[1:]]))
d['configuration']['defect_axis_file'] = [f'{w}/avaax_axis5']" read

reject S3 rj_bad_unit "unit" \
  "df = d['Defect']['species'][0]
df['frame'] = 'qubit'
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
reject S1a rj_type_nk_element "cluster.selection.max_per_order[1] must be an integer" "d['max_clusters_per_order'] = [0, '1', 0]"
reject S1a rj_type_bathadjust_element "bath_shift[0][2] must be a number" "d['configuration']['bath_shift'] = [[0.0, 0.0, '1.0']]"
reject S1a rj_type_method "method must be a string" "d['method'] = 1"

echo
echo "=== a key nobody reads, close to one that took its default"
# Until 2026-10-02 (D24) such a key was ignored without a word: "nstates" left nstate at 0, and
# the top-level "alphams" of the qubitfile form (which read qalphams) left the NV branch at +1; the same inside a block
# (Qubit.alphms for Qubit.alphams, A47; Qubit.alpha_stat for Qubit.alpha_state since A93, the default taken in its place).
reject S1a rj_typo_nstates "did you mean \"s_range\"" "d['configuration']['s_rnge'] = [1, 1]"
reject S1a rj_typo_qubitfile_alphams "did you mean \"alpha_state\"" "d['Qubit']['alpha_stat'] = d['Qubit'].pop('alpha_state')"

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
    if [ "$how" = full ] && [ ! -e /dev/full ]; then echo "  SKIP  $name (no /dev/full here: macOS)"; return; fi
    rm -rf "${outw:?}"; mkdir -p "$outw"
    python3 - "$CASE/ccein.json" "$outw/ccein.json" "$savemode" "$outw/out" <<'PY'
import json, sys; from input_layout import flat, nest
d = flat(json.load(open(sys.argv[1]))); d["savemode"] = sys.argv[3]; d["outfile"] = sys.argv[4]
json.dump(nest(d), open(sys.argv[2], "w"), indent=1)
PY
    [ "$how" = nodiv ] && mkdir "$outw/out_noDiv"
    [ "$how" = full ] && ln -s /dev/full "$outw/out_noDiv" && ln -s /dev/full "$outw/out_wiDiv"
    # relative paths in the case's input are read from the case's folder
    ( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$outw/ccein.json" 2>&1 < /dev/null ) > "$outw/proc"
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
    local name="$1" want="$2" edit="$3" request="$4"
    local w="$WORK/adaptive_work"; mkdir -p "$w"
    python3 - "$P1IN" "$w" "$name" "$request" <<PY
import json, sys; from input_layout import flat, nest
src, w, name, request = sys.argv[1:5]
d = flat(json.load(open(src)))
$edit
r = json.loads(request); r["response"] = f"{w}/{name}.response"
json.dump(r, open(f"{w}/{name}.request.json", "w"))
d["adaptive"] = {"work_file": f"{w}/{name}.request.json"}
d["outfile"] = f"{w}/{name}"
json.dump(nest(d), open(f"{w}/{name}.json", "w"), indent=1)
PY
    # relative paths in the input are the repository's
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    refused_with "$name" $rc "$want" "$w/proc_$name" || return
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
reject_work rj_factors_contributions "cannot be combined with output.cluster_contributions" "d['cluster_contributions'] = {'enabled': True}" "$F"
reject_work rj_factors_plan_only "probe, factors and except work require dynamics" "d['cluster_only'] = True" "$F"
reject_work rj_state_without_factors "state and step require factors mode" "" "{\"mode\": \"run\", $B, \"state\": 1, \"step\": 3}"
reject_work rj_except_gcce "factors/except work requires method cce" "d['method'] = 'gCCE'" "$E"
reject_work rj_except_order1 "factors/except work requires bundle order 2" "d['order'] = 1" "$E"
reject_work rj_except_probes "probes require probe mode" "" "${E%\}}, \"probes\": [[0]]}"
reject_work rj_except_empty "except mode requires a nonempty except list" "" "{\"mode\": \"except\", $B, \"except\": []}"
reject_work rj_except_twice "except lists a term twice" "" "{\"mode\": \"except\", $B, \"except\": [[0, 3], [0, 3]]}"
reject_work rj_except_not_a_term "an except entry is not a term of this partition's expansion" "" "{\"mode\": \"except\", $B, \"except\": [[0, 1]]}"
reject_work rj_except_plan_only "probe, factors and except work require dynamics" "d['cluster_only'] = True" "$E"
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
reject_case(){    # reject_case <case> <name> <expected message fragment> <python edit>
    local case="$1" name="$2" want="$3" edit="$4"
    local w="$WORK/evolution_update_refusals"; mkdir -p "$w"
    python3 - "$REPO/test/code_verification/check_numerics/cases/$case/ccein.json" "$w/$name.json" "$w/$name" <<PY
import json, sys; from input_layout import flat, nest
d = flat(json.load(open(sys.argv[1])))
d["outfile"] = sys.argv[3]
$edit
json.dump(nest(d), open(sys.argv[2], "w"), indent=1)
PY
    # relative paths in the case's input are read from the case's folder
    ( cd "$REPO/test/code_verification/check_numerics/cases/$case" &&
      "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    refused_with "$name" $rc "$want" "$w/proc_$name" || return
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
reject_case p1_cce2_single rj_update_value "evolution_update must be exact or recursive" "d['evolution_update'] = 'fast'"
reject_case p1_cce2_single rj_update_type "evolution_update must be exact or recursive" "d['evolution_update'] = 1"
reject_case p1_manual_hahn rj_update_gcce_expm "this run has method=gCCE with execution.propagator=expm" \
  "d['evolution_update'] = 'recursive'; d['evolution'] = 'matrix'; d['propagator'] = 'expm'"
reject_case p1_cce2_single rj_update_backend "this run has execution.backend=cuda" \
  "d['evolution_update'] = 'recursive'; d['backend'] = 'cuda'"

echo "=== gpu_options"
# The GPU settings are one block since 2026-10-12 (A26, docs/user/InputFile_Backend.md#tag_gpu_options): the flat keys
# before it are refused with the new form named, and so are an unknown key of the block, a block that is not an object,
# and a packed_max_dim other than 0, 2 or 4. packed_max_dim needs method=CCE and the GPU binary (this check runs
# bin/main.out), both refused before the bath is read. The accepted forms are in ../check_input_forms/gpu_options.sh.
reject_case p1_cce2_single rj_gpu_min_dim "gpu_min_dim is now execution.gpu_options.min_dim" "d['gpu_min_dim'] = 256"
reject_case p1_cce2_single rj_gpu_batch_steps "gpu_batch_steps is now execution.gpu_options.evolution_batch_steps" "d['gpu_batch_steps'] = 64"
reject_case p1_cce2_single rj_gpu_batch_clusters "gpu_batch_clusters is now execution.gpu_options.evolution_batch_size" "d['gpu_batch_clusters'] = 8"
reject_case p1_cce2_single rj_gpu_cluster_batch_mb "gpu_cluster_batch_mb is now execution.gpu_options.cluster_batch_mb" "d['gpu_cluster_batch_mb'] = 64"
reject_case p1_cce2_single rj_gpu_cache_mb "gpu_cache_mb is now execution.gpu_options.cache_mb" "d['gpu_cache_mb'] = 64"
reject_case p1_cce2_single rj_gpu_options_type "gpu_options must be an object" "d['gpu_options'] = 256"
reject_case p1_cce2_single rj_gpu_options_key "unknown key execution.gpu_options.packed_dim" "d['gpu_options'] = {'packed_dim': 4}"
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
reject_case p1_cce2_single rj_save_nodiv "save_nodiv is now output.savemode.nodiv" "d['save_nodiv'] = False"
reject_case p1_cce2_single rj_savemode_key 'savemode has no key "nodvi" (mode, nodiv)' "d['savemode'] = {'mode': 'coherence', 'nodvi': False}"
reject_case p1_cce2_single rj_savemode_type "savemode must be a string" "d['savemode'] = 1"
reject_case p1_cce2_single rj_savemode_mode_type "savemode.mode must be a string" "d['savemode'] = {'mode': 1}"
reject_case p1_cce2_single rj_savemode_nodiv_type "savemode.nodiv must be a boolean" "d['savemode'] = {'mode': 'coherence', 'nodiv': 0}"
reject_case p1_cce2_single rj_savemode_unknown 'savemode "averaged" is not one of' "d['savemode'] = 'averaged'"

echo
echo "=== propagator"
# propagator (A2, 2026-10-12): expm is gCCE only; sparse runs for one sampled state of the cCCE and gCCE evolution=vector
# with one qubit, backend=cpu, evolution_update=exact, without knight, hfmedi and delta_A - everything else is refused at the
# start (the paths that run, and the malformed values, are in ../check_input_forms/propagator.sh).
SP="d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_cce2_single rj_prop_expm_cce "propagator=expm is implemented for method=gCCE only" "d['propagator'] = 'expm'"
reject_case p1_cce2_ensemble rj_prop_sparse_ensemble "this run has the ensemble (no configuration.s_range)" "$SP"
reject_case p1_manual_hahn rj_prop_sparse_matrix "this run has method=gCCE with execution.evolution=matrix" "$SP; d['evolution'] = 'matrix'"
reject_case p1_cce2_single rj_prop_sparse_cuda "this run has execution.backend=cuda" "$SP; d['backend'] = 'cuda'"
reject_case p1_manual_hahn rj_prop_sparse_auto "this run has execution.backend=auto" "$SP; d['backend'] = 'auto'"
reject_case p1_cce2_single rj_prop_sparse_recursive "this run has execution.evolution_update=recursive" "$SP; d['evolution_update'] = 'recursive'"
reject_case p1_manual_hahn rj_prop_sparse_recursive_gcce "this run has execution.evolution_update=recursive" "$SP; d['evolution_update'] = 'recursive'"
reject_case natab_knight_shift rj_prop_sparse_knight "this run has the knight term" "$SP; d['configuration']['nstate'] = 1"
reject_case hbn_hfmedi rj_prop_sparse_hfmedi "this run has the hfmedi term" "$SP; d['configuration']['nstate'] = 1"
reject_case natab_knight_shift rj_prop_sparse_delta_A "this run has the delta_A term" \
  "$SP; d['configuration']['nstate'] = 1; d['hamiltonian'] = {'cCCE': {'delta_A': True}}"
# a second NV 20 A away added to the one-qubit Qubit block (the check_numerics cases hold one qubit each; the old
# natab_gcce2_two_qubits held one too, through its qubitfile, and is natab_gcce2_bathadjust since 2026-10-13)
TWO="q = d['Qubit']; q2 = dict(q['qubit'][0], name='qubit2', xyz=[x + 20.0 for x in q['qubit'][0]['xyz']]); q['qubit'].append(q2); q['intmap'].append(dict(q['intmap'][0], between=['qubit2', 'qubit2']))
for k in ('alpha_state', 'beta_state'):
    if 'ms' in q.get(k, {}): q[k]['ms'] = q[k]['ms']*2"
reject_case p1_gcce2_zfs rj_prop_sparse_two_qubits "this run has nqubit=2" "$SP; $TWO"

echo
echo "=== column_evolution (A62)"
# the block of the column-block one-state cCCE evolution: an object, enabled a JSON boolean, min_dim and columns positive
# integers, no other key (docs/user/InputFile_Backend.md#tag_column_evolution)
reject_case p1_pcce2_sk2 rj_colevo_type "column_evolution must be an object" "d['column_evolution'] = True"
reject_case p1_pcce2_sk2 rj_colevo_enabled "column_evolution.enabled must be a boolean" "d['column_evolution'] = {'enabled': 1}"
reject_case p1_pcce2_sk2 rj_colevo_min_dim "column_evolution.min_dim must be a positive 32-bit integer" "d['column_evolution'] = {'min_dim': 0}"
reject_case p1_pcce2_sk2 rj_colevo_columns "column_evolution.columns must be a positive 32-bit integer" "d['column_evolution'] = {'columns': 1.5}"
reject_case p1_pcce2_sk2 rj_colevo_key "column_evolution.colums" "d['column_evolution'] = {'colums': 64}"

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
reject_case $HB rj_tensor_defect_vs_spin "spin S = 1 defect (provenance line) but the qubit spin is 1.5" "d['hamiltonian']['common']['hyperfine']['file'] = '$V3'; d['Qubit']['qubit'][0]['spin'] = 1.5; d['Qubit'].update(alpha_state={'ms': [1.5]}, beta_state={'ms': [0.5]})"
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
reject_case natab_cce2_hahn rj_old_qubitfile "${OLDK}qubitfile\" is now an entry of Qubit.qubit" "d['qubitfile'] = d['Qubit'].pop('qubit')[0]['file']"
reject_case natab_cce2_hahn rj_old_rbathcut "${OLDK}rbathcut\" is now the inner radius of \"rbath\": [inner, outer]" "d['rbathcut'] = 1.0"
reject_case natab_cce2_hahn rj_old_addsubclus "${OLDK}addsubclus\" is now \"cluster.selection.add_subclusters\"" "d['addsubclus'] = 0"
reject_case natab_cce2_hahn rj_old_nk "${OLDK}nk\" is now \"cluster.selection.max_per_order\"" "d['nk'] = [0, 0, 0]"
reject_case natab_cce2_hahn rj_old_statefile "${OLDK}statefile\" is now \"configuration.bath_state_file\"" "d['statefile'] = 'state'"
reject_case natab_cce2_hahn rj_old_vector_mean_field "${OLDK}vector_mean_field\" is now \"hamiltonian.common.bath_bath_mean_field\"" "d['vector_mean_field'] = True"
reject_case natab_cce2_hahn rj_old_knight "${OLDK}knight\" is now \"hamiltonian.cCCE.knight\"" "d['knight'] = True"
reject_case hbn_cce2_hf3_qd0 rj_old_hf_readmode "${OLDK}hf_readmode\" is now \"hamiltonian.common.hyperfine\"" "d['hf_readmode'] = 3"
reject_case hbn_cce2_hf3_qd0 rj_old_DefectTotSpin "${OLDK}DefectTotSpin\" was removed: the qubit spin is used" "d['DefectTotSpin'] = 1.0"
reject_case hbn_cce2_hf3_qd0 rj_old_CorrTotSpin "${OLDK}CorrTotSpin\" was removed" "d['CorrTotSpin'] = 0.0"
reject_case p1_pcce2_sk2 rj_old_sK "${OLDK}sK\" is now \"cluster.pcce.partition_size\"" "d['sK'] = d['pcce'].pop('partition_size')"
reject_case p1_pcce2_sk2 rj_old_method_pcce "method \"pCCE\" is now \"method\": \"cCCE\" with \"cluster\": {\"mode\": \"pcce\"}" "d['method'] = 'pCCE'; d.pop('cluster_mode')"
reject_case natab_cce2_hahn rj_old_method_cce "method \"CCE\" is now \"cCCE\"" "d['method'] = 'CCE'"
reject_case natab_cce2_hahn rj_old_cluster_mode "cluster.mode \"adaptive\" is now \"bundle\"" "d['cluster_mode'] = 'adaptive'"
reject_case natab_cce2_hahn rj_old_bundle_max_dim "${OLDK}bundle_max_dim\" is now \"cluster.bundle.max_dim\"" "d['bundle_max_dim'] = 8"
reject_case natab_cce2_hahn rj_old_gcce_refinement "\"execution.gcce_refinement\": true|false is now \"execution.gcce_refinement\": {\"enabled\": true|false}" "d['gcce_refinement'] = False"
reject_case hbn_gcce2_oh rj_old_nqubit "${OLDK}Qubit.nqubit\" was removed" "d['Qubit']['nqubit'] = 1"
reject_case hbn_gcce2_oh rj_old_psia "${OLDK}Qubit.psia\" is now \"Qubit.alpha_state\": {\"vector\": [...]}" "d['Qubit']['psia'] = d['Qubit'].pop('alpha_state')['vector']"
reject_case hbn_gcce2_oh rj_old_overhaus "${OLDK}Qubit.overhaus\" is now \"hamiltonian.common.qubit_bath_mean_field\"" "d['Qubit']['overhaus'] = True"
reject_case p1_cce2_single rj_old_defect_list "\"Defect\": [...] is now \"Defect\": {\"species\": [...]}" "d['Defect'] = d['Defect']['species']"
reject_case p1_cce2_single rj_old_dfname "${OLDK}Defect.species[0].dfname\" is now \"name\"" "d['Defect']['species'][0]['dfname'] = d['Defect']['species'][0].pop('name')"
reject_case p1_cce2_single rj_old_exstatefile "${OLDK}exstatefile\" is now \"configuration.defect_state_file\"" "d['exstatefile'] = 'x'"
# A65 (2026-10-15): the keys of how it is computed, of the clusters, of the output and of the frame are in the execution,
# cluster, output and frame blocks; at the old places they are refused with the new one (__keep_old__: input_layout.py
# leaves them there)
reject_case natab_cce2_hahn rj_a65_backend "${OLDK}backend\" is now \"execution.backend\"" "d['__keep_old__'] = ['backend']; d['backend'] = 'cpu'"
reject_case natab_cce2_hahn rj_a65_cluster_mode "${OLDK}cluster_mode\" is now \"cluster.mode\"" "d['__keep_old__'] = ['cluster_mode']; d['cluster_mode'] = 'standard'"
reject_case natab_cce2_hahn rj_a65_outfile "${OLDK}outfile\" is now \"output.file\"" "d['__keep_old__'] = ['outfile']"
reject_case natab_cce2_hahn rj_a65_savemode "${OLDK}savemode\" is now \"output.savemode\"" "d['__keep_old__'] = ['savemode']; d['savemode'] = 'coherence'"
reject_case natab_cce2_hahn rj_a65_rotation "${OLDK}coordinate_frame_rotation\" is now \"frame.rotation\"" "d['__keep_old__'] = ['coordinate_frame_rotation']; d['coordinate_frame_rotation'] = {'enabled': True, 'qubit_axis': [0, 0, 1]}"
reject_case natab_cce2_hahn rj_a65_defect_axis_reference "${OLDK}defect_axis_reference\" was removed" "d['defect_axis_reference'] = 'qubit_axis'"
reject_case natab_cce2_hahn rj_a65_max_clusters "${OLDK}max_clusters_per_order\" is now \"cluster.selection.max_per_order\"" "d['__keep_old__'] = ['max_clusters_per_order']; d['max_clusters_per_order'] = [0, 0, 0]"
reject_case natab_cce2_hahn rj_a65_add_subclusters "${OLDK}add_subclusters\" is now \"cluster.selection.add_subclusters\"" "d['__keep_old__'] = ['add_subclusters']; d['add_subclusters'] = True"
reject_case natab_cce2_hahn rj_a65_probe_memory "${OLDK}probe_memory_gb\" is now \"cluster.adaptive.probe_memory_gb\"" "d['__keep_old__'] = ['probe_memory_gb']; d['probe_memory_gb'] = 1"
reject_case natab_cce2_hahn rj_a65_propagator "${OLDK}propagator\" is now \"execution.propagator\"" "d['__keep_old__'] = ['propagator']; d['propagator'] = 'eigen'"
reject_case p1_cce2_single rj_a65_n_nuclear_spins "${OLDK}Defect.species[0].n_nuclear_spins\" is now \"n_onsite_spins\"" "d['__keep_old__'] = ['n_nuclear_spins']"
reject_case p1_cce2_single rj_a65_n_configurations "${OLDK}Defect.species[0].n_configurations\" is now \"n_axis_states\"" "d['__keep_old__'] = ['n_configurations']"
# output.path (A65): an absolute output.file beside it names two places (reject_case gives every case an absolute file)
reject_case natab_cce2_hahn rj_output_path_absolute_file "is an absolute path and output.path" "d['output'] = {'path': 'out'}"
reject_case natab_cce2_hahn rj_output_path_type "output.path must be a nonempty folder name" "d['output'] = {'path': ''}"
# cluster_mode pcce runs with gCCE, meCCE and megCCE too (internal pgCCE, pmeCCE, pmegCCE); the partition cache stays with
# cCCE + pcce, and the branch's method names (adding_MECCE) are refused with the new form
reject_case p1_pcce2_sk2 rj_pg_cache "enabled cache requires method cCCE with cluster.mode pcce" \
  "d['method']='gCCE'; d['pcce']['cache']={'mode':'readwrite','file':'x'}"
reject_case p1_pcce2_sk2 rj_old_pgcce "method \"pgCCE\" is now method \"gCCE\" with \"cluster.mode\": \"pcce\" in the input" "d['method']='pgCCE'"
# gCCE with the cluster shapes that run with cCCE only: not implemented yet (from the input)
reject_case natab_dsj2 rj_gcce_dsj "method gCCE with cluster.mode \"dsj\" is not implemented yet" "d['method'] = 'gCCE'"
reject_case natab_dsj2 rj_gcce_itb "method gCCE with cluster.mode \"itb\" is not implemented yet" "d['method'] = 'gCCE'; d['cluster_mode'] = 'itb'"
reject_case natab_dsj2 rj_gcce_dsjitb "method gCCE with cluster.mode \"dsjitb\" is not implemented yet" "d['method'] = 'gCCE'; d['cluster_mode'] = 'dsjitb'"
# a hamiltonian setting of the other method: hamiltonian.cCCE in a gCCE run (written so, or the method changed beside a
# cCCE input's hamiltonian.cCCE), hamiltonian.gCCE
reject_case p1_manual_hahn rj_hamiltonian_ccce_in_gcce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE, this run is gCCE" \
  "d.setdefault('hamiltonian', {})['cCCE'] = {'knight': True}"
reject_case natab_knight_shift rj_hamiltonian_ccce_cli_gcce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE, this run is gCCE" "d['method'] = 'gCCE'"
reject_case natab_cce2_hahn rj_hamiltonian_gcce_key "hamiltonian.gCCE has no settings yet" "d['hamiltonian'] = {'gCCE': {'knight': True}}"
# the reserved method noise (A49): refused as not implemented yet
reject_case natab_cce2_hahn rj_method_noise "method \"noise\" is not implemented yet" "d['method'] = 'noise'"
# the adaptive block belongs to cluster_mode pcce or bundle; the pcce block to pcce
reject_case natab_cce2_hahn rj_adaptive_standard "the cluster.adaptive block belongs to cluster.mode pcce or bundle (this input has cluster.mode standard)" \
  "d['adaptive'] = {'output': 'x', 'atol': 1e-8, 'rtol': 1e-6}"
reject_case natab_cce2_hahn rj_pcce_block_standard "the cluster.pcce block belongs to cluster.mode pcce" "d['pcce'] = {'partition_size': 2}"
# delta_A (A83; what it computes is in ../check_analytic/) needs a spin-1 qubit with a zero-field splitting and is
# not taken with hfmedi, which already holds it; bath_bath_mean_field (A47) takes off, z or xyz
reject_case natab_knight_shift rj_delta_A_no_zfs "hamiltonian.cCCE.delta_A is the second-order term of a spin-1 qubit with a zero-field splitting" \
  "del d['Qubit']['intmap']; d['hamiltonian'] = {'cCCE': {'delta_A': True}}"
reject_case natab_knight_shift rj_delta_A_hfmedi "hamiltonian.cCCE.delta_A is refused with hfmedi: hfmedi already includes the delta-A (same-nucleus second-order) term" \
  "d['hamiltonian'] = {'cCCE': {'delta_A': True, 'hfmedi': True}}"
# knight divides by the same Delta (A83 follow-up)
reject_case natab_knight_shift rj_knight_no_zfs "hamiltonian.cCCE.knight divides by the zero-field splitting of the qubit" \
  "del d['Qubit']['intmap']; d['hamiltonian'] = {'cCCE': {'knight': True}}"
reject_case natab_cce2_hahn rj_mean_field_value "hamiltonian.common.bath_bath_mean_field must be \"off\", \"z\" or \"xyz\"" \
  "d['hamiltonian'] = {'common': {'bath_bath_mean_field': 'vector'}}"
# one qubit state, one form
reject_case hbn_gcce2_oh rj_state_index_and_vector "Qubit.alpha_state must have exactly one of \"ms\", \"index\" and \"vector\" (it has more than one)" \
  "d['Qubit']['alpha_state']['index'] = 0"
reject_case hbn_gcce2_oh rj_state_neither "Qubit.beta_state must have exactly one of \"ms\", \"index\" and \"vector\" (it has none)" \
  "d['Qubit']['beta_state'] = {}"

echo
echo "=== DEER bath pulses"
# pulse.bath (DEER, 2026-10-13): pulses on the bath spins of one defect species. Only gCCE and megCCE have the bath in
# their matrices (cluster_mode standard or pcce), and only their CPU paths (gCCE with propagator eigen or expm, megCCE)
# expand a pulse beyond the qubit; the block's form is checked as it is read, and the merged qubit+bath sequence before
# the run (each pulse exactly once). The refusal names the method as the input writes it (meCCE + pcce is pmeCCE inside).
DEER="d['pulse']['bath'] = {'count': 1, 'defect': 'P1', 'energy_shift': [84.946496], 'sequence': [['0.5','X',180.0],['0.5','I',0.0]]}"
reject_case p1_cce2_single rj_deer_cce "bath pulses (pulse.bath.count > 0) are implemented for methods gCCE and megCCE" "$DEER"
reject_case p1_cce2_single rj_deer_mecce "implemented for methods gCCE and megCCE (cluster.mode standard or pcce), this run is method meCCE" \
  "$DEER; d['method']='meCCE'"
reject_case p1_pcce2_sk2 rj_deer_pmecce "implemented for methods gCCE and megCCE (cluster.mode standard or pcce), this run is method meCCE" \
  "$DEER; d['method']='meCCE'"
reject_case p1_pcce2_sk2 rj_deer_pcce_ccce "implemented for methods gCCE and megCCE (cluster.mode standard or pcce), this run is method cCCE" \
  "$DEER"
reject_case p1_manual_hahn rj_deer_cuda "bath pulses (pulse.bath.count > 0) run on the CPU only" "$DEER; d['backend']='cuda'"
reject_case p1_manual_hahn rj_deer_detuning "bath pulses need pulse.detuning_factor = 1" "$DEER; d['pulse']['detuning_factor']=0.5"
reject_case p1_manual_hahn rj_deer_pulseiter "bath pulses cannot be combined with pulse.per_qubit = true" "$DEER; d['pulse']['per_qubit']=True"
reject_case p1_manual_hahn rj_deer_sparse "bath pulses need execution.propagator eigen or expm" "$DEER; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_manual_hahn rj_deer_meg_sparse "bath pulses with method megCCE need execution.propagator eigen (execution.propagator sparse" \
  "$DEER; d['method']='megCCE'; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_manual_hahn rj_deer_warp "bath pulses cannot use execution.gpu_options.warp_max_dim > 0" "$DEER; d['gpu_options'] = {'warp_max_dim': 27}"
reject_case p1_manual_hahn rj_deer_bundle "bath pulses need cluster.mode standard or pcce (this run has cluster.mode bundle)" \
  "$DEER; d['cluster_mode']='bundle'"
reject_case p1_manual_hahn rj_deer_negative "pulse.bath.count must be >= 0" "d['pulse']['bath'] = {'count': -1}"
reject_case p1_manual_hahn rj_deer_rows "pulse.bath.sequence must have pulse.bath.count+1 rows" "$DEER; d['pulse']['bath']['sequence']=[['1','I',0.0]]"
reject_case p1_manual_hahn rj_deer_sum "the fractions of pulse.bath.sequence sum to" "$DEER; d['pulse']['bath']['sequence']=[['0.5','X',180.0],['0.4','I',0.0]]"
reject_case p1_manual_hahn rj_deer_defect "pulse.bath.defect is required when pulse.bath.count > 0" "$DEER; del d['pulse']['bath']['defect']"
reject_case p1_manual_hahn rj_deer_fraction "pulse.bath.sequence[0] fraction \"1/x\" must be a/b" "$DEER; d['pulse']['bath']['sequence'][0][0]='1/x'"
reject_case p1_manual_hahn rj_deer_tones "pulse.bath.energy_shift needs at least one finite number" "$DEER; d['pulse']['bath']['energy_shift']=[]"
reject_case p1_manual_hahn rj_deer_unknown "pulse.bath has no key \"tone\"" "$DEER; d['pulse']['bath']['tone']=1"
# the pulse of the last row would come after the last interval: never applied, so a non-identity one is refused
reject_case p1_manual_hahn rj_deer_last_row "the last row of pulse.bath.sequence must be identity" "$DEER; d['pulse']['bath']['sequence']=[['0.5','X',180.0],['0.5','X',180.0]]"
# two pulses of one channel within 1e-9 would share one edge of the merged sequence and one of them be lost; the chained
# case: two bath pulses 1.8e-9 apart, each within 1e-9 of the qubit pi at 0.5
reject_case p1_manual_hahn rj_deer_two_qubit_pulses "pulse.sequence has two qubit pulses at the same time, bath pulses need them apart" \
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
# the jump lists in variables: macOS's bash 3.2 brace-expands a list written inside "$( ... )" ({a,b} -> a b)
J_E_MINUS="[{'bath_name':'e','operator':'-','rate':1e-4}]"
J_E_Q="[{'bath_name':'e','operator':'q','rate':1e-4}]"
J_E_PM="[{'bath_name':'e','operator':'+-','rate':1e-4}]"
J_E_NEG="[{'bath_name':'e','operator':'z','rate':-1}]"
J_NONAME="[{'operator':'z','rate':1e-4}]"
J_QUBIT_Z="[{'bath_name':'qubit','operator':'z','rate':1e-4}]"
J_E_P1X="[{'bath_name':'e','operator':'-','rate':1e-4},{'bath_name':'P1x','operator':'+','rate':1e-4}]"
J_P1X_Z="[{'bath_name':'P1x','operator':'z','rate':1e-4}]"
J_P1_Z="[{'bath_name':'P1','operator':'z','rate':1e-4}]"
JUMP="$(MEJ meCCE "$J_E_MINUS")"
reject_case p1_manual_hahn rj_me_section_gcce "hamiltonian.meCCE (jump_operators) belongs to method meCCE, this run is gCCE" "$JUMP"
reject_case p1_cce2_single rj_me_section_other "hamiltonian.megCCE (jump_operators) belongs to method megCCE, this run is meCCE" \
  "d['method']='meCCE'; $(MEJ megCCE "$J_E_MINUS")"
# the run's method as the input writes it: the default cCCE, not the internal name CCE (an old name the reader refuses)
reject_case p1_cce2_single rj_me_section_ccce "hamiltonian.meCCE (jump_operators) belongs to method meCCE, this run is cCCE" \
  "d.pop('method'); $JUMP"
reject_case p1_cce2_single rj_me_section_key "hamiltonian.meCCE has no key \"knight\" (jump_operators)" \
  "d['method']='meCCE'; d['hamiltonian'] = {'meCCE': {'knight': True}}"
reject_case p1_cce2_single rj_old_jump_operators "${OLDK}jump_operators\" is now \"hamiltonian.meCCE.jump_operators\"" \
  "d['jump_operators'] = [{'bath_name':'e','operator':'-','rate':1e-4}]"
reject_case p1_cce2_single rj_me_op "jump_operators[0].operator must be one of - + z x y" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_E_Q")"
reject_case p1_cce2_single rj_me_old_notation "the old names +- and -+ are now - and +" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_E_PM")"
reject_case p1_cce2_single rj_me_rate "jump_operators[0].rate must be a finite number >= 0" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_E_NEG")"
reject_case p1_cce2_single rj_me_name "jump_operators[0] needs bath_name, operator and rate" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_NONAME")"
reject_case natab_manual_hahn rj_me_manual "method meCCE takes its pulses from pulse.count" \
  "d['method']='meCCE'; d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',180],['1/2','I',0]]}"
reject_case p1_cce2_single rj_me_recursive "this run has method=meCCE. Use execution.evolution_update=exact" \
  "d['method']='meCCE'; d['evolution_update']='recursive'"
reject_case p1_cce2_single rj_me_qubit_jump "method meCCE does not take qubit jump operators" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_QUBIT_Z")"
# meCCE with cluster_mode pcce (inside pmeCCE) projects the qubit out as well; the message names the method as written
reject_case p1_pcce2_sk2 rj_pme_qubit_jump "method meCCE does not take qubit jump operators" \
  "d['method']='meCCE'; $(MEJ meCCE "$J_QUBIT_Z")"
# the method as the input writes it in every refusal a gCCE / meCCE / megCCE + pcce run (inside pgCCE, pmeCCE, pmegCCE) can
# meet: the internal names must not reach the messages
reject_case p1_pcce2_sk2 rj_pme_recursive "this run has method=meCCE" "d['method']='meCCE'; d['evolution_update']='recursive'"
reject_case p1_pcce2_sk2 rj_pme_vector "evolution=vector is implemented for method=gCCE only (this run has method=meCCE)" \
  "d['method']='meCCE'; d['evolution']='vector'"
reject_case p1_pcce2_sk2 rj_pmeg_expm "propagator=expm is implemented for method=gCCE only (this run has method=megCCE)" \
  "d['method']='megCCE'; d['propagator']='expm'"
reject_case p1_pcce2_sk2 rj_pme_pulse "method meCCE calculates every pulse as an ideal pi flip" \
  "d['method']='meCCE'; d['pulse']={'count':1,'name':'Manual','sequence':[['1/2','X',90],['1/2','I',0]]}"
reject_case p1_pcce2_sk2 rj_pme_rdip_ladder "rdip_ladder needs method cCCE (cluster.mode standard or pcce) or gCCE (cluster.mode standard), this run is method meCCE" \
  "d['method']='meCCE'; d['rdip_ladder']=[600]"
reject_case p1_pcce2_sk2 rj_pg_packed "execution.gpu_options.packed_max_dim works with method=CCE (this run: gCCE, cluster.mode pcce)" \
  "d['method']='gCCE'; d['gpu_options']={'packed_max_dim': 4}"
reject_case p1_pcce2_sk2 rj_pg_warp "execution.gpu_options.warp_max_dim works with method=CCE or gCCE (this run: gCCE, cluster.mode pcce)" \
  "d['method']='gCCE'; d['gpu_options']={'warp_max_dim': 27}"
for name in rj_pme_qubit_jump rj_pme_recursive rj_pme_vector rj_pmeg_expm rj_pme_pulse rj_pme_rdip_ladder rj_pg_packed rj_pg_warp \
            rj_deer_pmecce; do
    if grep -i "error" "$WORK/evolution_update_refusals/proc_$name" | grep -qiE "p(g|me|meg)cce"; then
        echo "  FAIL  $name : an internal method name in the message"; FAIL=$((FAIL+1))
    else
        echo "  PASS  $name names no internal method"; PASS=$((PASS+1))
    fi
done
reject_case p1_cce2_single rj_me_sparse "method meCCE needs execution.propagator eigen" \
  "d['method']='meCCE'; d['propagator'] = {'method': 'sparse', 'min_dim': 1}"
reject_case p1_cce2_single rj_me_shape "method meCCE with cluster.mode \"dsj\" is not implemented" \
  "d['method']='meCCE'; d['cluster_mode']='dsj'"
reject_case p1_manual_hahn rj_meg_pulseiter "method megCCE needs pulse.per_qubit = false" "d['method']='megCCE'; d['pulse']['per_qubit']=True"
reject_case p1_manual_hahn rj_meg_detuning "method megCCE needs pulse.detuning_factor = 1" "d['method']='megCCE'; d['pulse']['detuning_factor']=0.5"
reject_case p1_manual_hahn rj_meg_density_matrix "method megCCE does not write output.savemode cluster_density_matrix" \
  "d['method']='megCCE'; d['savemode']='cluster_density_matrix'"
reject_case p1_manual_hahn rj_meg_ccce "hamiltonian.cCCE (knight, hfmedi, delta_A) belongs to method cCCE, this run is megCCE" \
  "d['method']='megCCE'; d.setdefault('hamiltonian', {})['cCCE'] = {'knight': True}"
# quantity dm is not implemented for the ME methods
reject_case p1_cce2_single rj_me_dm "method meCCE computes quantity coherence only" "d['method']='meCCE'; d['quantity']='dm'"
reject_case p1_manual_hahn rj_meg_dm_cli "method megCCE computes quantity coherence only" "d['method']='megCCE'; d['quantity']='dm'"
# the pcce block in an ME run names the run's cluster_mode (standard), not the method
reject_case natab_cce2_hahn rj_me_pcce_block "the cluster.pcce block belongs to cluster.mode pcce (this input has cluster.mode standard)" \
  "d['method']='meCCE'; d['pcce'] = {'partition_size': 2}"
# rdip_ladder assembles cCCE, gCCE and pCCE products only; refused for ME before the bath is read (the ladder's own check
# comes after the clusters)
reject_case p1_cce2_single rj_me_rdip_ladder "rdip_ladder needs method cCCE (cluster.mode standard or pcce) or gCCE (cluster.mode standard), this run is method meCCE" \
  "d['method']='meCCE'; d['rdip_ladder']=[600]"
# guarded already by v2's own whitelists: checked here so a later change cannot open them for ME silently
reject_case p1_manual_hahn rj_meg_cuda "execution.backend=cuda supports method cCCE (cluster.mode standard or pcce), or method gCCE with execution.propagator=eigen" \
  "d['method']='megCCE'; d['backend']='cuda'"
reject_case p1_manual_hahn rj_meg_vector "evolution=vector is implemented for method=gCCE only" "d['method']='megCCE'; d['evolution']='vector'"
reject_case p1_manual_hahn rj_meg_expm "propagator=expm is implemented for method=gCCE only (this run has method=megCCE)" \
  "d['method']='megCCE'; d['propagator']='expm'"
reject_case p1_cce2_single rj_me_plan "output.cluster_only supports cCCE, gCCE and cluster.mode pcce plans" "d['method']='meCCE'; d['cluster_only']=True"
# the qubit count (two NV, S5): meCCE takes one qubit, and a qubit jump operator needs one qubit
reject S5 rj_me_two_qubits "method meCCE supports one qubit" "d['method']='meCCE'"
reject S5 rj_me_qubit_jump_two "a qubit jump operator needs one qubit" \
  "d['method']='megCCE'; d['hamiltonian'] = {'megCCE': {'jump_operators': [{'bath_name':'qubit','operator':'z','rate':1e-4}]}}"
# a jump whose bath_name is neither the qubit nor a bath spin of the run (bath file and Defect sub-spins) would run at
# rate 0 without a word: refused once the bath is read, with the names the bath has. meCCE on check_master_equation's
# S=1/2 electron-spin system (spins named e; the first entry matches, the second does not), megCCE on the P1 bath
# (nstep and rbath small, so that a build without the check runs these in seconds)
MEI="$REPO/test/code_verification/check_master_equation/inputs"
reject_case p1_cce2_single rj_me_jump_name 'hamiltonian.meCCE.jump_operators[1].bath_name "P1x" names no bath spin of this run (bath spins, compared ignoring case: e)' \
  "d = {'method': 'meCCE', 'quantity': 'coherence', 'Qubit': {'file': '$MEI/qubit_espin', 'spin': 0.5, 'alpha_state': {'ms': [0.5]}, 'beta_state': {'ms': [-0.5]}}, 'configuration': {'gyrofile': '$MEI/gyro_espin', 'bathfile': ['$MEI/bath_espin'], 'nstate': 0}, 'order': 2, 'bfield': 300, 'rbath': 500, 'rdip': 300, 'deltat': 0.04, 'nstep': 5, 'pulse': {'count': 1}}; $(MEJ meCCE "$J_E_P1X")"
reject_case p1_manual_hahn rj_meg_jump_name 'hamiltonian.megCCE.jump_operators[0].bath_name "P1x" names no bath spin of this run (bath spins, compared ignoring case: P1)' \
  "d['method']='megCCE'; d['nstep']=5; d['rbath']=200; $(MEJ megCCE "$J_P1X_Z")"
# a bath_name that is the qubit's name and a bath spin's: the megCCE kernel takes it as a qubit jump and drops the bath
# spins, so it is refused once the bath is read (the qubit of p1_gcce2_zfs named P1, as its bath spins)
reject_case p1_gcce2_zfs rj_meg_jump_both_names 'hamiltonian.megCCE.jump_operators[0].bath_name "P1" names both the qubit and bath spins of this run' \
  "d['method']='megCCE'; d['nstep']=5; d['Qubit']['qubit'][0]['name']='P1'; d['Qubit']['intmap'][0]['between']=['P1','P1']; $(MEJ megCCE "$J_P1_Z")"
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
echo "=== -f as getopt reads it, keys without regard to case (review 20261013-180049 I01, I03)"
# I01: the input file was found by a search for an argument equal to "-f", so "-fFILE" (and "-vfFILE") skipped the adaptive
# block: a bundle input with one ran as a plain bundle run (curve off by up to 2.84), and adaptive with cluster_mode
# standard was not refused. The file is found by getopt now; the refusal in both attached forms, then the merge-grow
# regression input with an adaptive block run spaced and attached: both refine, to the same final curve.
reject_attached(){    # reject_attached <option word before the file name> <name>: rj_adaptive_standard's input, the file attached
    local w="$WORK/evolution_update_refusals" rc
    ( cd "$REPO/test/code_verification/check_numerics/cases/natab_cce2_hahn" &&
      "$MPIRUN" -n 1 "$BIN_DEV" "$1$w/rj_adaptive_standard.json" > "$w/proc_$2" 2>&1 < /dev/null )
    rc=$?
    if [ $rc -ne 0 ] && grep -qF "the cluster.adaptive block belongs to cluster.mode pcce or bundle (this input has cluster.mode standard)" "$w/proc_$2"
    then echo "  PASS  $2 (exit=$rc)"; PASS=$((PASS+1))
    else echo "  FAIL  $2 : exit=$rc, see $w/proc_$2"; FAIL=$((FAIL+1)); fi
}
reject_attached -f rj_adaptive_standard_attached
reject_attached -vf rj_adaptive_standard_combined
adaptive_attached(){
    local w="$WORK/adaptive_attached" f ok=1; rm -rf "${w:?}"; mkdir -p "$w"
    python3 - "$REPO/test/cluster_verification/check_adaptive/regression/inputs/merge-grow-0.json" "$w" <<'PY'
import json, sys; from input_layout import flat, nest
src, w = sys.argv[1:3]
d = flat(json.load(open(src)))
for form in ("spaced", "attached", "upper"):
    # upper: the adaptive keys in another case (I03), read as the lower-case ones
    keys = ("Workdir", "ATOL", "Rtol") if form == "upper" else ("workdir", "atol", "rtol")   # workdir: output until A94
    d["adaptive"] = dict(zip(keys, (f"{w}/{form}_adaptive", 1e-8, 2e-6)))
    d["outfile"] = f"{w}/{form}"
    json.dump(nest(d), open(f"{w}/{form}.json", "w"), indent=1)
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
    if [ $ok -eq 1 ]; then echo "  PASS  adaptive_attached (-f FILE, -fFILE and Workdir/ATOL/Rtol refine alike)"; PASS=$((PASS+1))
    else echo "  FAIL  adaptive_attached : see $w"; FAIL=$((FAIL+1)); fi
}
adaptive_attached
# accept_same <case> <name> <python edit> <python edit of the reference>: the edited input runs, and its raw result is the
# reference's, byte for byte
accept_same(){
    local case="$1" name="$2" edit="$3" ref="$4"
    local w="$WORK/accepted_forms" c="$REPO/test/code_verification/check_numerics/cases/$case" f ok=1; mkdir -p "$w"
    python3 - "$c/ccein.json" "$w/$name" <<PY
import json, sys; from input_layout import flat, nest
d = flat(json.load(open(sys.argv[1]))); d["outfile"] = sys.argv[2]; $edit
json.dump(nest(d), open(sys.argv[2] + ".json", "w"), indent=1)
d = flat(json.load(open(sys.argv[1]))); d["outfile"] = sys.argv[2] + "_ref"; $ref
json.dump(nest(d), open(sys.argv[2] + "_ref.json", "w"), indent=1)
PY
    rm -f "$w/$name"_*Div "$w/$name"_state* "$w/$name"_ref_*
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/${name}_ref.json" > "$w/proc_${name}_ref" 2>&1 < /dev/null ) || ok=0
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" > "$w/proc_$name" 2>&1 < /dev/null ) || ok=0
    local n=0
    for f in "$w/${name}_ref"_*; do
        n=$((n+1)); cmp -s "$f" "$w/$name${f#$w/${name}_ref}" || ok=0
    done
    [ $n -gt 0 ] || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  $name"; PASS=$((PASS+1))
    else echo "  FAIL  $name : see $w/proc_$name"; FAIL=$((FAIL+1)); fi
}
# (I02, the -m option over the file's method, went with -m in A63: ac_gcce_pcce_cli_ccce, ac_gcce_pcce_cli_ccce_first and
# rj_file_method_under_cli; the method is the input's alone)
# I03: the keys of the savemode and propagator blocks in another case
accept_same natab_cce2_hahn ac_savemode_keys_case "d['savemode'] = {'Mode': 'coherence', 'NoDiv': True}" \
  "d['savemode'] = {'mode': 'coherence', 'nodiv': True}"
accept_same natab_cce2_hahn ac_propagator_keys_case "d['propagator'] = {'Method': 'eigen'}" "d['propagator'] = {'method': 'eigen'}"
# A89: max_per_order entries 1 and 2 cut nothing (only orders 3 and up are cut): the same result as all 0, and one warning
# (A94: the one wording of a key read and not used)
accept_same p1_cce2_single ac_max_per_order_low "d['max_clusters_per_order'] = [0, 3, 5]" "d['max_clusters_per_order'] = [0, 0, 0]"
if [ "$(grep -c 'Warning: input key "cluster.selection.max_per_order" is not used by this run: entries 1, 2 are ignored, only orders 3 and up are cut' \
        "$WORK/accepted_forms/proc_ac_max_per_order_low")" = 1 ] \
   && ! grep -q 'input key "cluster.selection.max_per_order" is not used' "$WORK/accepted_forms/proc_ac_max_per_order_low_ref"; then
    echo "  PASS  ac_max_per_order_low warning (once, and not for all 0)"; PASS=$((PASS+1))
else echo "  FAIL  ac_max_per_order_low warning : see $WORK/accepted_forms/proc_ac_max_per_order_low"; FAIL=$((FAIL+1)); fi

echo "=== A67 : units of a key, the memory check"
# a unit not in the key's row (case counts), an object that is not {"value", "unit"}, a value that is not a number or a list
reject S1a rj_unit_unknown 'deltat.unit "sec" is not known, the unit is one of s, ms, us, ns' "d['deltat'] = {'value': 1, 'unit': 'sec'}"
reject S1a rj_unit_case 'bfield.unit "MT" is not known, the unit is one of G, mT, T' "d['bfield'] = {'value': 50, 'unit': 'MT'}"
reject S1a rj_unit_values_key 'rdip must be a number (in A) or {"value": number or list of numbers, "unit": "..."} with the unit one of A, nm' \
  "d['rdip'] = {'values': [1, 10], 'unit': 'A'}"
reject S1a rj_unit_value_string 'rbath must be a number (in A) or {"value"' "d['rbath'] = {'value': '30', 'unit': 'A'}"
reject S1a rj_unit_extra_key 'deltat must be a number (in ms) or {"value"' "d['deltat'] = {'value': 1, 'unit': 'us', 'note': 'x'}"
# the memory check (no key): a Slurm allocation of 1 MB per node is less than any run needs - refused before the first
# cluster, with the dimension, the two sizes and what to do
SLURM_MEM_PER_NODE=1 reject S1a rj_memory_not_enough "Error: not enough memory: the largest cluster (dimension" "pass" read
if grep -q "Main calculation" "$WORK/S1a/proc_rj_memory_not_enough"; then
    echo "  FAIL  rj_memory_not_enough : refused after the calculation started"; FAIL=$((FAIL+1)); fi

echo "=== A94 : the input key review ((internal review))"
P1D=../../inputs/p1_1ppm
# configuration.s_range [first, last] in place of nstate (an old nstate refused with the new form), its forms
reject_case p1_cce2_single rj_a94_nstate 'input key "configuration.nstate" is now "configuration.s_range"' \
  "d['configuration'].pop('s_range'); d['configuration']['nstate'] = 2; d['__keep_old__'] = ['nstate']"
reject_case p1_cce2_single rj_s_range_number 'configuration.s_range must be [first, last], two integers of 1 or more with first <= last (both included), a list also for one ([n, n])' \
  "d['configuration']['s_range'] = 2"
reject_case p1_cce2_single rj_s_range_reversed 'configuration.s_range must be [first, last]' "d['configuration']['s_range'] = [2, 1]"
reject_case p1_cce2_single rj_s_range_zero 'configuration.s_range must be [first, last]' "d['configuration']['s_range'] = [0, 2]"
# configuration.seed: a number (the first configuration's) or one per configuration
reject_case natab_cce2_hahn rj_seed_list_length 'configuration.seed as a list has one seed per configuration, 1 here' \
  "d['configuration']['seed'] = [19, 20]"
reject_case natab_cce2_hahn rj_seed_list_entry 'configuration.seed[0] must be an integer' "d['configuration']['seed'] = ['19']"
# a file the input names and the run reads that does not open: refused (until A94 a warning, and states or axes drawn
# from the seed, or a default)
reject_case p1_cce2_single rj_missing_state_file "configuration.bath_state_file[0]: cannot read \"$P1D/state_3\"" \
  "d['configuration']['s_range'] = [1, 3]"
reject_case p1_cce2_single rj_missing_state_name_s "configuration.bath_state_file[0]: cannot read \"$P1D/st2ate\"" \
  "d['configuration']['bath_state_file'] = ['$P1D/st{s}ate']; d['configuration']['s_range'] = [2, 2]"
reject_case p1_cce2_single rj_missing_axis_file "configuration.defect_axis_file[0]: cannot read \"$P1D/no_axis\"" \
  "d['configuration']['defect_axis_file'] = ['$P1D/no_axis']"
reject_case p1_cce2_single rj_missing_defect_state "configuration.defect_state_file[0]: cannot read \"$P1D/no_onsite1\"" \
  "d['configuration']['defect_state_file'] = ['$P1D/no_onsite']"
reject_case natab_cce2_hahn rj_missing_bathfile 'configuration.bathfile[0]: cannot read "no_bath"' "d['configuration']['bathfile'] = ['no_bath']"
reject_case natab_cce2_hahn rj_missing_gyrofile 'configuration.gyrofile: cannot read "no_gyro"' "d['configuration']['gyrofile'] = 'no_gyro'"
reject_case hbn_cce2_hf3_qd2 rj_missing_hyperfine_file 'hamiltonian.common.hyperfine.file: cannot read "no_hf"' \
  "d['hamiltonian']['common']['hyperfine']['file'] = 'no_hf'"
reject_case hbn_cce2_hf3_qd2 rj_missing_quadrupole_file 'hamiltonian.common.quadrupole.file: cannot read "no_qd"' \
  "d['hamiltonian']['common']['quadrupole']['file'] = 'no_qd'"
reject_case p1_pcce2_sk2 rj_missing_pcce_cache 'cluster.pcce.cache.file: cannot read "no_cache"' \
  "d['pcce']['cache'] = {'mode': 'read', 'file': 'no_cache'}"
# the ensemble reads no state files: reported, not checked
accept_same natab_cce2_hahn ac_ensemble_state_file "d['configuration']['bath_state_file'] = ['no_state_']" "pass"
grep -qF 'Warning: input key "configuration.bath_state_file" is not used by this run: the ensemble' \
     "$WORK/accepted_forms/proc_ac_ensemble_state_file" \
  && { echo "  PASS  ac_ensemble_state_file warning"; PASS=$((PASS+1)); } \
  || { echo "  FAIL  ac_ensemble_state_file warning : see $WORK/accepted_forms/proc_ac_ensemble_state_file"; FAIL=$((FAIL+1)); }
# the Qubit.file form is a Qubit.qubit entry (G12); a qubit's detuning is a number (G21); a gyro takes a unit (G20)
reject_case natab_cce2_hahn rj_a94_qubit_file 'input key "Qubit.file" is now an entry of Qubit.qubit' \
  "q = d['Qubit']; q['file'] = q.pop('qubit')[0]['file']; d['__keep_old__'] = ['Qubit.file']"
reject_case natab_cce2_hahn rj_a94_qubit_spin 'input key "Qubit.spin" is now the "spin" of the Qubit.qubit entry' "d['Qubit']['spin'] = 1.0"
reject_case natab_cce2_hahn rj_a94_qubit_name 'input key "Qubit.name" is now the "name" of the Qubit.qubit entry' "d['Qubit']['name'] = 'NV'"
reject_case natab_cce2_hahn rj_qubit_detuning_list 'Qubit.qubit[0].detuning is a number (kHz) or {"value": ..., "unit": ...}, not a list' \
  "d['Qubit']['qubit'][0]['detuning'] = [0.0]"
reject_case natab_cce2_hahn rj_qubit_gyro_unit 'Qubit.qubit[0].gyro.unit = "kHz/T" is not supported' \
  "d['Qubit']['qubit'][0]['gyro'] = {'value': -17608.59705, 'unit': 'kHz/T'}"
reject_case natab_gcce2_hahn rj_state_vector_component 'Qubit.alpha_state.vector[0] must be a number or [re, im]' \
  "d['Qubit']['alpha_state'] = {'vector': [[1, 0, 0], 0, 0]}"
# pulse.name: five names, and the count a name fixes (R1); finite pulses are gCCE's (R2)
reject_case natab_cce2_hahn rj_pulse_name_unknown 'pulse.name must be Ramsey, HahnEcho, CPMG, WAHUHA or Manual (not "Hahn")' \
  "d['pulse'] = {'count': 1, 'name': 'Hahn'}"
reject_case natab_cce2_hahn rj_pulse_name_cpmg0 'pulse.name CPMG needs pulse.count 2 or more (this input has 0)' \
  "d['pulse'] = {'count': 0, 'name': 'CPMG'}"
reject_case natab_cce2_hahn rj_pulse_name_ramsey1 'pulse.name Ramsey needs pulse.count 0 (this input has 1)' \
  "d['pulse'] = {'count': 1, 'name': 'Ramsey'}"
reject_case natab_cce2_hahn rj_pulse_name_hahn3 'pulse.name HahnEcho needs pulse.count 1 (this input has 3)' \
  "d['pulse'] = {'count': 3, 'name': 'HahnEcho'}"
reject_case natab_cce2_hahn rj_ccce_duration 'method cCCE needs pulse.detuning_factor = 1 and pulse.duration = 0' \
  "d['pulse'] = {'count': 1, 'duration': 30}"
reject_case natab_cce2_hahn rj_ccce_detuning_factor 'method cCCE needs pulse.detuning_factor = 1 and pulse.duration = 0' \
  "d['pulse'] = {'count': 1, 'detuning_factor': 2.0}"
reject_case p1_cce2_single rj_mecce_duration 'method meCCE needs pulse.duration = 0' \
  "d['method'] = 'meCCE'; d['pulse'] = {'count': 1, 'duration': 30}"
# pulse.bath: spin is the driven spins' own, alphams/betams are levels (G6)
reject_case p1_cce2_single rj_bath_pulse_spin 'pulse.bath.spin was removed' \
  "d['pulse']['bath'] = {'count': 0, 'spin': 0.5}; d['__keep_old__'] = ['pulse.bath']"
reject_case p1_cce2_single rj_bath_pulse_alphams 'pulse.bath.alphams and betams are now pulse.bath.levels' \
  "d['pulse']['bath'] = {'count': 0, 'alphams': 0.5}; d['__keep_old__'] = ['pulse.bath']"
reject_case p1_cce2_single rj_bath_pulse_levels 'pulse.bath.levels must be [m_a, m_b]' "d['pulse']['bath'] = {'count': 0, 'levels': 0.5}"
# cluster.pcce.print_trials is the number printed (G7)
reject_case p1_pcce2_sk2 rj_print_limit 'cluster.pcce.print_limit is now cluster.pcce.print_trials' "d['pcce']['print_limit'] = 5"
reject_case p1_pcce2_sk2 rj_print_trials_bool 'cluster.pcce.print_trials is a number' "d['pcce']['print_trials'] = True"
# bundle connections checked can be written; it takes no rdip (G8)
reject_case natab_cce2_hahn rj_connections_checked_rdip 'cluster.bundle.connections checked takes no rdip' \
  "d['cluster_mode'] = 'bundle'; d['bundle'] = {'max_dim': 8, 'connections': 'checked'}"
# Defect: the axis belongs to zfs (G5), value or values (R3), units with their case (R3)
reject_case p1_cce2_single rj_defect_axis '"Defect.species[0].axis" is now "zfs.axis"' "d['Defect']['species'][0]['axis'] = [1, 1, 1]; d['__keep_old__'] = ['axis']"
reject_case p1_cce2_single rj_defect_zfs_no_axis 'Defect.species[P1].zfs uses the D form but zfs.axis is missing' \
  "d['Defect']['species'][0]['zfs'] = {'D': 0.0, 'unit': 'MHz'}"
reject_case p1_cce2_single rj_defect_value_and_values 'Defect.species[P1].detuning has both "value" and "values"' \
  "r = d['Defect']['species'][0]['detuning']; d['Defect']['species'][0]['detuning'] = {'values': r, 'value': r, 'unit': 'MHz'}"
reject_case p1_cce2_single rj_defect_unit_case 'Defect.species[P1].detuning.unit "mhz" is written "MHz"' \
  "r = d['Defect']['species'][0]['detuning']; d['Defect']['species'][0]['detuning'] = {'values': r, 'unit': 'mhz'}"
# the hamiltonian sections stay named by method (G1 reverted, 2026-10-17): A94's corrections and dissipation are refused
reject_case natab_cce2_hahn rj_hamiltonian_corrections_section 'hamiltonian.corrections is now "hamiltonian.cCCE"' \
  "d['hamiltonian'] = {'corrections': {'knight': False}}"
reject_case p1_cce2_single rj_hamiltonian_dissipation_section 'hamiltonian.dissipation is now "hamiltonian.meCCE" or "hamiltonian.megCCE"' \
  "d['method']='meCCE'; d['hamiltonian'] = {'dissipation': {'jump_operators': []}}"
# cluster.adaptive: workdir (G10), keys with _ (G23)
ADB="d['cluster_mode'] = 'bundle'; d['bundle'] = {'max_dim': 8}; d.pop('rdip')"
reject_case natab_cce2_hahn rj_adaptive_output 'cluster.adaptive.output is now cluster.adaptive.workdir' \
  "$ADB; d['adaptive'] = {'output': 'x', 'atol': 1e-3, 'rtol': 0}"
reject_case natab_cce2_hahn rj_adaptive_dash_key 'cluster.adaptive.max-rounds is spelled "max_rounds" in the input' \
  "$ADB; d['adaptive'] = {'workdir': 'x', 'atol': 1e-3, 'rtol': 0, 'max-rounds': 2}"

# accept_runs <case> <name> <python body>: the body writes inputs with put(tag, d) (outputs <w>/<tag>_...) and names output
# pairs with same(a, b) (two output names, byte for byte the same); every input runs, every pair is compared
accept_runs(){
    local case="$1" name="$2" body="$3"
    local w="$WORK/a94_forms/$name" c="$REPO/test/code_verification/check_numerics/cases/$case" ok=1 a b n=0 f
    rm -rf "${w:?}"; mkdir -p "$w"
    python3 - "$c/ccein.json" "$w" <<PY || ok=0
import copy, json, os, sys; from input_layout import flat, nest
src, w = sys.argv[1:3]
base = flat(json.load(open(src)))
pairs = []
def put(tag, d):
    d = copy.deepcopy(d); d["outfile"] = f"{w}/{tag}"
    json.dump(nest(d), open(f"{w}/{tag}.json", "w"), indent=1)
def same(a, b): pairs.append(f"{a} {b}")
d = copy.deepcopy(base)
$body
open(f"{w}/pairs", "w").write("\n".join(pairs) + "\n")
PY
    for f in "$w"/*.json; do
        ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$f" > "${f%.json}.proc" 2>&1 < /dev/null ) || { ok=0; echo "        $(basename "$f") exited $?"; }
    done
    while read -r a b; do
        [ -n "$a" ] || continue
        for f in wiDiv noDiv; do
            n=$((n+1)); cmp -s "$w/${a}_$f" "$w/${b}_$f" || { ok=0; echo "        ${a}_$f != ${b}_$f"; }
        done
    done < "$w/pairs"
    [ $n -gt 0 ] || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  $name"; PASS=$((PASS+1)); else echo "  FAIL  $name : see $w"; FAIL=$((FAIL+1)); fi
}
# states drawn from the seed (no state files): s_range [2, 3] calculates states 2 and 3 of a run of [1, 3]
accept_runs p1_cce2_single ac_s_range_tail "
d['nstep'] = 40; d['configuration'].pop('bath_state_file')
d['configuration']['s_range'] = [1, 3]; put('all', d)
d['configuration']['s_range'] = [2, 3]; put('tail', d)
same('all_state2', 'tail_state2'); same('all_state3', 'tail_state3')"
# {s} in a state file name: the same files as the number appended
accept_runs p1_cce2_single ac_state_name_s "
d['nstep'] = 40; put('append', d)
d['configuration']['bath_state_file'] = ['$P1D/state_{s}']; put('braces', d)
same('append_state1', 'braces_state1'); same('append_state2', 'braces_state2')"
# configuration.seed with k_range: a number is the first configuration's, k takes seed + (k - first); a list one per k
accept_runs p1_cce2_single ac_seed_of_a_range "
d['nstep'] = 40; c = d['configuration']; c.pop('bath_state_file'); c['s_range'] = [1, 1]
for k in (1, 2): os.symlink(os.path.abspath(c['bathfile'][0]) if os.path.isabs(c['bathfile'][0]) else
                            os.path.normpath(os.path.join(os.path.dirname(src), c['bathfile'][0])), f'{w}/bath_{k}')
c['bathfile'] = [f'{w}/bath_{{k}}']; c['k_range'] = [1, 2]; c['seed'] = 19; put('number', d)
c['seed'] = [19, 20]; put('list', d)
c['seed'] = [5, 20]; put('list5', d)
c.pop('k_range'); c['bathfile'] = [f'{w}/bath_2']; c['seed'] = 20; put('alone', d)
same('number_k2_state1', 'alone_state1'); same('list_k1_state1', 'number_k1_state1'); same('list_k2_state1', 'alone_state1')
same('list5_k2_state1', 'alone_state1')"
# the other forms that mean the same: value for values (R3), a gyro with its unit (G20), a complex vector with im 0
# (G13), connections checked written (G8), a Defect species with detunings only needs no frame with the rotation on (G4)
accept_same p1_cce2_single ac_defect_value \
  "r = d['Defect']['species'][0]['detuning']; d['Defect']['species'][0]['detuning'] = {'value': r, 'unit': 'MHz'}" \
  "r = d['Defect']['species'][0]['detuning']; d['Defect']['species'][0]['detuning'] = {'values': r, 'unit': 'MHz'}"
accept_same natab_cce2_hahn ac_qubit_gyro_unit "d['Qubit']['qubit'][0]['gyro'] = {'value': -17608.597050, 'unit': 'radkHz/G'}" \
  "d['Qubit']['qubit'][0]['gyro'] = -17608.597050"
accept_same natab_gcce2_hahn ac_state_vector_complex "d['Qubit']['alpha_state'] = {'vector': [[1, 0], [0, 0], 0]}" \
  "d['Qubit']['alpha_state'] = {'vector': [1, 0, 0]}"
accept_same natab_cce2_hahn ac_connections_checked \
  "d['cluster_mode'] = 'bundle'; d['bundle'] = {'max_dim': 64, 'connections': 'checked'}; d.pop('rdip'); d['nstep'] = 20; d['rbath'] = 10" \
  "d['cluster_mode'] = 'bundle'; d['bundle'] = {'max_dim': 64}; d.pop('rdip'); d['nstep'] = 20; d['rbath'] = 10"
# A96: connections checked with bundles of more than one site (max_dim 8, merge_cutoff 1) ran its audit through calculate()
# without clusters and stopped with a segmentation fault; it runs to the end and reports the audit
checked_merged(){
    local w="$WORK/checked_merged" ok=1; rm -rf "${w:?}"; mkdir -p "$w"
    python3 - "$REPO/test/cluster_verification/check_bundle/checked_ensemble.json" "$w" <<'PY'
import json, sys
src, w = sys.argv[1:3]
d = json.load(open(src))
d["cluster"]["bundle"] = {"max_dim": 8, "merge_cutoff": 1.0}
d["output"]["file"] = f"{w}/checked"
json.dump(d, open(f"{w}/checked.json", "w"), indent=1)
PY
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/checked.json" > "$w/proc" 2>&1 < /dev/null ) || ok=0
    grep -q '^\[adaptive-check\] status=' "$w/proc" || ok=0
    [ -s "$w/checked_wiDiv" ] || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  checked_merged (connections checked, bundles of several sites)"; PASS=$((PASS+1))
    else echo "  FAIL  checked_merged : see $w/proc"; FAIL=$((FAIL+1)); fi
}
checked_merged
# pulse.bath_drive (continuous bath drive, docs/methods/bath-drive.md): what it does not run is refused, not ignored.
# S1a is an NV in a natural 13C bath (cCCE-2, the ensemble); the run-level checks come after the bath is read ("read").
DRIVE="d['pulse']['bath_drive'] = {'tones': [{'lines': [0.0], 'rabi': 1.0}]}"
reject S1a rj_drive_gcce "pulse.bath_drive runs with method cCCE" "$DRIVE; d['method'] = 'gCCE'" read
reject S1a rj_drive_mean_field "pulse.bath_drive.mean_field must be pfull, projected, z or none" "$DRIVE; d['pulse']['bath_drive']['mean_field'] = 'hartree'" read
reject S1a rj_drive_nuclear "drives electron-like spin-1/2 bath spins only" "$DRIVE" read
reject S1a rj_drive_key 'pulse.bath_drive has no key "rabi"' "$DRIVE; d['pulse']['bath_drive']['rabi'] = 1.0"
reject S1a rj_drive_rabi "pulse.bath_drive.tones[0].rabi must be >= 0" "d['pulse']['bath_drive'] = {'tones': [{'lines': [0.0], 'rabi': -1.0}]}"

accept_same p1_cce2_single ac_defect_frame_not_needed \
  "d['nstep'] = 40; d['coordinate_frame_rotation'] = {'enabled': True, 'qubit_axis': [0.0, 0.0, 1.0]}" "d['nstep'] = 40"

echo
echo "==========================================="
echo " refused combinations : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
