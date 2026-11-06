#!/usr/bin/env bash
#
# Bath-coordinate rotation : test suite
# ---------------------------------------------------------------------------
# Two halves.
#
#   test_rotation_unit.cpp   the geometry alone -- the matrix, the fixed point,
#                            what the transform preserves. Linked against the
#                            objects the main build already produced.
#
#   this script              everything that needs a real run: that a disabled
#                            rotation changes nothing, that a rotated run keeps
#                            the bath order, labels and every distance, that all
#                            bath files get the SAME transform, that the bath
#                            files on disk are untouched, and that the refused
#                            combinations really are refused.
#
# The rotated positions are not compared against a hard-coded table. The checker
# rebuilds R from the convention in the spec (ez = q, ex = q x z, ey = ez x ex,
# as ROWS) and requires CCEX to agree -- so a sign flip anywhere in the C++ shows
# up here rather than being copied into the expected values.
#
# Usage:  bash test/code_verification/check_rotation/run.sh
# Needs:  bin/main.out already built  (bash do_compile.sh wsl-gnu)
# ---------------------------------------------------------------------------

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# One run at a time per folder: a second run would delete the first one's work files (D29).
# A script this run starts (check_refusals runs legacy.sh) inherits the lock and does not take it again.
if [ "${CCEX_RUN_LOCK:-}" != "$HERE" ]; then
    exec 9>"$HERE/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $HERE (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$HERE"; export CCEX_RUN_LOCK
fi
cd "$HERE"
REPO="$(cd "$HERE/../../.." && pwd)"

CCEX_BIN="${CCEX_BIN:-$REPO/bin/main.out}"   # the binary to test; a CMake build keeps it in its own folder

# mpirun, in order of preference: CCEX_MPI, the Intel MPI the loaded module points at,
# whatever is on PATH, and finally the WSL location this suite was written against. The
# absolute path used to be the only choice, which made the suite WSL-only.
MPIRUN="${CCEX_MPI:-}"
if [ -z "$MPIRUN" ] || [ ! -x "$MPIRUN" ]; then
    for __c in "${I_MPI_ROOT:+$I_MPI_ROOT/bin/mpirun}" "$(command -v mpirun 2>/dev/null)" \
               /opt/intel/oneapi/mpi/2021.18/bin/mpirun; do
        [ -n "$__c" ] && [ -x "$__c" ] && { MPIRUN="$__c"; break; }
    done
fi
export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"

WORK="$HERE/_work"

# Same toolchain AND the same flags as the main build. The flags are not cosmetic: obj/*.o is built with
# -O3 -march=native, which on an AVX-512 box raises Eigen's EIGEN_MAX_ALIGN_BYTES. Compiling this file with
# a different -march gives it a different layout for the very same Eigen types, and the first MatrixXcd
# returned across the boundary corrupts the heap ("double free or corruption"). Keep them in step.
#
# The build says how (D26, 2026-10-02): test/unit_build_env.sh brings it up to date and prints its compiler,
# flags, includes, libraries and objects -- the Makefile of do_compile.sh (unit-env) or a CMake build of any
# generator (CCEX_BUILD_DIR=<its build folder>, unit-build.env). Until then this script read them out of the
# Makefile, CMakeCache.txt and flags.make itself and compared every object with its sources.
UNIT_ENV_FILE="$(mktemp)"
bash "$REPO/test/unit_build_env.sh" > "$UNIT_ENV_FILE" || { rm -f "$UNIT_ENV_FILE"; exit 1; }
unit_value(){ sed -n "s/^$1='\(.*\)'\$/\1/p" "$UNIT_ENV_FILE"; }
CXX="$(unit_value CCEX_UNIT_CXX)"
CXXFLAGS="$(unit_value CCEX_UNIT_CXXFLAGS) $(unit_value CCEX_UNIT_INCLUDES)"
UNIT_LIBS="$(unit_value CCEX_UNIT_LIBS)"
UNIT_OBJS="$(unit_value CCEX_UNIT_OBJS)"
for kv in $(unit_value CCEX_UNIT_ENV); do export "$kv"; done
rm -f "$UNIT_ENV_FILE"
[ -n "$CXX" ] && [ -n "$UNIT_OBJS" ] || { echo "ERROR: test/unit_build_env.sh gave no compiler or no objects" >&2; exit 1; }

npass=0
nfail=0
pass(){ echo "  ok    $1"; npass=$((npass+1)); }
fail(){ echo "  FAIL  $1"; nfail=$((nfail+1)); }

[ -x "$CCEX_BIN" ] || { echo "ERROR: no binary at $CCEX_BIN -- build it: cd $REPO && bash do_compile.sh wsl-gnu" >&2; exit 1; }
[ -x "$MPIRUN" ]   || { echo "ERROR: mpirun not found: $MPIRUN" >&2; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK"

# --- input files must come out of this untouched (test 12) -------------------
MD5_BEFORE="$(md5sum inputfiles/* | sort)"

###############################################################################
echo
echo "############################################################"
echo "# unit : geometry"
echo "############################################################"
###############################################################################
UNIT_BIN="$WORK/test_rotation_unit.out"

"$CXX" test_rotation_unit.cpp $UNIT_OBJS -o "$UNIT_BIN" \
    $CXXFLAGS -w \
    -I"$REPO/external/eigen" -I"$REPO/external/uthash/include/" -I"$REPO/include/" \
    $UNIT_LIBS \
    || { echo "ERROR: could not build the unit test" >&2; exit 1; }

"$UNIT_BIN"
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

###############################################################################
echo
echo "############################################################"
echo "# integration : runs"
echo "############################################################"
###############################################################################

# write_ccein <dest> <python dict literal of overrides>
write_ccein(){
    python3 - "$1" "$2" <<'PY'
import ast, json, sys
dest, extra = sys.argv[1], sys.argv[2]
cfg = {
    "method"   : "cCCE",
    "quantity" : "coherence",
    "Qubit"    : {"file": "./inputfiles/qubit_offcenter"},
    "gyrofile" : "./inputfiles/gyro_13C",
    "bathfile" : ["./inputfiles/bath_a"],
    "order"    : 2,
    "bfield"   : 300,
    "rbath"    : 10,
    "rdip"     : 10,
    "deltat"   : 0.1,
    "nstep"    : 2,
    "pulse"    : {"count": 1},
    "seed"     : 1,
}
# literal_eval, not eval: the overrides are plain dict/list/bool literals written
# in this file, and there is no reason for the generator to execute anything.
cfg.update(ast.literal_eval(extra))
from configuration_block import configuration_block   # the A64 block (configuration_block.py)
configuration_block(cfg)
json.dump(cfg, open(dest, "w"), indent=4)
PY
}

# run_case <name> <overrides> -> $WORK/<name>.log , exit code in RC
run_case(){
    local name="$1" extra="$2"
    write_ccein "$WORK/$name.json" "{'outfile': '$WORK/$name', $extra}"
    "$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/$name.json" -v > "$WORK/$name.log" 2>&1
    RC=$?
}

ROT111="'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0]}"

echo
echo "--- [1,2,3] a disabled or identity rotation must change nothing ---"
run_case norot    ""
[ $RC -eq 0 ] || fail "norot run failed (rc=$RC, see $WORK/norot.log)"
run_case disabled "'coordinate_frame_rotation': {'enabled': False, 'qubit_axis': [1.0,1.0,1.0]}"
[ $RC -eq 0 ] || fail "disabled run failed (rc=$RC)"
run_case identity "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [0.0,0.0,1.0]}"
[ $RC -eq 0 ] || fail "identity run failed (rc=$RC)"

for tag in disabled identity; do
    same=1
    for suf in _noDiv _wiDiv; do
        cmp -s "$WORK/norot$suf" "$WORK/$tag$suf" || same=0
    done
    if [ $same -eq 1 ]; then
        pass "[$tag] coherence output is byte-identical to the run with no rotation key"
    else
        fail "[$tag] coherence output DIFFERS from the run with no rotation key"
    fi
done

# The rotated run must actually differ, or the checks above prove nothing.
echo
echo "--- [4,5] the enabled rotation is applied ---"
run_case rot111 "$ROT111"
[ $RC -eq 0 ] || fail "rot111 run failed (rc=$RC, see $WORK/rot111.log)"

if cmp -s "$WORK/norot_noDiv" "$WORK/rot111_noDiv"; then
    fail "[rot111] output is identical to the unrotated run -- the rotation did nothing"
else
    pass "[rot111] output differs from the unrotated run"
fi

grep -q "R \* normalize(qubit_axis)" "$WORK/rot111.log" \
    && pass "[rot111] the rotation diagnostic is printed" \
    || fail "[rot111] no rotation diagnostic in the log"

nrotmsg=$(grep -c "Coordinate frame rotation \.\.\." "$WORK/rot111.log")
[ "$nrotmsg" -eq 1 ] \
    && pass "[rot111] the diagnostic is printed exactly once, not once per bath spin" \
    || fail "[rot111] diagnostic printed $nrotmsg times (expected 1)"

echo
echo "--- [6,7,8,9] geometry and bookkeeping across the rotation ---"
python3 check_geometry.py "$WORK/norot.log" "$WORK/rot111.log" \
    --qubit-axis 1 1 1 --r0 10 20 30
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "--- [9] the sampled bath states are unchanged ---"
# Same seed, same spin order, same number of spins -> the RNG must hand out exactly
# the same states. This is what makes a rotated run comparable to an unrotated one.
run_case norot_st1  "'s_range': [1, 1]"
[ $RC -eq 0 ] || fail "norot_st1 run failed (rc=$RC)"
run_case rot111_st1 "'s_range': [1, 1], $ROT111"
[ $RC -eq 0 ] || fail "rot111_st1 run failed (rc=$RC)"

if diff -q <(grep -E "^\s+bath\[[0-9]+\]\.state" "$WORK/norot_st1.log") \
           <(grep -E "^\s+bath\[[0-9]+\]\.state" "$WORK/rot111_st1.log") > /dev/null; then
    nst=$(grep -cE "^\s+bath\[[0-9]+\]\.state" "$WORK/norot_st1.log")
    pass "bath-spin states are identical across the rotation ($nst lines)"
else
    fail "bath-spin states DIFFER across the rotation"
fi

echo
echo "--- [P1] a defect bath : avaax assignments and per-axis detuning ---"
# A P1 electron bath with n_nuclear_spins = 0, n_configurations = 12, scalar per-axis detuning and a
# FIXED Defect.axis_file (no RNG anywhere). The Defect section is read as written -- it is
# expected to be in the qubit-aligned frame already -- so nothing in it may move.
#
# Two comparisons, because they prove different things:
#   identity rotation vs no rotation  -> byte-identical coherence. The geometry does
#                                        not move, so this isolates the feature itself:
#                                        the whole defect pipeline is untouched by it.
#   [111] rotation   vs no rotation   -> the geometry DOES move, so the coherence must
#                                        differ; what may not differ is paxes[i] and
#                                        the detuning[iax] table they index into.
P1DEF_SPECIES="'name': 'P1', 'n_nuclear_spins': 0, 'n_configurations': 12,
                   'detuning': [[1,'e',1.0],[2,'e',2.0],[3,'e',3.0],[4,'e',4.0],
                                [5,'e',5.0],[6,'e',6.0],[7,'e',7.0],[8,'e',8.0],
                                [9,'e',9.0],[10,'e',10.0],[11,'e',11.0],[12,'e',12.0]]"
P1DEF_NOFRAME="'gyrofile': './inputfiles/gyro_P1', 'bathfile': ['./inputfiles/bath_p1'],
       'Defect': {'axis_file': './inputfiles/avaax_p1', 'species': [{$P1DEF_SPECIES}]}"
# frame "qubit": the species is read as written (A90 asks for it with the rotation on; it was the hidden default)
P1DEF="'gyrofile': './inputfiles/gyro_P1', 'bathfile': ['./inputfiles/bath_p1'],
       'Defect': {'axis_file': './inputfiles/avaax_p1', 'species': [{$P1DEF_SPECIES, 'frame': 'qubit'}]}"

run_case p1_norot    "$P1DEF"
[ $RC -eq 0 ] || fail "p1_norot run failed (rc=$RC, see $WORK/p1_norot.log)"
run_case p1_identity "$P1DEF, 'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [0.0,0.0,1.0]}"
[ $RC -eq 0 ] || fail "p1_identity run failed (rc=$RC)"
run_case p1_rot111   "$P1DEF, $ROT111"
[ $RC -eq 0 ] || fail "p1_rot111 run failed (rc=$RC, see $WORK/p1_rot111.log)"

grep -q "Read AvaaxFile" "$WORK/p1_norot.log" \
    && pass "[P1] the fixed avaaxfile was read (assignments are not random)" \
    || fail "[P1] the avaaxfile was NOT read -- paxes fell back to the RNG"

same=1
for suf in _noDiv _wiDiv; do cmp -s "$WORK/p1_norot$suf" "$WORK/p1_identity$suf" || same=0; done
[ $same -eq 1 ] \
    && pass "[P1] identity rotation leaves the defect run byte-identical" \
    || fail "[P1] identity rotation CHANGED the defect run"

for tag in identity rot111; do
    if diff -q <(grep -E "^\s+paxes\[[0-9]+\]" "$WORK/p1_norot.log") \
               <(grep -E "^\s+paxes\[[0-9]+\]" "$WORK/p1_$tag.log") > /dev/null; then
        npx=$(grep -cE "^\s+paxes\[[0-9]+\]" "$WORK/p1_norot.log")
        pass "[P1/$tag] avaax assignments paxes[i] are unchanged ($npx spins)"
    else
        fail "[P1/$tag] avaax assignments paxes[i] CHANGED"
    fi

    if diff -q <(grep -E "^\s+detuning\[[0-9]+\]" "$WORK/p1_norot.log") \
               <(grep -E "^\s+detuning\[[0-9]+\]" "$WORK/p1_$tag.log") > /dev/null; then
        nde=$(grep -cE "^\s+detuning\[[0-9]+\]" "$WORK/p1_norot.log")
        pass "[P1/$tag] per-axis detuning values are unchanged ($nde axes)"
    else
        fail "[P1/$tag] per-axis detuning values CHANGED"
    fi
done

# The rotated defect run must still move the bath, or the two checks above are vacuous.
if cmp -s "$WORK/p1_norot_noDiv" "$WORK/p1_rot111_noDiv"; then
    fail "[P1/rot111] output is identical to the unrotated run -- the rotation did nothing"
else
    pass "[P1/rot111] the bath really was rotated (coherence differs)"
fi

python3 check_geometry.py "$WORK/p1_norot.log" "$WORK/p1_rot111.log" \
    --qubit-axis 1 1 1 --r0 10 20 30
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# defect_axis_reference only named the rotation's qubit axis in the log; it was removed on 2026-10-15 (A65) and is refused,
# on the defect system as anywhere
run_case p1_rot111_defaxis "$P1DEF, $ROT111, 'defect_axis_reference': 'qubit_axis'"
if [ $RC -ne 0 ] && grep -q 'input key "defect_axis_reference" was removed' "$WORK/p1_rot111_defaxis.log"; then
    pass "[P1] defect_axis_reference is refused as removed (A65)"
else
    fail "[P1] defect_axis_reference was not refused as removed (rc=$RC)"
fi

echo
echo "--- [NV ZFS] shared D/E and tensor forms, units, frames, and legacy compatibility ---"
# n_configurations is three frozen configurations here, selected deterministically as 1,2,3,... .
# It is not three NV orientations: the one physical NV_0 axis belongs to the Defect.
NVBASE="'gyrofile': './inputfiles/gyro_NV', 'bathfile': ['./inputfiles/bath_nv']"
NVAX="'axis_file': './inputfiles/avaax_nv'"
NVDET="'detuning': [[1,'e',-2.14],[2,'e',0.0],[3,'e',2.14]]"
NVDET_MHZ_OBJ="'detuning': {'values': [[1,'e',-2.14],[2,'e',0.0],[3,'e',2.14]],
                            'unit': 'MHz'}"
NVDET_GHZ_OBJ="'detuning': {'values': [[1,'e',-0.00214],[2,'e',0.0],[3,'e',0.00214]],
                            'unit': 'GHz'}"
NVTENSOR="[[0.0,956.6666666666666,956.6666666666666],
           [956.6666666666666,0.0,956.6666666666666],
           [956.6666666666666,956.6666666666666,0.0]]"
NVTENSOR_FLAT="[0.0,956.6666666666666,956.6666666666666,
                956.6666666666666,0.0,956.6666666666666,
                956.6666666666666,956.6666666666666,0.0]"

NV_DE_BATH="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[1.0,1.0,1.0], 'frame':'file',
    'zfs':{'D':2870.0,'unit':'MHz'}, $NVDET}]}"
NV_TENSOR_BATH="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[1.0,1.0,1.0], 'frame':'file',
    'zfs':{'tensor':$NVTENSOR,'unit':'MHz'}, $NVDET}]}"
NV_LEGACY_BATH="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'frame':'file',
    'zfs':[[1,'e',$NVTENSOR_FLAT],[2,'e',$NVTENSOR_FLAT],[3,'e',$NVTENSOR_FLAT]], $NVDET}]}"
NV_DE_GHZ="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[1.0,1.0,1.0], 'frame':'file',
    'zfs':{'D':2.87,'unit':'GHz'}, $NVDET}]}"
NV_DE_QUBIT="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[0.0,0.0,1.0], 'frame':'qubit',
    'zfs':{'D':2870.0,'unit':'MHz'}, $NVDET}]}"
NV_DE_DET_MHZ="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[1.0,1.0,1.0], 'frame':'file',
    'zfs':{'D':2870.0,'unit':'MHz'}, $NVDET_MHZ_OBJ}]}"
NV_DE_DET_GHZ="$NVBASE, 'Defect': {$NVAX, 'species': [{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
    'axis':[1.0,1.0,1.0], 'frame':'file',
    'zfs':{'D':2870.0,'unit':'MHz'}, $NVDET_GHZ_OBJ}]}"

run_case nv_de_bath     "$NV_DE_BATH, $ROT111"
[ $RC -eq 0 ] || fail "[NV ZFS] D/E bath-frame run failed (rc=$RC, see $WORK/nv_de_bath.log)"
run_case nv_tensor_bath "$NV_TENSOR_BATH, $ROT111"
[ $RC -eq 0 ] || fail "[NV ZFS] shared tensor bath-frame run failed (rc=$RC)"
run_case nv_legacy_bath "$NV_LEGACY_BATH, $ROT111"
[ $RC -eq 0 ] || fail "[NV ZFS] legacy tensor bath-frame run failed (rc=$RC)"
run_case nv_de_ghz      "$NV_DE_GHZ, $ROT111"
[ $RC -eq 0 ] || fail "[NV ZFS] GHz D/E run failed (rc=$RC)"
run_case nv_de_qubit    "$NV_DE_QUBIT, $ROT111"
[ $RC -eq 0 ] || fail "[NV ZFS] qubit-frame D/E run failed (rc=$RC)"
run_case nv_det_mhz     "$NV_DE_DET_MHZ, $ROT111"
[ $RC -eq 0 ] || fail "[NV detuning] unit-aware MHz run failed (rc=$RC)"
run_case nv_det_ghz     "$NV_DE_DET_GHZ, $ROT111"
[ $RC -eq 0 ] || fail "[NV detuning] unit-aware GHz run failed (rc=$RC)"

# compare_numeric_outputs <left tag> <right tag> <description>
compare_numeric_outputs(){
    local left="$1" right="$2" description="$3"
    python3 - "$WORK/${left}_noDiv" "$WORK/${right}_noDiv" \
              "$WORK/${left}_wiDiv" "$WORK/${right}_wiDiv" <<'PY'
import math, re, sys
number = re.compile(r'(?<![A-Za-z_])[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?')
for a, b in zip(sys.argv[1::2], sys.argv[2::2]):
    av = [float(x) for x in number.findall(open(a).read())]
    bv = [float(x) for x in number.findall(open(b).read())]
    if len(av) != len(bv) or any(not math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-12)
                                 for x,y in zip(av,bv)):
        raise SystemExit(1)
PY
    if [ $? -eq 0 ]; then
        pass "$description"
    else
        fail "$description -- numerical outputs differ"
    fi
}

compare_numeric_outputs nv_de_bath nv_tensor_bath "[NV ZFS] D/E and shared tensor forms are numerically equivalent"
compare_numeric_outputs nv_de_bath nv_legacy_bath "[NV ZFS] shared and legacy per-index tensors are numerically equivalent"
compare_numeric_outputs nv_de_bath nv_de_ghz "[NV ZFS] 2870 MHz and 2.87 GHz are numerically equivalent"
compare_numeric_outputs nv_de_bath nv_de_qubit "[NV ZFS] bath-frame [111] and qubit-frame [001] describe the same NV"
compare_numeric_outputs nv_de_bath nv_det_mhz "[NV ZFS] legacy and unit-aware MHz detuning are numerically equivalent"
compare_numeric_outputs nv_de_bath nv_det_ghz "[NV ZFS] -2.14 MHz and -0.00214 GHz detuning are numerically equivalent"

grep -q "Defect templates.*1 transformed from file to qubit frame" "$WORK/nv_de_bath.log" \
    && pass "[NV ZFS] bath-frame Defect template is transformed" \
    || fail "[NV ZFS] no transformed Defect template in the diagnostic"
grep -q "zfs\[shared\]" "$WORK/nv_de_bath.log" \
    && pass "[NV ZFS] shared ZFS is reported once, not as three independent inputs" \
    || fail "[NV ZFS] shared ZFS report is missing"

echo
echo "--- [Defect units] hyperfine, efg, gyros, and quadrupole_moments object forms ---"
UNITBASE="'gyrofile':'./inputfiles/gyro_unit_defect',
          'bathfile':['./inputfiles/bath_unit_defect'], 'order':1"
UNITAX="'axis_file':'./inputfiles/avaax_unit_defect'"
UNITPREFIX="$UNITBASE, 'Defect':{$UNITAX, 'species':[{'name':'DUT','n_nuclear_spins':1,
    'types':['14N'],'spins':[1.0],'n_configurations':1,'frame':'qubit',
    'nuclear_positions':[[1,'14N',[0.1,0.0,0.0]]],"
HYPF_MHZ="'hyperfine':{'values':[[1,'14N',[[81.312,0.0,0.0],[0.0,81.312,0.0],[0.0,0.0,114.0264]]]],'unit':'MHz'}"
HYPF_GHZ="'hyperfine':{'values':[[1,'14N',[[0.081312,0.0,0.0],[0.0,0.081312,0.0],[0.0,0.0,0.1140264]]]],'unit':'GHz'}"
EFG_AU="'efg':{'values':[[1,'14N',[[0.521910,0.0,0.0],[0.0,0.522771,-0.000086],[0.0,-0.000086,-1.044681]]]],'unit':'Hartree/Bohr^2'}"
EFG_VA2="'efg':{'values':[[1,'14N',[[50.715887885787069,0.0,0.0],[0.0,50.79955437899406,-0.0083569319579576713],[0.0,-0.0083569319579576713,-101.51544226478113]]]],'unit':'V/angstrom^2'}"
GYRO_INTERNAL="'gyros':{'values':[1.9337792],'unit':'radkHz/G'}"
GYRO_MHZT="'gyros':{'values':[3.0777051852829089],'unit':'MHz/T'}"
EQ_INTERNAL="'quadrupole_moments':{'values':[2.044],'unit':'1e-30 m^2'}"
EQ_MBARN="'quadrupole_moments':{'values':[20.44],'unit':'mbarn'}"

run_case unit_legacy "$UNITPREFIX
    'gyros':[1.9337792], 'quadrupole_moments':[2.044],
    'hyperfine':[[1,'14N',[81.312,0.0,0.0,0.0,81.312,0.0,0.0,0.0,114.0264]]],
    'efg':[[1,'14N',[0.521910,0.0,0.0,0.0,0.522771,-0.000086,0.0,-0.000086,-1.044681]]]}]}"
[ $RC -eq 0 ] || fail "[Defect units] legacy run failed (rc=$RC)"
run_case unit_objects "$UNITPREFIX $GYRO_INTERNAL, $EQ_INTERNAL, $HYPF_MHZ, $EFG_AU}]}"
[ $RC -eq 0 ] || fail "[Defect units] canonical object run failed (rc=$RC)"
run_case unit_hypf_ghz "$UNITPREFIX $GYRO_INTERNAL, $EQ_INTERNAL, $HYPF_GHZ, $EFG_AU}]}"
[ $RC -eq 0 ] || fail "[Defect units] hypf GHz run failed (rc=$RC)"
run_case unit_gyro_mhzt "$UNITPREFIX $GYRO_MHZT, $EQ_INTERNAL, $HYPF_MHZ, $EFG_AU}]}"
[ $RC -eq 0 ] || fail "[Defect units] gyro MHz/T run failed (rc=$RC)"
run_case unit_eq_mbarn "$UNITPREFIX $GYRO_INTERNAL, $EQ_MBARN, $HYPF_MHZ, $EFG_AU}]}"
[ $RC -eq 0 ] || fail "[Defect units] eqs mbarn run failed (rc=$RC)"
run_case unit_efg_va2 "$UNITPREFIX $GYRO_INTERNAL, $EQ_INTERNAL, $HYPF_MHZ, $EFG_VA2}]}"
[ $RC -eq 0 ] || fail "[Defect units] EFG V/angstrom^2 run failed (rc=$RC)"

compare_numeric_outputs unit_legacy unit_objects "[Defect units] canonical object forms preserve legacy output"
compare_numeric_outputs unit_objects unit_hypf_ghz "[Defect units] 81.312 MHz and 0.081312 GHz hypf are equivalent"
compare_numeric_outputs unit_objects unit_gyro_mhzt "[Defect units] radkHz/G and MHz/T gyros are equivalent"
compare_numeric_outputs unit_objects unit_eq_mbarn "[Defect units] 2.044 x 1e-30 m^2 and 20.44 mbarn are equivalent"
compare_numeric_outputs unit_objects unit_efg_va2 "[Defect units] atomic-unit and V/angstrom^2 EFG are equivalent"

echo
echo "--- [11] every bath file gets the same transform ---"
run_case norot_two  "'bathfile': ['./inputfiles/bath_a','./inputfiles/bath_b']"
[ $RC -eq 0 ] || fail "norot_two run failed (rc=$RC)"
run_case rot111_two "'bathfile': ['./inputfiles/bath_a','./inputfiles/bath_b'], $ROT111"
[ $RC -eq 0 ] || fail "rot111_two run failed (rc=$RC)"

python3 check_geometry.py "$WORK/norot_two.log" "$WORK/rot111_two.log" \
    --qubit-axis 1 1 1 --r0 10 20 30 --expect-nspin 9
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "--- [13,14,15] external HF / QD tensors move with the frame ---"
# Real DFT tensors, from the h-BN V_B set committed under test/code_verification/.
# They are non-diagonal, which is the whole point: a wrong convention shows up in the
# off-diagonal elements and nowhere else.
HEX="$REPO/test/code_verification/bath_data/hexagonal_data"

# hex_case <name> <hf_readmode> <qd_readmode> <extra overrides>
hex_case(){
    local name="$1" hf="$2" qd="$3" extra="$4"
    python3 - "$WORK/$name.json" "$HEX" "$WORK/$name" "$hf" "$qd" "$extra" <<'PYEOF'
import ast, json, sys
dest, hexdir, outfile, hf, qd, extra = sys.argv[1:7]
cfg = {
    "method": "cCCE", "quantity": "coherence",
    "Qubit": {"file": hexdir + "/defect", "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
    "gyrofile" : hexdir + "/../ensemble_data/Gyro_h_BN",
    "bathfile" : [hexdir + "/bath_1"],
    "order": 1, "bfield": 30000, "rbath": 8, "rdip": 8,
    "deltat": 0.002, "nstep": 2, "pulse": {"count": 1}, "seed": 1,
    "outfile": outfile,
}
# the readmodes of the code: hyperfine 0 point/no contact, 1 point+contact, 2 dft, 3 dft+contact; quadrupole 0 off,
# 2 plain, 3 strained, 4 bubble (the hamiltonian block, A47)
common = {}
if int(hf) != 0:
    common["hyperfine"] = {"dipole": "dft" if int(hf) >= 2 else "point", "contact": int(hf) in (1, 3),
                           "file": hexdir + "/Afile_vertex", "missing_row": "stop"}
if int(qd) != 0:
    common["quadrupole"] = {"structure": {2: "plain", 3: "strained", 4: "bubble"}[int(qd)],
                            "file": hexdir + "/Qfile_vertex"}
    if int(qd) != 2:
        common["quadrupole"]["file_without_defect"] = hexdir + "/Qfile_vertex"
if common:
    cfg["hamiltonian"] = {"common": common}
def merge(into, extra):     # a block of the overrides changes the keys it names
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(into.get(k), dict): merge(into[k], v)
        else: into[k] = v
if extra.strip():
    # the overrides arrive as a bare "key: value, ..." fragment, like run_case's
    merge(cfg, ast.literal_eval("{" + extra + "}"))
from configuration_block import configuration_block   # the A64 block (configuration_block.py)
configuration_block(cfg)
json.dump(cfg, open(dest, "w"), indent=4)
PYEOF
    "$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/$name.json" -v > "$WORK/$name.log" 2>&1
    RC=$?
}

# The h-BN defect sits at [23.75635, 25.71802, 23.38568]; the rotation is taken about it.
HEXROT="'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0]}"
HEXROT_HF="$HEXROT, 'hamiltonian': {'common': {'hyperfine': {'frame': 'file'}}}"
HEXROT_QD="$HEXROT, 'hamiltonian': {'common': {'quadrupole': {'frame': 'file'}}}"
HEXROT_HFQ="$HEXROT, 'hamiltonian': {'common': {'hyperfine': {'frame': 'qubit'}}}"

hex_case hex_hf2_norot 2 0 ""            ; [ $RC -eq 0 ] || fail "hex_hf2_norot failed (rc=$RC)"
hex_case hex_hf2_rot   2 0 "$HEXROT_HF"  ; [ $RC -eq 0 ] || fail "hex_hf2_rot failed (rc=$RC, see $WORK/hex_hf2_rot.log)"
hex_case hex_hf2_rotq  2 0 "$HEXROT_HFQ" ; [ $RC -eq 0 ] || fail "hex_hf2_rotq failed (rc=$RC)"
hex_case hex_qd2_norot 0 2 ""            ; [ $RC -eq 0 ] || fail "hex_qd2_norot failed (rc=$RC)"
hex_case hex_qd2_rot   0 2 "$HEXROT_QD"  ; [ $RC -eq 0 ] || fail "hex_qd2_rot failed (rc=$RC, see $WORK/hex_qd2_rot.log)"

# quadrupole strained/bubble (qd_readmode 3/4) need both the defect and no-defect files. Use one bath spin that
# matches the defect file so this remains a fast reader test, and move every vertex
# in the no-defect copy by +100 A along x. The two v1 reports then prove that the
# second boundary was read from quadrupole.file_without_defect rather than from quadrupole.file
# a second time.
python3 - "$WORK/qd34_bath" "$HEX/Qfile_vertex" "$WORK/Qfile_wo_shifted" <<'PYEOF'
import re, sys
bath, source, shifted = sys.argv[1:4]
with open(bath, "w") as f:
    f.write("1.00000 0.00000 0.00000 0.00000\n")
    f.write("22.45382 26.47004 23.38568 11B\n")
lines = []
for line in open(source):
    m = re.match(r'^(v[1-8]\s*:\s*)([-+0-9.eE]+)(.*)$', line)
    if m:
        line = "%s%.15g%s\n" % (m.group(1), float(m.group(2)) + 100.0,
                                  m.group(3).rstrip("\n"))
    lines.append(line)
with open(shifted, "w") as f:
    f.writelines(lines)
PYEOF

for qdmode in 3 4; do
    hex_case "hex_qd${qdmode}_files" 0 "$qdmode" \
        "'bathfile': ['$WORK/qd34_bath'],
         'hamiltonian': {'common': {'quadrupole': {'file_without_defect': '$WORK/Qfile_wo_shifted',
                                                   'cell': [100.0,100.0,100.0]}}}"
    if [ $RC -ne 0 ]; then
        fail "hex_qd${qdmode}_files failed (rc=$RC, see $WORK/hex_qd${qdmode}_files.log)"
        continue
    fi
    python3 - "$WORK/hex_qd${qdmode}_files.log" <<'PYEOF'
import math, re, sys
text = open(sys.argv[1]).read()
xs = [float(x) for x in re.findall(r'^\s*v1\s*:\s*([-+0-9.eE]+)', text, re.M)]
# printHfInfo_BD emits vertex, face-center and normal blocks (three v1 lines)
# for each file, so the first woDef vertex is xs[3].
ok = len(xs) >= 4 and math.isclose(xs[3] - xs[0], 100.0, abs_tol=2e-6)
sys.exit(0 if ok else 1)
PYEOF
    if [ $? -eq 0 ]; then
        pass "[qd_readmode=$qdmode] woDef boundary comes from quadrupole.file_without_defect"
    else
        fail "[qd_readmode=$qdmode] woDef boundary did not come from its own file"
    fi
done

# strained and bubble (qd_readmode 3, 4): a spin outside the defect file's range is folded by the cell into the
# range of quadrupole.file_without_defect and takes that file's tensor. One 11B 100 A from the defect along x, cell
# [100, 0, 0], and a no-defect copy whose tensors are x10: its quadrupole tensor must be 10 times that of the run
# with the plain copy. Both range forms of the files: MinDif/MaxDif (Qfile) and the eight vertices (Qfile_vertex).
# Before 2.1.0 strained looked the folded spin up in the defect file again, and a vertex file never ended the fold.
python3 - "$WORK" "$HEX" <<'PYEOF'
import sys
work, hexdir = sys.argv[1:3]
for src in ("Qfile", "Qfile_vertex"):
    out = []
    for line in open(hexdir + "/" + src).read().split("\n"):
        f = line.split()
        if len(f) == 13 and f[0] in ("B", "N", "C"):
            v = [float(x) for x in f[1:]]
            line = "\t".join([f[0]] + ["%.5f" % x for x in v[:3]] + ["%.5f" % (10 * x) for x in v[3:]])
        out.append(line)
    open(work + "/" + src + "_x10", "w").write("\n".join(out))
d = [float(x) for x in open(hexdir + "/defect").read().split()[:3]]
open(work + "/qd_far_bath", "w").write("1.0 0.0 0.0 0.0\n%.5f %.5f %.5f 11B\n" % (d[0] - 3.75537 + 100, d[1] + 0.72281, d[2]))
PYEOF
for qdmode in 3 4; do
    for src in Qfile Qfile_vertex; do
        for wo in plain x10; do
            WOFILE="$HEX/$src"; [ $wo = x10 ] && WOFILE="$WORK/${src}_x10"
            hex_case "hex_qd${qdmode}_fold_${src}_$wo" 0 "$qdmode" \
                "'bathfile': ['$WORK/qd_far_bath'], 'rbath': 200,
                 'hamiltonian': {'common': {'quadrupole': {'file': '$HEX/$src', 'file_without_defect': '$WOFILE',
                                                           'cell': [100.0,0.0,0.0]}}}"
            [ $RC -eq 0 ] || fail "hex_qd${qdmode}_fold_${src}_$wo failed (rc=$RC)"
        done
        python3 - "$WORK/hex_qd${qdmode}_fold_${src}_plain.log" "$WORK/hex_qd${qdmode}_fold_${src}_x10.log" <<'PYEOF'
import re, sys
def q(f):
    m = re.search(r'quad\[0\]\s*:\s*\[\s*([-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?)', open(f).read())
    return float(m.group(1)) if m else None
a, b = q(sys.argv[1]), q(sys.argv[2])
sys.exit(0 if a and b and abs(b / a - 10.0) < 1e-3 else 1)
PYEOF
        if [ $? -eq 0 ]; then
            pass "[qd_readmode=$qdmode, $src] a folded spin takes the tensor of quadrupole.file_without_defect"
        else
            fail "[qd_readmode=$qdmode, $src] a folded spin did not take the tensor of quadrupole.file_without_defect"
        fi
    done
done

echo "  [hf_readmode=2, hyperfine.frame=file]"
python3 check_tensors.py "$WORK/hex_hf2_norot.log" "$WORK/hex_hf2_rot.log" \
    --kind hypf --qubit-axis 1 1 1 --expect rotated
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo "  [qd_readmode=2, quadrupole.frame=file]"
python3 check_tensors.py "$WORK/hex_qd2_norot.log" "$WORK/hex_qd2_rot.log" \
    --kind quad --qubit-axis 1 1 1 --expect rotated
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# hyperfine.frame = "qubit" is NOT "leave every tensor alone". The file's tensors are
# already in the computational basis and are kept; the point-dipole fallback that
# readHftensorfile builds for a spin the file does not cover was computed from the
# source geometry and still has to be transformed. Both live in the same BathArray.
echo "  [hyperfine.frame=qubit : file tensors kept, point-dipole fallbacks transformed]"
python3 check_tensors.py "$WORK/hex_hf2_norot.log" "$WORK/hex_hf2_rotq.log" \
    --kind hypf --qubit-axis 1 1 1 --expect mixed
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# ... while the POSITIONS moved in that run all the same, which is what separates
# "the frame flag was honoured" from "the rotation never ran".
python3 check_geometry.py "$WORK/hex_hf2_norot.log" "$WORK/hex_hf2_rotq.log" \
    --qubit-axis 1 1 1 --r0 23.75635 25.71802 23.38568 > "$WORK/hex_rotq_geom.txt" 2>&1
if [ $? -eq 0 ]; then
    pass "[hyperfine.frame=qubit] the bath positions were still rotated"
else
    fail "[hyperfine.frame=qubit] the positions did NOT move -- see $WORK/hex_rotq_geom.txt"
fi

echo
echo "--- [10] invalid and unsupported input is refused ---"
# expect_fail <name> <overrides> <string the message must contain>
expect_fail(){
    local name="$1" extra="$2" needle="$3"
    run_case "$name" "$extra"
    if [ $RC -eq 0 ]; then
        fail "[$name] run SUCCEEDED, expected a fatal error"
    elif grep -q "$needle" "$WORK/$name.log"; then
        pass "[$name] refused (rc=$RC) : \"$needle\""
    else
        fail "[$name] failed (rc=$RC) but without the expected message \"$needle\""
    fi
}

expect_fail defect_zfs_nonzero_E \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[1,1,1],'frame':'file','zfs':{'D':2870,'E':5,'unit':'MHz'}}]}" \
    "zfs.E was removed (only 0 was taken)"
expect_fail defect_zfs_missing_axis \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'frame':'file','zfs':{'D':2870,'unit':'MHz'}}]}" \
    "axis is missing"
# A90: with the rotation on every species states its frame (until then an absent one was "qubit"); with it off the frame
# is not used, so it is not asked for
expect_fail defect_zfs_missing_frame \
    "$NVBASE, $ROT111, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[1,1,1],'zfs':{'D':2870,'unit':'MHz'}}]}" \
    "frame.rotation is on, so Defect.species\[NV_0\].frame must say"
# A94 (G4): a species with no nuclear_positions, hyperfine, efg or zfs (a P1 of detunings) has nothing a frame turns, so
# none is asked for with the rotation on
run_case defect_nodirection_noframe "$P1DEF_NOFRAME, $ROT111"
[ $RC -eq 0 ] \
    && pass "[defect, rotation on] a species of detunings alone runs without a frame (A94)" \
    || fail "[defect, rotation on] a species of detunings alone was refused without a frame (rc=$RC)"
expect_fail defect_old_name \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[1,1,1],'coordinate_frame':'bath','zfs':{'D':2870,'unit':'MHz'}}]}" \
    '"Defect.species\[0\].coordinate_frame" is now "frame"'
expect_fail defect_old_value \
    "$NVBASE, $ROT111, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[1,1,1],'frame':'bath','zfs':{'D':2870,'unit':'MHz'}}]}" \
    'Defect.species\[NV_0\].frame "bath" is now "file"'
run_case defect_axis_norot_noframe "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[0,0,1],'zfs':{'D':2870,'unit':'MHz'}}]}"
[ $RC -eq 0 ] \
    && pass "[defect, rotation off] axis and object-form zfs without a frame run (the frame is not used)" \
    || fail "[defect, rotation off] a species without frame was refused with the rotation off (rc=$RC)"
expect_fail defect_zfs_bad_unit \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'axis':[1,1,1],'frame':'file','zfs':{'D':2870,'unit':'cm-1'}}]}" \
    "Supported frequency units are Hz, kHz, MHz and GHz"
expect_fail defect_zfs_bad_tensor \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'frame':'file','zfs':{'tensor':[[1,0],[0,1]],'unit':'MHz'}}]}" \
    "must be a flat 9-number or nested 3x3 array"
expect_fail defect_detuning_bad_unit \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'detuning':{'values':[[1,'e',1.0]],'unit':'cm-1'}}]}" \
    "Defect.species\[NV_0\].detuning.unit"
expect_fail defect_detuning_missing_values \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'detuning':{'unit':'MHz'}}]}" \
    "detuning.values must be an indexed array"
expect_fail defect_detuning_bad_index \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'detuning':{'values':[[4,'e',1.0]],'unit':'MHz'}}]}" \
    "axis state must be an integer in 1..3"
expect_fail defect_detuning_bad_spin \
    "$NVBASE, 'Defect':{$NVAX, 'species':[{'name':'NV_0','n_nuclear_spins':0,'n_configurations':3,
      'detuning':{'values':[[1,'n',1.0]],'unit':'MHz'}}]}" \
    "spin label must be .e."
expect_fail defect_gyro_bad_unit \
    "$UNITPREFIX 'gyros':{'values':[1.0],'unit':'rpm/G'}}]}" \
    "gyros.unit"
expect_fail defect_gyro_bad_count \
    "$UNITPREFIX 'gyros':{'values':[1.0,2.0],'unit':'radkHz/G'}}]}" \
    "gyros.values must be an array of exactly 1 numbers"
expect_fail defect_eq_bad_unit \
    "$UNITPREFIX $GYRO_INTERNAL, 'quadrupole_moments':{'values':[2.0],'unit':'C m^2'}}]}" \
    "quadrupole_moments.unit"
expect_fail defect_hypf_bad_unit \
    "$UNITPREFIX $GYRO_INTERNAL, 'hyperfine':{'values':[[1,'14N',[1,0,0,0,1,0,0,0,1]]],'unit':'Tesla'}}]}" \
    "hyperfine.unit"
expect_fail defect_efg_bad_unit \
    "$UNITPREFIX $GYRO_INTERNAL, 'efg':{'values':[[1,'14N',[1,0,0,0,1,0,0,0,1]]],'unit':'MHz'}}]}" \
    "efg.unit"
expect_fail defect_hypf_missing_frame \
    "$UNITBASE, $ROT111, 'Defect':{$UNITAX, 'species':[{'name':'DUT','n_nuclear_spins':1,'types':['14N'],'spins':[1.0],
      'gyros':[1.9337792],'n_configurations':1,
      'hyperfine':{'values':[[1,'14N',[1,0,0,0,1,0,0,0,1]]],'unit':'MHz'}}]}" \
    "frame.rotation is on, so Defect.species\[DUT\].frame must say"

expect_fail zeronorm \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [0.0,0.0,0.0]}" \
    "qubit_axis has zero norm"
# A90: bath_axis (it only turned x and y about qubit_axis), reference_qubit (the origin, a translation) and
# qubit_position_frame (one value) were removed, and are refused whatever their value
expect_fail bath_axis_removed \
    "'coordinate_frame_rotation': {'enabled': True, 'bath_axis': [1.0,1.0,0.0], 'qubit_axis': [1.0,1.0,1.0]}" \
    '"frame.rotation.bath_axis" was removed'
expect_fail bath_axis_default_removed \
    "'coordinate_frame_rotation': {'enabled': False, 'bath_axis': [0.0,0.0,1.0], 'qubit_axis': [1.0,1.0,1.0]}" \
    '"frame.rotation.bath_axis" was removed'
expect_fail antiparallel \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [0.0,0.0,-1.0]}" \
    "antiparallel to z"
expect_fail shortaxis \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0]}" \
    "qubit_axis must be an array of exactly 3 numbers"
expect_fail longaxis \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0,0.0]}" \
    "qubit_axis must be an array of exactly 3 numbers"
expect_fail notarray \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': 1.0}" \
    "qubit_axis must be an array of exactly 3 numbers"
# A non-number element is the dangerous one: it used to read as 0.0, turning the axis
# into a DIFFERENT but perfectly valid-looking direction, and the run went on silently.
expect_fail strelement \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,'1',1.0]}" \
    "qubit_axis\[1\] is not a number"
expect_fail nullelement \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,None,1.0]}" \
    "qubit_axis\[1\] is not a number"
expect_fail boolelement \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,True,1.0]}" \
    "qubit_axis\[1\] is not a number"
expect_fail notobject \
    "'coordinate_frame_rotation': True" \
    "frame.rotation must be an object"
expect_fail notbool \
    "'coordinate_frame_rotation': {'enabled': 'true', 'qubit_axis': [1.0,1.0,1.0]}" \
    "must be a boolean"
# The tensor frame has no default: not saying it is the one mistake that cannot be
# caught downstream, because a missed or doubled rotation still looks like a tensor.
expect_fail hfmode_noframe \
    "'hamiltonian': {'common': {'hyperfine': {'dipole': 'dft', 'contact': False, 'file': './inputfiles/gyro_13C'}}}, $ROT111" \
    "frame.rotation is on, so hamiltonian.common.hyperfine.frame must say"
expect_fail qdmode_noframe \
    "'hamiltonian': {'common': {'quadrupole': {'structure': 'plain', 'file': './inputfiles/gyro_13C'}}}, $ROT111" \
    "frame.rotation is on, so hamiltonian.common.quadrupole.frame must say"
expect_fail hfmode_badframe \
    "'hamiltonian': {'common': {'hyperfine': {'dipole': 'dft', 'contact': False, 'file': './inputfiles/gyro_13C', 'frame': 'crystal'}}}, $ROT111" \
    "hamiltonian.common.hyperfine.frame must be .file."
expect_fail qdmode_badframe \
    "'hamiltonian': {'common': {'quadrupole': {'structure': 'plain', 'file': './inputfiles/gyro_13C', 'frame': 'lab'}}}, $ROT111" \
    "hamiltonian.common.quadrupole.frame must be .file."
# the old value and the old place, refused with the new ones (A90)
expect_fail hfmode_oldvalue \
    "'hamiltonian': {'common': {'hyperfine': {'dipole': 'dft', 'contact': False, 'file': './inputfiles/gyro_13C', 'frame': 'bath'}}}, $ROT111" \
    'hamiltonian.common.hyperfine.frame "bath" is now "file"'
expect_fail hfmode_oldplace \
    "'hamiltonian': {'common': {'hyperfine': {'dipole': 'dft', 'contact': False, 'file': './inputfiles/gyro_13C'}}}, 'coordinate_frame_rotation': {'qubit_axis': [1.0,1.0,1.0], 'hf_tensor_frame': 'bath'}" \
    '"frame.rotation.hf_tensor_frame" is now "hamiltonian.common.hyperfine.frame"'
expect_fail qdmode_oldplace \
    "'hamiltonian': {'common': {'quadrupole': {'structure': 'plain', 'file': './inputfiles/gyro_13C'}}}, 'coordinate_frame_rotation': {'qubit_axis': [1.0,1.0,1.0], 'qd_tensor_frame': 'qubit'}" \
    '"frame.rotation.qd_tensor_frame" is now "hamiltonian.common.quadrupole.frame"'
# with the rotation off a written frame is not used, and said so once
run_case hf_frame_norot "'hamiltonian': {'common': {'hyperfine': {'dipole': 'point', 'frame': 'file'}}}"
[ $RC -eq 0 ] && [ "$(grep -c 'input key "hamiltonian.common.hyperfine.frame" is not used by this run: frame.rotation is off' "$WORK/hf_frame_norot.log")" -eq 1 ] \
    && pass "[frame, rotation off] a written hyperfine.frame runs and is warned about once" \
    || fail "[frame, rotation off] no single warning for an unused hyperfine.frame (rc=$RC)"

# An unknown hyperfine.dipole is caught where it is read, before any file is opened. (The out-of-range
# hf_readmode of the old keys cannot be written since A47; the HYPF_UNSET sweep at the end of readHftensorfile stays
# as an internal invariant and is not reachable from input -- which is the point of checking here.)
expect_fail hf_dipole_bad \
    "'hamiltonian': {'common': {'hyperfine': {'dipole': 'x', 'file': './inputfiles/gyro_13C'}}}" \
    "hamiltonian.common.hyperfine.dipole must be"

# An unknown quadrupole structure or value used to slip through (as qd_readmode 9) and leave the run with no
# quadrupole interaction at all, silently, so both of these must stop at option-parse time and never reach
# readQdtensorfile. (The removed qd_readmode 1 cannot be written since A47; an unknown string takes its place.)
expect_fail qd_structure_bad \
    "'hamiltonian': {'common': {'quadrupole': {'structure': 'x', 'file': './inputfiles/gyro_13C'}}}" \
    "hamiltonian.common.quadrupole.structure must be"
grep -q "Read the Quadrupole interaction" "$WORK/qd_structure_bad.log" \
    && fail "[qd_structure_bad] the run reached readQdtensorfile before failing" \
    || pass "[qd_structure_bad] rejected before readQdtensorfile ran"

expect_fail qd_value_bad \
    "'hamiltonian': {'common': {'quadrupole': 'on'}}" \
    "hamiltonian.common.quadrupole must be \"off\" or an object"
grep -q "Read the Quadrupole interaction" "$WORK/qd_value_bad.log" \
    && fail "[qd_value_bad] the run reached readQdtensorfile before failing" \
    || pass "[qd_value_bad] an unknown quadrupole value rejected before readQdtensorfile ran"

echo
echo "--- [S1] option surface : new name, alias, removed keys ---"
run_case newname "$ROT111"
[ $RC -eq 0 ] || fail "newname run failed (rc=$RC)"
ALIAS111="'bath_coordinate_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0]}"
# the old names (bath_coordinate_rotation, accepted with a warning until A65; coordinate_frame_rotation) are refused with
# the new place frame.rotation
run_case aliased "$ALIAS111"
if [ $RC -ne 0 ] && grep -q '"bath_coordinate_rotation" is now "frame.rotation"' "$WORK/aliased.log"; then
    pass "[alias] bath_coordinate_rotation is refused with the new place frame.rotation"
else
    fail "[alias] bath_coordinate_rotation was not refused with frame.rotation (rc=$RC)"
fi

# frame.rotation is on when the block is given (A65): without "enabled" it is the run with "enabled": true (newname);
# "enabled": false is the run without the block (the [disabled] case above)
run_case noenabled "'coordinate_frame_rotation': {'qubit_axis': [1.0,1.0,1.0]}"
same=1; for suf in _noDiv _wiDiv; do cmp -s "$WORK/newname$suf" "$WORK/noenabled$suf" || same=0; done
[ $RC -eq 0 ] && [ $same -eq 1 ] \
    && pass "[frame.rotation] a block without \"enabled\" is on: byte-identical to \"enabled\": true" \
    || fail "[frame.rotation] a block without \"enabled\" is not the rotation on (rc=$RC)"

grep -q "is deprecated" "$WORK/newname.log" \
    && fail "[newname] the new name wrongly warns about deprecation" \
    || pass "[newname] the new name warns about nothing"

expect_fail bothnames \
    "$ROT111, 'bath_coordinate_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0]}" \
    "\"bath_coordinate_rotation\" is now \"frame.rotation\""
expect_fail refq_removed \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0], 'reference_qubit': 'q0'}" \
    '"frame.rotation.reference_qubit" was removed'
expect_fail posframe_removed \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0], 'qubit_position_frame': 'bath'}" \
    '"frame.rotation.qubit_position_frame" was removed'

echo
echo "--- [S2] bfield with the rotation on ---"
# With the rotation on, the computational z axis IS the qubit axis, and "bfield" is read
# in that computational frame. So [0,0,500] is 500 G along +qubit_axis, [0,0,-500] is the
# same magnitude the other way, and anything with a transverse component is a tilted
# field, which this version does not support alongside the rotation.
for spec in "pos:[0.0,0.0,500.0]" "neg:[0.0,0.0,-500.0]" "zero:[0.0,0.0,0.0]"; do
    tag=${spec%%:*}; vec=${spec#*:}
    run_case "bf_$tag" "$ROT111, 'bfield': $vec"
    [ $RC -eq 0 ] \
        && pass "[bfield $tag] rotation + bfield=$vec accepted" \
        || fail "[bfield $tag] rotation + bfield=$vec was REFUSED (rc=$RC)"
done

run_case bf_scalar "$ROT111, 'bfield': 500"
[ $RC -eq 0 ] \
    && pass "[bfield scalar] rotation + scalar bfield=500 accepted" \
    || fail "[bfield scalar] rotation + scalar bfield=500 was refused (rc=$RC)"

same=1; for suf in _noDiv _wiDiv; do cmp -s "$WORK/bf_pos$suf" "$WORK/bf_scalar$suf" || same=0; done
[ $same -eq 1 ] \
    && pass "[bfield] scalar 500 and [0,0,500] give byte-identical output" \
    || fail "[bfield] scalar 500 and [0,0,500] DIFFER"

# +z and -z are both along the qubit-axis line, so both are accepted -- and they are
# genuinely different physics, so they must not give the same answer.
if cmp -s "$WORK/bf_pos_noDiv" "$WORK/bf_neg_noDiv"; then
    fail "[bfield] +500 and -500 gave identical output"
else
    pass "[bfield] +500 and -500 are both accepted and give different results"
fi

# A scalar bfield is Bz in the computational frame, so it composes with the rotation instead of conflicting with it.
# (Until A63 this was checked with -B 500 over "bfield": [0,0,0]; -B was removed and its value is the scalar bfield.)
python3 check_srcfield.py "$WORK/bf_scalar.log" --qubit-axis 1 1 1 --bz 500 | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# R^T takes the computational field back to the source frame: R^T*[0,0,B] = B*q_hat.
python3 check_srcfield.py "$WORK/bf_pos.log" --qubit-axis 1 1 1 --bz 500
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

expect_fail bf_tilted_x "$ROT111, 'bfield': [1.0,0.0,500.0]" \
    "requires bfield to be aligned"
expect_fail bf_tilted_y "$ROT111, 'bfield': [0.0,2.0,500.0]" \
    "requires bfield to be aligned"
expect_fail bf_transverse "$ROT111, 'bfield': [1.0,1.0,0.0]" \
    "requires bfield to be aligned"

expect_fail bf_string_element "$ROT111, 'bfield': [0.0,0.0,'500']" \
    "must be a number"
expect_fail bf_short_array "$ROT111, 'bfield': [0.0,500.0]" \
    "bfield is too short"

# Without the rotation the check is a no-op and an arbitrary field vector still works.
run_case bf_tilted_norot "'bfield': [1.0,0.0,500.0]"
[ $RC -eq 0 ] \
    && pass "[bfield] a tilted field is still fine without the rotation" \
    || fail "[bfield] a tilted field was refused even without the rotation (rc=$RC)"

expect_fail bf_removed_option \
    "$ROT111, 'magnetic_field_axis': {'axis': [1.0,1.0,1.0]}" \
    "was removed"

echo
echo "--- [S3] defect_axis_reference ---"
# removed on 2026-10-15 (A65): it only named the rotation's qubit axis in the log; refused whatever its value
expect_fail defaxis           "$ROT111, 'defect_axis_reference': 'qubit_axis'" \
    'input key "defect_axis_reference" was removed'
expect_fail defaxis_badstring "$ROT111, 'defect_axis_reference': 'bfield'" \
    'input key "defect_axis_reference" was removed'
expect_fail defaxis_norot     "'defect_axis_reference': 'qubit_axis'" \
    'input key "defect_axis_reference" was removed'

echo
echo "--- [MQ] multi-qubit rotation ---"
# Two spin-1 qubits, gCCE (cCCE is single-qubit only and stays that way). NV1 is placed
# off every symmetry direction so that a wrong R or a wrong rotation centre shows up.
mq_case(){
    local name="$1" extra="$2"
    python3 - "$WORK/$name.json" "$WORK/$name" "$extra" <<'PYEOF'
import ast, json, sys
dest, outfile, extra = sys.argv[1:4]
cfg = {
    "method": "gCCE", "quantity": "coherence",
    "gyrofile": "./inputfiles/gyro_13C", "bathfile": ["./inputfiles/bath_a"],
    "order": 1, "bfield": [0.0, 0.0, 500.0], "rbath": 10, "rdip": 10,
    "deltat": 0.1, "nstep": 2, "pulse": {"count": 1}, "s_range": [1, 1], "seed": 1,
    "Qubit": {"qubit": [
        {"name": "NV0", "spin": 1.0, "gyro": -17608.597050,
         "xyz": [10.0, 20.0, 30.0]},
        {"name": "NV1", "spin": 1.0, "gyro": -17608.597050,
         "xyz": [17.0, 23.0, 41.0]}]},
    "outfile": outfile,
}
if extra.strip():
    cfg.update(ast.literal_eval("{" + extra + "}"))
from configuration_block import configuration_block   # the A64 block (configuration_block.py)
configuration_block(cfg)
json.dump(cfg, open(dest, "w"), indent=4)
PYEOF
    "$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/$name.json" -v > "$WORK/$name.log" 2>&1
    RC=$?
}

# the rotation is about the first qubit, NV0 (A90; reference_qubit chose it until then)
MQROT="'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0]}"
MQID="'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [0.0,0.0,1.0]}"

echo "  [A] two-qubit identity rotation changes nothing"
mq_case mq_norot ""      ; [ $RC -eq 0 ] || fail "mq_norot run failed (rc=$RC, see $WORK/mq_norot.log)"
mq_case mq_ident "$MQID" ; [ $RC -eq 0 ] || fail "mq_ident run failed (rc=$RC, see $WORK/mq_ident.log)"

same=1; for suf in _state1_noDiv _state1_wiDiv; do cmp -s "$WORK/mq_norot$suf" "$WORK/mq_ident$suf" || same=0; done
[ $same -eq 1 ] \
    && pass "[MQ/A] identity rotation gives byte-identical gCCE coherence" \
    || fail "[MQ/A] identity rotation CHANGED the gCCE coherence"

# The rotated run prints its tensors again after the rotation, so the two logs do not
# have the same number of lines. The checkers key on the index and take the last
# occurrence, which is the state that actually reaches the Hamiltonian.
for what in intmap hypf; do
    python3 check_tensors.py "$WORK/mq_norot.log" "$WORK/mq_ident.log" \
        --kind $what --qubit-axis 0 0 1 --expect unchanged | sed 's/^/  /'
    if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi
done

echo "  [B] two-qubit [001] -> [111], about the first qubit NV0"
mq_case mq_rot "$MQROT" ; [ $RC -eq 0 ] || fail "mq_rot run failed (rc=$RC, see $WORK/mq_rot.log)"
python3 check_qubits.py "$WORK/mq_rot.log" --qubit-axis 1 1 1 \
    --origin NV0 --expect-nqubit 2
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

python3 check_geometry.py "$WORK/mq_norot.log" "$WORK/mq_rot.log" \
    --qubit-axis 1 1 1 --r0 10 20 30 --skip-qubit-check
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

python3 check_tensors.py "$WORK/mq_norot.log" "$WORK/mq_rot.log" \
    --kind hypf --qubit-axis 1 1 1 --expect rotated
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo "  [E,F] explicit and automatic intmap tensors"
# NV0 gets a ZFS written in the source frame; NV1 gets a DIFFERENT one, standing for an
# NV with another crystallographic orientation. Both are transformed by the same global
# R -- that is what the common computational frame means -- and the off-diagonal
# elements that appear are correct, not something to diagonalize away.
ZFS0="[[2870.0,0.0,0.0],[0.0,2870.0,0.0],[0.0,0.0,-5740.0]]"
ZFS1="[[1900.0,120.0,-80.0],[120.0,2100.0,45.0],[-80.0,45.0,-4000.0]]"
PAIR="[[12.0,3.0,-5.0],[3.0,-7.0,2.0],[-5.0,2.0,-5.0]]"

MQ_QUBITS_PREFIX="'Qubit': {'qubit': [
    {'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},
    {'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}]"

IM_BATH="'intmap': [{'between':['NV0','NV0'],'frame':'file','tensor':$ZFS0},
                    {'between':['NV1','NV1'],'frame':'file','tensor':$ZFS1},
                    {'between':['NV0','NV1'],'frame':'file','tensor':$PAIR}]"
IM_QUBIT="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','tensor':$ZFS0},
                     {'between':['NV1','NV1'],'frame':'qubit','tensor':$ZFS1},
                     {'between':['NV0','NV1'],'frame':'qubit','tensor':$PAIR}]"

mq_case mq_im_norot   "'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}], $IM_BATH}"
[ $RC -eq 0 ] || fail "mq_im_norot run failed (rc=$RC)"
mq_case mq_im_bath    "$MQROT, 'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}], $IM_BATH}"
[ $RC -eq 0 ] || fail "mq_im_bath run failed (rc=$RC, see $WORK/mq_im_bath.log)"
mq_case mq_im_qubit   "$MQROT, 'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}], $IM_QUBIT}"
[ $RC -eq 0 ] || fail "mq_im_qubit run failed (rc=$RC, see $WORK/mq_im_qubit.log)"

# Qubit.intmap historically means kHz when unit is absent. Exercise every supported
# explicit unit against that legacy path using the same non-zero self tensor.
IM_ONE_LEGACY="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','tensor':$ZFS0}]"
IM_ONE_KHZ="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':'kHz','tensor':$ZFS0}]"
IM_ONE_MHZ="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':'MHz',
    'tensor':[[2.87,0.0,0.0],[0.0,2.87,0.0],[0.0,0.0,-5.74]]}]"
IM_ONE_GHZ="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':'GHz',
    'tensor':[[0.00287,0.0,0.0],[0.0,0.00287,0.0],[0.0,0.0,-0.00574]]}]"
IM_ONE_HZ="'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':'Hz',
    'tensor':[[2870000.0,0.0,0.0],[0.0,2870000.0,0.0],[0.0,0.0,-5740000.0]]}]"

for spec in legacy kHz MHz GHz Hz; do
    case "$spec" in
        legacy) im="$IM_ONE_LEGACY" ;;
        kHz)    im="$IM_ONE_KHZ" ;;
        MHz)    im="$IM_ONE_MHZ" ;;
        GHz)    im="$IM_ONE_GHZ" ;;
        Hz)     im="$IM_ONE_HZ" ;;
    esac
    mq_case "mq_im_unit_$spec" "${MQ_QUBITS_PREFIX}, ${im}}"
    [ $RC -eq 0 ] || fail "[intmap unit/$spec] run failed (rc=$RC)"
done

for spec in kHz MHz GHz Hz; do
    python3 check_tensors.py "$WORK/mq_im_unit_legacy.log" "$WORK/mq_im_unit_$spec.log" \
        --kind intmap --qubit-axis 0 0 1 --expect unchanged --tol 1e-6 \
        | sed 's/^/  /'
    if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi
done

echo "    frame=file : ZFS of both qubits and the explicit pair tensor"
python3 check_tensors.py "$WORK/mq_im_norot.log" "$WORK/mq_im_bath.log" \
    --kind intmap --qubit-axis 1 1 1 --expect rotated --tol 0.5 | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo "    frame=qubit : already computational, left alone"
python3 check_tensors.py "$WORK/mq_im_norot.log" "$WORK/mq_im_qubit.log" \
    --kind intmap --qubit-axis 1 1 1 --expect unchanged --tol 0.5 | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# NV1's ZFS is not axially symmetric about the new z, so it must still carry off-diagonal
# elements after the transform -- nothing diagonalizes it.
python3 - "$WORK/mq_im_bath.log" <<'PYEOF'
import re, sys
pat = re.compile(r"^\s*intmap\[1\]\[1\]\s*:\s*\[(.+)\]\s*$")
cplx = re.compile(r"(-?\d+\.\d+)[-+]\d+\.\d+j")
vals = None
for line in open(sys.argv[1]):
    m = pat.match(line)
    if m:
        vals = [float(v) for v in cplx.findall(m.group(1))]
off = max(abs(vals[i]) for i in (1, 2, 3, 5, 6, 7)) if vals else 0.0
ok = vals is not None and off > 1.0
print(("  ok    " if ok else "  FAIL  ") +
      "NV1 ZFS keeps its off-diagonal elements after R*T*R^T (max |off-diag| = %.3g)" % off)
sys.exit(0 if ok else 1)
PYEOF
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# The automatic pair tensor: no explicit intmap at all, so intmap[0][1] is the
# point-dipole tensor setIntmap_dipAuto built from the source qubit positions.
python3 check_tensors.py "$WORK/mq_norot.log" "$WORK/mq_rot.log" \
    --kind intmap --qubit-axis 1 1 1 --expect mixed | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo "  [D,G,H] validation"
MQ2="'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}]"

# mq_expect_fail <name> <extra> <needle>
mq_expect_fail(){
    local name="$1" extra="$2" needle="$3"
    mq_case "$name" "$extra"
    if [ $RC -eq 0 ]; then
        fail "[$name] run SUCCEEDED, expected a fatal error"
    elif grep -q "$needle" "$WORK/$name.log"; then
        # The point of the early validator: no reader may have run.
        if grep -qE "Read the (Hyperfine|Quadrupole) interaction from DFT inputfile" "$WORK/$name.log"; then
            fail "[$name] refused, but only after a tensor reader had already started"
        else
            pass "[$name] refused before any reader ran (rc=$RC) : \"$needle\""
        fi
    else
        fail "[$name] failed (rc=$RC) but without the expected message \"$needle\""
    fi
}

mq_expect_fail mq_xyz_string \
    "$MQROT, 'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,'41']}]}" \
    "must be a number"

mq_expect_fail mq_refq_removed \
    "'coordinate_frame_rotation': {'enabled': True, 'qubit_axis': [1.0,1.0,1.0], 'reference_qubit': 'NV1'}" \
    '"frame.rotation.reference_qubit" was removed'
mq_expect_fail mq_dupname \
    "$MQROT, 'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}]}" \
    "share the name"

mq_expect_fail mq_im_noframe \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'tensor':$ZFS0}]}" \
    "frame.rotation is on, so Qubit.intmap\[\"NV0\",\"NV0\"\].frame must say"
mq_expect_fail mq_im_badframe \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'frame':'lab','tensor':$ZFS0}]}" \
    "\].frame must be .file."
mq_expect_fail mq_im_frame_notstring \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'frame':7,'tensor':$ZFS0}]}" \
    "\].frame must be .file."
mq_expect_fail mq_im_old_name \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'tensor_frame':'bath','tensor':$ZFS0}]}" \
    '"Qubit.intmap\[0\].tensor_frame" is now "frame"'

mq_expect_fail mq_im_badunit \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':'radkHz','tensor':$ZFS0}]}" \
    "Supported frequency units are Hz, kHz, MHz and GHz"
mq_expect_fail mq_im_unit_notstring \
    "$MQROT, $MQ2, 'intmap': [{'between':['NV0','NV0'],'frame':'qubit','unit':7,'tensor':$ZFS0}]}" \
    ".unit must be a string"

mq_expect_fail mq_hfmode \
    "$MQROT, 'hamiltonian': {'common': {'hyperfine': {'dipole': 'dft', 'contact': False, 'file': './inputfiles/gyro_13C', 'frame': 'file'}}}" \
    "a hyperfine tensor file (hamiltonian.common.hyperfine) with 2 qubits is not supported"
mq_expect_fail mq_qdmode \
    "$MQROT, 'hamiltonian': {'common': {'quadrupole': {'structure': 'plain', 'file': './inputfiles/gyro_13C', 'frame': 'file'}}}" \
    "a quadrupole tensor file (hamiltonian.common.quadrupole) with 2 qubits is not supported"

echo
echo "--- [A/B] hand-transformed input vs automatic rotation ---"
# The strongest check available: run A is the SAME physical system written out by hand
# in the qubit-aligned frame, with no rotation option at all; run B is the original
# input with the rotation on. Every component check above verifies one transformation
# in isolation -- this one would catch a frame-dependent quantity that nothing
# transforms, because A never had one to begin with.
#
# Byte identity is not the criterion. The two take different floating-point routes to
# the same numbers (A parses coordinates that B computes), so the test is numerical;
# whether the files also came out byte-identical is reported as a bonus.
mkdir -p "$WORK/ab"

echo "  [A/B-1] point-dipole system (hf_readmode=0)"
python3 make_manual_frame.py --qubit-axis 1 1 1 --r0 10 20 30 \
    --bath inputfiles/bath_a "$WORK/ab/bath_a_qframe" | sed 's/^/  /'

# A : pre-rotated bath, no rotation, field written directly in the computational frame
write_ccein "$WORK/ab/A1.json" "{'outfile': '$WORK/ab/A1',
    'bathfile': ['$WORK/ab/bath_a_qframe'], 'bfield': [0.0,0.0,500.0]}"
"$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/ab/A1.json" -v > "$WORK/ab/A1.log" 2>&1
rcA=$?
# B : original bath, rotation on, field as a magnitude along the physical [111]
write_ccein "$WORK/ab/B1.json" "{'outfile': '$WORK/ab/B1', $ROT111, 'bfield': [0.0,0.0,500.0]}"
"$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/ab/B1.json" -v > "$WORK/ab/B1.log" 2>&1
rcB=$?
[ $rcA -eq 0 ] && [ $rcB -eq 0 ] || fail "[A/B-1] a run failed (A=$rcA B=$rcB)"

python3 check_ab.py "$WORK/ab/A1.log" "$WORK/ab/B1.log" \
    --a-out "$WORK/ab/A1" --b-out "$WORK/ab/B1" --kind hypf
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "  [A/B-2] DFT hyperfine and quadrupole tensors (hf_readmode=2, qd_readmode=2)"
QR0="23.75635 25.71802 23.38568"
python3 make_manual_frame.py --qubit-axis 1 1 1 --r0 $QR0 \
    --bath "$HEX/bath_1" "$WORK/ab/bath_1_qframe" \
    --tensor 13 file "$HEX/Afile_vertex" "$WORK/ab/Afile_qframe" \
    --tensor 12 file "$HEX/Qfile_vertex" "$WORK/ab/Qfile_qframe" | sed 's/^/  /'

hex_case ab/A2 2 2 "'bathfile': ['$WORK/ab/bath_1_qframe'],
    'hamiltonian': {'common': {'hyperfine': {'file': '$WORK/ab/Afile_qframe'},
                               'quadrupole': {'file': '$WORK/ab/Qfile_qframe'}}}"
rcA=$RC
hex_case ab/B2 2 2 "$HEXROT, 'hamiltonian': {'common': {'hyperfine': {'frame': 'file'}, 'quadrupole': {'frame': 'file'}}}"
rcB=$RC
[ $rcA -eq 0 ] && [ $rcB -eq 0 ] || fail "[A/B-2] a run failed (A=$rcA B=$rcB)"

python3 check_ab.py "$WORK/ab/A2.log" "$WORK/ab/B2.log" \
    --a-out "$WORK/ab/A2" --b-out "$WORK/ab/B2" --kind hypf --kind quad --pos-tol 2e-3
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "  [A/B-3] hf_readmode=1 : the file supplies only the isotropic Fermi contact"
# readmode 1 throws the file's anisotropic tensor away and rebuilds it as a point-dipole
# tensor from the SOURCE geometry, keeping only the scalar fc. So the stored tensor is
# file-frame whatever hyperfine.frame says, and both spellings below must give the same
# answer as the hand-transformed input.
python3 make_manual_frame.py --qubit-axis 1 1 1 --r0 $QR0 \
    --tensor 13 qubit "$HEX/Afile_vertex" "$WORK/ab/Afile_poskept" | sed 's/^/  /'

hex_case ab/A3 1 0 "'bathfile': ['$WORK/ab/bath_1_qframe'], 'hamiltonian': {'common': {'hyperfine': {'file': '$WORK/ab/Afile_qframe'}}}"
rcA=$RC
hex_case ab/B3 1 0 "$HEXROT"
rcB=$RC
[ $rcA -eq 0 ] && [ $rcB -eq 0 ] || fail "[A/B-3] a run failed (A=$rcA B=$rcB)"

grep -q "hyperfine.frame is not used with dipole" "$WORK/ab/B3.log" \
    && fail "[A/B-3] hyperfine.frame was reported as given, but it was omitted" \
    || pass "[A/B-3] hyperfine.frame is not required at hf_readmode = 1"

echo "    hyperfine.frame omitted:"
python3 check_ab.py "$WORK/ab/A3.log" "$WORK/ab/B3.log" \
    --a-out "$WORK/ab/A3" --b-out "$WORK/ab/B3" --kind hypf | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

# Spelling it out must change nothing, and must say so.
hex_case ab/B3q 1 0 "$HEXROT_HFQ"
[ $RC -eq 0 ] || fail "[A/B-3] the hyperfine.frame=qubit run failed (rc=$RC)"

grep -q "hyperfine.frame is not used with dipole" "$WORK/ab/B3q.log" \
    && pass "[A/B-3] hyperfine.frame=qubit at readmode 1 is reported as not used" \
    || fail "[A/B-3] no note that hyperfine.frame does not apply at readmode 1"

echo "    hyperfine.frame=qubit:"
python3 check_ab.py "$WORK/ab/A3.log" "$WORK/ab/B3q.log" \
    --a-out "$WORK/ab/A3" --b-out "$WORK/ab/B3q" --kind hypf | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "  [A/B-4] hf_readmode=2 and 3 with hyperfine.frame=qubit"
# The A side declares the same thing the B side does: the file's components are already
# computational, so they are copied across unrotated -- only the position columns and
# the boundary vertices move, because A's bath has moved.
for M in 2 3; do
    hex_case "ab/A4m$M" $M 0 "'bathfile': ['$WORK/ab/bath_1_qframe'], 'hamiltonian': {'common': {'hyperfine': {'file': '$WORK/ab/Afile_poskept'}}}"
    rcA=$RC
    hex_case "ab/B4m$M" $M 0 "$HEXROT_HFQ"
    rcB=$RC
    [ $rcA -eq 0 ] && [ $rcB -eq 0 ] || fail "[A/B-4] hf_readmode=$M failed (A=$rcA B=$rcB)"

    echo "    hf_readmode=$M :"
    python3 check_ab.py "$WORK/ab/A4m$M.log" "$WORK/ab/B4m$M.log" \
        --a-out "$WORK/ab/A4m$M" --b-out "$WORK/ab/B4m$M" --kind hypf | sed 's/^/  /'
    if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

    # ... and the split really is a split, for readmode 3 as well as 2.
    hex_case "ab/N4m$M" $M 0 ""
    python3 check_tensors.py "$WORK/ab/N4m$M.log" "$WORK/ab/B4m$M.log" \
        --kind hypf --qubit-axis 1 1 1 --expect mixed | sed 's/^/  /'
    if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi
done

echo
echo "  [A/B-5] multi-qubit : hand-transformed input vs automatic rotation"
# A : qubit positions, bath positions and intmap tensors all transformed by hand into
#     the computational frame, run with no rotation option.
# B : the same system left in the source frame, with the rotation doing the work.
python3 make_manual_frame.py --qubit-axis 1 1 1 --r0 10 20 30 \
    --bath inputfiles/bath_a "$WORK/ab/bath_a_mq" | sed 's/^/  /'

python3 - "$WORK/ab/A5.json" "$WORK/ab/A5" "$WORK/ab/bath_a_mq" <<'PYEOF'
import json, math, sys
dest, outfile, bathfile = sys.argv[1:4]

def norm(v):
    n = math.sqrt(sum(c*c for c in v)); return [c/n for c in v]
def cross(u, v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
q, b = norm([1,1,1]), norm([0,0,1])
ez = q; ex = norm(cross(q, b)); ey = cross(ez, ex)
R = [ex, ey, ez]
r0 = [10.0, 20.0, 30.0]

def rot_pos(p):
    d = [p[k]-r0[k] for k in range(3)]
    return [r0[k] + sum(R[k][j]*d[j] for j in range(3)) for k in range(3)]
def rot_T(T):
    return [[sum(R[i][k]*T[k][l]*R[j][l] for k in range(3) for l in range(3))
             for j in range(3)] for i in range(3)]

ZFS0 = [[2870.0,0.0,0.0],[0.0,2870.0,0.0],[0.0,0.0,-5740.0]]
ZFS1 = [[1900.0,120.0,-80.0],[120.0,2100.0,45.0],[-80.0,45.0,-4000.0]]
PAIR = [[12.0,3.0,-5.0],[3.0,-7.0,2.0],[-5.0,2.0,-5.0]]

cfg = {
    "method": "gCCE", "quantity": "coherence",
    "gyrofile": "./inputfiles/gyro_13C", "bathfile": [bathfile],
    "order": 1, "bfield": [0.0, 0.0, 500.0], "rbath": 10, "rdip": 10,
    "deltat": 0.1, "nstep": 2, "pulse": {"count": 1}, "s_range": [1, 1], "seed": 1,
    "Qubit": {"qubit": [
        {"name": "NV0", "spin": 1.0, "gyro": -17608.597050,
         "xyz": rot_pos([10.0,20.0,30.0])},
        {"name": "NV1", "spin": 1.0, "gyro": -17608.597050,
         "xyz": rot_pos([17.0,23.0,41.0])}],
        "intmap": [
            {"between": ["NV0","NV0"], "tensor": rot_T(ZFS0)},
            {"between": ["NV1","NV1"], "tensor": rot_T(ZFS1)},
            {"between": ["NV0","NV1"], "tensor": rot_T(PAIR)}]},
    "outfile": outfile,
}
from configuration_block import configuration_block   # the A64 block (configuration_block.py)
configuration_block(cfg)
json.dump(cfg, open(dest, "w"), indent=4)
PYEOF
"$MPIRUN" -n 1 "$CCEX_BIN" -f "$WORK/ab/A5.json" -v > "$WORK/ab/A5.log" 2>&1
rcA=$?
mq_case ab/B5 "$MQROT, 'Qubit': {'qubit': [{'name':'NV0','spin':1.0,'gyro':-17608.597050,'xyz':[10.0,20.0,30.0]},{'name':'NV1','spin':1.0,'gyro':-17608.597050,'xyz':[17.0,23.0,41.0]}], $IM_BATH}"
rcB=$RC
[ $rcA -eq 0 ] && [ $rcB -eq 0 ] || fail "[A/B-5] a run failed (A=$rcA B=$rcB)"

python3 check_ab.py "$WORK/ab/A5.log" "$WORK/ab/B5.log" \
    --a-out "$WORK/ab/A5_state1" --b-out "$WORK/ab/B5_state1" \
    --kind hypf --kind intmap --nqubit 2 --tensor-tol 0.5 | sed 's/^/  /'
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "--- [Z] gCCE with E != 0 turned about z, the qubit state turned too (A91) ---"
# The rotation about the qubit axis is a symmetry once psi0 turns with the system; with the m_s states kept the run is
# another preparation phase and differs (physically, check_analytic has the exact answer). check_zrot_gcce.py says how.
python3 check_zrot_gcce.py "$CCEX_BIN" "$MPIRUN" "$WORK/zrot"
if [ $? -eq 0 ]; then npass=$((npass+1)); else nfail=$((nfail+1)); fi

echo
echo "--- [12] the bath files on disk are untouched ---"
if [ "$MD5_BEFORE" = "$(md5sum inputfiles/* | sort)" ]; then
    pass "every file in inputfiles/ has the same md5 as before the suite ran"
else
    fail "an input file CHANGED on disk"
    diff <(echo "$MD5_BEFORE") <(md5sum inputfiles/* | sort)
fi

###############################################################################
echo
echo "############################################################"
echo "# $npass passed, $nfail failed"
echo "############################################################"
echo
[ $nfail -eq 0 ] || exit 1
exit 0
