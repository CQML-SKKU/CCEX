#!/usr/bin/env bash
#
# The command line of A63 (2026-10-15): -f FILE, -v, -h/--help, -o DIR and -e.
#   X1  -e writes ccex_example.jsonc alone, and with its bathfile pointed at 00_quickstart's bath (a sed copy, run.jsonc,
#       the path written with a '//' inside the string) it runs to the end (2 ranks): the // comments are read, a '//'
#       in a string is kept, and the outputs are run_* (output.file's default without .jsonc)
#   X5  the example has no trailing comments (each explanation is on comment lines above its key, indented by a multiple
#       of 4), and the ':' of the keys that start a line (active or commented out) are in one column per object:
#       "key" : value, one space before the ':' of the longest key and one after (objects whose '{' ends its line)
#   X2  -e -o DIR writes the same file into DIR (made when missing)
#   X3  the example names every input key the readers accept (example_keys.py: the contract on the source)
#   X6  the example's auto_convergence lines, uncommented, are accepted by the runner, for every choice of the alternative
#       lines (a key written more than once: the line that finds a value and the one that fixes it; example_autoconv.py)
#   O1  no output.file, no -o: the outputs are <input name without .json> in the working directory
#   O2  no output.file, -o DIR (two levels, missing): DIR/<input name>, byte for byte O1's curves
#   O3  a relative output.file with -o DIR: DIR/<output.file>; -o before -f is the same
#   O4  an absolute output.file without -o: that path
#   O5  output.cluster_only true with -o DIR: DIR/<input name>.clusters.json and no curves (the plan of the old -P)
#   O6  output.path (A65) of two missing levels: made, the outputs are output.path/<input name>, O1's curves
#   O7  -o DIR beside output.path: DIR wins, output.path is not made
# The refusals (removed options, an absolute outfile with -o, a file in DIR's place, -e over a file or with -f) are in
# ../check_refusals.
#
# Usage:  bash test/code_verification/check_input_forms/cli.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/cli"
rm -rf "$W"; mkdir -p "$W/e1"
ok(){ echo "  PASS  $1"; PASS=$((PASS+1)); }
no(){ echo "  FAIL  $1"; FAIL=$((FAIL+1)); }
same_curves(){ cmp -s "$1_noDiv" "$2_noDiv" && cmp -s "$1_wiDiv" "$2_wiDiv"; }

( cd "$W/e1" && "$MPIRUN" -n 1 "$BIN_DEV" -e > e.log 2>&1 < /dev/null ) && [ -s "$W/e1/ccex_example.jsonc" ] \
    && [ "$(ls "$W/e1" | tr '\n' ' ')" = "ccex_example.jsonc e.log " ] || no "X1 -e did not write ccex_example.jsonc alone (see $W/e1/e.log)"
sed 's|"bathfile" *: \["my_bath.txt"\]|"bathfile": ["'"$REPO"'/example/tutorial/basic/01_first_run/bath//bath_Diamond_NV_natab_1"]|' \
    "$W/e1/ccex_example.jsonc" > "$W/e1/run.jsonc"
grep -q 'bath//bath_Diamond' "$W/e1/run.jsonc" || no "X1 the example has no \"bathfile\" : [\"my_bath.txt\"] line to point at a bath"
( cd "$W/e1" && "$MPIRUN" -n 2 "$BIN_DEV" -f run.jsonc > run.log 2>&1 < /dev/null ); rc=$?
lines=$(grep -c . "$W/e1/run_wiDiv" 2>/dev/null || echo 0)
nstep=$(python3 -c "import sys; sys.path.insert(0, '$REPO/example/tools/InputKeys'); import json; from convert_input_keys import strip_comments; print(json.loads(strip_comments(open(sys.argv[1]).read()))['nstep'])" "$W/e1/run.jsonc")
if [ $rc -eq 0 ] && [ "$lines" -eq "$nstep" ] && [ -s "$W/e1/run_noDiv" ] && ! grep -q "Warning: input key" "$W/e1/run.log"
then ok "X1 -e writes the template, and with a bath it runs ($lines time points, no unread key, outputs run_*)"
else no "X1 the template with a bath did not run: exit $rc, $lines lines for nstep $nstep (see $W/e1/run.log)"; fi

# X5: a trailing comment is a '//' outside a string after the code of an active line or of a commented-out key line
# ('//' and two or more spaces, then a JSON character); there must be none. The ':' check follows the brackets.
bad=$(python3 - "$W/e1/ccex_example.jsonc" <<'PY'
import re, sys
lines = open(sys.argv[1]).read().splitlines()
body = lines.index('{')
bad, keys, stack, nid, block = [], {}, [], 0, set()
for n, line in enumerate(lines, 1):
    if re.match(r'\s*//', line) and not re.match(r'//\s{2,}["{}\[\]]', line):     # a whole-line comment
        if (len(line) - len(line.lstrip())) % 4: bad.append(n)
        continue
    code, ins, esc = line, False, False
    for i in range(2 if line.startswith('//') else 0, len(line) - 1):
        c = line[i]
        if ins:
            if esc: esc = False
            elif c == '\\': esc = True
            elif c == '"': ins = False
        elif c == '"': ins = True
        elif line.startswith('//', i):
            bad.append(n)
            code = line[:i].rstrip()
            break
    if n <= body: continue
    code = '  ' + code[2:] if code.startswith('//') else code
    if re.match(r'\s*"', code) and stack and stack[-1][0] == '{':          # a key starting the line
        m = re.match(r'\s*"(?:[^"\\]|\\.)*"( *): (?! )', code)
        keys.setdefault(stack[-1][1], []).append((n, len(m.group(1)) if m else -1, m.end(1) if m else -1))
    ins = esc = False
    for i, c in enumerate(code):
        if ins:
            if esc: esc = False
            elif c == '\\': esc = True
            elif c == '"': ins = False
        elif c == '"': ins = True
        elif c in '{[':
            nid += 1; stack.append((c, nid))
            if c == '{' and not code[i + 1:].strip(): block.add(nid)
        elif c in '}]': stack.pop()
colon = sorted(k[0] for o, ks in keys.items() if o in block
               and (len({k[2] for k in ks}) != 1 or min(k[1] for k in ks) != 1) for k in ks)
print(' '.join(map(str, bad)) + '|' + ' '.join(map(str, colon)))
PY
) || bad="(the check did not run)|"
[ "$bad" = "|" ] && [ -s "$W/e1/ccex_example.jsonc" ] \
    && ok "X5 the template has no trailing comments, comment lines on 4-space indents, the ':' of each object in one column" \
    || no "X5 trailing or misindented comments on lines: ${bad%%|*}; ':' not aligned on lines: ${bad#*|}"

( cd "$W" && "$MPIRUN" -n 1 "$BIN_DEV" -e -o e2/sub > e2.log 2>&1 < /dev/null )
if cmp -s "$W/e1/ccex_example.jsonc" "$W/e2/sub/ccex_example.jsonc" && [ "$(ls "$W/e2/sub")" = "ccex_example.jsonc" ]
then ok "X2 -e -o DIR writes the same file into DIR"; else no "X2 -e -o DIR (see $W/e2.log)"; fi

if python3 "$HERE/example_keys.py" "$W/e1/ccex_example.jsonc" "$REPO" | sed 's/^/    /'; [ "${PIPESTATUS[0]}" -eq 0 ]
then ok "X3 the example names every input key the readers accept"; else no "X3 the example misses input keys (above)"; fi

# A94 (review 20261016 R6): the template's auto_convergence lines, uncommented, are accepted by the runner (every choice of
# the alternative lines; example_autoconv.py)
if python3 "$HERE/example_autoconv.py" "$W/e1/ccex_example.jsonc" "$REPO" | sed 's/^/    /'; [ "${PIPESTATUS[0]}" -eq 0 ]
then ok "X6 the template's auto_convergence lines are accepted by the runner"; else no "X6 the runner refuses the template's auto_convergence lines (above)"; fi

# O1-O5 on the template with its bath (run.jsonc of X1)
E="$W/e1"
cp "$E/run.jsonc" "$E/plain.json"
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f plain.json > o1.log 2>&1 < /dev/null )
[ -s "$E/plain_noDiv" ] && [ -s "$E/plain_wiDiv" ] && ok "O1 no outfile, no -o: <input name> in the working directory" \
    || no "O1 no plain_noDiv/_wiDiv in $E (see $E/o1.log)"

( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f plain.json -o "$W/o2/a/b" > o2.log 2>&1 < /dev/null )
same_curves "$E/plain" "$W/o2/a/b/plain" && ok "O2 no outfile, -o DIR (made): DIR/<input name>, the same curves" \
    || no "O2 -o DIR without outfile (see $E/o2.log)"

python3 - "$E/plain.json" "$E/rel.json" "$E/abs.json" "$E/plan.json" "$W/o4/abs" "$E/path.json" "$E/pathb.json" <<'PY'
import re, sys
src, rel, absolute, plan, abspath, path, pathb = sys.argv[1:8]
text = open(src).read()
def put(text, line):   # a key line before the active "pulse" block, which ends the active part
    return re.sub(r'(?m)^    "pulse" *: \{', line + '\n    "pulse": {', text, count=1)
open(rel, "w").write(put(text, '    "output": {"file": "named"},'))
open(absolute, "w").write(put(text, f'    "output": {{"file": "{abspath}"}},'))
open(plan, "w").write(put(text, '    "output": {"cluster_only": true},'))
open(path, "w").write(put(text, '    "output": {"path": "pp/q"},'))     # A65: the outputs' folder in the input
open(pathb, "w").write(put(text, '    "output": {"path": "ppb"},'))
PY
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f rel.json -o "$W/o3" > o3.log 2>&1 < /dev/null )
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -o "$W/o3b" -f rel.json > o3b.log 2>&1 < /dev/null )
same_curves "$E/plain" "$W/o3/named" && same_curves "$E/plain" "$W/o3b/named" \
    && ok "O3 a relative outfile with -o DIR: DIR/<outfile>, -o before or after -f" || no "O3 (see $E/o3.log, $E/o3b.log)"

mkdir -p "$W/o4"
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f abs.json > o4.log 2>&1 < /dev/null )
same_curves "$E/plain" "$W/o4/abs" && ok "O4 an absolute outfile without -o: that path" || no "O4 (see $E/o4.log)"

( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f plan.json -o "$W/o5" > o5.log 2>&1 < /dev/null ); rc=$?
if [ $rc -eq 0 ] && [ -s "$W/o5/plan.clusters.json" ] && ! ls "$W/o5"/*Div > /dev/null 2>&1
then ok "O5 cluster_only true: DIR/<input name>.clusters.json, no curves"; else no "O5 cluster_only (exit $rc, see $E/o5.log)"; fi

# O6, O7 output.path (A65): the folder of the outputs, made when missing; -o DIR takes its place
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f path.json > o6.log 2>&1 < /dev/null )
same_curves "$E/plain" "$E/pp/q/path" && ok "O6 output.path (two levels, made): output.path/<input name>, the same curves" \
    || no "O6 output.path (see $E/o6.log)"
( cd "$E" && "$MPIRUN" -n 1 "$BIN_DEV" -f pathb.json -o "$W/o7" > o7.log 2>&1 < /dev/null )
same_curves "$E/plain" "$W/o7/pathb" && [ ! -e "$E/ppb" ] && ok "O7 -o DIR wins over output.path (not made)" \
    || no "O7 -o DIR with output.path (see $E/o7.log)"

echo "cli: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
