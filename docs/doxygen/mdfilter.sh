#!/bin/sh
# Doxygen input filter for the Markdown pages (docs/user, docs/methods and the tutorial
# README.md files of example/tutorial/{basic,advanced}/NN_*/).
# Only POSIX sh, sed and awk features: the CI action runs it with busybox.
#  - A link to a file that is not part of the site (docs/validation, docs/development)
#    becomes its text followed by the repository path in code, so the site has no dead
#    link while GitHub keeps the real one. A link ../../<path> from docs/ into the
#    rest of the repository is treated the same way.
#  - A link written as [text](Page.md#anchor) or [text](../user/Page.md#anchor) works on
#    GitHub, but Doxygen cannot turn it into a page link; it can resolve the bare anchor,
#    which is global. So the path and page part are dropped before Doxygen reads the file.
#  - A link to a tutorial README (../../example/tutorial/LEVEL/NN_name/README.md from docs/,
#    ../NN_name/README.md or ../../LEVEL/NN_name/README.md from another tutorial) becomes a link to its page,
#    tutorial_bNN (basic) or tutorial_aNN (advanced).
#  - An image written as ![text](../images/<dir>/file.png) is found by GitHub through
#    that path, but Doxygen copies an image into the site only when it is named bare
#    and found in IMAGE_PATH. So the directory part is dropped too.
#  - In a tutorial README: the title gets the page label tutorial_bNN or tutorial_aNN; every
#    figure is called expected.png in its calculation folder, so ![..](calc/expected.png)
#    becomes ![..](bNN_calc.png) (and ![..](expected.png) bNN.png), the names of the links in
#    docs/images/tutorials/; a line "Input: [..](calc/file.jsonc)" or "Inputs: [..](..),
#    [..](..)" is replaced by the input files themselves, one code block each.
#  - Math is written for GitHub, $...$ inline and $$...$$ displayed (on one line or over
#    several). Doxygen reads \f$...\f$ and \f[...\f], so the delimiters are rewritten,
#    outside code blocks and code spans only. A $ that does not open math (no closing $
#    on the line, a space after it) is left alone, as GitHub does.
#  - A [TOC] line goes under the title of every page, for the table of contents beside it.
nn=""
case "$1" in
    */example/tutorial/basic/[0-9][0-9]_*/README.md)
        nn=$(basename "$(dirname "$1")"); nn=b${nn%%_*} ;;
    */example/tutorial/advanced/[0-9][0-9]_*/README.md)
        nn=$(basename "$(dirname "$1")"); nn=a${nn%%_*} ;;
esac
dir=$(dirname "$1")
tutorial() {
    if [ -n "$nn" ]; then
        sed -E \
          -e "1s/^# .*/& {#tutorial_$nn}/" \
          -e "s#\]\(([A-Za-z0-9_]+)/expected\.png\)#](${nn}_\1.png)#g" \
          -e "s#\]\(expected\.png\)#](${nn}.png)#g" \
          -e "s#\]\(\.\./([0-9][0-9])_[A-Za-z0-9_]*/(README\.md)?\)#](@ref tutorial_${nn%??}\1)#g" \
          -e 's#\]\(\.\./\.\./basic/([0-9][0-9])_[A-Za-z0-9_]*/(README\.md)?\)#](@ref tutorial_b\1)#g' \
          -e 's#\]\(\.\./\.\./advanced/([0-9][0-9])_[A-Za-z0-9_]*/(README\.md)?\)#](@ref tutorial_a\1)#g' \
          -e 's#\]\((\.\./)+docs/(user|methods)/([A-Za-z0-9_-]+\.md)?\##](\##g' \
          "$1"
    else
        cat "$1"
    fi
}
tutorial "$1" | sed -E -e '1{/^# /a\
[TOC]
}' \
  -e 's#\[([^]]*)\]\((\.\./)+(docs/)?(validation|development)/([^)]*)\)#\1 (`docs/\4/\5`)#g' \
  -e 's#\]\((\.\./)+example/tutorial/basic/([0-9][0-9])_[A-Za-z0-9_]*/(README\.md)?\)#](@ref tutorial_b\2)#g' \
  -e 's#\]\((\.\./)+example/tutorial/advanced/([0-9][0-9])_[A-Za-z0-9_]*/(README\.md)?\)#](@ref tutorial_a\2)#g' \
  -e 's#\[([^]]*)\]\(\.\./\.\./([^.)][^)]*)\)#\1 (`\2`)#g' \
  -e 's#\]\((\.\./(user|methods)/)?[A-Za-z0-9_-]+\.md\##](\##g' \
  -e 's#\]\(\.\./images/([A-Za-z0-9_]+/)*#](#g' \
  | awk -v dir="$dir" '
BEGIN { fence = ""; disp = 0 }
function conv(s,    out, i, n, c, j, k, run, rl, body) {
    out = ""; i = 1; n = length(s)
    while (i <= n) {
        c = substr(s, i, 1)
        if (c == "`") {                          # code span: copy through the matching run
            j = i; while (j <= n && substr(s, j, 1) == "`") j++
            run = substr(s, i, j - i); rl = j - i
            k = index(substr(s, j), run)
            if (k > 0) { out = out substr(s, i, rl + k - 1 + rl); i = j + k - 1 + rl }
            else { out = out run; i = j }
            continue
        }
        if (disp) {                              # inside a displayed formula
            if (substr(s, i, 2) == "$$") { out = out "\\f]"; disp = 0; i += 2 }
            else { out = out c; i++ }
            continue
        }
        if (c == "\\" && substr(s, i + 1, 1) == "$") { out = out "\\$"; i += 2; continue }
        if (substr(s, i, 2) == "$$") {
            k = index(substr(s, i + 2), "$$")
            if (k > 0) { out = out "\\f[" substr(s, i + 2, k - 1) "\\f]"; i = i + 2 + k - 1 + 2 }
            else { out = out "\\f["; disp = 1; i += 2 }
            continue
        }
        if (c == "$" && i < n && substr(s, i + 1, 1) != " ") {
            j = i + 1; k = 0
            while (j <= n) {
                if (substr(s, j, 1) == "\\") { j += 2; continue }
                if (substr(s, j, 1) == "$") {
                    if (substr(s, j - 1, 1) != " " && substr(s, j + 1, 1) !~ /[0-9]/) k = j
                    break
                }
                j++
            }
            if (k > 0) { body = substr(s, i + 1, k - i - 1); out = out "\\f$" body "\\f$"; i = k + 1; continue }
        }
        out = out c; i++
    }
    return out
}
function inputs(s,    p, f, name, l, n) {   # Input: [..](calc/file.jsonc), ... -> the files
    while (match(s, /\]\([^)]*\.jsonc\)/)) {
        p = substr(s, RSTART + 2, RLENGTH - 3); s = substr(s, RSTART + RLENGTH)
        f = dir "/" p; n = 0
        print "`" p "`:"; print ""; print "```js"
        while ((getline l < f) > 0) { print l; n++ }
        close(f)
        if (n == 0) print "mdfilter.sh: input file " f " not found or empty" > "/dev/stderr"
        print "```"; print ""
    }
}
{
    line = $0
    if (fence == "" && line ~ /^Inputs?: .*\.jsonc\)/) { inputs(line); next }
    if (fence == "" && (line ~ /^[ ]*```[^`]*$/ || line ~ /^[ ]*~~~/)) { fence = (line ~ /^[ ]*```/) ? "```" : "~~~"; print; next }
    if (fence != "") { if (index(line, fence) > 0 && line ~ /^[ ]*(```|~~~)[ ]*$/) fence = ""; print; next }
    print conv(line)
}'
