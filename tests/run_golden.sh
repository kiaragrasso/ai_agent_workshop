#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on data/.
#
# Usage: ./tests/run_golden.sh
#        MYTOOLS=/path/to/other/mytools ./tests/run_golden.sh
#
# To test the harness rather than the code:  MYTOOLS=$(command -v bedtools) ./tests/run_golden.sh
# That must report every case passing. If it does not, the bug is in here.
#
# mytools is resolved from the repo this script lives in, never from PATH. There is a
# ~/.local/bin/mytools symlink pointing at one particular clone, and several copies of
# this repo exist when agents work in parallel -- trusting bare `mytools` grades a
# different branch than the one you are sitting in, and exits 0 while doing it.
# See CLAUDE.md, Testing.

set -uo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(dirname "$here")
MYTOOLS=${MYTOOLS:-$root/mytools}
DATA=$root/data

command -v bedtools >/dev/null 2>&1 || {
  echo "bedtools not found -- the golden tests need the oracle" >&2; exit 2; }
[[ -x $MYTOOLS ]] || { echo "not executable: $MYTOOLS" >&2; exit 2; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# Compare what we captured. Row order is part of the answer, so never sort first.
report() {
  local name=$1 got_rc=$2 want_rc=$3
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return 0
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
  return 0
}

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", compares exit codes then stdout.
#   stderr is deliberately not compared: bedtools' wording is its own.
check() {
  local name=$1; shift 2
  "$MYTOOLS" "$@" >"$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools    "$@" >"$tmp/want" 2>/dev/null;      local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# check_stdin <name> <file> -- <args...>
#   same, but feeds <file> on stdin to both. Use with `-i -` or `-a -`.
check_stdin() {
  local name=$1 infile=$2; shift 3
  "$MYTOOLS" "$@" <"$infile" >"$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools    "$@" <"$infile" >"$tmp/want" 2>/dev/null;      local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# merge needs sorted input, and a.bed/b.bed are deliberately unsorted.
bedtools sort -i "$DATA/a.bed" > "$tmp/a.sorted.bed"
bedtools sort -i "$DATA/b.bed" > "$tmp/b.sorted.bed"
: > "$tmp/empty.bed"

# --- sort (#5) ----------------------------------------------------------------
check "sort a.bed"                     -- sort -i "$DATA/a.bed"
check "sort b.bed"                     -- sort -i "$DATA/b.bed"
check "sort empty"                     -- sort -i "$tmp/empty.bed"
check_stdin "sort stdin" "$DATA/a.bed" -- sort -i -

# --- merge (#6) ---------------------------------------------------------------
check "merge sorted a"                          -- merge -i "$tmp/a.sorted.bed"
check "merge -d 10 sorted a"                    -- merge -d 10 -i "$tmp/a.sorted.bed"
check "merge sorted b"                          -- merge -i "$tmp/b.sorted.bed"
check_stdin "merge stdin" "$tmp/a.sorted.bed"   -- merge -i -

# --- intersect (#7) -----------------------------------------------------------
check "intersect a b"      -- intersect     -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -u a b"   -- intersect -u  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -v a b"   -- intersect -v  -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wa a b"  -- intersect -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
check_stdin "intersect stdin" "$DATA/a.bed" -- intersect -a - -b "$DATA/b.bed"

# No swapped-argument case (-a b.bed -b a.bed). bedtools cannot use a.bed as -b: the
# zero-length interval at coordinate 0 (a12, "chr2 0 0") makes its tree build die with
#   ERROR: Received illegal bin number -1 from getBin call.
# Verified on this VM: removing that one line makes it exit 0, and that line alone as
# -b reproduces it. The other zero-length intervals (chr1 500 500, chr2 300 300) are
# fine, so it is specifically zero-length AT position 0.
# This is a bedtools crash, not behaviour worth requiring mytools to reproduce, so the
# asymmetry check is left out rather than encoded.

# --- subtract (#8) ------------------------------------------------------------
check "subtract a b"       -- subtract -a "$DATA/a.bed" -b "$DATA/b.bed"
check_stdin "subtract stdin" "$DATA/a.bed" -- subtract -a - -b "$DATA/b.bed"
# Swapped arguments omitted here too, same bedtools crash -- see the note above.

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
