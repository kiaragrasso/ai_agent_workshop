#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on data/.
#
# Usage: ./tests/run_golden.sh
#        MYTOOLS=/path/to/other/mytools ./tests/run_golden.sh
#
# To test the harness rather than the code:  MYTOOLS=$(command -v bedtools) ./tests/run_golden.sh
# Every check()/check_stdin() case must pass under that, since it compares bedtools
# with itself. If one does not, the bug is in here. The check_rc() cases are the
# exception: they are graded against SPEC.md rather than the oracle and bedtools fails
# them by design -- it exits 1 where SPEC.md §7 says 2, and accepts two flag
# combinations SPEC.md §5 forbids. See SPEC.md §8.
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

# check_rc <name> <want_rc> -- <args...>
#   asserts mytools' own exit code, with no oracle comparison. For cases where SPEC.md
#   and bedtools genuinely disagree: SPEC.md §7 makes every caller error a 2, while
#   bedtools exits 1 (or 0, where it accepts a combination we reject). See SPEC.md §8.
#   stdout must also be empty -- errors go to stderr, stdout is data (CLAUDE.md).
check_rc() {
  local name=$1 want_rc=$2; shift 3
  "$MYTOOLS" "$@" >"$tmp/got" 2>"$tmp/got.err"; local got_rc=$?
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, wanted $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return 0
  fi
  if [[ -s $tmp/got ]]; then
    echo "FAIL $name (exit $got_rc as wanted, but stdout was not empty)"
    sed 's/^/      /' "$tmp/got" | head -3
    (( fail++ )); return 0
  fi
  echo "ok   $name"; (( pass++ )); return 0
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

# Mutually exclusive flags. Not run through check(): SPEC.md §5 makes all three
# exclusive and §7 makes that a 2, whereas bedtools rejects only -u with -v (exit 1)
# and quietly accepts -u -wa and -v -wa. Graded against the spec, not the oracle.
check_rc "intersect -u -v rejected"   2 -- intersect -u -v  -a "$DATA/a.bed" -b "$DATA/b.bed"
check_rc "intersect -u -wa rejected"  2 -- intersect -u -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
check_rc "intersect -v -wa rejected"  2 -- intersect -v -wa -a "$DATA/a.bed" -b "$DATA/b.bed"
check_rc "intersect missing -b"       2 -- intersect -a "$DATA/a.bed"
check_rc "intersect missing -a"       2 -- intersect -b "$DATA/b.bed"
check_rc "intersect -b does not exist" 2 -- intersect -a "$DATA/a.bed" -b "$DATA/nope.bed"

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
