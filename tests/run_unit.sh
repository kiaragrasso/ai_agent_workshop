#!/usr/bin/env bash
# Unit tests: every tests/test_*.R, run with Rscript.
#
# These need no bedtools and no network -- that is the point of them. They run in
# milliseconds and, when one fails, it names the behaviour rather than just telling you
# that output differs. See tests/README.md, "Two kinds of test".
#
# Usage: ./tests/run_unit.sh

set -uo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
failed=0

shopt -s nullglob
files=("$here"/test_*.R)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "no test_*.R files found in $here" >&2
  exit 2
fi

for f in "${files[@]}"; do
  echo "== $(basename "$f")"
  Rscript "$f" || failed=1
done

echo "======"
if [[ $failed -ne 0 ]]; then
  echo "unit tests FAILED"
  exit 1
fi
echo "unit tests passed"
