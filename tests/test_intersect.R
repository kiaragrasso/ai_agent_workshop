#!/usr/bin/env Rscript
#
# Unit tests for R/intersect.R. No bedtools, no network (tests/README.md) -- where a
# case encodes surprising oracle behaviour, the comment says so and records the
# bedtools run it came from.
#
# One test per edge case that needed thinking about: bookended, nested, identical,
# position 0, zero-length, a chromosome on one side only, and the -u/-wa multiplicity
# split. Coordinates are taken from data/a.bed and data/b.bed so a failure here points
# at a named fixture row.

here <- dirname(normalizePath(sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
root <- dirname(here)
source(file.path(root, "R", "bed.R"))
source(file.path(root, "R", "intersect.R"))
source(file.path(here, "assert.R"))

# --- helpers ------------------------------------------------------------------

bed_file <- function(lines) {
  path <- tempfile(fileext = ".bed")
  writeLines(lines, path)
  path
}

# Run the subcommand end to end and return its stdout as a character vector.
run_intersect <- function(a_lines, b_lines, flags = character(0)) {
  a <- bed_file(a_lines)
  b <- bed_file(b_lines)
  on.exit(unlink(c(a, b)))
  capture.output(cmd_main(c(flags, "-a", a, "-b", b)))
}

# The per-chromosome index entry, for testing find_overlaps() directly.
index_for <- function(b_lines, chrom) {
  b <- bed_file(b_lines)
  on.exit(unlink(b))
  idx <- build_b_index(b)
  if (!exists(chrom, envir = idx, inherits = FALSE)) return(NULL)
  get(chrom, envir = idx, inherits = FALSE)
}

hits <- function(b_lines, chrom, a_start, a_end) {
  length(find_overlaps(index_for(b_lines, chrom), a_start, a_end))
}

# ==============================================================================
# The overlap predicate, as intersect uses it
#
# R/bed.R has the predicate's own tests. These are the four shapes the issue calls
# out, exercised through the index so a bad binary-search bound fails here too.
# ==============================================================================

B <- "chr1\t100\t200\tb\t0\t+"

test_eq("predicate via index: plain overlap",   hits(B, "chr1", 150, 250), 1L)
test_eq("predicate via index: nested",          hits(B, "chr1", 120, 180), 1L)
test_eq("predicate via index: containing",      hits(B, "chr1",  50, 250), 1L)
test_eq("predicate via index: identical",       hits(B, "chr1", 100, 200), 1L)
test_eq("predicate via index: disjoint",        hits(B, "chr1", 300, 400), 0L)

# Bookended: a01 ends at 100 where a02 begins. Strict < on both sides means they do
# not overlap (CLAUDE.md, interval semantics) -- the single commonest off-by-one.
test_eq("bookended a.end == b.start does not overlap", hits(B, "chr1", 50, 100), 0L)
test_eq("bookended b.end == a.start does not overlap", hits(B, "chr1", 200, 300), 0L)

# ==============================================================================
# Fixture shapes
# ==============================================================================

# Nested: b08 (chr1 750 760) sits wholly inside a09 (chr1 700 800). The intersected
# region is the inner interval, and it must survive the cummax_end lower bound --
# a bound taken on b.start alone would skip a nested b whose neighbours end earlier.
test_eq("nested b inside a: one hit, region is the inner interval",
        run_intersect("chr1\t700\t800\ta09\t60\t+",
                      c("chr1\t500\t600\tearlier\t0\t+",
                        "chr1\t750\t760\tb08\t16\t+")),
        "chr1\t750\t760\ta09\t60\t+")

# Identical coordinates: a09 and a10 are the same interval twice. Both are reported,
# independently -- there is no dedup across -a features.
test_eq("identical -a features are each reported",
        run_intersect(c("chr1\t700\t800\ta09\t60\t+",
                        "chr1\t700\t800\ta10\t60\t+"),
                      "chr1\t750\t760\tb08\t16\t+"),
        c("chr1\t750\t760\ta09\t60\t+",
          "chr1\t750\t760\ta10\t60\t+"))

# Position 0: a01/b01 both start at 0, and b16 is chrX 0 5. The left edge is where an
# off-by-one hides, because a `-1` or a falsy-zero guard looks harmless there.
test_eq("interval at position 0 intersects correctly",
        run_intersect("chr1\t0\t100\ta01\t10\t+",
                      "chr1\t0\t50\tb01\t11\t+"),
        "chr1\t0\t50\ta01\t10\t+")
test_eq("no negative coordinate is ever emitted at the left edge",
        run_intersect("chrX\t10\t20\ta20\t85\t+",
                      "chrX\t0\t5\tb16\t24\t+"),
        character(0))

# ==============================================================================
# Zero-length intervals -- ORACLE BEHAVIOUR
#
# These encode what real bedtools prints, not what the half-open model predicts. The
# plain predicate is strict on both sides, so a zero-length interval would never
# overlap anything; bedtools inflates it to [start-1, end+1] first. Established by
# running the oracle directly (see the probe table in R/intersect.R). CLAUDE.md:
# encode bedtools' behaviour, with a comment saying why -- do not "fix" these.
# ==============================================================================

test_eq("zero-length -b (b02, chr1 100 100) is inflated: hits a01, region 99..100",
        run_intersect("chr1\t0\t100\ta01\t10\t+",
                      "chr1\t100\t100\tb02\t0\t-"),
        "chr1\t99\t100\ta01\t10\t+")

test_eq("the same b02 also hits a02 from the other side, region 100..101",
        run_intersect("chr1\t100\t200\ta02\t20\t-",
                      "chr1\t100\t100\tb02\t0\t-"),
        "chr1\t100\t101\ta02\t20\t-")

test_eq("zero-length -b two bases clear of -a does not hit",
        run_intersect("chr1\t100\t200\ta02\t20\t-",
                      "chr1\t99\t99\tb\t0\t-"),
        character(0))

# a07 (chr1 500 500) against b07 (chr1 500 500): both zero-length. bedtools prints
# a07's own coordinates -- the region is clipped to -a, which for a zero-length -a is
# the point itself, whichever side of it the -b sits.
test_eq("zero-length -a and -b at the same point: prints -a's coordinates",
        run_intersect("chr1\t500\t500\ta07\t0\t+",
                      "chr1\t500\t500\tb07\t0\t-"),
        "chr1\t500\t500\ta07\t0\t+")

test_eq("zero-length -a against the -b that merely touches it (a08's side, 500..600)",
        run_intersect("chr1\t500\t600\ta08\t35\t-",
                      "chr1\t500\t500\tb07\t0\t-"),
        "chr1\t500\t501\ta08\t35\t-")

# a12 is "chr2 0 0" -- zero-length AT coordinate 0, the row that makes bedtools' own
# tree build die when a.bed is used as -b (see tests/run_golden.sh). As -a it is fine
# and must still report its hit against b10 (chr2 0 10).
test_eq("zero-length -a at coordinate 0 (a12) still hits, and prints 0 0",
        run_intersect("chr2\t0\t0\ta12\t0\t+",
                      "chr2\t0\t10\tb10\t18\t-"),
        "chr2\t0\t0\ta12\t0\t+")

test_eq("zero-length -a at 0 does not reach a -b starting at 1",
        run_intersect("chr2\t0\t0\ta12\t0\t+",
                      "chr2\t1\t10\tb\t0\t-"),
        character(0))

# a16 (chr2 300 300) sits at the end of b12 (chr2 200 300): inflated, it reaches back
# one base and hits.
test_eq("zero-length -a at the far end of a -b feature (a16/b12) hits",
        run_intersect("chr2\t300\t300\ta16\t0\t+",
                      "chr2\t200\t300\tb12\t20\t+"),
        "chr2\t300\t300\ta16\t0\t+")

# ==============================================================================
# Chromosomes present on one side only
# ==============================================================================

# b19 is "chr3 100 200" and no -a feature is on chr3. Nothing may be emitted for it
# and nothing may crash -- the index simply holds a chromosome nobody looks up.
test_eq("a chromosome in -b but not -a emits nothing",
        run_intersect("chr1\t100\t200\ta\t0\t+",
                      c("chr1\t100\t200\tb\t0\t+",
                        "chr3\t100\t200\tb19\t27\t+")),
        "chr1\t100\t200\ta\t0\t+")

# The reverse: a chromosome in -a with no entry in the index at all. Every feature on
# it is a -v hit, and find_overlaps must handle the missing entry rather than error.
test_eq("a chromosome in -a but not -b yields no hits",
        run_intersect("chrY\t100\t200\ta\t0\t+",
                      "chr1\t100\t200\tb\t0\t+"),
        character(0))
test_eq("...and every feature on it is a -v hit",
        run_intersect(c("chrY\t100\t200\ta\t0\t+",
                        "chrY\t300\t400\ta2\t0\t+"),
                      "chr1\t100\t200\tb\t0\t+", "-v"),
        c("chrY\t100\t200\ta\t0\t+",
          "chrY\t300\t400\ta2\t0\t+"))
test_eq("find_overlaps on an absent chromosome returns nothing",
        hits("chr1\t100\t200\tb\t0\t+", "chrY", 100, 200), 0L)

# ==============================================================================
# Multiplicity: -u once, -wa once per overlapping -b
# ==============================================================================

THREE_B <- c("chr1\t110\t120\tb1\t0\t+",
             "chr1\t130\t140\tb2\t0\t+",
             "chr1\t150\t160\tb3\t0\t+")
ONE_A   <- "chr1\t100\t200\ta\t0\t+"

test_eq("-u emits the -a feature once though three -b features overlap it",
        run_intersect(ONE_A, THREE_B, "-u"), ONE_A)

test_eq("-wa emits it three times in that same case",
        run_intersect(ONE_A, THREE_B, "-wa"), rep(ONE_A, 3L))

test_eq("default emits the three intersected regions",
        run_intersect(ONE_A, THREE_B),
        c("chr1\t110\t120\ta\t0\t+",
          "chr1\t130\t140\ta\t0\t+",
          "chr1\t150\t160\ta\t0\t+"))

test_eq("-v emits nothing when there is an overlap",
        run_intersect(ONE_A, THREE_B, "-v"), character(0))

# ==============================================================================
# Hit order and input order
# ==============================================================================

# Hits come out in -b FILE order, not coordinate order. bedtools walks its bin tree in
# insertion order and every fixture interval lands in one bin, so b.bed's order is the
# answer. a02 is the row that shows it: b03 (start 180) is printed before b02 (start
# 100) because b03 is the first line of b.bed.
test_eq("hits follow -b file order, not start order",
        run_intersect("chr1\t100\t200\ta02\t20\t-",
                      c("chr1\t180\t220\tb03\t12\t+",
                        "chr1\t0\t50\tb01\t11\t+",
                        "chr1\t100\t100\tb02\t0\t-")),
        c("chr1\t180\t200\ta02\t20\t-",
          "chr1\t100\t101\ta02\t20\t-"))

# -a is never sorted for us and must never be sorted by us (SPEC.md §5).
test_eq("-a input order is preserved, unsorted",
        run_intersect(c("chr1\t300\t400\tsecond\t0\t+",
                        "chr1\t100\t200\tfirst\t0\t+"),
                      "chr1\t0\t1000\tb\t0\t+", "-wa"),
        c("chr1\t300\t400\tsecond\t0\t+",
          "chr1\t100\t200\tfirst\t0\t+"))

# ==============================================================================
# Trailing columns and BED3
# ==============================================================================

test_eq("default carries -a's trailing columns onto the intersected region",
        run_intersect("chr1\t100\t200\tname\t42\t-\textra",
                      "chr1\t150\t250\tb\t0\t+"),
        "chr1\t150\t200\tname\t42\t-\textra")

test_eq("BED3 -a stays BED3",
        run_intersect("chr1\t100\t200", "chr1\t150\t250"),
        "chr1\t150\t200")

# ==============================================================================
# Argument errors -- SPEC.md §7 makes every one of these a 2
# ==============================================================================

test_raises("-u and -v together are rejected",
            parse_args(c("-u", "-v", "-a", "x", "-b", "y")),
            status = EXIT_USAGE, pattern = "mutually exclusive")
test_raises("-u and -wa together are rejected too",
            parse_args(c("-u", "-wa", "-a", "x", "-b", "y")),
            status = EXIT_USAGE, pattern = "mutually exclusive")
test_raises("-v and -wa together are rejected too",
            parse_args(c("-v", "-wa", "-a", "x", "-b", "y")),
            status = EXIT_USAGE, pattern = "mutually exclusive")
test_raises("missing -a", parse_args(c("-b", "y")),
            status = EXIT_USAGE, pattern = "-a is required")
test_raises("missing -b", parse_args(c("-a", "x")),
            status = EXIT_USAGE, pattern = "-b is required")
test_raises("-a without a value", parse_args(c("-a")),
            status = EXIT_USAGE, pattern = "requires a value")
test_raises("-b may not be stdin: it is read fully into memory",
            parse_args(c("-a", "x", "-b", "-")),
            status = EXIT_USAGE, pattern = "must be a file")
test_raises("unknown flag", parse_args(c("-q", "-a", "x", "-b", "y")),
            status = EXIT_USAGE, pattern = "unrecognised")
test_raises("a nonexistent -b file is a usage error, not a data error",
            build_b_index(file.path(tempdir(), "definitely-not-here.bed")),
            status = EXIT_USAGE, pattern = "no such file")

# Flags parse to the mode they name.
test_eq("default mode", parse_args(c("-a", "x", "-b", "y"))$mode, "default")
test_eq("-u mode",  parse_args(c("-u",  "-a", "x", "-b", "y"))$mode, "u")
test_eq("-v mode",  parse_args(c("-v",  "-a", "x", "-b", "y"))$mode, "v")
test_eq("-wa mode", parse_args(c("-wa", "-a", "x", "-b", "y"))$mode, "wa")

test_summary()
