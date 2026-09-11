#!/usr/bin/env Rscript
#
# Unit tests for R/bed.R -- the overlap predicate, the BED reader, and the exit-code
# contract. These run without bedtools, in milliseconds (tests/README.md).
#
# One test per edge case that needed thinking about. When a golden test fails and gets
# fixed, the unit test that would have caught it belongs here (CLAUDE.md, Testing).

here <- dirname(normalizePath(sub("^--file=", "",
                grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])))
root <- dirname(here)
source(file.path(root, "R", "bed.R"))
source(file.path(here, "assert.R"))

# ==============================================================================
# The overlap predicate
#
# BED is 0-based half-open, so overlap is a.start < b.end AND b.start < a.end, with
# STRICT < on both sides. CLAUDE.md requires this to have its own test.
# ==============================================================================

test_ok("predicate: plain overlap",            bed_overlaps(100, 200, 150, 250))
test_ok("predicate: overlap is symmetric",     bed_overlaps(150, 250, 100, 200))
test_ok("predicate: nested counts as overlap", bed_overlaps(100, 200, 120, 180))
test_ok("predicate: containing counts too",    bed_overlaps(120, 180, 100, 200))
test_ok("predicate: identical intervals",      bed_overlaps(100, 200, 100, 200))
test_ok("predicate: single shared base",       bed_overlaps(100, 101, 100, 101))

test_eq("predicate: bookended a.end == b.start does NOT overlap",
        bed_overlaps(100, 200, 200, 300), FALSE)
test_eq("predicate: bookended the other way round does not either",
        bed_overlaps(200, 300, 100, 200), FALSE)
test_eq("predicate: disjoint with a gap",
        bed_overlaps(100, 200, 300, 400), FALSE)

# Position 0: the left edge is where an off-by-one hides, because 0 is falsy in many
# languages and a '<= 0' guard looks harmless.
test_ok("predicate: interval starting at 0 overlaps one above it",
        bed_overlaps(0, 100, 50, 150))
test_eq("predicate: interval at 0 is bookended with one starting at 100",
        bed_overlaps(0, 100, 100, 200), FALSE)

# Zero-length intervals, worked out from the predicate rather than from intuition -- the
# first version of these tests asserted that a zero-length interval overlaps nothing,
# reasoning that "a.start < a.end is false". The predicate never tests that, and the
# assertion was simply wrong.
#
# A zero-length interval at p overlaps b iff b.start < p < b.end: strictly INTERIOR.
# Touching either edge of b fails, because one half of the comparison becomes p < p.
#
# This pins the PREDICATE as SPEC.md §4 defines it, and nothing more. What bedtools does
# with zero-length features inside intersect, merge and subtract is oracle behaviour that
# does not follow from this predicate -- merge expands chr1 500 500 into 499 600 -- and
# belongs in those subcommands' own tests (#6, #7, #8).
test_ok("predicate: zero-length point strictly inside an interval overlaps it",
        bed_overlaps(500, 500, 400, 600))
test_eq("predicate: zero-length point on b's start edge does not overlap",
        bed_overlaps(400, 400, 400, 600), FALSE)
test_eq("predicate: zero-length point on b's end edge does not overlap",
        bed_overlaps(600, 600, 400, 600), FALSE)
test_eq("predicate: zero-length against an identical zero-length",
        bed_overlaps(500, 500, 500, 500), FALSE)
test_eq("predicate: zero-length at position 0 against an interval starting at 0",
        bed_overlaps(0, 0, 0, 100), FALSE)

# Vectorised, because intersect tests one -a feature against a whole -b index at once.
test_eq("predicate: vectorises over b",
        bed_overlaps(100, 200, c(50, 200, 150), c(150, 300, 160)),
        c(TRUE, FALSE, TRUE))

# ==============================================================================
# Gaps -- why bookended features merge at -d 0 while not overlapping
# ==============================================================================

test_eq("gap: bookended features have gap 0",      bed_gap(200, 200), 0)
test_eq("gap: ten bases apart",                    bed_gap(200, 210), 10)
test_eq("gap: overlapping gives a negative gap",   bed_gap(200, 150), -50)

# ==============================================================================
# Parsing
# ==============================================================================

p <- function(line) bed_parse_line(line, "t.bed", 1L)

r <- p("chr1\t100\t200\tname\t40\t+")
test_eq("parse: chrom",  r$chrom, "chr1")
test_eq("parse: start",  r$start, 100)
test_eq("parse: end",    r$end,   200)
test_eq("parse: all six columns retained", length(r$fields), 6L)
test_eq("parse: raw line retained verbatim", r$line, "chr1\t100\t200\tname\t40\t+")

test_eq("parse: BED3 is accepted", length(p("chr1\t100\t200")$fields), 3L)
test_eq("parse: BED4 is accepted", length(p("chr1\t100\t200\tn")$fields), 4L)
test_eq("parse: BED5 is accepted", length(p("chr1\t100\t200\tn\t0")$fields), 5L)

test_eq("parse: zero-length interval is legal", p("chr1\t500\t500")$start, 500)
test_eq("parse: zero-length start equals end",
        local({ z <- p("chr1\t500\t500"); z$end - z$start }), 0)
test_eq("parse: position 0 parses as 0, not 1", p("chr1\t0\t100")$start, 0)
test_eq("parse: zero-length at position 0", p("chr2\t0\t0")$end, 0)

# Coordinates beyond 2^31 must not silently become NA, which is what as.integer() does.
test_eq("parse: coordinate above 2^31 survives",
        p("chr1\t0\t3000000000")$end, 3e9)

# Error paths. SPEC.md §7: data problems are exit 1, and the message names file and line.
test_raises("parse: start > end is exit 1, message names the coordinates",
            p("chr1\t500\t400"), status = EXIT_DATA, pattern = "start > end \\(500 > 400\\)")
test_raises("parse: error message carries file and line number",
            p("chr1\t500\t400"), pattern = "^t\\.bed:1:")
test_raises("parse: non-integer start is exit 1",
            p("chr1\tfoo\t200"), status = EXIT_DATA)
test_raises("parse: negative coordinate is rejected",
            p("chr1\t-5\t200"), status = EXIT_DATA)
test_raises("parse: float coordinate is rejected",
            p("chr1\t1.5\t200"), status = EXIT_DATA)
test_raises("parse: fewer than three fields is exit 1",
            p("chr1\t100"), status = EXIT_DATA)
test_raises("parse: space-separated is not tab-separated",
            p("chr1 100 200"), status = EXIT_DATA)

# ==============================================================================
# Skippable lines
# ==============================================================================

test_ok("skip: comment",        bed_is_skippable("# a comment"))
test_ok("skip: track line",     bed_is_skippable("track name=x"))
test_ok("skip: browser line",   bed_is_skippable("browser position chr1"))
test_ok("skip: blank",          bed_is_skippable(""))
test_ok("skip: whitespace only", bed_is_skippable("   "))
test_eq("skip: a data line is not skipped",
        bed_is_skippable("chr1\t100\t200"), FALSE)

# ==============================================================================
# Reading a file
# ==============================================================================

tmp <- tempfile(fileext = ".bed")
writeLines(c("# header comment",
             "track name=test",
             "",
             "chr1\t100\t200\ta\t1\t+",
             "chr1\t0\t100\tb\t2\t-",
             "chr1\t500\t500\tc\t3\t+"), tmp)

recs <- bed_read_all(tmp)
test_eq("read: comments, track and blank lines skipped", length(recs), 3L)
test_eq("read: first data record is the first non-skipped line", recs[[1]]$chrom, "chr1")
test_eq("read: input order preserved", vapply(recs, function(r) r$start, numeric(1)),
        c(100, 0, 500))
test_eq("read: line numbers count skipped lines, so errors point at the real line",
        recs[[1]]$lineno, 4L)

# Reading the same content through the streaming reader must agree with bed_read_all.
rd <- bed_reader(tmp)
streamed <- list()
repeat {
  rec <- rd$next_record()
  if (is.null(rec)) break
  streamed[[length(streamed) + 1L]] <- rec
}
rd$close()
test_eq("read: streaming reader agrees with bed_read_all",
        vapply(streamed, function(r) r$start, numeric(1)),
        vapply(recs,     function(r) r$start, numeric(1)))

empty <- tempfile(fileext = ".bed")
invisible(file.create(empty))
test_eq("read: empty file yields no records", length(bed_read_all(empty)), 0L)

test_raises("read: missing file is exit 2, not exit 1",
            bed_read_all(file.path(tempdir(), "definitely-not-here.bed")),
            status = EXIT_USAGE)

unlink(c(tmp, empty))

# ==============================================================================
# Fixture-derived cases
#
# The real fixtures, so these tests fail if data/a.bed is ever "tidied" (CLAUDE.md says
# not to). Named after the features in the tests/README.md table.
# ==============================================================================

a <- bed_read_all(file.path(root, "data", "a.bed"))
by_name <- function(recs, nm) {
  hit <- Filter(function(r) length(r$fields) >= 4 && r$fields[4] == nm, recs)
  if (length(hit) == 0L) NULL else hit[[1]]
}

a01 <- by_name(a, "a01"); a02 <- by_name(a, "a02")
a05 <- by_name(a, "a05"); a06 <- by_name(a, "a06")
a07 <- by_name(a, "a07"); a09 <- by_name(a, "a09")
a10 <- by_name(a, "a10"); a12 <- by_name(a, "a12")

test_eq("fixture: a01 starts at position 0", a01$start, 0)
test_eq("fixture: a07 is zero-length at 500", a07$end - a07$start, 0)
test_eq("fixture: a12 is zero-length at position 0", c(a12$start, a12$end), c(0, 0))

test_eq("fixture: a01 and a02 are bookended at 100, so they do not overlap",
        bed_overlaps(a01$start, a01$end, a02$start, a02$end), FALSE)
test_eq("fixture: a01 and a02 merge at -d 0 because their gap is 0",
        bed_gap(a01$end, a02$start), 0)
test_ok("fixture: a06 is nested inside a05",
        bed_overlaps(a05$start, a05$end, a06$start, a06$end) &&
          a06$start >= a05$start && a06$end <= a05$end)
test_ok("fixture: a09 and a10 are identical and overlap",
        a09$start == a10$start && a09$end == a10$end &&
          bed_overlaps(a09$start, a09$end, a10$start, a10$end))

test_summary()
