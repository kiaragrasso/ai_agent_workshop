# mytools intersect -- report -a features that overlap -b features.
#
#   mytools intersect -a <file|-> -b <file> [-u | -v | -wa]
#
# -b is read fully into memory and indexed; -a streams (SPEC.md §6).
#
# Base R only (CLAUDE.md). Real bedtools is the oracle: where its behaviour is
# surprising, the comments below record what it does and how that was established.

# --- zero-length intervals ----------------------------------------------------
#
# BED is 0-based half-open and the overlap predicate is strict on both sides
# (R/bed.R, bed_overlaps), so a zero-length interval [p,p] can never overlap
# anything under the plain rule. bedtools nevertheless reports hits for them, and
# the fixtures carry five (a07, a12, a16, b02, b07) precisely so we find out.
#
# What bedtools actually does, established by probing it directly:
#
#   A=[100,200] B=[100,100] -> c 100 101      A=[500,500] B=[500,600] -> c 500 500
#   A=[100,200] B=[200,200] -> c 199 200      A=[500,500] B=[400,500] -> c 500 500
#   A=[100,200] B=[150,150] -> c 149 151      A=[500,500] B=[501,501] -> c 500 500
#   A=[100,200] B=[99,99]   -> (none)         A=[500,500] B=[498,498] -> (none)
#   A=[0,0]     B=[0,10]    -> c 0 0          A=[500,500] B=[502,502] -> (none)
#
# Both columns are explained by one rule: a zero-length interval is inflated by one
# base on each side, to [start-1, end+1], before the overlap test. That is oracle
# behaviour, not a derivation from the half-open model -- do not "simplify" it back
# to the plain predicate (CLAUDE.md: encode what bedtools does, not what it should).
#
# start - 1 can be -1 for a zero-length interval at coordinate 0 (a12, "chr2 0 0").
# That is fine: the inflated bounds are only ever compared, never printed.

inflate_start <- function(start, end) ifelse(start == end, start - 1, start)
inflate_end   <- function(start, end) ifelse(start == end, end   + 1, end)

# --- the -b index -------------------------------------------------------------
#
# Grouped by chromosome, sorted by (inflated) start, so a query is two binary
# searches rather than a scan (SPEC.md §6). A scan per -a feature is quadratic: it
# passes every golden test on these 22-line fixtures and falls over at 10^6.
#
# Sorting by start alone only bounds the search from the right. cummax_end -- the
# running maximum of the ends -- is non-decreasing and so is binary-searchable too,
# and it bounds it from the left: every b before the first cummax_end > a.start has
# every end at or below a.start and cannot overlap. Nested intervals make a plain
# lower bound on end wrong; this handles them.
build_b_index <- function(path) {
  recs <- bed_read_all(path)

  chrom <- vapply(recs, function(r) r$chrom, character(1))
  start <- vapply(recs, function(r) r$start, numeric(1))
  end   <- vapply(recs, function(r) r$end,   numeric(1))

  eff_start <- inflate_start(start, end)
  eff_end   <- inflate_end(start, end)

  index <- new.env(parent = emptyenv())
  for (chr in unique(chrom)) {
    on_chr <- which(chrom == chr)
    # Ties broken by file position, so hit order stays deterministic.
    ord <- on_chr[order(eff_start[on_chr], on_chr)]
    assign(chr, list(
      eff_start  = eff_start[ord],
      eff_end    = eff_end[ord],
      cummax_end = cummax(eff_end[ord]),
      # Position in b as written, which is the order hits are emitted in -- see
      # find_overlaps().
      file_pos   = ord
    ), envir = index)
  }
  index
}

# Positions within the chromosome's index arrays, in the order bedtools emits them,
# of every -b feature overlapping one -a feature. Empty integer(0) when there are
# none.
#
# Hits come back in *b-file* order, not start order. bedtools stores -b in a UCSC
# bin tree and walks each bin in insertion order, which for these fixtures (every
# interval below 128kb, so all in one bin) is exactly the order of b.bed. It shows
# up on a02: b03 (start 180) is emitted before b02 (start 100) because b03 comes
# first in the file. Sorting hits by coordinate gets that pair backwards.
find_overlaps <- function(b, a_start, a_end) {
  if (is.null(b)) return(integer(0))

  a_eff_start <- inflate_start(a_start, a_end)
  a_eff_end   <- inflate_end(a_start, a_end)

  # Right bound: need b.start < a.end. Coordinates are whole numbers held as
  # doubles, so that is b.start <= a.end - 1, which findInterval answers directly.
  hi <- findInterval(a_eff_end - 1, b$eff_start)
  if (hi == 0L) return(integer(0))

  # Left bound: need b.end > a.start. Everything up to and including `lo` has
  # cummax_end <= a.start and is excluded.
  lo <- findInterval(a_eff_start, b$cummax_end)
  if (lo >= hi) return(integer(0))

  cand <- (lo + 1L):hi
  hit  <- cand[bed_overlaps(a_eff_start, a_eff_end,
                            b$eff_start[cand], b$eff_end[cand])]
  if (length(hit) == 0L) return(integer(0))

  # Back into b-file order. The index arrays are in start order, so sorting the
  # positions themselves would emit coordinate order -- sort on file_pos instead.
  hit[order(b$file_pos[hit])]
}

# --- output -------------------------------------------------------------------

# The intersected region for default mode, carrying -a's trailing columns.
#
# Clipped to -a, using -a's *original* coordinates against -b's inflated ones. A
# zero-length -a feature therefore prints its own coordinates unchanged, which is
# what bedtools does: a07 prints "chr1 500 500" whichever of b07's bounds it met.
write_intersection <- function(a, b_eff_start, b_eff_end) {
  start <- max(a$start, b_eff_start)
  end   <- min(a$end,   b_eff_end)
  bed_write_record(a$chrom, start, end, a$fields[-(1:3)])
}

# --- argument parsing ---------------------------------------------------------

parse_args <- function(args) {
  opts <- list(a = NULL, b = NULL, mode = "default")
  modes <- character(0)

  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "-a" || arg == "-b") {
      if (i == length(args)) usage_error(sprintf("%s requires a value", arg))
      key <- substring(arg, 2L)
      if (!is.null(opts[[key]])) usage_error(sprintf("%s given more than once", arg))
      opts[[key]] <- args[i + 1L]
      i <- i + 2L
    } else if (arg %in% c("-u", "-v", "-wa")) {
      modes <- c(modes, arg)
      i <- i + 1L
    } else {
      usage_error(sprintf("intersect: unrecognised argument: %s", arg))
    }
  }

  # SPEC.md §5 makes all three mutually exclusive and §7 makes a caller error a 2.
  # bedtools differs -- it rejects only -u with -v, and with exit 1 -- so this case
  # cannot be graded against the oracle. See tests/run_golden.sh and SPEC.md §8.
  if (length(modes) > 1L) {
    usage_error(sprintf("intersect: %s are mutually exclusive",
                        paste(unique(modes), collapse = " and ")))
  }
  if (length(modes) == 1L) opts$mode <- substring(modes, 2L)

  if (is.null(opts$a)) usage_error("intersect: -a is required")
  if (is.null(opts$b)) usage_error("intersect: -b is required")
  # -b is read fully into memory, so it has to be seekable-ish: a real file
  # (SPEC.md §2, §6). Only -a may be stdin.
  if (identical(opts$b, "-")) usage_error("intersect: -b must be a file, not '-'")

  opts
}

# --- main ---------------------------------------------------------------------

cmd_main <- function(args) {
  opts  <- parse_args(args)
  index <- build_b_index(opts$b)

  reader <- bed_reader(opts$a)
  on.exit(reader$close())

  repeat {
    a <- reader$next_record()
    if (is.null(a)) break

    # A chromosome in -a but not in -b (chrX-style) has no entry at all: no hits,
    # no crash. The reverse (chr3, in -b only) is simply never looked up.
    b    <- if (exists(a$chrom, envir = index, inherits = FALSE)) {
              get(a$chrom, envir = index, inherits = FALSE)
            } else NULL
    hits <- find_overlaps(b, a$start, a$end)

    if (opts$mode == "v") {
      # -v: no overlap at all. Input order preserved, unsorted as -a was written.
      if (length(hits) == 0L) bed_write_line(a$line)
    } else if (length(hits) == 0L) {
      next
    } else if (opts$mode == "u") {
      bed_write_line(a$line)              # at most once, however many -b it met
    } else if (opts$mode == "wa") {
      for (h in hits) bed_write_line(a$line)   # once per overlapping -b feature
    } else {
      for (pos in hits) {
        write_intersection(a, b$eff_start[pos], b$eff_end[pos])
      }
    }
  }

  invisible(NULL)
}
