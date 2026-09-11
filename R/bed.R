# Shared BED I/O and interval predicates for mytools.
#
# No subcommand logic lives here -- sort, merge, intersect and subtract each own their
# own file. This is the plumbing they share: reading BED, the overlap predicate, and the
# error/exit contract from SPEC.md §7.
#
# Base R only (CLAUDE.md): no data.table, no readr, no stringr.

EXIT_OK    <- 0L   # success, including empty output
EXIT_DATA  <- 1L   # bad input data
EXIT_USAGE <- 2L   # caller error: bad flags, missing file

# --- errors -------------------------------------------------------------------
#
# Errors are signalled as conditions rather than printed-and-quit, so that unit tests
# can provoke them without killing the test process. mytools catches them at the top
# level and maps them to stderr + an exit code.

bed_error <- function(msg, status = EXIT_DATA) {
  stop(structure(
    class = c("bed_error", "error", "condition"),
    list(message = msg, call = NULL, status = status)
  ))
}

# Data problem at a known location. SPEC.md §7 wants file:line in the message.
bed_data_error <- function(source, lineno, msg) {
  bed_error(sprintf("%s:%d: %s", source, lineno, msg), EXIT_DATA)
}

usage_error <- function(msg) bed_error(msg, EXIT_USAGE)

# --- interval semantics -------------------------------------------------------
#
# BED is 0-based half-open: chr1 100 200 covers bases 100..199.
#
# Two intervals overlap iff a.start < b.end AND b.start < a.end. Note the STRICT < on
# both sides: bookended intervals (a.end == b.start) do NOT overlap. Every off-by-one in
# this project lives in this comparison (CLAUDE.md).
#
# Vectorised over all four arguments, so it works on a whole -b index at once.

bed_overlaps <- function(a_start, a_end, b_start, b_end) {
  a_start < b_end & b_start < a_end
}

# Gap between two intervals on the same chromosome, assuming a precedes b.
# Bookended intervals have gap 0, which is why they merge at -d 0 while not overlapping.
bed_gap <- function(a_end, b_start) b_start - a_end

# --- coordinates --------------------------------------------------------------
#
# Coordinates are held as doubles, which represent every integer up to 2^53 exactly --
# well past any chromosome length -- and avoid NA-on-overflow from as.integer(). They
# must never be printed by R's default formatting, which would give "1e+08".

fmt_coord <- function(x) sprintf("%.0f", x)

COORD_RE <- "^[0-9]+$"

# --- reading ------------------------------------------------------------------

# "-" means stdin (SPEC.md §2). Anything else is a file, which must exist.
bed_connection <- function(path) {
  if (identical(path, "-")) return(file("stdin", open = "r"))
  if (!file.exists(path)) usage_error(sprintf("no such file: %s", path))
  file(path, open = "r")
}

bed_source_name <- function(path) if (identical(path, "-")) "stdin" else path

# Lines that are skipped silently (SPEC.md §3): blank, and track/browser/# headers.
bed_is_skippable <- function(line) {
  if (!nzchar(line)) return(TRUE)
  if (grepl("^[[:space:]]*$", line)) return(TRUE)
  substr(line, 1, 1) == "#" ||
    startsWith(line, "track") ||
    startsWith(line, "browser")
}

# Parse one BED3-BED6 line into a record.
#
# Returns list(chrom, start, end, fields, line, lineno), where `fields` is every column
# as a character vector (so trailing columns survive untouched, whatever their number)
# and `line` is the raw input, which sort re-emits verbatim.
bed_parse_line <- function(line, source, lineno) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]

  if (length(fields) < 3) {
    bed_data_error(source, lineno,
                   sprintf("expected at least 3 tab-separated fields, got %d",
                           length(fields)))
  }

  start_raw <- fields[2]
  end_raw   <- fields[3]

  if (!grepl(COORD_RE, start_raw)) {
    bed_data_error(source, lineno, sprintf("start is not a non-negative integer (%s)",
                                           start_raw))
  }
  if (!grepl(COORD_RE, end_raw)) {
    bed_data_error(source, lineno, sprintf("end is not a non-negative integer (%s)",
                                           end_raw))
  }

  start <- as.numeric(start_raw)
  end   <- as.numeric(end_raw)

  # start == end is legal: zero-length intervals are in the fixtures on purpose.
  if (start > end) {
    bed_data_error(source, lineno,
                   sprintf("start > end (%s > %s)", fmt_coord(start), fmt_coord(end)))
  }

  list(chrom = fields[1], start = start, end = end,
       fields = fields, line = line, lineno = lineno)
}

# Streaming reader. Returns a list of closures; next_record() yields one record at a
# time or NULL at EOF, skipping comments and blanks without disturbing line numbering.
bed_reader <- function(path) {
  con    <- bed_connection(path)
  source <- bed_source_name(path)
  lineno <- 0L

  next_record <- function() {
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) return(NULL)
      lineno <<- lineno + 1L
      if (bed_is_skippable(line)) next
      return(bed_parse_line(line, source, lineno))
    }
  }

  list(next_record = next_record,
       close       = function() close(con),
       source      = source)
}

# Read every record. sort needs the whole input (SPEC.md §6), and intersect/subtract
# need all of -b resident.
bed_read_all <- function(path) {
  reader <- bed_reader(path)
  on.exit(reader$close())
  out <- list()
  repeat {
    rec <- reader$next_record()
    if (is.null(rec)) break
    out[[length(out) + 1L]] <- rec
  }
  out
}

# --- writing ------------------------------------------------------------------
#
# stdout carries data only; errors go to stderr (CLAUDE.md).

bed_write_line <- function(line) cat(line, "\n", sep = "")

# Emit chrom/start/end plus any trailing columns, tab-separated.
bed_write_record <- function(chrom, start, end, rest = character(0)) {
  bed_write_line(paste(c(chrom, fmt_coord(start), fmt_coord(end), rest), collapse = "\t"))
}

# --- argument helpers ---------------------------------------------------------

# Value of a flag taking one argument, or NULL if absent.
arg_value <- function(args, flag) {
  i <- which(args == flag)
  if (length(i) == 0L) return(NULL)
  if (length(i) > 1L) usage_error(sprintf("%s given more than once", flag))
  if (i == length(args)) usage_error(sprintf("%s requires a value", flag))
  args[i + 1L]
}

has_flag <- function(args, flag) flag %in% args

# Exactly one input may be stdin (SPEC.md §2).
check_one_stdin <- function(...) {
  paths <- unlist(list(...))
  if (sum(paths == "-") > 1L) {
    usage_error("at most one input may be '-' (stdin)")
  }
  invisible(TRUE)
}

# --- top level ----------------------------------------------------------------

# Run a subcommand, mapping bed_error conditions onto stderr plus an exit code.
# Returns the status so the caller can quit() with it.
run_guarded <- function(fn) {
  tryCatch({
    fn()
    EXIT_OK
  }, bed_error = function(e) {
    cat(conditionMessage(e), "\n", sep = "", file = stderr())
    e$status
  })
}
