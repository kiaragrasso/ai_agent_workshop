# Minimal unit-test harness. Base R only, so the unit tests run without bedtools and
# without any third-party package (CLAUDE.md: standard library only).
#
# Each test prints one line. A failing test names what it expected, so a broken test
# points at the bug rather than just reporting a number.

.tests <- new.env(parent = emptyenv())
.tests$pass <- 0L
.tests$fail <- 0L

test_ok <- function(desc, actual) {
  if (isTRUE(actual)) {
    cat("ok   ", desc, "\n", sep = "")
    .tests$pass <- .tests$pass + 1L
  } else {
    cat("FAIL ", desc, "\n       expected TRUE, got ", deparse(actual), "\n", sep = "")
    .tests$fail <- .tests$fail + 1L
  }
  invisible(NULL)
}

test_eq <- function(desc, actual, expected) {
  if (identical(actual, expected) ||
      (is.numeric(actual) && is.numeric(expected) &&
       length(actual) == length(expected) && all(actual == expected))) {
    cat("ok   ", desc, "\n", sep = "")
    .tests$pass <- .tests$pass + 1L
  } else {
    cat("FAIL ", desc, "\n",
        "       expected: ", paste(deparse(expected), collapse = " "), "\n",
        "       got:      ", paste(deparse(actual), collapse = " "), "\n", sep = "")
    .tests$fail <- .tests$fail + 1L
  }
  invisible(NULL)
}

# Assert that `expr` raises a bed_error, optionally with `status` and a message matching
# `pattern`. This is why R/bed.R signals conditions instead of calling quit().
test_raises <- function(desc, expr, status = NULL, pattern = NULL) {
  got <- tryCatch({
    force(expr)
    NULL
  }, bed_error = function(e) e)

  if (is.null(got)) {
    cat("FAIL ", desc, "\n       expected a bed_error, none raised\n", sep = "")
    .tests$fail <- .tests$fail + 1L
    return(invisible(NULL))
  }
  if (!is.null(status) && !identical(got$status, status)) {
    cat("FAIL ", desc, "\n       expected exit status ", status,
        ", got ", got$status, "\n", sep = "")
    .tests$fail <- .tests$fail + 1L
    return(invisible(NULL))
  }
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(got))) {
    cat("FAIL ", desc, "\n       message did not match /", pattern, "/\n",
        "       got: ", conditionMessage(got), "\n", sep = "")
    .tests$fail <- .tests$fail + 1L
    return(invisible(NULL))
  }
  cat("ok   ", desc, "\n", sep = "")
  .tests$pass <- .tests$pass + 1L
  invisible(NULL)
}

# Call at the end of a test file. Exits non-zero if anything failed.
test_summary <- function() {
  cat("---\n", .tests$pass, " passed, ", .tests$fail, " failed\n", sep = "")
  quit(save = "no", status = if (.tests$fail > 0L) 1L else 0L)
}
