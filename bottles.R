#!/usr/bin/env Rscript

# "99 Bottles of Beer" — warm-up exercise.
# The bottom of the song is all edge cases: the plural stops being plural at 1,
# and the count becomes a word ("no more") at 0.

# How many bottles, as it appears mid-sentence: 2 -> "2 bottles", 1 -> "1 bottle",
# 0 -> "no more bottles".
bottles <- function(n) {
  if (n == 0) {
    "no more bottles"
  } else if (n == 1) {
    "1 bottle"
  } else {
    paste0(n, " bottles")
  }
}

# Same thing at the start of a sentence, so "no more" needs a capital.
bottles_sentence_start <- function(n) {
  s <- bottles(n)
  sub("^no", "No", s)
}

start <- 99

for (n in start:1) {
  cat(sprintf(
    "%s of beer on the wall, %s of beer.\n",
    bottles_sentence_start(n), bottles(n)
  ))
  # "Take it down" when there is only one left, "Take one down" otherwise.
  take <- if (n == 1) "Take it down" else "Take one down"
  cat(sprintf(
    "%s and pass it around, %s of beer on the wall.\n\n",
    take, bottles(n - 1)
  ))
}

# Final verse: out of beer, go restock.
cat(sprintf(
  "%s of beer on the wall, %s of beer.\n",
  bottles_sentence_start(0), bottles(0)
))
cat(sprintf(
  "Go to the store and buy some more, %s of beer on the wall.\n",
  bottles(start)
))
