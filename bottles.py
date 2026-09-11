#!/usr/bin/env python3
"""Print the full lyrics of "99 Bottles of Beer" to stdout.

The interesting part is the bottom of the song, where three things change at once:
the plural drops at 1, the count becomes the word "no more" at 0, and the verse
for the last bottle says "Take it down" rather than "Take one down".
"""

START = 99


def count(n):
    """The bottle count as it appears mid-sentence."""
    if n == 0:
        return "no more bottles"
    if n == 1:
        return "1 bottle"
    return f"{n} bottles"


def count_capitalised(n):
    """The same count where it opens a sentence: "no more" -> "No more"."""
    phrase = count(n)
    return phrase[0].upper() + phrase[1:]


lines = []

for n in range(START, 0, -1):
    lines.append(f"{count_capitalised(n)} of beer on the wall, {count(n)} of beer.")
    take = "Take it down" if n == 1 else "Take one down"
    lines.append(f"{take} and pass it around, {count(n - 1)} of beer on the wall.")
    lines.append("")

# Out of beer: the closing verse restocks the wall.
lines.append(f"{count_capitalised(0)} of beer on the wall, {count(0)} of beer.")
lines.append(f"Go to the store and buy some more, {count(START)} of beer on the wall.")

print("\n".join(lines))
