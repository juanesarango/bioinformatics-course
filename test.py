"""File to run the code that will be input into the
embedded code editor of the Coursera Online Classroom.
Inputs come via `stdin` from `text_input.txt`.
Outputs must be `stdout` by printing them.

To run it from bash:

    $ `python test.py < test_input.txt`
"""

import sys
import pandas as pd


BASE_PATTERN = 'ACGT'


def hamming_distance(text_p, text_q):
    """Returns the hamming distance between 2 segments.
    This distance is the total sum of mismatches for eachj position
    `i` if `p[i] != q[i]`.
    """
    if len(text_p) != len(text_q):
        return -1
    return sum([text_p[i] != text_q[i] for i in range(len(text_p))])


def neighbors(pattern, d):
    """Generates the d-neighborhood, the set of all k-mers whose Hamming
    distance from Pattern does not exceed d. Version 3.0
    """
    if not pattern or not d:
        return pattern
    if len(pattern) == 1:
        return list(BASE_PATTERN)
    neighborhood = []
    sufix_neighbors = neighbors(pattern[1:], d)
    for neighbor in sufix_neighbors:
        if hamming_distance(pattern[1:], neighbor) < d:
            for base in BASE_PATTERN:
                neighborhood.append(base + neighbor)
        else:
            neighborhood.append(pattern[0] + neighbor)
    return neighborhood


pattern, d = sys.stdin.read().splitlines()
output = neighbors(pattern, d)
for o in output:
    print(o)
