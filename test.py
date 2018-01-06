"""File to run the code that will be input into the
embedded code editor of the Coursera Online Classroom.
Inputs come via `stdin` from `text_input.txt`.
Outputs must be `stdout` by printing them.

To run it from bash:

    $ `python test.py < test_input.txt`
"""


import sys 

BASE_PATTERN = 'ACGT'

def number_to_symbol(number):
    """Convert weights to nucleotides 0123 -> ACGT."""
    return BASE_PATTERN[number]


def number_to_pattern(number, k):
    """Recursively the pattern can be calculated by taking dividing
    the number by the weight of that position. The remainer is used to
    get the correct symbol, and the integer division goes recursively
    to find the symbol of the next position.
    """
    if k == 1:
        return number_to_symbol(number)
    prefix_number = int(number/4)
    remainer = number - prefix_number * 4
    symbol = number_to_symbol(remainer)
    prefix_pattern = number_to_pattern(prefix_number, k - 1)
    return prefix_pattern + symbol

number, k = sys.stdin.read().splitlines()

print(number_to_pattern(int(number), int(k)))