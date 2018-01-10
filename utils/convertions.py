BASE_PATTERN = 'ACGT'

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


def get_complementary(text):
    """Get complementary nucleotide chain."""
    return ''.join(reversed([COMPLEMENT[n] for n in text]))


def symbol_to_number(symbol):
    """Convert to nucleotides to weights ACGT -> 0123."""
    return BASE_PATTERN.find(symbol)


def number_to_symbol(number):
    """Convert weights to nucleotides 0123 -> ACGT."""
    return BASE_PATTERN[number]


# Normal Version
def pattern_to_number_1(pattern):
    """Having a 4^k array with every k-mer combination sorted lexicographically,
    the functions returns the position of a given pattern in the array."""
    position_weigths = [4**i for i in reversed(range(0, len(pattern)))]
    number = 0
    for index, letter in enumerate(pattern):
        number += position_weigths[index] * BASE_PATTERN.find(letter)
    return number


def number_to_pattern_1(number, pattern_length):
    """Having a 4^k array with every k-mer combination sorted lexicographically,
    the functions return the pattern given the position number in the array."""
    position_weigths = [4**i for i in reversed(range(0, pattern_length))]
    pattern = ''
    for weight in position_weigths:
        pattern_weight = int(number/weight)
        pattern += BASE_PATTERN[pattern_weight]
        number -= weight * pattern_weight
    return pattern


# Recursive Version
def pattern_to_number(pattern):
    """Having a 4^k array with every k-mer combination sorted lexicographically,
    the functions returns the position of a given pattern in the array.

    This takes into account that if I remove a letter of a k-mer. The length of the
    array is 4 times the array of the (k-1)-mer.
    """
    if not pattern:
        return 0
    symbol = pattern[-1:]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol)


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


def complementary_number(number, k):
    pattern = number_to_pattern(number, k)
    pattern_c = get_complementary(pattern)
    return pattern_to_number(pattern_c)
