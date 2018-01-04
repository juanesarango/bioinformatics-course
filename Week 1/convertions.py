
BASE_PATTERN = 'ACGT'

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}

# Get complement nucleotide
def get_complementary(text):
    return ''.join(reversed([COMPLEMENT[n] for n in text]))

# Convert to weights ACGT <-> 0123 
def symbol_to_number(symbol):
    return BASE_PATTERN.find(symbol)

def number_to_symbol(number):
    return BASE_PATTERN[number]

# Normal Version
def pattern_to_number_1(pattern):
    position_weigths = [4**i for i in reversed(range(0, len(pattern)))]
    number = 0
    for index, letter in enumerate(pattern):
        number += position_weigths[index] * BASE_PATTERN.find(letter)
    return number

def number_to_pattern_1(number, pattern_length):
    position_weigths = [4**i for i in reversed(range(0, pattern_length))]
    pattern = ''
    for weight in position_weigths:
        pattern_weight = int(number/weight)
        pattern += BASE_PATTERN[pattern_weight]
        number -= weight * pattern_weight
    return pattern

# Recursive Version
def pattern_to_number(pattern):
    if not pattern:
        return 0
    symbol = pattern[-1:]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + symbol_to_number(symbol) 

def number_to_pattern(number, k):
    if k == 1:
        return number_to_symbol(number)
    prefix_number = int(number/4)
    remainer = number - prefix_number
    symbol = number_to_symbol(remainer)
    prefix_pattern = number_to_pattern(prefix_number, k - 1)
    return prefix_pattern + symbol
