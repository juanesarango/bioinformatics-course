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

COMPLEMENT = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}


def get_complementary(text):
    """Get complementary nucleotide chain."""
    return ''.join(reversed([COMPLEMENT[n] for n in text]))


def number_to_symbol(number):
    """Convert weights to nucleotides 0123 -> ACGT."""
    return BASE_PATTERN[number]


def symbol_to_number(symbol):
    """Convert to nucleotides to weights ACGT -> 0123."""
    return BASE_PATTERN.find(symbol)


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


def complementary_number(number, k):
    pattern = number_to_pattern(number, k)
    pattern_c = get_complementary(pattern)
    return pattern_to_number(pattern_c)


def find_array_positions(num, array):
    """Returns and array with all positions where the number is
    found within the array.
    """
    return [index for index, item in enumerate(array) if item == num]


def hamming_distance(text_p, text_q):
    """Returns the hamming distance between 2 segments.
    This distance is the total sum of mismatches for eachj position
    `i` if `p[i] != q[i]`.
    """
    if len(text_p) != len(text_q):
        return -1
    return sum([text_p[i] != text_q[i] for i in range(len(text_p))])


def find_kmers(pattern, text, d):
    """Returns the positions where the k-mer matchs the pattern
    with a hamming distance of d or less.
    """
    return [i for i in range(0, len(text) - len(pattern) + 1)
            if hamming_distance(pattern, text[i:i + len(pattern)]) <= d]


def aprox_pattern_count(pattern, text, d):
    """Implementation of Count_d(text, pattern) as the number of times
    a k-mer matchs the pattern with a hamming distance of d or less.
    """
    return len(find_kmers(pattern, text, d))


def frequent_kmers(text, k, d):
    """Find the most frequent k-mers with mismatches <= d in a string
    """
    freq_array = [
        aprox_pattern_count(number_to_pattern(kmer_index, k), text, d)
        for kmer_index in range(0, 4**k)]
    return [
        number_to_pattern(kmer_index, k)
        for kmer_index in find_array_positions(max(freq_array), freq_array)]


def frequent_kmers_complements(text, k, d):
    """Find the most frequent k-mers (with mismatches and reverse
    complements) in a string. Returns all k-mers pattern maximizing the sum
    Count_d(Text, Pattern)+ Count_d(Text, Pattern_rc) over all possible k-mers.
    """
    freq_array = [
        aprox_pattern_count(number_to_pattern(kmer_index, k), text, d)
        for kmer_index in range(0, 4**k)]
    comp_sum = [
        freq + freq_array[complementary_number(kmer_index, k)]
        for kmer_index, freq in enumerate(freq_array)]
    return [
        number_to_pattern(kmer_index, k)
        for kmer_index in find_array_positions(max(comp_sum), comp_sum)]


def ls2str(ls):
    """A lot of the algorithms returns an array or list.
    But to introduce the result in the classroom a string is required.
    """
    return ' '.join([str(i) for i in ls])


def skew(genome):
    """Calculate the skew for each nucleotide position
    by the difference of G-C accross the genome.

        ex:
         C  A  T  G  G  G  C  A  T  C  G  G  C  C  A  T  A  C  G  C  C
        0 -1 -1 -1  0  1  2  1  1  1  0  1  2  1  0  0  0  0 -1  0 -1 -2
    """
    skew_list = [0]
    for ntd in genome:
        skew_change = 1 if ntd == 'G' else -1 if ntd == 'C' else 0
        skew_list.append(skew_list[-1] + skew_change)
    return skew_list


def find_array_positions(num, array):
    """Returns and array with all positions where the number is
    found within the array.
    """
    return [index for index, item in enumerate(array) if item == num]


# genome, = sys.stdin.read().splitlines()
data = pd.read_csv('salmonella_enterica.txt', sep='\n', header=None)
genome = ' '.join([i[0] for i in data.values.tolist()])
k = 9
d = 1
skew_se = skew(genome)
ori = find_array_positions(min(skew_se), skew_se)
subgenome = genome[ori - 50:ori + 50]
output = frequent_kmers_complements(subgenome, k, d)
print(ori)
