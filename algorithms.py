import requests

from utils.convertions import *


# Week 1
def pattern_count(text, pattern):
    """Counts the number of times a pattern is in text. """
    return len([i
                for i in range(0, len(text) - len(pattern) + 1)
                if text[i:i + len(pattern)] == pattern])


def frequent_words(text, k):
    """Returns the patterns(k-mers) that are more frequent. """
    frequent_patterns = []
    count = {}
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        count[i] = pattern_count(text, pattern)
    max_count = max(count.values()) if count else 0
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        if count[i] == max_count and pattern not in frequent_patterns:
            frequent_patterns.append(text[i:i+k])
    return frequent_patterns


def frequent_words_t(text, k, t):
    """Returns the patterns(k-mers) whose frequency of
    occurrence (count) is greater than t.
    """
    frequent_patterns = []
    count = {}
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        count[i] = pattern_count(text, pattern)
        if count[i] >= t and pattern not in frequent_patterns:
            frequent_patterns.append(text[i:i+k])
    return frequent_patterns


def find_text_positions(pattern, text):
    """Returns and array with all positions where the pattern is
    found within the text.
    """
    positions = []
    i = 0
    while text[i:] and text[i:].find(pattern) != -1:
        position = text[i:].find(pattern) + i
        positions.append(position)
        i = position + 1
    return positions


def find_clumps(text, k, L, t):
    """ It slides a window L and it finds the most common
    patterns (k-mers) whose value is greater than t"""
    clumps = []
    k_mers = frequent_words_t(text, k, t)
    for k_mer in k_mers:
        positions = find_text_positions(k_mer, text)
        for position in positions:
            subtext = text[position:position + L]
            count = pattern_count(subtext, k_mer)
            if count >= t and k_mer not in clumps:
                clumps.append(k_mer)
    return clumps


def compute_freq(text, k):
    """A k-mer can be arranged in a 4^k ordered array.
    This function returns an array of the frequency of each of the k-mers
    in the text. The position of the array can be matched with the
    pattern using pattern_to_number and number_to_pattern.
    """
    freq_array = [0 for i in range(0, 4**k)]
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i + k]
        j = pattern_to_number(pattern)
        freq_array[j] += 1
    # return ' '.join([str(i) for i in freq_array])
    return freq_array


def faster_frequent_words(text, k):
    """A better version to find frequent k-mers
    using the compute_freq function.
    """
    frequent_patterns = []
    freq_array = compute_freq(text, k)
    max_count = max(freq_array)
    for i in range(0, len(text)-k+1):
        if freq_array[i] == max_count:
            pattern = number_to_pattern(i, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


def frequent_words_by_sorting(text, k):
    """An alternative method of finding frequent words."""
    frequent_patterns = []
    index = []
    count = []
    for i in range(0, len(text) - k + 1):
        pattern = text[i:i + k]
        index[i] = pattern_to_number(pattern)
        count[i] = 1
    sorted_index = sorted(index)
    for i in range(0, len(text) - k + 1):
        if sorted_index[i] == sorted_index[i-1]:
            count[i] = count[i -1] + 1
    max_count = max(count)
    for i in range(0, len(text) - k + 1):
        if count[i] == max_count:
            pattern = number_to_pattern(sorted_index[i], k)
            frequent_patterns.append(pattern)
    return frequent_patterns


def clumps_finding(text, k, t, L):
    """Same that find_clumps but using the improved functions
    to find frequent patterns
    """
    frequent_patterns = []
    clumps = [0 for i in range(0, 4**k)]
    for i in range(0, len(text) - L + 1):
        subtext = text[i:i + L]
        freq_array = compute_freq(subtext, k)
        for index, freq in enumerate(freq_array):
            if freq >= t:
                clumps[index] = 1
    for index, clump in enumerate(clumps):
        if clump == 1:
            pattern = number_to_pattern(index, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


def better_clumps_finding(text, k, t, L):
    """Same that clumps_finding but instead of calling the
    compute_freq for each L is only done once. And updates
    with the first and last k-mer.
    """
    frequent_patterns = []
    clumps = [0 for i in range(0, 4**k)]
    first_subtext = text[:L]
    freq_array = compute_freq(first_subtext, k)
    for index, freq in enumerate(freq_array):
        if freq >= t:
            clumps[index] = 1
    for i in range(1, len(text) - L + 1):
        old_kmer = text[i - 1:i - 1 + k]
        old_kmer_number = pattern_to_number(old_kmer)
        freq_array[old_kmer_number] -= 1
        new_kmer = text[i + L:i + L + k]
        new_kmer_number = pattern_to_number(new_kmer)
        freq_array[new_kmer_number] += 1
        if freq_array[new_kmer_number] >= t:
            clumps[new_kmer_number] = 1
    for index, clump in enumerate(clumps):
        if clump == 1:
            pattern = number_to_pattern(index, k)
            frequent_patterns.append(pattern)
    return frequent_patterns


# Week 2

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


def ls2str(ls):
    """A lot of the algorithms returns an array or list.
    But to introduce the result in the classroom a string is required.
    """
    return ' '.join([str(i) for i in ls])


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


def find_neighbors(pattern, d):
    """
    """
    pass


def frequent_kmers(text, k, d):
    """Find the most frequent k-mers with mismatches <= d in a string.
    """
    freq_kmers = []
    max_count = 0
    for kmer_index in range(0, 4**k):
        kmer = number_to_pattern(kmer_index, k)
        count = aprox_pattern_count(kmer, text, d)
        if count > max_count:
            max_count = count
            freq_kmers = []
        if count == max_count:
            freq_kmers.append(kmer)
    return freq_kmers


def frequent_kmers_2(text, k, d):
    """Find the most frequent k-mers with mismatches <= d in a string.
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
