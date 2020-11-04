import math
import random
import requests

import numpy as np

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
    kmers = frequent_words_t(text, k, t)
    for kmer in kmers:
        positions = find_text_positions(kmer, text)
        for position in positions:
            subtext = text[position:position + L]
            count = pattern_count(subtext, kmer)
            if count >= t and kmer not in clumps:
                clumps.append(kmer)
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


def neighbors_1(pattern):
    """Generates the d-neighborhood, the set of all k-mers whose Hamming
    distance from Pattern does not exceed d. Version 1.0
    """
    neighborhood = []
    for i in range(len(pattern)):
        symbol = pattern[i]
        for base in BASE_PATTERN:
            neighbor = pattern[:i] + base + pattern[i+1:]
            if base != symbol:
                neighborhood.append(neighbor)
    return neighborhood


def neighbors_2(pattern):
    """Generates the d-neighborhood, the set of all k-mers whose Hamming
    distance from Pattern does not exceed d. Version 2.0
    """
    return [
        pattern[:index] + base + pattern[index + 1:]
        for index in range(len(pattern))
        for base in BASE_PATTERN
        if pattern[index] != base]


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


def iterative_neighbors(pattern, d):
    """Iterative """
    neighborhood = list(BASE_PATTERN)
    for i in range(1, d + 1):
        for neighbor in neighborhood:
            neighborhood.append([])
    pass

    """
        IterativeNeighbors(Pattern, d)
        Neighborhood ← set consisting of single string Pattern
        for j = 1 to d
            for each string Pattern’ in Neighborhood
                add ImmediateNeighbors(Pattern') to Neighborhood
                remove duplicates from Neighborhood
        return Neighborhood


        FrequentWordsWithMismatches(Text, k, d)
        FrequentPatterns ← an empty set
        Neighborhoods ← an empty list
        for i ← 0 to |Text| − k
            add Neighbors(Text(i, k), d) to Neighborhoods
        form an array NeighborhoodArray holding all strings in Neighborhoods
        for i ← 0 to |Neighborhoods| − 1
            Pattern ← NeighborhoodArray(i)
            Index(i) ← PatternToNumber(Pattern)
            Count(i) ← 1
        SortedIndex ← Sort(Index)
        for i ← 0 to |Neighborhoods| − 2
            if SortedIndex(i) = SortedIndex(i + 1)
                Count(i + 1) ← Count(i) + 1
       maxCount ← maximum value in array Count
       for i ← 0 to |Neighborhoods| − 1
           if Count(i) = maxCount
               Pattern ← NumberToPattern(SortedIndex(i), k)
               add Pattern to FrequentPatterns
       return FrequentPatterns
    """


def odds_of_kmer_in_strings(k, L, N):
    """Return the expected number of k-mers that can be found
    in N strings of lenght L
    """
    odds_kmer = 1/(4**k)
    kmer_in_string = L - k + 1
    return odds_kmer * kmer_in_string * N


# Week 3
def find_d_string(pattern, text, d):
    """Returns the positions where all the patterns are found
    with a maximum of d mismatches.
    """
    return [index
            for index in range(len(text) - len(pattern) + 1)
            if hamming_distance(text[index:index+len(pattern)], pattern) <= d]


def find_motifs(patterns, k, d):
    """Find hidden k-mer motifs with at least a d hamming-distance
    across the different list of patterns.
    These are called (k-d)-motifs.
    """
    motifs = []
    for i in range(4**k):
        kmer = number_to_pattern(i, k)
        if all([find_d_string(kmer, pattern, d) for pattern in patterns]) \
                and kmer not in motifs:
            motifs.append(kmer)
    return motifs


def entropy_score(list_motifs):
    """The entropy of a motif matrix is defined as the sum of the
    entropies of its columns.
    The entropy of each column is a measurement of the uncertainty
    of a probability distribution p. And it's defined as:

                          N
        H(p[1],…,p[N]) = −∑ (p[i]·log2(p[i]))
                         i=1
    """
    rows = len(list_motifs)
    columns = len(list_motifs[0])
    entropy = 0
    for i in range(0, columns):
        count_column = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }
        column_entropy = 0
        for j in range(0, rows):
            nucleotide = list_motifs[j][i]
            count_column[nucleotide] += 1
        for count in count_column.values():
            prob_dist = count/rows
            if prob_dist > 0:
                column_entropy -= prob_dist * math.log2(prob_dist)
        entropy += column_entropy
    return entropy


def median_string(dna_list, k):
    """Given a list of DNA strings, the method finds the k-mer
    that minimizes the hamming distance across all the list.
    This means this median string k-mer is a motif.
    """
    median = ''
    min_distance = len(dna_list) * len(dna_list[0])
    for i in range(4**k):
        kmer = number_to_pattern(i, k)
        list_distance = d_dna_list(kmer, dna_list)
        if list_distance < min_distance:
            min_distance = list_distance
            median = kmer
    return median


def d_dna_list(pattern, dna_list):
    """Minimum hamming distance of the list of DNA strings."""
    return sum([d_dna_string(pattern, dna_string)[0]
                for dna_string in dna_list])


def d_dna_string(pattern, dna_string):
    """Minimum hamming distance of the pattern with a
    k-mer found in the DNA string.
    """
    min_kmer = ''
    min_distance = len(dna_string)
    for i in range(len(dna_string) - len(pattern) + 1):
        kmer = dna_string[i:i + len(pattern)]
        distance = hamming_distance(kmer, pattern)
        if distance < min_distance:
            min_distance = distance
            min_kmer = kmer
    return min_distance, min_kmer


def string_probability(string, profile):
    """" Given a motif matrix profile, the probability of a get an
    specific string is the probability of getting each nucleotide at
    each of the positions.

    profile = {
        'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
        'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
        'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
    }

    Create profile matrix from stdin:
    ```
    BASE_PATTERN = 'ACGT'
    string, k, *list_profile = sys.stdin.read().splitlines()
    profile = {BASE_PATTERN[i]: [float(n) for n in list_profile[i].split(' ')]
            for i in range(len(ntds))}
    ```
    """
    return np.prod([profile[string[i]][i] for i in range(len(string))])


def most_probable_kmer(string, k, profile):
    """Given a motif matrix profile, this method returns the most probable
    k-mer that can be obtain from the string.
    """
    max_prob = -1
    best_kmer = ''
    for i in range(len(string) - k + 1):
        kmer = string[i:i + k]
        kmer_prob = string_probability(kmer, profile)
        if kmer_prob > max_prob:
            max_prob = kmer_prob
            best_kmer = kmer
    return best_kmer


def create_profile(motifs):
    """This method creates a motif profile, by using a list of k-mers
    Profile is a dictionary for each base where its value is a list
    of lenght k.
    """
    k = len(motifs[0])
    # Initialize in 1 to apply Laplace's rule of sucession to avoid 0's.
    profile = {base: [1] * k for base in BASE_PATTERN}
    for index in range(k):
        for motif in motifs:
            profile[motif[index]][index] += 1
    return {base: [value/(len(motifs) + len(BASE_PATTERN)) for value in items]
            for base, items in profile.items()}


def greedy_motif_search(dna_list, k):
    """Given a list of DNA strings. It creates a profile by updating the
    matrix with each k-mer of the list. Then it finds the most probable
    k-mer to find the best motifs.
    """
    best_motifs = [dna_string[:k] for dna_string in dna_list]
    first_dna = dna_list[0]
    best_score = k * len(dna_list)
    for i in range(len(first_dna) - k + 1):
        motifs = [first_dna[i:i + k]]
        for dna_string in dna_list[1:]:
            profile = create_profile(motifs)
            motifs.append(most_probable_kmer(dna_string, k, profile))
        score = entropy_score(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


# Week 4
def find_consensus(motifs):
    """Returns the first string consensus from a given
    list of motifs.
    """
    k = len(motifs[0])
    base_len = len(BASE_PATTERN)
    profile = create_profile(motifs)
    profile = np.array([bij for _, bi in profile.items()
                       for bij in bi]).reshape(base_len, k).T
    return ''.join([BASE_PATTERN[np.argmax(array_in_position)]
                    for array_in_position in profile])


def get_score(motifs):
    """Returns the sum distance from all strings away from
    the consensus motif.
    """
    pattern = find_consensus(motifs)
    return d_dna_list(pattern, motifs)


def get_random_motifs(dna_list, k):
    """This method takes a random k-mer from each of the DNA
    strings of the list. It returns a list of random k-mers.
    """
    motifs = []
    for dna_string in dna_list:
        random_index = random.randint(0, len(dna_list[0]))
        motifs.append(dna_string[random_index:random_index + k])
    return motifs


def get_motifs(profile, dna_list):
    """Given a profile and a DNA list, it returns the most probable
    k-mer for each string of the list.
    """
    k = len(profile['A'])
    return [most_probable_kmer(dna_string, k, profile)
            for dna_string in dna_list]


def randomized_motif_search(dna_list, k):
    """Given a list of DNA strings. It creates a profile with the results of the
    previous iteration. The first set of k-mer motifs is random.
    """
    motifs = get_random_motifs(dna_list, k)
    best_motifs = motifs
    best_score = get_score(best_motifs)
    while True:
        profile = create_profile(motifs)
        motifs = get_motifs(profile, dna_list)
        score = get_score(motifs)
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs
