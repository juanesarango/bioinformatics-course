"""File to run the code that will be input into the
embedded code editor of the Coursera Online Classroom.
Inputs come via `stdin` from `text_input.txt`.
Outputs must be `stdout` by printing them.

To run it from bash:

    $ `python test.py < test_input.txt`
"""

import sys
import math
import numpy as np
import random

BASE_PATTERN = 'ACGT'


def hamming_distance(text_p, text_q):
    """Returns the hamming distance between 2 segments.
    This distance is the total sum of mismatches for eachj position
    `i` if `p[i] != q[i]`.
    """
    if len(text_p) != len(text_q):
        return -1
    return sum([text_p[i] != text_q[i] for i in range(len(text_p))])


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


def find_consensus(motifs):
    """Returns the first string consensus from a given
    list of motifs.
    """
    k = len(motifs[0])
    base_len = len(BASE_PATTERN)
    profile = create_profile(motifs)
    # /print(profile)
    profile = np.array([bij for _, bi in profile.items()
                       for bij in bi]).reshape(base_len, k).T
    return ''.join([BASE_PATTERN[np.argmax(array_in_position)]
                    for array_in_position in profile])


def get_score(motifs):
    """Returns the sum distance from all strings away from
    the consensus motif.
    """
    pattern = find_consensus(motifs)
    # for motif in motifs:
    #     print(f'{motif} Motif')
    # print('--------------------')
    # print(f'{pattern} Consensus')
    return d_dna_list(pattern, motifs)


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


def get_random_motifs(dna_list, k):
    """This method takes a random k-mer from each of the DNA
    strings of the list. It returns a list of random k-mers.
    """
    motifs = []
    for dna_string in dna_list:
        random_index = random.randint(0, len(dna_list[0]) - k)
        motifs.append(dna_string[random_index:random_index + k])
    return motifs


def get_motifs(profile, dna_list):
    """Given a profile and a DNA list, it returns the most probable
    k-mer for each string of the list.
    """
    k = len(profile['A'])
    return [most_probable_kmer(dna_string, k, profile)
            for dna_string in dna_list]


def randomized_motif_search_r(dna_list, k):
    """Given a list of DNA strings. It creates a profile with the results of the
    previous iteration. The first set of k-mer motifs is random.
    """
    motifs = get_random_motifs(dna_list, k)
    best_motifs = motifs
    best_score = get_score(best_motifs)
    while True:
        # print(f'Score0 {best_score}')
        profile = create_profile(motifs)
        motifs = get_motifs(profile, dna_list)
        score = get_score(motifs)
        # print(f'Score1 {score}')
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            return best_motifs


def randomized_motif_search(dna_list, k):
    """Given a list of DNA strings. It creates a profile with the results of the
    previous iteration. The first set of k-mer motifs is random.
    """
    motifs = []
    for dna_string in dna_list:
        random_index = random.randint(0, len(dna_list[0]) - k)
        motifs.append(dna_string[random_index:random_index + k])
    best_motifs = motifs
    best_score = get_score(best_motifs)
    n = 0
    while True:
        n += 1
        # print(f'Score0 {best_score} {n}')
        profile = create_profile(motifs)
        motifs = get_motifs(profile, dna_list)
        score = get_score(motifs)
        # print(f'Score1 {score}')
        if score < best_score:
            best_motifs = motifs
            best_score = score
        else:
            # print(f'Iterations {n}')
            return best_motifs


def run_n_times(n, dna_list, k):
    best_motifs = []
    best_score = 100
    for i in range(n):
        motifs = randomized_motif_search(dna_list, k)
        score = get_score(motifs)
        # print(f'{motifs} -> {score}')
        if score <= best_score:
            best_motifs = motifs
            best_score = score
    return best_motifs

k_t, *dna_list = sys.stdin.read().splitlines()
k = int(k_t.split(' ')[0])

output = run_n_times(1000, dna_list, k)
for o in output:
    print(o)
print(f'{entropy_score(output)}')