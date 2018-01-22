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

BASE_PATTERN = 'ACGT'


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
    profile = {base: [0] * k for base in BASE_PATTERN}
    for index in range(k):
        for motif in motifs:
            profile[motif[index]][index] += 1
    return {base: [value/(len(motifs) + 0) for value in items]
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
            best_kmer = most_probable_kmer(dna_string, k, profile)
            motifs.append(best_kmer)
        score = entropy_score(motifs)
        print(f'motifs -> {motifs[0]} score -> {score}')
        if score < best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


k_t, *dna_list = sys.stdin.read().splitlines()
k = int(k_t.split(' ')[0])
output = greedy_motif_search(dna_list, k)
for o in output:
    print(o)
