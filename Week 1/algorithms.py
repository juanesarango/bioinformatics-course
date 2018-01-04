from convertions import *
import requests
import numpy as np


def pattern_count(text, pattern):
    """Counts the number of times a pattern is in text. """
    return len([i 
        for i in range(0, len(text)-len(pattern)+1) 
        if text[i:i+len(pattern)]== pattern])

def frequent_words(text, k):
    """Returns the patterns(k-mers) that are more frequent. """    
    frequent_patterns = []
    count={}
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        count[i] = pattern_count(text, pattern)
    max_count = max(count.values()) if len(count) else 0
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
    count={}
    for i in range(0, len(text)-k+1):
        pattern = text[i:i+k]
        count[i] = pattern_count(text, pattern)
        if count[i] >= t and pattern not in frequent_patterns:
            frequent_patterns.append(text[i:i+k])
    return frequent_patterns

def find_position(pattern, text):
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
        positions = find_position(k_mer, text)
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
