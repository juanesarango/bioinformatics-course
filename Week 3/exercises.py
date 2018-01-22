from algorithms import *
import sys


def exercise_test():
    k, *dna = sys.stdin.read().splitlines()
    output = median_string(dna, k)
    return output

print(exercise_test())

def exercise_graded():
    k_t, *dna_list = sys.stdin.read().splitlines()
    k = int(k_t.split(' ')[0])
    output = greedy_motif_search(dna_list, k)
    for o in output:
        print(o)