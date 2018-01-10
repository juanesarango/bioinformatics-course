from algorithms import *
import requests
import sys


def exercise_test():
    data = requests.get('http://bioinformaticsalgorithms.com/data/extradatasets/replication/minimum_skew.txt')
    args = data.text.split('\r\n')
    genome = args[1]
    expected_output = args[3]

    skew_array = skew(genome)

    output = find_array_positions(min(skew_array), skew_array)
    print(output)
    assert ls2str(output) == expected_output
    return output


def exercise_1():
    genome = sys.stdin.read().splitlines()
    skew_array = skew(genome)

    output = find_array_positions(min(skew_array), skew_array)

    return output
