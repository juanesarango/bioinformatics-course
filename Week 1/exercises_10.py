from algorithms import *
import sys


def exercise_1():
    text, pattern = sys.stdin.read().splitlines()
    
    print(pattern_count(text, pattern))

def exercise_2():
    text, k = sys.stdin.read().splitlines()
    
    for word in frequent_words(text, k):
        print(word)

def exercise_3():
    text, = sys.stdin.read().splitlines()

    print(get_complementary(text))

def exercise_4():
    pattern, text = sys.stdin.read().splitlines()

    print(' '.join([str(i) for i in find_position(pattern, text)]))

def exercise_5():
    text, params = sys.stdin.read().splitlines()
    k, L, t = [int(p) for p in params.split(' ')]
    
    for clump in find_clumps(text, k, L, t):
        print(clump)

