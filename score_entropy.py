import sys
from algorithms import *

motifs = sys.stdin.read().splitlines()
print(entropy_score(motifs))