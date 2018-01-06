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
    return ' '.join([str(i) for i in skew_list])
