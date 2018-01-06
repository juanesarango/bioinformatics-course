from algorithms import *


def first_week_exercises():
    # # data = requests.get('http://bioinformaticsalgorithms.com/data/realdatasets/Replication/Vibrio_cholerae.txt')
    # # text = data.text
    # # pattern_count(text, 'GCGCGGCG')
    # # frequent_words(text, 'GCGCGGCG')

    data = requests.get('http://bioinformaticsalgorithms.com/data/realdatasets/Rearrangements/E_coli.txt')
    find_clumps(data.text, 9, 500, 3)

    # file_data = open( '/Users/arangooj/Downloads/dataset_4_5.txt')
    # input = list(file_data)
    # genome = input[0].replace('\n', '')
    # params = input[1].replace('\n', '').split(' ')
    # k = int(params[0])
    # L = int(params[1])
    # t = int(params[2])
    # find_clumps(genome, k, L, t)

    # data = requests.get('https://stepik.org/media/attachments/bioinformatics/FrequencyArray.txt')
    # params = data.text.split('\r\n')
    # text = params[1]
    # k = int(params[2])
    # expected_output = params[4]
    # obtained_output = compute_freq(text, k)
    # assert expected_output == obtained_output

    # pattern_to_number(pattern, k) = 4 * pattern_to_number(pattern[:-1], k) + BASE_PATTERN.find(pattern[-1:])
