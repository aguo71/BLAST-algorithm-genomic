import sys
import numpy as np

# calculates ungapped local alignment between 2 sequences
def align(str1, str2, scoring_matrix, letter_map):
    gap_penalty = float('-inf')

    M = len(str1)
    N = len(str2)

    edit_graph = [[float('-inf')]*(N+1) for _ in range(M+1)] # M+1 x N+1
    edit_graph[0][0] = 0
    for i in range(1, M+1):
        edit_graph[i][0] = 0
    for j in range(1, N+1):
        edit_graph[0][j] = 0

    optimal_score = 0
    for j in range(1, N+1):
        for i in range(1, M+1):
            unit_score = scoring_matrix[letter_map[str1[i-1]]][letter_map[str2[j-1]]] 
            edit_graph[i][j] = max(0, edit_graph[i-1][j-1] + unit_score, edit_graph[i-1][j] + gap_penalty, 
                                edit_graph[i][j-1] + gap_penalty)
            if edit_graph[i][j] >= optimal_score:
                optimal_score = edit_graph[i][j]
    
    return optimal_score

# generates a random sequence given values with probabilities
def generate_random_seq(freq_list, min_len, max_len, val_count):
    values = []
    freqs = []
    for entry in freq_list:
        values.append(entry[0])
        freqs.append(entry[1]/float(val_count))

    length = np.random.randint(min_len, max_len+1)
    seq = np.random.choice(values, size=length, p=freqs)
    return ''.join(seq)

# counts occurrence of each unique value in database in map, and returns total number of values
def get_freqs(database):
    freq_map = {}
    total_count = 0
    for entry in database:
        for value in entry:
            total_count += 1
            if value in freq_map:
                freq_map[value] += 1
            else:
                freq_map[value] = 1

    return (freq_map, total_count)

# analyzes statistics of dataset
def statistics():
    try:
        database_file = sys.argv[1]
        matrix_file = sys.argv[2]
    except:
        print("ERROR: Should be 2 inputs to statistics script.")
        return
    
    with open(database_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        database = [] # list of database entries (strings)
        max_len = float('-inf')
        min_len = float('inf')
        while entry != '':
            database.append(entry)
            if len(entry) > max_len:
                max_len = len(entry)
            if len(entry) < min_len:
                min_len = len(entry)
            entry = file.readline().replace('\n', "").upper()

    letter_map = {}

    with open(matrix_file, 'r') as file:
        labels = file.readline().replace('\n', "").split()
        scoring_matrix = [[0]*(len(labels)-1) for i in range(len(labels)-1)]
        for i in range(1, len(labels)):
            letter_map[labels[i]] = i-1
        for i in range(len(labels)-1):
            scores = file.readline().replace('\n', "").split()
            for j in range(1, len(scores)):
                scoring_matrix[i][j-1] = int(scores[j])

    (freq_map, val_count) = get_freqs(database)
    freq_list = list(freq_map.items())
    freq_list = sorted(freq_list, key=lambda x: x[1], reverse=True)
    for key in freq_list:
        print(key[0] + ": " + str(key[1]) + '/' + str(val_count) + ' = ' + str(float(key[1])/val_count))

    random_sequences = []
    for _ in range(100):
        random_sequences.append(generate_random_seq(freq_list, min_len, max_len, val_count))
    with open('random.txt', 'w') as f:
        f.write('\n'.join(random_sequences))
    print(random_sequences)
    
    alignments = []
    for i in range(0, 100, 2):
        alignments.append(align(random_sequences[i], random_sequences[i+1], scoring_matrix, letter_map))
    print(sorted(alignments))

if __name__ == "__main__":
    statistics()