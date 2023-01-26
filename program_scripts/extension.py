import sys

# computes extended version of given seed
def extend_seed(seed, scoring_matrix, letter_map, X, S, database, query, k):
    database_hit = database[seed[0]][seed[1]:seed[1]+k]
    query_hit = query[seed[2]:seed[2]+k]
    max_score_right = get_score(scoring_matrix, letter_map, database_hit, query_hit)
    max_counter_right = 0
    counter = 1
    while not (seed[1]+k+counter > len(database[seed[0]]) or seed[2]+k+counter > len(query)):
        extended_database = database[seed[0]][seed[1]:seed[1]+k+counter]
        extended_query = query[seed[2]:seed[2]+k+counter]
        new_score = get_score(scoring_matrix, letter_map, extended_database, extended_query)
        if max_score_right - new_score > X:
            break
        if new_score >= max_score_right:
            max_score_right = new_score
            max_counter_right = counter
        counter += 1

    max_score_left = get_score(scoring_matrix, letter_map, database_hit, query_hit)
    max_counter_left = 0
    counter = 1
    while not (seed[1]-counter < 0 or seed[2] - counter < 0):
        extended_database = database[seed[0]][seed[1]-counter:seed[1]+k]
        extended_query = query[seed[2]-counter:seed[2]+k]
        new_score = get_score(scoring_matrix, letter_map, extended_database, extended_query)
        if max_score_left - new_score > X:
            break
        if new_score >= max_score_left:
            max_score_left = new_score
            max_counter_left = counter
        counter += 1

    best_database = database[seed[0]][seed[1]-max_counter_left:seed[1]+k+max_counter_right] 
    best_query = query[seed[2]-max_counter_left:seed[2]+k+max_counter_right] 
    best_score = get_score(scoring_matrix, letter_map, best_database, best_query)
    
    if best_score > S:
        # returns (alignment score, length of alignment, database sequence number, pos #, query-pos #)
        return [(best_score, max_counter_right + k + max_counter_left, seed[0], seed[1]-max_counter_left, seed[2]-max_counter_left)]
    else:
        return []

# calculates ungapped alignment score between two same length sequences
def get_score(scoring_matrix, letter_map, seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        score += scoring_matrix[letter_map[seq1[i]]][letter_map[seq2[i]]]
    return score

# computes extension phase of BLAST algorithm
def BLAST_extension():
    try:
        database_file = sys.argv[1]
        matrix_file = sys.argv[2]
        seed_file = sys.argv[3]
        X = int(sys.argv[4])
        S = int(sys.argv[5])
    except:
        print("ERROR: Should be 5 inputs to BLAST extension algorithm.")
        return
    
    with open(database_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        database = [] # list of database entries (strings)
        while entry != '':
            database.append(entry)
            entry = file.readline().replace('\n', "").upper()
    
    letter_map = {} # maps protein/nucleotides with number representations to help index into scoring_matrix
    with open(matrix_file, 'r') as file:
        labels = file.readline().replace('\n', "").split()
        scoring_matrix = [[0]*(len(labels)-1) for i in range(len(labels)-1)]
        for i in range(1, len(labels)):
            letter_map[labels[i]] = i-1
        for i in range(len(labels)-1):
            scores = file.readline().replace('\n', "").split()
            for j in range(1, len(scores)):
                scoring_matrix[i][j-1] = int(scores[j])

    with open(seed_file, 'r') as file:
        query = file.readline().replace('\n', '').upper()
        k = int(file.readline().replace('\n', ''))
        T = int(file.readline().replace('\n', ''))
        file.readline() # skips total seeds found
        seeds = [] # list of tuples in form (seq #, pos #, query-pos #)
        
        entry = file.readline().replace('\n', "").split()
        while len(entry) != 0:
            seeds.append((int(entry[1]), int(entry[3]), int(entry[5])))
            entry = file.readline().replace('\n', "").split()

    if X <= 0:
        print("ERROR: X value should be positive.")
        return

    if S < 0:
        print("ERROR: S value should be positive.")
        return

    extended_seeds = set()
    for seed in seeds:
        # returns (alignment score, length of alignment, database sequence number, pos #, query-pos #)
        extended_seeds = extended_seeds.union(extend_seed(seed, scoring_matrix, letter_map, X, S, database, query, k))

    # sort
    sorted_seeds = sorted(extended_seeds, key=lambda seed: (seed[0], seed[1], -seed[2], -seed[3], -seed[4]), reverse=True)

    print(query)
    print(k)
    print(T)
    print(len(sorted_seeds))
    print('-----------------------------------')
    for seed in sorted_seeds:
        print("Sequence " + str(seed[2]) + " Position " + str(seed[3]) + " Q-index " + str(seed[4]))
        print(query[seed[4]:seed[4]+seed[1]])
        print(database[seed[2]][seed[3]:seed[3]+seed[1]])
        print(seed[0])
        print('-----------------------------------')

if __name__ == "__main__":
    BLAST_extension()