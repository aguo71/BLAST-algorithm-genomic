import sys

# calculates seeds from alphabet which score greater than T against any seed from the query
def get_seeds(alphabet_kmers, query_kmers, scoring_matrix, letter_map, T):
    seeds = [] # list of tuples in form (seq #, pos #, query-pos #)
    for a in alphabet_kmers:
        for q in query_kmers:
            if get_score(scoring_matrix, letter_map, a[0], q[0]) > T:
                seeds.append((a[1], a[2], q[1]))
    return seeds

# calculates ungapped alignment score between two same length sequences
def get_score(scoring_matrix, letter_map, seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        score += scoring_matrix[letter_map[seq1[i]]][letter_map[seq2[i]]]
    return score

# calculates all kmers in one entry of database
def get_database_kmers(entry, k, i):
    kmers = set()
    start = 0
    end = k
    while end <= len(entry):
        kmer = (entry[start: end], i, start) # saves tuple of kmer, database sequence #, index position #
        if kmer not in kmers:
            kmers.add(kmer)
        start += 1
        end += 1

    return kmers

# calculates all kmers from query
def get_query_kmers(query, k):
    kmers = set()
    start = 0
    end = k
    while end <= len(query):
        kmer = (query[start: end], start) #saves tuple of kmer, index position #
        if kmer not in kmers:
            kmers.add(kmer)
        start += 1
        end += 1

    return kmers

# computes seed phase of BLAST algorithm
def BLAST_seeding():
    try:
        database_file = sys.argv[1]
        query_file = sys.argv[2]
        matrix_file = sys.argv[3]
        k = int(sys.argv[4])
        T = int(sys.argv[5])
    except:
        print("ERROR: Should be 5 inputs to BLAST seeding algorithm.")
        return
    
    with open(database_file, 'r') as file:
        entry = file.readline().replace('\n', "").upper()
        database = [] # list of database entries (strings)
        while entry != '':
            database.append(entry)
            entry = file.readline().replace('\n', "").upper()

    with open(query_file, 'r') as file:
        query = file.readline().replace('\n', "").upper()
    
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

    if k <= 0:
        print("ERROR: k value should be positive.")
        return

    query_kmers = get_query_kmers(query, k)
    alphabet_kmers = set()
    for i in range(len(database)):
        alphabet_kmers = alphabet_kmers.union(get_database_kmers(database[i], k, i))
    seeds = get_seeds(alphabet_kmers, query_kmers, scoring_matrix, letter_map, T)
    seeds = sorted(seeds)

    print(query)
    print(k)
    print(T)
    print(len(seeds))
    for seed in seeds:
        print("Sequence " + str(seed[0]) + " Position " + str(seed[1]) + " Q-index " + str(seed[2]))

if __name__ == "__main__":
    BLAST_seeding()