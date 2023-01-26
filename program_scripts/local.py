import sys

# implements Smith-Waterman local alignment algorithm, compatible with nucleotide or protein scoring matrices
def localAlignment():
    try:
        seq_file = sys.argv[1]
        matrix_file = sys.argv[2]
        gap_penalty = sys.argv[3]
    except:
        print("ERROR: Should be 3 inputs to local alignment algorithm.")
        return
    
    with open(seq_file, 'r') as file:
        str1 = file.readline().replace('\n', "").upper() #len = M rows
        str2 = file.readline().replace('\n', "").upper() #len = N cols

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

    if gap_penalty == 'negInf':
        gap_penalty = float('-inf')
    else:
        gap_penalty = int(gap_penalty)
        if gap_penalty > 0:
            print("ERROR: Gap Penalty should be non-negative.")

    M = len(str1)
    N = len(str2)

    edit_graph = [[float('-inf')]*(N+1) for _ in range(M+1)] # M+1 x N+1
    edit_graph[0][0] = 0
    for i in range(1, M+1):
        edit_graph[i][0] = 0
    for j in range(1, N+1):
        edit_graph[0][j] = 0

    optimal_score = 0
    optimal_i = M
    optimal_j = N
    for j in range(1, N+1):
        for i in range(1, M+1):
            unit_score = scoring_matrix[letter_map[str1[i-1]]][letter_map[str2[j-1]]] 
            edit_graph[i][j] = max(0, edit_graph[i-1][j-1] + unit_score, edit_graph[i-1][j] + gap_penalty, 
                                edit_graph[i][j-1] + gap_penalty)
            if edit_graph[i][j] >= optimal_score:
                optimal_score = edit_graph[i][j]
                optimal_i = i
                optimal_j = j

    alignment1 = ""
    alignment2 = ""
    i = optimal_i
    j = optimal_j
    while i > 0 or j > 0:
        if edit_graph[i][j] == 0:
            break
        if i > 0 and (edit_graph[i][j] - gap_penalty) == edit_graph[i-1][j]:
            alignment1 = str1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        elif j > 0 and (edit_graph[i][j] - gap_penalty) == edit_graph[i][j-1]:
            alignment1 = '-' + alignment1
            alignment2 = str2[j-1] + alignment2
            j -= 1
        else:
            alignment1 = str1[i-1] + alignment1
            alignment2 = str2[j-1] + alignment2
            i -= 1
            j -= 1

    print(alignment1)
    print(alignment2)
    print(optimal_score)

if __name__ == "__main__":
    localAlignment()