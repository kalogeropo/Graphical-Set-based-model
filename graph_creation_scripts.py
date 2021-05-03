import numpy


##############################################################################################################################
### PART 1: NIKOLAOS SKAMNELOS usefull graph functions:

#CREATES a graph using a fixed sized window, by spliting the original textfile into sub-textfiles
def splitFileConstantWindow(file, window):
    # Open the file and split it into words
    inputFile = open(file, 'r').read().split()
    outputFile = []

    # Join words according to window
    for i in range(0, len(inputFile), window):
        outputFile.append(' '.join(inputFile[i:i + window]))

    # print(outputData)
    return outputFile


def CreateAdjMatrixFromInvIndexWithWindow(terms, file, window_size):
    # print("Adj_Matrix = %d * %d " % (len(terms), len(tf)))
    # print(terms)
    adj_matrix = numpy.zeros(shape=(len(terms), len(terms)))
    split_file = splitFileConstantWindow(file, window_size)
    for subfile in split_file:
        window_terms = subfile.split()
        for term in window_terms:
            # print("\n")
            # print(term)
            row_index = terms.index(term)
            # print("TERM:",row_index)
            for x in range(0, len(window_terms)):
                col_index = terms.index(window_terms[x])
                # print("Y TERM:",col_index)
                adj_matrix[row_index][col_index] += 1
            adj_matrix[row_index][row_index] -= 1
    # print(adj_matrix)
    # fullsize = rows.size + row.size + col.size + adj_matrix.size
    # print(fullsize / 1024 / 1024)
    return (adj_matrix)

### This function covers the constant window function, but if window == 0 then a percentage of length of the text
### is considered as the window size. Althoigh the size is lower bound with size 5.
def splitFilePercentageWindow(file, window):
    # Open the file and split it into words
    inputFile = open(file, 'r').read().split()
    num_of_words = len(inputFile)
    outputFile = []
    perc = 0.01055
    # If window is equal to zero get window according to length
    if window == 0:
        window = int(num_of_words * perc) + 1 #0.01055 IS THE PERCENTAGE LENGTH OF WINDOW
        # print("Window Size: ", window)
        if window < 5:
            window = 5

    # Join words according to window
    for i in range(0, num_of_words, window):
        outputFile.append(' '.join(inputFile[i:i + window]))

    # print(outputData)
    return outputFile

### OVERLOADING the percentage variable into the arguments for automation
def splitFilePercentageWindow(file, window, per_window):
    # Open the file and split it into words
    inputFile = open(file, 'r').read().split()
    num_of_words = len(inputFile)
    outputFile = []

    # If window is equal to zero get window according to length or if percentage window flag is true
    if window == 0:
        # print(per_window)
        window = int(num_of_words * per_window) + 1
        # print("Window Size: ", window)
        if window < 8:
            window = 8

    # Join words according to window
    for i in range(0, num_of_words, window):
        outputFile.append(' '.join(inputFile[i:i + window]))

    # print(outputData)
    return outputFile


#here the graph is created using a overlapping sliding window as Graph of word dictates
def CreateAdjMatrix_Vazirgiannis_implementation(terms, file, window_size):
    # print("Adj_Matrix = %d * %d " % (len(terms), len(tf)))
    #print(terms)
    adj_matrix = numpy.zeros(shape=(len(terms),len(terms)))
    split_file = open(file, 'r').read().split() #splitFileConstantWindow(file, window_size)
    counter = 0
    for term in split_file:
        row_index = terms.index(term)
        for x in range(0,window_size):
            try:
                col_index = terms.index(split_file[counter + x])
                adj_matrix[row_index][col_index]+=1
            except IndexError:
                break
        counter+=1
        adj_matrix[row_index][row_index]-=1
    return (adj_matrix)

#PARAGRAPH and SENTANCE importance differation implementation
def CreateAdjMatrixFromInvIndexWithSenParWindow(terms, file, sen_window_size, par_window_size):
    matrix_size = len(terms)

    # Create the matrices
    sen_adj_matrix = numpy.zeros(shape=(matrix_size, matrix_size))
    par_adj_matrix = numpy.zeros(shape=(matrix_size, matrix_size))

    # Get the adjacency matrix for each window
    sen_adj_matrix = CreateAdjMatrixFromInvIndexWithWindow(terms, file, sen_window_size)
    par_adj_matrix = CreateAdjMatrixFromInvIndexWithWindow(terms, file, par_window_size)

    # Create the final Matrix
    final_adj_matrix = numpy.zeros(shape=(matrix_size, matrix_size))

    # Create coefficients a and b
    a = 1.0
    b = 0.05

    # Add the two matrices
    final_adj_matrix = [[a * sen_adj_matrix[r][c] + b * par_adj_matrix[r][c] for c in range(len(sen_adj_matrix[0]))] for
                        r in range(matrix_size)]
    # print(final_adj_matrix)

    return final_adj_matrix
########################################################################################################################
######### PART 2:KALOGEROPOULOS graph creation proccess and usefull graph functions


def createInvertedIndexFromFile(file, postingl):
    with open(file, 'r') as fd:
        # list containing every word in text document
        text = fd.read().split()
        uninque_terms = []
        termFreq = []
        for term in text:
            if term not in uninque_terms:
                uninque_terms.append(term)
                termFreq.append(text.count(term))
            if term not in postingl:
                postingl.append(term)
                postingl.append([file, text.count(term)])
            else:
                existingtermindex = postingl.index(term)
                if file not in postingl[existingtermindex + 1]:
                    postingl[existingtermindex + 1].extend([file, text.count(term)])
    # print(len(uninque_terms))
    # print(termFreq)
    return (uninque_terms, termFreq, postingl, len(text))

    ###############################lemmas################################
    # Weight_of_edge(i,j) = No.occurencies_of_i * No.occurencies_of_j   #
    #####################################################################


# using as an input the terms and the term frequency it creates the adjacency matrix of the graph
# in the main diagon we have the Win of each node of the graph and by the sum of each colume
# except the element of the diagon  is the  Wout of each node
# For more info see LEMMA 1 and LEMMA 2 of P: A graph based extension for the Set-Based Model, A: Doukas-Makris
def CreateAdjMatrixFromInvIndex(terms, tf):
    # print("Adj_Matrix = %d * %d " % (len(terms), len(tf)))
    rows = numpy.array(tf)
    row = numpy.transpose(rows.reshape(1, len(rows)))
    col = numpy.transpose(rows.reshape(len(rows), 1))
    adj_matrix = numpy.array(numpy.dot(row, col))
    #fullsize = rows.size + row.size + col.size + adj_matrix.size
    #print(fullsize / 1024 / 1024)
    for i in range(len(adj_matrix)):
        for j in range(len(adj_matrix)):
            if i == j:
                adj_matrix[i][j] = rows[i] * (rows[i] + 1) * 0.5
    # print(adj_matrix)
    del row, rows, col
    return (adj_matrix)

    ################################################################################################
    # For each node we calculate the sum of the elements of the respective row or colum of its index#
    # as its degree                                                                                 #
    ################################################################################################


# computes the degree of every node using adj matrix
def Woutdegree(mat):
    list_of_degrees = numpy.sum(mat, axis=0)
    list_of_degrees = numpy.asarray(list_of_degrees)
    id = []
    # print(list_of_degrees)
    # print(numpy.size(list_of_degrees))
    for k in range(numpy.size(list_of_degrees)):
        id.append(k)
        list_of_degrees[k] -= mat[k][k]
    list_of_degrees.tolist()
    return list_of_degrees, id


def sortByDegree(val):
    return val[0]

