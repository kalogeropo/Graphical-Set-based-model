import math
from math import log

import numpy
import sys

import collections
import matplotlib.pyplot as plt
import networkx as nx

import string

import nltk
from nltk.corpus import stopwords

translator = str.maketrans('', '', string.punctuation)

from operator import itemgetter

from networkx import core_number, k_core
import csv

postinglist = []
docinfo = []
docs_without_main_core = []

# buggy na xrisimopoii8ei i proepeksergasia pou exei dimiourgi8ei san ksexwristo .py

def preproccess(file):
    with open(file, 'r') as fd:
        text = fd.read().split()
    fd.close()
    with open(file, 'w'):
        pass
    with open(file, 'a') as fd:
        fd.write('')
        for term in text:
            term = term.translate(translator)
            term = term.upper()
            if len(term) != 1:
                fd.write("%s \n" % term)
    fd.close()
    return 1


# todo: more efficient way to calculate max length of path it doesnt work on realistic scale (CANT BE DONE BECAUSE THE COMPLEXITY)
# finds the maximum distance which exists in the graph
def findMaxDistance(gr, adjmatrix):
    maxlist = []
    for adi in range(len(adjmatrix)):

        for adj in range(len(adjmatrix)):

            if adj < adi:
                # cut down the number of calculated paths path(i,j) == path(j,i) so we need only the upper
                # or lower  tri of adj_matrix
                path = list(nx.shortest_simple_paths(gr, adj, adi, weight='weight'))
                print('Longest Path for (%d,%d) is: ' % (adj, adi))
                print(path[-1])
                i = 0
                weightsum = 0
                for item in range(len(path[-1]) - 1):
                    indexI = path[-1][i]
                    indexIpp = path[-1][i + 1]
                    i += 1
                    weightsum += adjmatrix[indexI][indexIpp]
                maxlist.append(weightsum)
    print(len(maxlist))
    return max(maxlist)


# computes the cosine similarity of 2 vectors
def cos_sim(a, b):
    a = numpy.asarray(a)
    b = numpy.asarray(b)
    if numpy.all(a == 0) or numpy.all(b == 0):
        ret = 0
    else:
        dot_p = numpy.dot(a, b)
        normA = numpy.linalg.norm(a)
        normB = numpy.linalg.norm(b)
        ret = dot_p / (normA * normB)
    return ret


# finds the maximum and the minimum similarity between the nodes of the graph
def node_simi(adjmatrix):
    max = 0
    min = 1
    for adi in range(len(adjmatrix)):
        for adj in range(len(adjmatrix)):
            if adj < adi:
                temp = cos_sim(adjmatrix[adi], adjmatrix[adj])
                if temp > max:
                    max = temp
                if temp < min:
                    min = temp
    # print(max, min)
    return max, min


# deletes by re drawing the graph edges of the graph given a minimum similarity
def pruneGraphbySimilarity(aMatrix, pers, minsim, termName):
    g = nx.Graph()
    for i in range(len(aMatrix)):
        for j in range(len(aMatrix)):
            if i > j:
                if cos_sim(aMatrix[i], aMatrix[j]) > minsim - ((minsim * pers)):
                    g.add_node(termName[i], id=i)
                    g.add_node(termName[j], id=j)
                    g.add_edge(i, j, weight=aMatrix[i][j])
    Matrix = nx.adjacency_matrix(g)
    # print(Matrix)
    return g, Matrix.todense()


def graphUsingAdjMatrix(adjmatrix, termlist, *args, **kwargs):
    gr = nx.Graph()
    filename = kwargs.get('filename', None)
    if not filename:
        filename = 'Name not found!'  # used when i want to visualize graphs with name

    for i in range(0, len(adjmatrix)):
        gr.add_node(i, term=termlist[i])
        for j in range(len(adjmatrix)):
            if i > j:
                gr.add_edge(i, j, weight=adjmatrix[i][j])
    # graphToPng(gr,filename = filename)
    return gr


# ------------------Graph visualization---------------

def getGraphStats(graph, filename, graphPng, degreePng):
    if nx.is_connected(graph):
        print("IT IS CONNECTED")
    name = filename[10:]
    if graphPng:
        graphToPng(graph=graph, filename=str(name))
    if degreePng:
        plot_degree_dist(graph=graph, filename=str(name))


def graphToPng(graph, *args, **kwargs):
    options = {
        'node_color': 'yellow',
        'node_size': 50,
        'linewidths': 0,
        'width': 0.1,
        'font_size': 8,
    }
    filename = kwargs.get('filename', None)
    if not filename:
        filename = 'Union graph'
    plt.figure(filename, figsize=(17, 8))
    plt.suptitle(filename)
    pos_nodes = nx.circular_layout(graph)
    nx.draw(graph, pos_nodes, with_labels=True, **options)
    pos_atrs = {}
    for node, coords in pos_nodes.items():
        pos_atrs[node] = (coords[0], coords[1] + 0.01)

    node_attrs = nx.get_node_attributes(graph, 'term')
    cus_node_att = {}
    for node, attr in node_attrs.items():
        cus_node_att[node] = attr
    nx.draw_networkx_labels(graph, pos_atrs, labels=cus_node_att, font_color='red', font_size=8)

    labels = nx.get_edge_attributes(graph, 'weight')
    nx.draw_networkx_edge_labels(graph, pos_nodes, edge_labels=labels)
    # plt.show()
    plt.savefig('figures/allq/' + str(filename) + '.png', format="PNG", dpi=600)


def plot_degree_dist(graph, *args, **kwargs):
    filename = kwargs.get('filename', None)
    degree_sequence = sorted([d for n, d in graph.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color="b")

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    plt.setp(ax.get_xticklabels(), rotation=90, horizontalalignment='right', fontsize=3)
    ax.set_xticklabels(deg)

    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = graph.subgraph(sorted(nx.connected_components(graph), key=len, reverse=True)[0])
    pos = nx.spring_layout(graph)
    plt.axis("off")
    nx.draw_networkx_nodes(graph, pos, node_size=20)
    nx.draw_networkx_edges(graph, pos, alpha=0.4)

    plt.savefig('figures/allq/' + str(filename) + '_degree.png', format="PNG", dpi=600)
#TODO: paths
def stopwordsStats(kcore,term_list,file):

    stopword = stopwords.words('english')
    stopword_list = [x.upper() for x in stopword]
    print(stopword)
    stopword_count = 0
    stopwords_in_file = 0
    print(len(kcore.nodes))
    print(len(term_list))

    for i in kcore.nodes:
        #print(term_list[i])
        if term_list[i] in stopword_list:
            stopword_count += 1
            #print(stopword_count)
    for i in term_list:
        if i in stopword_list:
            stopwords_in_file += 1
            #print(i)
    print(stopword_count)
    print(stopwords_in_file)
    stopwords_in_file_per = float(stopword_count/stopwords_in_file)
    stopwords_per = float(stopword_count/len(kcore.nodes))
    #print(stopwords_per)
    fw=open('stopwords_stats.txt','a')
    string_to_write = "File " + str(file) + " stopwords in kcore percentage : " + str(stopwords_per) + " and stopwords percentage in file: "+ str(stopwords_in_file_per) + "\n"
    fw.write(string_to_write)
    fw.close()
    fw=open('stopwords_kcore_stats.txt','a')
    fw.write(str(stopwords_per))
    fw.close()
    fw=open('stopwords_file_stats.txt','a')
    fw.write(str(stopwords_in_file_per))
    fw.close()


# -----------Union Graph to inverted index-------------
def graphToIndex(id, terms, calc_term_w, plist, *args, **kwargs):
    filename = kwargs.get('filename', None)
    if not filename:
        filename = 'inverted index.dat'
    f = open(filename, "a+")
    data = ','.join(
        [str(i) for i in plist])  # join list to a string so we can write it in the inv index and load it with ease
    f.write('%d;%s;%f;%s;\n' % (id, terms, calc_term_w, data))
    f.close()
    return 1



# calculating the weight and write the inverted index file using graphToIndex method
# NO USE
def w_and_write_to_file(listofdeg, Umatrix, collection_terms, union_graph_termlist_id, collection_term_freq):
    print('here')
    for i in range(len(listofdeg[0])):
        Wout = listofdeg[0][i]
        Win = collection_term_freq[i]
        nbrs = numpy.count_nonzero(Umatrix[i])
        VarA = 1
        VarB = 10
        Alog = 1 + VarA * ((Wout / (nbrs + 1)) / (Win + 1))
        Blog = 1 + VarB * (1 / (nbrs + 1))
        temp = log(Alog) * log(Blog)
        print(temp)

        indexofw = postinglist.index(collection_terms[i])  # maybe not the best way of implementing the
        graphToIndex(union_graph_termlist_id[i], collection_terms[i], temp, postinglist[indexofw + 1])
    return 1


def w_and_write_to_filev2(wout, collection_terms, union_graph_termlist_id, collection_term_freq, postinglist, file):
    # wout |[[term][wout][neibours]]\
    print('Calculating weights and create inveted index')
    print(file)
    for entry in collection_terms:
        term = entry
        id = collection_terms[entry]
        win = collection_term_freq[id]
        for sublist in wout:
            if term in sublist[0]:
                # print('here')
                Wo = sublist[1]
                nbrs = sublist[2]
                break
        VarA = 1
        VarB = 10
        Alog = 1 + VarA * ((Wo / (nbrs + 1)) / (win + 1))
        Blog = 1 + VarB * (1 / (nbrs + 1))
        temp = log(Alog) * log(Blog)
        indexofwordinlist = 2 * id + 1
        graphToIndex(id, term, temp, postinglist[indexofwordinlist], filename=file)
    return 1


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def union(lista, listb):
    c = []
    for i in lista + listb:
        if i not in c:
            c.append(i)
    return c


# generate new k+1 itemsets !!!reference to: Fast Algorithms for Mining Association Rules by Argawal
# a nice change would be the closed sets idea => less sets than frequent termsets
def apriori(l1, minfreq):
    final_list = []
    final_list.append(l1)
    k = 2
    l = l1
    #print('=========Generating frequent sets ============')
    while (l != []):
        c = apriori_gen(l)
        l = apriori_prune(c, minfreq)
        # print('Frequent  %d-termset is: %s'%(k,l))
        # print(len(l))
        final_list.append(l)
        # print('====for k = %d the l list is' %k )
        # print(l)
        k += 1
    return final_list


def apriori_gen(itemset):
    candidate = []
    ck = []
    texts = []
    length = len(itemset)
    for i in range(length):
        #[[leksi,emfanisi],[.... , ...],[... , ....]]
        ele = itemset[i][0]
        for j in range(i + 1, length):
            ele1 = itemset[j][0]
            # print(ele, ele1)
            if ele[0:len(ele) - 1] == ele1[0:len(ele1) - 1]:  # and ele1 != ele:
                texts.append(intersection(itemset[i][1], itemset[j][1]))
                candidate.append([union(ele, ele1), intersection(itemset[i][1], itemset[j][1])])
    return candidate


def apriori_prune(termsets_list, min_support):
    prunedlist = []
    for j in termsets_list:
        if len(j[1]) > min_support:
            # print('-----------')
            # print(j[0],len(j[1]))
            # print('-----------')
            prunedlist.append([j[0], j[1]])

    return prunedlist


def printmenu():
    # menu implementation
    print("1.create index file seperate graphs and  union graph")
    print("2.load index file and then quering \n \n")

    # x = input('Insert option: ')
    hargs = int(sys.argv[1])
    print(hargs)
    S = float(sys.argv[2])
    print(S)
    x = int(sys.argv[3])
    print(x)

    return x, hargs, S


def doc_rep(doc_vec, idf_vec, *args, **kwargs):
    args = list(args)
    if not args:
        nw = []
        for i in range(len(idf_vec)):
            nw.append(1)
    else:
        nw = args[0]
    # print(nw)
    test = numpy.zeros((len(doc_vec), len(idf_vec)))
    for i in range(len(doc_vec)):
        # print(docs[i])
        for j in range(len(idf_vec)):
            # print(doc_vec[i][j])
            if doc_vec[i][j] > 0:
                test[i][j] = (1 + log(doc_vec[i][j])) * idf_vec[j] * float(nw[j])
            else:
                test[i][j] = 0
    # with open('debuglog.dat', 'a') as fd:
    #    fd.write('doc representa \n')
    #   for doci in test:
    #        fd.write('%s \n' %str(len(doci)))
    # fd.close()
    return test


def load_inv_index(*args):
    arg = list(args)
    if not arg:
        invindex = 'inverted index.dat'
    else:
        invindex = arg[0]
    ids = []
    trms = []
    W = []
    plist = []
    with open(invindex, 'r') as csvf:
        reader = csv.reader(csvf, delimiter=";")
        for row in reader:
            if row[0] not in ids:
                ids.append(row[0])
                trms.append(row[1])
                W.append(row[2])
                plist.append(row[3].split(','))
    csvf.close()
    # print(len(ids))
    return ids, trms, W, plist


def load_doc_info(*args):
    args = list(args)
    if not args:
        docinfofile = "docinfo.dat"
    else:
        docinfofile = args[0]
    info = []
    with open(docinfofile, "r") as fh:
        lines = [line.split() for line in fh]
        for line in lines:
            if line not in info:
                info.append(line)
    return info


# input the Query Q as a list of words consisting the Query
def one_termsets(Q, trms, plist, minfreq):
    termsets = []
    One_termsets = []
    for word in Q:
        if word in trms:
            i = trms.index(word)
            doc = plist[(i + 1)]
            doc = doc[::2]
            word = [''.join(word)]
            if len(doc) > minfreq:
                One_termsets.append([word, doc])
        else:
            print('word %s has not required support or it already exists:' % word)
    return One_termsets


def fij_calculation(docinfo, final_list, plist, trms):
    doc_vectors = []
    docs = []
    # counting the number of apperances of each termset in each doc in which the set exists

    for document in docinfo:
        # print('============doc name =============')
        print(document[0])
        docs.append([document[0]])
        print(final_list)# not needed
        # k = 1
        # test is a temp list which contains the termfreq of TSets for the current doc, then we append that list to create
        # a matrix of [Docs X Termsets].
        test = []
        for itemsets in final_list:
            # print('------------------------------%d- termset is:--------------------------------'%k)
            print(itemsets)
            # k+=1
            # print(len(itemsets))
            for i in range(len(itemsets)):
                # calculating termset frequency
                if document[0] in itemsets[i][1]:
                    # option 2: na xrisimopoisw ti lista pou exw dimiourgisei apo to inverted file
                    # print('............')
                    # print(itemsets[i])
                    sum = 0
                    for term in itemsets[i][0]:
                        if term in trms:
                            termindx = trms.index(term)
                            #print(termindx)
                            #print(ids[termindx])
                            #print(trms[termindx])
                            #print(plist[termindx])
                            if document[0] in plist[termindx]:
                                docindex = plist[termindx].index(document[0])
                                print('edw:%s'%plist[termindx][docindex+1])
                                sum += int(plist[termindx][
                                               docindex + 1])  # to amesws epomeno stoixeio antistoixei ston ari8mo emfanisis tou orou TERM(i) sto Document(j)
                    test.append(sum)
                else:
                    test.append(0)
        #if len(test) == 1213:
        #    print(len(test))
        #    print('problem')
        #    exit(-2)
        doc_vectors.append(test)
    with open('debuglog.dat', 'a') as fd:
        fd.write('doc vectors \n')
        for doci in doc_vectors:
            fd.write('%s \n' % str(len(doci)))
    fd.close()
    return docs, doc_vectors


def calculate_idf(termsetsL, numofdocs):
    #print('=====calculating idfs============')
    idf_vector = []
    for ts in termsetsL:  # iterate based on the number of terms in termset
        for item in ts:  # iterate all termsets with the same number of terms in set
            Nt = len(item[1])
            N = numofdocs
            if Nt != 0:
                idf = log(1 + (N / Nt))
                idf_vector.append(idf)
            else:
                idf_vector.append(0)
                #print(item[1], len(item[1]))
    return idf_vector


# doukas weight
def calculate_termset_W(termsetsL, W, terms):
    #print("=======Calculating W ======")
    termset_W_vector = []
    for ts in termsetsL:  # iterate based on the number of terms in termset
        for item in ts:  # iterate all termsets with the same number of terms in set
            product = 1
            for term in item[0]:
                if term in terms:
                    tindx = terms.index(term)
                    weight = W[tindx]
                    product *= float(weight)
            termset_W_vector.append(product)

    return termset_W_vector


def q_D_similarities(q, docmatrix, docs):
    ret_list = []
    cnt = 0
    # debug
    # with open('debuglog.dat', 'a') as fd:
    #    fd.write('doc matrix \n')
    #   for doci in docmatrix:
    #       fd.write('%s \n' %str(len(doci)))
    # fd.close()
    for doci in docmatrix:
        # the 0 array issue is fixed inside the cosine similarity function
        try:
            temp = cos_sim(q, doci)
            ret_list.append([docs[cnt], temp])
            cnt += 1
        except ValueError:
            print(doci)
            print(len(doci))
            print(cnt)
            exit(-1)
    return ret_list


def Woutusinggraph(inputgraph):
    nd = inputgraph.nodes()
    woutlist = []
    for n in nd:
        # print(n,inputgraph.degree(n,weight= 'weight'))
        woutlist.append([n, inputgraph.degree(n, weight='weight'), len(list(inputgraph.neighbors(n)))])

    print('success')
    return woutlist


def density(A_graph):
    graph_edges = A_graph.number_of_edges()
    # print(graph_edges)
    graph_nodes = len(list(A_graph.nodes))
    # print(graph_nodes)
    dens = graph_edges / (graph_nodes * (graph_nodes - 1))
    return dens


# given points A and B it caluclates the distance of a point P from the line AB
def distance_to_line(starting_point, end_point, point):
    dist = -9999
    # spoint = (x1,y1)
    x1 = starting_point[0]
    y1 = starting_point[1]
    # end point = (x2,y2)
    x2 = end_point[0]
    y2 = end_point[1]
    # point = (x0,y0)
    print(point)
    x0 = point[0]
    y0 = point[1]
    dist = (abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y1 * x1)) / (math.sqrt(((y2 - y1) ** 2) + ((x2 - x1) ** 2)))

    return dist


def elbow(listofpoints):
    # at first we need to create a line between first and last element
    if len(listofpoints) == 1:
        bestindex = 0
    elif len(listofpoints) == 2:
        if listofpoints[0] > listofpoints[1]:
            bestindex = 0
        else:
            bestindex = 1
    elif len(listofpoints) > 2:
        # p1 the starting point of line and p2 the last point of line
        # using that we will calulate the distance of each point of our starting list
        # from the line using the known forumla
        p1 = numpy.array([listofpoints[0], 0])
        p2 = numpy.array([listofpoints[-1], (len(listofpoints) - 1)])
        distance = []
        # print(p1,p2)
        # print(listofpoints)
        pnt = []
        for point in listofpoints:
            pnt.append(point)
            pnt.append(listofpoints.index(point) + 1)
            print(pnt)
            distance.append(distance_to_line(p1, p2, pnt))
            pnt = []
        bestdistance = max(distance)
        bestindex = distance.index(bestdistance)
    return bestindex
# test method
