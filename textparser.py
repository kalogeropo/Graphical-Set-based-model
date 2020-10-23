import nltk
#from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords

#nltk.data.path.append("E:/python-projects/nltk")

import os
import string
from nltk.corpus import stopwords

path = "/home/nik/Desktop/texts/"
dest_path = "/home/nik/Desktop/texts/parsed texts"

def List_to_write(filename,title_list,abstract_list):
    print("Writing on file ",filename)
    print("TITLE  ",title_list)
    print("ABSTRACT ",abstract_list)
    #TI and AB and EX removal
    firstfilter=["TI","AB","EX"]

    list_of_words = title_list + abstract_list
    start_size = len(list_of_words)
    if start_size==0:
        print("POTENTIAL ERROR HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    #print(start_size)
    # remove punctation
    table = str.maketrans('', '', string.punctuation)
    stripped = [w.translate(table) for w in list_of_words]
    #print(len(stripped))

    #words = [word for word in stripped if word.isalpha()]
    # stopword filtering
    stop_words = set(stopwords.words('english'))
    #print(stop_words)
    #print("the stopwords are %d" % len(stop_words))
    list_of_words = [w for w in stripped if not w in firstfilter]
    list_of_words = [word.upper() for word in list_of_words]
    list_of_words_size = len(list_of_words)
    #print(start_size - list_of_words_size)
    rewrite_file(filename,list_of_words)
    return list_of_words

def rewrite_file(filename,list_of_words_to_write):
    #print('Re-writing data on %s......../texts/parsed texts/' %filename)
    filename="/home/nik/Desktop/texts/parsed texts/"+filename #change that!!!!!!!
    fd = open(filename,'w')
    for w in list_of_words_to_write:
        fd.write("%s\n"% w)
    #print("done")
    return 0
print("start")
print(path)

filecount = 0
count = 0
for file in os.listdir(path):
    filecount+=1
    full_path = path + file
    if os.path.isfile(full_path):
        #print(file)
        fd = open(full_path, 'rt')
        text = fd.readline()
        #print("PN: ",text)
        RN= fd.readline()
        rn=RN.split()
        #print("RN: ",RN)
        filename = rn[1]
        #print(filename)
        line=fd.readline()
        sline=line.split()
        while(line!="END"):
            if(line!="\n" and line!=""):
                #print(filename)
                #handling AN multilined
                if len(sline)>0 and sline[0] =="AN":
                    #print(sline)
                    line = fd.readline()
                    sline = line.split()
                while len(sline)>0 and sline[0] != "AU":
                    #print(sline)
                    line =  fd.readline()
                    sline=line.split()

                # handling AU multilined
                if len(sline)>0 and sline[0]=="AU":
                    #print(sline)
                    line = fd.readline()
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "TI":
                    #print(sline)
                    line = fd.readline()
                    sline = line.split()
                title = []
                    # handling TI multilined
                if len(sline)>0 and sline[0] == "TI":
                    #print(sline)
                    title.extend(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "SO":
                    #print(sline)
                    title.extend(sline)
                    line = fd.readline(500)
                    sline = line.split()
                #print(title)
                # handling SO multilined
                if len(sline)>0 and sline[0] == "SO":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "MJ":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                # handling SO multilined
                if len(sline)>0 and sline[0] == "MJ":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "MN":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()

                if len(sline)>0 and sline[0] == "MN":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and (sline[0] != "AB" and sline[0] != "EX" ):
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                abstract = []
                if len(sline)>0 and (sline[0] == "AB" or sline[0] == "EX"):
                    abstract.extend(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "RF" :
                    abstract.extend(sline)
                    line = fd.readline(500)
                    sline = line.split()
                    #abstract.extend(sline)

                #print(abstract)
                if len(sline)>0 and sline[0] == "RF":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0 and sline[0] != "CT":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                if len(sline)>0 and sline[0] == "CT":
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                while len(sline) > 0:
                    #print(sline)
                    line = fd.readline(500)
                    sline = line.split()
                List_to_write(filename,title,abstract)

            else:
                count+=1
                text = fd.readline()
                #print("PN: ", text)
                RN = fd.readline()
                #print("RN: ", RN)
                if len(RN)>0:
                    rn = RN.split()
                    filename = rn[1]
                    #print(filename)
                line = fd.readline()
                sline = line.split()
                if len(sline)==0 or sline[0]=="END":
                    print(count)
                    break
        print("END")
print(count)