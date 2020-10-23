import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords

nltk.data.path.append("E:/python-projects/nltk")

import os
import string
from nltk.corpus import stopwords


def load_read_textfile(filename):
    print("loading %s ...." %filename)
    fd = open(filename,'rt')
    text =fd.read()
    fd.close()
    return text
def rewrite_file(filename,list_of_words_to_write):
    print('Re-writing data on %s......' %filename)
    fd = open(filename,'w')
    for w in list_of_words_to_write:
        fd.write("%s\n"% w)
    return 0
#debug perpose
def write_preproccess_stats(filename , size_before ,size_after):
    debugfile = 'preproccessed_files\debug.dat'
    if os.path.exists(debugfile):
        with open(debugfile ,"a") as df:
            df.write("***********************************\n")
            df.write("At file %s the no of words before preproccess was %d and after is %d WITH SIZE DIFFERENCE = %d \n" %(filename,size_before,size_after, size_before-size_after) )
    else:
        with open(debugfile ,"w") as df:
            df.write("***********************************\n")
            df.write("At file %s the no of words before preproccess was %d and after is %d WITH SIZE DIFFERENCE = %d \n" %(filename,size_before,size_after, size_before-size_after) )
    return 0

pathfile = "E:/python-projects/k_core_decomp/preproccessed_files/"
for file in os.listdir(pathfile):
    if file.isdigit():
        print(file)
        full_path = pathfile+file
        a = load_read_textfile(full_path)
        #Word spliting using NLTK
        tokens = word_tokenize(a)
        size = len(tokens)
        #remove punctation
        table = str.maketrans('', '', string.punctuation)
        stripped = [w.translate(table) for w in tokens]
        #remove numbers
        words = [word for word in stripped if word.isalpha()]

        test1 = len(words)
        #stopword filtering
        stop_words = set(stopwords.words('english'))
        print("the stopwords are %d" %len(stop_words))
        words = [w for w in words if not w in stop_words]
        words = [word.upper() for word in words]
        final_word_list_size = len(words)
        print(test1 - final_word_list_size)
        #rewrite
        rewrite_file(full_path,words)
        #debug
        write_preproccess_stats(full_path,size,final_word_list_size)





