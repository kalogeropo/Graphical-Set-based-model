import os
import shutil
import xlsxwriter

""""IMPORTANT NOTICE: depending on the testing and the option the main function arguments needs to be configured properly
THATS depends on the testing which at the moment is on board. FOR EXAMPLE:
os.system('python main.py 700 0.2 1 1 3 200 0.091') would result in an error because the current main function is configured on"""

#creates the initial structure of the model
core_path = "C:/Users/nrk_pavilion/PycharmProjects/Graphical-Set-based-model"
print(core_path)

core_path = os.getcwd()

print(core_path)
#exit(10)

def init_function(test_path):
    #create txtfiles: the collection storage directory
    temp = [test_path, "/txtfiles"]
    path = "".join(temp)
    try:
        os.mkdir(path)
    except:
        print("EXISTS")
    #create figures: the results files
    temp =[test_path, "/figures"]
    path = "".join(temp)
    try:
        os.mkdir(path)
    except:
        print("EXISTS")
    #create backup: usefull directory to store items collection not code important though
    temp=[test_path, "/backup"]
    path = "".join(temp)
    try:
        os.mkdir(path)
    except:
        print("EXISTS")

    return 0

if "txtfiles" not in os.listdir(core_path):
    init_function(core_path)
    print("Init_process_starts")
else:
    print("TESTING!")

def PrepareForNextRound(counter,test_str):
    path = "figures"
    namelist = [str(counter)," ", test_str]
    name = "".join(namelist)
    print(name)
    dirpath = os.path.join(path, name)
    print(dirpath)
    try:
        os.mkdir("figures/allq")
    except:
        print("EXISTS")
    try:
        os.mkdir(dirpath)
    except:
        print("EXISTS")

    targets = ["DotSplit.dat","NegMain.dat","invertedindex.dat","PerSplit.dat",
               "docinfo.dat","ConstantWindFile.dat","SenParConWind.dat",
                "densfile.dat","newfile.dat","debuglog.dat","invertedindex.dat","CoreRankfile.dat","new_res.xlsx"]
    for file in os.listdir(core_path):
        if file in targets:
            try:
                shutil.copy2(file, dirpath)
            except FileNotFoundError:
                print(file + "NOT EXists")
            try:
                os.remove(file)
            except FileNotFoundError:
                print(file + "NOT EXists to delete")

    workbook = xlsxwriter.Workbook('new_res.xlsx')
    workbook.close()
    print("Copy paste your PARSED collection in TXTFILES and start over")
    exit(0)


counter = 0

#name of main: sys.args[1] =  h - sys.args[2] =  p - sys.args[3] = menu 1 - sys.args[4] = choice of menu 1 -
#sys.args[5] = constant window size, sys.args[6] = paragraph size - sys.args[7] = percentage window.
#example:  main_v2.py 700 0.2 1 1 3 200 0.091')

#os.system('python main_v3.py 700 0.2 1 1 20 200 0.091')
#os.system('python main_v3.py 700 0.2 1 2 50 200 0.091')
#os.system('python main_v3.py 700 0.2 1 3 0 200 0.091')
#os.system('python main_v3.py 700 0.2 1 4 3 200 0.091')
#os.system('python main_v3.py 700 0.2 1 5 3 200 0.091')
#os.system('python main_v3.py 700 0.2 1 6 3 200 0.091')


os.system('python main_v3.py 700 0.2 2')

#name for file in figures directory
testname = "test3"


PrepareForNextRound(counter,testname)
counter+=1

