import os
import shutil
""""IMPORTANT NOTICE: depending on the testing and the option the main function arguments needs to be configured properly
THATS depends on the testing which at the moment is on board. FOR EXAMPLE:
os.system('python main.py 700 0.2 1 1 3 200 0.091') would result in an error because the current main function is configured on"""
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
    #os.rename("figures/allq",dirpath)
    #os.chmod(dirpath, 0o444)
    #os.mkdir("figures/allq")
    targets = ["DotSplit.dat","NegMain.dat","invertedindex.dat","PerSplit.dat",
               "docinfo.dat","ConstantWindFile.dat","SenParConWind.dat",
                "densfile.dat","newfile.dat","debuglog.dat","invertedindex.dat","CoreRankfile.dat","new_res.xlsx"]
    for file in os.listdir("C:/Users/nrk_pavilion/PycharmProjects/Graphical-Set-based-model"):
        if file in targets:
            try:
                shutil.copy2(file, dirpath)
            except FileNotFoundError:
                print(file + "NOT EXists")
            try:
                os.remove(file)
            except FileNotFoundError:
                print(file + "NOT EXists to delete")


counter = 0

"""os.system('python main.py 700 0.2 1 1 3 200 0.091')
os.system('python main.py 700 0.2 1 2 3 200 0.091')
os.system('python main.py 700 0.2 1 3 0 200 0.091')
os.system('python main.py 700 0.2 1 4 3 200 0.091')
os.system('python main.py 700 0.2 2')

testname = "700 0.2 3 200 0.091 noprune"

PrepareForNextRound(counter,"testname")
counter+=1

"""
#name of main: sys.args[1] =  h - sys.args[2] =  p - sys.args[3] = menu 1 - sys.args[4] = choice of menu 1 -
#sys.args[5] = constant window size, sys.args[6] = paragraph size - sys.args[7] = percentage window.
#main_v2.py 700 0.2 1 1 3 200 0.091')

#os.system('python main_v2.py 700 0.2 1 1 3 200 0.091')
#os.system('python main_v2.py 700 0.2 1 2 3 200 0.091')
#os.system('python main_v2.py 700 0.2 1 3 0 200 0.091')
#os.system('python main_v2.py 700 0.2 1 4 3 200 0.091')
#os.system('python main_v2.py 700 0.2 1 5 3 200 0.091')

#os.system('python main_v2.py 700 0.2 2')

testname = "testname"


PrepareForNextRound(counter,testname)
counter+=1

