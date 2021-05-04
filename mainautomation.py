import os
import shutil
""""IMPORTANT NOTICE: depending on the testing and the option the main function arguments needs to be configured properly
THATS depends on the testing which at the moment is on board. FOR EXAMPLE:
os.system('python main.py 700 0.2 1 1 3 200 0.091') would result in an error because the current main function is configured on"""
def PrepareForNextRound(counter,test_str):
    path = "figures"
    namelist = [str(counter)," ", test_str]
    name = "".join(namelist)
    dirpath = os.path.join(path, name)
    os.rename("figures/allq",dirpath)
    os.chmod(dirpath, 0o444)
    os.mkdir("figures/allq")
    shutil.copy2("C:/Users/nrk_pavilion/PycharmProjects/Graphical-Set-based-model/new_res.xlsx", dirpath)
    targets = ["DotSplit.dat","NegMain.dat","invertedindex.dat","PerSplit.dat","docinfo.dat","ConstantWindFile.dat","SenParConWind.dat"]
    for file in os.listdir("C:/Users/nrk_pavilion/PycharmProjects/Graphical-Set-based-model"):
        if file in targets:
            shutil.copy2(file, dirpath)

counter = 0

#os.system('python main.py 700 0.2 1 1 3 200 0.091')
#os.system('python main.py 700 0.2 1 2 3 200 0.091')
#os.system('python main.py 700 0.2 1 3 0 200 0.091')
#os.system('python main.py 700 0.2 1 4 3 200 0.091')
os.system('python main.py 700 0.2 2')

PrepareForNextRound(counter,"700 0.2 3 200 0.091 noprune")
counter+=1
#os.remove("DotSplit.dat")
#os.remove("NegMain.dat")
#os.remove("invertedindex.dat")
#os.remove("PerSplit.dat")
#os.remove("ConstantWindFile.dat")
#os.remove("SenParConWind.dat")
#os.remove("docinfo.dat")

