import os
import shutil



#h2000ee12avn
#executes the command python3| main.py| h| S| menu_one answer| menu_two_anwser
#where main.py has been modified to get the anwsers from arguments and not as a keyboard input!


#for windows 10 os: os.system('python main.py H P MENU_OPTION_1 MENU_OPTION_2')
#at quering proccess menu_option_1 =2 and there is no menu_option_2

#os.system('python main.py 2000 1.2 1 3')
#os.system('python main.py 2000 1.2 1 4')
#os.system('python main.py 2000 1.2 2')

#os.system('python main.py 700 0.2 1 1')
#os.system('python main.py 700 0.2 1 2')
#os.system('python main.py 700 0.2 1 3')
#os.system('python main.py 700 0.2 1 4')
os.system('python main.py 700 0.2 2')


#for linux os: os.system('python3 main.py H P MENU_OPTION_1 MENU_OPTION_2')
#at quering proccess menu_option_1 =2 and there is no menu_option_2

#os.system('python3 main.py 700 0.3 1 1')
#os.system('python3 main.py 700 0.3 1 2')
#os.system('python3 main.py 700 0.3 1 3')
#os.system('python3 main.py 700 0.3 1 4')
#os.system('python3 main.py 700 0.3 2')

#os.system('python3 main.py 700 0.4 1 1')
#os.system('python3 main.py 700 0.4 1 2')
#os.system('python3 main.py 700 0.4 1 3')
#os.system('python3 main.py 700 0.4 1 4')
#os.system('python3 main.py 700 0.4 2')

#os.system('python3 main.py 700 0.6 1 1')
#os.system('python3 main.py 700 0.6 1 2')
#os.system('python3 main.py 700 0.6 1 3')
#os.system('python3 main.py 700 0.6 1 4')
#os.system('python3 main.py 700 0.6 2')


#os.system('python3 main.py 700 1.2 1 1')
#os.system('python3 main.py 700 1.2 1 2')
#os.system('python3 main.py 700 1.2 1 3')
#os.system('python3 main.py 700 1.2 1 4')
#os.system('python3 main.py 700 1.2 2')

#-----------------------------------------------------

#os.system('python3 main.py 1500 0.2 1 1')
#os.system('python3 main.py 1500 0.2 1 2')
#os.system('python3 main.py 1500 0.2 1 3')
#os.system('python3 main.py 1500 0.2 1 4')
#os.system('python3 main.py 1500 0.2 2')


#os.system('python3 main.py 1500 0.3 1 1')
#os.system('python3 main.py 1500 0.3 1 2')
#os.system('python3 main.py 1500 0.3 1 3')
#os.system('python3 main.py 1500 0.3 1 4')
#os.system('python3 main.py 1500 0.3 2')

#os.system('python3 main.py 1500 0.4 1 1')
#os.system('python3 main.py 1500 0.4 1 2')
#os.system('python3 main.py 1500 0.4 1 3')
#os.system('python3 main.py 1500 0.4 1 4')
#os.system('python3 main.py 1500 0.4 2')

#os.system('python3 main.py 1500 0.6 1 1')
#os.system('python3 main.py 1500 0.6 1 2')
#os.system('python3 main.py 1500 0.6 1 3')
#os.system('python3 main.py 1500 0.6 1 4')
#os.system('python3 main.py 1500 0.6 2')


#os.system('python3 main.py 1500 1.2 1 1')
#os.system('python3 main.py 1500 1.2 1 2')
#os.system('python3 main.py 1500 1.2 1 3')
#os.system('python3 main.py 1500 1.2 1 4')
#os.system('python3 main.py 1500 1.2 2')

#------------------------des kai ena h 5000 ee 03 av

#os.system('python3 main.py 5000 0.3 1 1')
#os.system('python3 main.py 5000 0.3 1 2')
#os.system('python3 main.py 5000 0.3 1 3')
#os.system('python3 main.py 5000 0.3 1 4')
#os.system('python3 main.py 5000 0.3 2')

#------------------------------- h 1000

#os.system('python3 main.py 1000 0.3 1 1')
#os.system('python3 main.py 1000 0.3 1 2')
#os.system('python3 main.py 1000 0.3 1 3')
#os.system('python3 main.py 1000 0.3 1 4')
#os.system('python3 main.py 1000 0.3 2')

#---------------------------h 800 ->h200 (0.3 0.2) ->h100
os.system('python3 main.py 100 0.3 1 1')
os.system('python3 main.py 100 0.3 1 2')
os.system('python3 main.py 100 0.3 1 3')
os.system('python3 main.py 100 0.3 1 4')
os.system('python3 main.py 100 0.3 2')


h=100
S= "03"
src = "figures/allq"
path = "figures"

namelist = ['h',str(h),'ee',str(S),"avn"]
name = "".join(namelist)
dirpath = os.path.join(path, name)
print(dirpath)

shutil.move(src, dirpath)
targets = ["CoreRankfile.dat","newfile.dat","invertedindex.dat","densfile.dat","docinfo.dat"]

os.mkdir(src)
for file in os.listdir("C:/Users/nrkal/PycharmProjects/Graphical-Set-based-model"):
    if file in targets:
        shutil.move(file, dirpath)

