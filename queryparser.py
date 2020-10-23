relevant = []
q=[]
sizes = []
#files = 'my_q','new_Qs.txt' IMPORTANT:use the file my_q
with open('my_q', 'r') as fd:
    line =fd.readline()
    sline = line.split()
    print(line,sline)
    while line and line!='QN' and line !="END":
        line = fd.readline()
        sline = line.split()
        reltemp=[]
        punctuation = ['(', ')', '?', ':', ';', ',', '.', '!', '/', '"', "'"]
        if len(sline)>0 and sline[0]=='QU':
            print("for Query %s ||||"%line)
            temp = line.replace('?',"")
            temp = temp.replace('(', "")
            temp = temp.replace(')', "")
            temp = temp.replace('\n','')
            temp = temp.replace('QU', '')
            for item in punctuation:
                temp = temp.replace(item, '')
            q.append(temp)
        elif len(sline)>0 and sline[0]=='RD':#read all relevant lines

            while line!='\n'and line!='END':
                print(line)
                if sline[0]=='RD':
                    reltemp += sline[1::2]
                else:
                    reltemp+=sline[::2]
                line = fd.readline()
                sline=line.split()

        if len(reltemp)!=0:
            relevant.append(reltemp)
            sizes.append(len(reltemp))
print('Num Of Queries =  ',len(q))
print(q)
print(relevant)
print(sizes)

