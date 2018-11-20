from conkit.applications import CCMpredCommandline
import subprocess
import random
import numpy as np

distance = np.empty((71,71), dtype=object)

for i in range(len(distance)):
    for j in range(len(distance)):
        distance[i,j]=np.random.uniform(3, 15)

#Parsing aln file et Reindexing
with open("1atzA.aln", "r") as aln_file:
    query = aln_file.readline().split('\n')[0]

    #repport index of gaps
    i_gap =[i for i, e in enumerate(query) if e == "-"]

    #Calculs des score de co-ev puis generation de la matrice top_ouputs
    ccmpred_cline = CCMpredCommandline(
        cmd ='./bin/CCMpred/bin/ccmpred', alnfile= "1atzA.aln", matfile= "contact.mat"
    )
    ccmpred_cline()
    #default_value for the tops ~ 30
    subprocess.call(
        "./bin/CCMpred/scripts/top_couplings.py contact.mat > top_output.mat",
        shell=True
    )
#parsing top_coupling
index = 1
dic_top_score= {}
with open("top_output.mat", "r") as file:
    file.readline()
    for line in file:
        aa1 = float(line[0:2])
        aa2 = float(line[3:5])
        conf = float(line[6:19])
        if aa1 not in i_gap and aa2 not in i_gap:
            dic_top_score[index] = [aa1-len([gap for gap in i_gap if gap < aa1]),aa2-len([gap for gap in i_gap if gap < aa1]), conf]
            index +=1

TP = 0
for i, val in dic_top_score.items():
    pos_aa1 = int(val[0])
    pos_aa2 = int(val[1])
        if distance[pos_aa1,pos_aa2] < 8.5:
            TP +=1
score_co_evolution = TP/len(dic_top_score)

print(score_co_evolt)





distance = np.empty((query_size, query_size), dtype=object)


    

