score_co_evolt = 0
#Parsing aln file et Reindexing
with open(self.aln, "r") as aln_file:
    query = aln_file.readline().split('\n')[0]

    #Reindexing
    list_index_pos_nongaps = []
    list_index_pos_nongaps.append([i for i, e in enumerate(query) if e != "-"]) #ajoute les position "ali" des res

    #Calculs des score de co-ev puis generation de la matrice top_ouputs
    ccmpred_cline = CCMpredCommandline(
        cmd ='./bin/CCMpred/bin/ccmpred', alnfile= self.aln, matfile= "contact.mat"
    )
    ccmpred_cline()
    #default_value for the tops ~ 30
    subprocess.call(
        "./bin/CCMpred/scripts/top_couplings.py contact.mat > top_output.mat",
        shell=True
    )
#Lecture top_outputs et some verifs
with open("top_output.mat", "r") as top_output_file:
    lines = top_output_file.read()
    for line in lines.split('\n'):
        #verif si les tops ne correspondent pas des gaps_positions
        if line[0:2] in list_index_pos_nongaps and line[3:5] in list_index_pos_nongaps:
            #checker dans la matrice de distance si les positions < 8angströms
                pos_i = list_index_pos_nongaps.index(line[0:2])
                pos_j = list_index_pos_nongaps.index(line[3:5])
                #incrementer score_co_evolution si verifié
                if distance_matrix[pos_i][pos_j] < 8:
                    score_co_evolt +=1
        ligne = top_output_file.readline()
print(score_co_evolt)
return score_co_evolt
