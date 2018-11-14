#parser fasta -> give fasta sequence
sequences = []
with open('A.fasta', 'r') as fin:
    sequence = ''
    for line in fin:
        if line.startswith('>'):
            sequences.append(sequence)
            sequence = ''
        else:
            sequence += line.strip()
len(sequence)

#parser query from aln
with open ('1atzA.aln','r') as file:
    data = file.read()
    query = data.split('\n', 1)[0]
#list of comprehention = repport indice of -
i_gap =[i for i, e in enumerate(query) if e == "-"]



for aa in range(len i_
        if i > i_gap[pos_gap]
        i = i
    print(i_gap[0])


#parser le top_coupling
dic_top_score = {}
with open("top_output.mat", "r") as file:
     index = 1
     file.readline()
     for line in file:
         dic_top_score[index] = [float(line[0:2]), float(line[3:5]), float(line[6:19])]
         index += 1
print(dic_top_score)
#dic_top_score = dic_top_score.pop("0", None) 

for i, val in dic_top_score.items():
    while val[1] < i_gap[0]:
        val[1] = val[1]
    for pos_gap in range(len(i_gap)-1):
        while val[1] > i_gap[pos_gap] and val[1] < i_gap[pos_gap+1]:
            val[1]= val[1]+1
            
 




with open("top_output.mat", "r") as file:
	for line in file :
		
	data = file.read()
	data2 = data[16:]
	values = data. split('\n')

data = data.replace('\n',',')
data =data.replace('\t',',')

data2.split()
	
	
	


for top in range(1,30):
    dic["top"=top]=



#remplir dico [res1+res2] =

energy[i, j] = dope_dict[row_res.name + col_res.name][interval_index]

    dope_dict = {}
    with open(, "r") as file:
        for line in file:
            # get the line with C-alpha for both amino acids
            if line[4:6] == "CA" and line[11:13] == "CA":
                res_1 = seq1(line[0:3])
                res_2 = seq1(line[7:10])
                dope_dict[res_1+res_2] = np.fromstring(line[14:-1], dtype=float, sep=" ")
return dope_dict

matrice distance =

query = ""
with open("1l2y.fasta", "r") as file:
    for line in file:
        if line != "^\>" + "\n"
            query += str(line)
            print(query)
query = "".join(query)
r
query.count('-')
len(re.findall("a", my_string))

header = re.search(r'^>\w+', line)
