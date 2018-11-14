#! coding: utf-8
# Meet-U

#Script to compute fonction for getting from an query seq an top couplings co-ev scores

#Librairies imported
import Bio
import subprocess

from conkit.applications import CCMpredCommandline


def calculate_co_evolution_score(aln_file,distance_matrix):
    """
        Calculates co-evolution score of the query based on the MSA alignment
        Co-evolution score measures co-occurence of a pair of amino acid in
        ortholog sequences. Two amino acid have co-evoluated if the occurence
        of one of this amino never occur whithout the other.

        Args:
            aln_file:clustal file in .aln
            distance_matrix: distance matrix of the query against itself
        Returns:
            score_co-evolution
    """
    #Parsing aln file in order to get query + gaps
    with open(aln_file,'r') as file:
        data = file.read()
        query = data.split('\n', 1)[0]
        #repport index of gaps
    i_gap =[i for i, e in enumerate(query) if e == "-"]
    #Predict contacts
    ccmpred_cline = CCMpredCommandline(
        cmd ='bin/CCMpred/bin/ccmpred', alnfile= aln_file, matfile= "contact.mat"
    )
    ccmpred_cline()
    #Extract 30 top coupling
    subprocess.call(
        "./bin/CCMpred/scripts/top_couplings.py contact.mat > top_output.mat",
        shell=True
    )
    #parsing top_coupling
    dic_top_score= {}
    with open("top_output.mat", "r") as file:
        index = 1
        file.readline()
        for line in file:
            if float(line[0:2]) not in i_gap and float(line[3:5]) not in i_gap:
                dic_top_score[index] = [float(line[0:2]), float(line[3:5]),
                float(line[6:19])]
                index += 1
    #re-indexing
    for i, val in dic_top_score.items():
        for pos_gap in range(len(i_gap)-1):
            if val[0] > i_gap[pos_gap] and val[0] < i_gap[pos_gap+1]:
                val[0] = val[0]+1
            if val[1] > i_gap[pos_gap] and val[1] < i_gap[pos_gap+1]:
                val[1] = val[1]+1
    #compare top 30/ distance_matrix
    for i, val in dic_top_score.items():
        if distance_matrix[val[0],val[1]] < 7:
            score_co-evo = score_co-evo +1
        else:
            score_co-evo = score_co-evo -1
    return score_co-evo
