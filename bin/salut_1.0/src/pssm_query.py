#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation de la PSSM de la query
Ce script fait appel a des fonctions du script pssm.py 
D'abord il lit un fichier de sortie MUSCLE 
Puis la pssm est calculee
"""

import sys
import re
from pssm import *

def read_multi(namefile, name) :
    '''
    Lit un fichier de sorti de multi-alignement MUSCLE et donne
    une liste des sequences alignees ainsi que la seq d'interet
    '''
    file = open(namefile, "r")
    #line_one = file.readline() #Premiere ligne avec le nom de la query
    #namematrix = line_one.split(" ")[0][1:-1] #Nom de la query
    begin = re.compile("^>")
    seq_list = []
    cpt_seq = 0
    flag = 0
    for line in file : 
        if begin.search(line) :
            seq_list.append('') #on cree une nouvelle seq
            cpt_seq += 1
            if name in line : #Il s'agit de la seq d'interet
                flag = 1
            else :
                flag = 0 #si jamais on a passe la seq
        else :
            seq_list[cpt_seq-1] += line[:-1] #ajout a la seq actuelle
        if flag == 1 : #s'il s'agit de la seq d'interet
            seq_principale = seq_list[cpt_seq-1]
    file.close()
    return([seq_list, seq_principale])

def create_quer_pssm(namefile, namematrix, seq_list, aa, bg_freq, beta, seq, AA) :
    '''
    Cree la PSSM de la query en faisant appel aux fonctions
    du script pssm.py
    '''
    poids_list = poids_seq(namefile, sys.argv[1])
    if len(poids_list) != 0 :
        print("PSSM calculée")
        return(freq_matrix(namematrix, seq_list, aa, bg_freq, beta, poids_list, "data/queries/", seq, AA))
    else :
        print("PSSM non calculée")




def main() :
    namefile = sys.argv[1]+".mfasta"
    query_name = sys.argv[2]
    #Muscle multi ali
    read_out = read_multi(namefile, query_name)
    seq_list = read_out[0]
    seq = read_out[1] #sequence d'interet
    freq_matrice = create_quer_pssm(namefile, query_name, seq_list, aa, bg_freq, BETA, seq, AA)



    
if __name__=='__main__':
    main()