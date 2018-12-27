#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# """
# Creation des PSSM pour la base de donnees templates/
# Ce script fait appel a des fonctions du script pssm.py 
# D'abord il lit un fichier map et le transforme en mfasta qu'il met dans 
# un le dossier tmp_template
# Puis la pssm est calculee
# Dans les fichiers map, la sequence d'interet n'est pas forcement la 1ere,
# ce script va donc chercher la sequence d'interet correspondant au
# code pdb du fichier METAFOLD.txt et l'utiliser comme sequence de reference.
# Si jamais le code pdb n'existe pas ou alors il n'existe mais n'est pas
# trouve dans le fichier .map, alors ce script mettra par defaut la premiere 
# sequence du .map comme sequence de reference.
# """

import sys
import re
from pssm import *

def pdb_code(template_name):
    """
    Recupere le code pdb du template voulu
    """
    file = open("data/metafold.list", "r")
    flag = 0
    line = "123"
    code_pdb = 0
    while flag == 0 and line != "":#tant qu'on ne l'a pas trouve on continue
        line = file.readline()
        if template_name == line.split(" ")[0] :
            code_pdb = (line.split(" ")[1][:-1]).split(".")[0] #On obtient le code pdb
            flag = 1
    file.close()
    return(code_pdb)


def multi_from_map(seq_list, name):
    """
    Cree un fichier texte mfasta a partir d'une liste de sequences
    """
    file = open("data/templates/"+name+"/"+name+".mfasta", "w")
    for seq in seq_list :
        file.write(">Sequence\n")
        file.write(str(seq))
        file.write("\n")
    file.close()


def read_map(namefile, namematrix):
    """
    Lit un fichier map de la database templates/ et donne le nom de
    la proteine  et une liste des sequences alignees
    """
    file = open(namefile, "r")
    #namematrix = (namefile.split("/")[-1]).split(".")[0] #Nom de la sequence
    begin = re.compile("^>")
    info = re.compile("^sequence")
    struct = re.compile("^structure")
    seq_list = []
    code_pdb = pdb_code(namematrix)
    flag = 0
    flag2 = 0
    cpt_seq = 0
    for line in file :
        if begin.search(line) :#on est sur une nouvelle seq
            seq_list.append("") #on cree une nouvelle seq
            cpt_seq += 1
            if code_pdb in line : #il s'agit de la seq d'interet
                flag = 1
                flag2 = 1 #on l'a trouve dans le fichier map
            else :
                flag = 0 #si jamais on a passe la seq
        elif not (info.search(line) or struct.search(line)) : #c'est un bout de seq
            #if len(seq) > 0 :#pour eviter e rajouter des seq quand on lit
            #deux lignes de suite sans seq
            if line[-2] == "*" :
                seq_list[cpt_seq-1] += line[:-2] #ajout a la seq actuelle sans l'etoile
            else :
                seq_list[cpt_seq-1] += line[:-1] #ajout a la seq actuelle
        if flag == 1 : #s'il s'agit de la seq d'interet
            seq_principale = seq_list[cpt_seq-1]
    file.close()
    multi_from_map(seq_list, namematrix)
    if code_pdb == 0 or flag2 == 0 : #soit on l'a pas dans le metafold, soit pas dans le map
        seq_principale = seq_list[0]
    return([seq_list, seq_principale])


def create_hom_pssm(template_name, seq_listhom, aa, bg_freq, beta, seq, AA):
    """
    Cree la PSSM de la seq template en faisant appel aux fonctions
    du script pssm.py
    """
    poids_list = poids_seq(sys.argv[1]+".mfasta", sys.argv[1])
    #on utilise le fichier mfasta cree
    if len(poids_list) != 0 :
        return(freq_matrix(template_name, seq_listhom, aa, bg_freq, beta, poids_list, "data/templates/", seq, AA))
    else :
        print("non calculé : {}".format(template_name))

def main():
    map_path = sys.argv[1]+".map"
    #Map multi ali
    template_name = sys.argv[2]
    read_out = read_map(map_path, template_name)
    seq_listhom = read_out[0]
    seq = read_out[1]
    freq_matricehom = create_hom_pssm(template_name, seq_listhom, aa, bg_freq, BETA, seq, AA)

if __name__== "__main__":
    main()
