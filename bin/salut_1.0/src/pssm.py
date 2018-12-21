#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Creation des matrices de frequences.
Pour la query :
- Lit le fichier de sortie de MUSCLE
- Recupere le nom de la query
- Cree une liste de sequences alignees
Pour les templates :
- Lit le fichier .map des templates
- Cree une liste de sequences alignees
- Cree un mfasta depuis cette liste (<name>.mfasta)
Pour les deux :
- Calcule le poids associe a chaque sequences (liste de poids)
- Calcule la frequence brute (somme des poids)
- Lisse les frequences avec des frequences background
- Retire les gaps sur la sequence principale

Ce script cree des matrices PSSM avec pour nom <nom>.aamtx
"""

import re
import numpy as np
import os
from random import randint

#Listes dont a besoin les scripts de creation des pssm
aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-', 'X', 'B', 'Z', 'J'] #X est n'importe quel aa
bg_freq=[0.074, 0.052, 0.045, 0.053, 0.025, 0.034, 0.054, 0.074, 
0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051,
0.013, 0.032, 0.073, 0.0, 0.0]
AA={"J":["L","I"],"B":["N","D"],"Z":["E","Q"]}

BETA = 1 #Parametre pour la correction des frequences, pourrait etre change


def poids_seq(namefile, fold) :
    """
    Utilise le script perl compute_weight_sequence_position.pl et recupere les poids
    associes a chaque sequence. Cree un fichier poids.out et le lit
    """
    os.system("perl bin/salut_1.0/src/compute_weight_sequence_position.pl {} > {}.poids.out".format(namefile, fold))
    file = open(fold+".poids.out", "r")
    poids_list = [float(line[18:-1]) for line in file if re.search("^Normalized_weight", line)]
    file.close()
    return(poids_list)


def gap_out(freq_matrice, seq_list, seq) :
    """
    Enleve les positions ou il a que des gap
    C'est a dire les positions avec que des 0'
    Et enleve les positions dans les sequences de seq_list aussi
    """
    cpt = 0
    for aa in seq :
        if aa == "-" : #il s'agit d'un gap
            #On supprime la position cpt de la matric
            freq_matrice = np.delete(freq_matrice, (cpt), axis = 0) 
            seq = seq[:cpt] + seq[(cpt+1):] #On enleve le gap dans la seq d'interet
            for i in range(len(seq_list)) : #Il faut raccourcir chaque sequence
                #On concatene autour de la position
                seq_list[i] = seq_list[i][:cpt] + seq_list[i][(cpt+1):] 
            #Si c'etait un gap, alors on a raccourci la liste et la matrice, 
            #donc pas besoin d'avancer d'un cran
        else :
            cpt += 1 #Si ce n'etait pas un gap, il faut avancer d'un cran
    return ([freq_matrice, seq_list, seq])




def freq_correction(freq_matrice, bg_freq, beta) :
    """
    Corrige les frequences brutes avec les frequences background
    """
    for i in range(freq_matrice.shape[0]) : #i est la position
        for j in range(freq_matrice.shape[1]) : #j est l'aa
            freq_matrice[i,j]=(freq_matrice[i,j]+(beta*bg_freq[j]))/(1+beta)
    return(freq_matrice)


def freq_correctionbis(cpt_matrice, beta) :
    """
    Corrige les counts avec equiprobabilites
    """
    for i in range(cpt_matrice.shape[0]) :
        for j in range(cpt_matrice.shape[1]) :
            cpt_matrice[i,j] = ((cpt_matrice[i,j] + (1/20) )/(1+beta))
    return(cpt_matrice)


def freq_matrix(namefile, namematrix, seq_list, aa, bg_freq, beta, poids_list, fold, seq, AA) :
    """
    Calcule les frequences des aa a chaque position en commencant par les frequences brutes
    Puis en les corrigeant avec la fonction freq_correction
    Ecrit un fichier aamtx
    """
    cpt_matrice = np.zeros((len(seq_list[0]), 20)) #Matrice initialisee a 0 avec 
    if len(poids_list) == 0 :
        return(0)
    else :
        # autant de positions que de aa dans la premiere seq et 20 aa
        for j in range(len(seq_list)) : #On parcourt toutes les seq une a une
            for i in range(len(seq_list[j])) : #Dans une seq, on regarde les aa un par un
                num_aa = aa.index(seq_list[j][i]) #On recupere l'emplacement dans l'alphabet
                if num_aa <20 : #l'aa n'est pas un gap ni un X
                    cpt_matrice[i, num_aa] += poids_list[j] #On ajoute le poids de la seq
                elif num_aa > 21 :#Ce n'est pas un gap ni un X
                    num_aa = aa.index(AA[seq_list[j][i]][randint(0,1)])
                    cpt_matrice[i, num_aa] += poids_list[j]
        #freq_matrice = cpt_matrice / len(seq_list)
        #On enleve les Gap
        #freq_matrice = gap_out(freq_matrice, seq_list)[0]
        sans_gap = gap_out(cpt_matrice, seq_list, seq)
        cpt_matrice = sans_gap[0]
        seq_list = sans_gap[1]
        seq = sans_gap[2]
        #ATTENTION IL FAUT CHOISIR LA CORRECTION EN ENLEVANT LE DIESE ET EN METTANT L'AUTRE EN COMMENTAIRE
        #Correction avec freq background blosum62
        freq_matrice = freq_correction(cpt_matrice, bg_freq, beta)
        #Correction avec freq background 1/20
        #freq_matrice = freq_correctionbis(cpt_matrice, beta)
        #Ecriture du fichier aamtx
        matrix = open(fold+"/"+namematrix+"/"+namematrix+".aamtx", "w")
        matrix.write(">"+namematrix+"\n")
        matrix.write(seq+"\n") #La premiere seq est la query ou template
        for i in range(freq_matrice.shape[0]) :
            for j in range(freq_matrice.shape[1]) :
                matrix.write("{:.7f} ".format(freq_matrice[i,j]))
            matrix.write("\n")
        matrix.close()
        
        return(freq_matrice)
