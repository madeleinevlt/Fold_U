#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

def addInDict(dico, key, value) :
    '''
    Ajoute dans un dictionnaire la valeur de la clé en argument si elle existe.
    Crée une nouvelle clé si elle n'existe pas et ajoute sa valeur.
    '''
    if key in dico.keys() :
        dico[key] += value
    else :
        dico[key] = value
    return dico

#################### MAIN PROGRAM #####################

def main():

    #Choix du round
    n_round = sys.argv[2]

    #dictionnaire :
        #clé = code Uniprot
        #valeur = [start, stop, sequence]
    dict_sbjct = {}

    #Expressions régulières pour sélectionner des lignes
    re_score = re.compile("^ Score")
    re_sbjct = re.compile("^Sbjct")
    re_query = re.compile("^Query")
    re_end = re.compile("^  Database|^Searching")  #Pour arrêter la lecture

    with open (sys.argv[1]+'.blast_out', 'r') as file1 :
        line = file1.readline()
        
        #Lit jusqu'à trouver les résultats du round choisi
        while not re.search('Results from round ' + str(n_round), line) :
            line = file1.readline()

        while not re_end.match(line) :
            
            #Pour les lignes commençant par >
            if re.search('^>', line) :
                #décompose la ligne de l'ID
                ref = re.match(
                    '((>|> )UniRef[0-9]*_)([A-Z0-9]*)([A-Za-z0-9 .,_-]*)', line)
                start = ''
                stop = ''   
                seq = '' 
                #si nscore > 1 le résultat ne sera pas pris en compte     
                nscore = 0
                line = file1.readline()

                #Tant qu'il ne trouve pas un autre ID ou l'indicateur de fin
                while (not re.search('>', line)) and (not re_end.match(line)):

                    line = file1.readline()

                    if re_score.match(line) :
                        #For one ID, add 1 for each alignment
                        #Pour 1 ID, ajoute 1 pour chaque alignement
                        nscore += 1

                    #Si c'est la 1ere ligne d'une subjct et que c'est le 1er
                    #alignement (meilleur score pour l'ID)
                    if re_sbjct.match(line) and start == '' and nscore < 2 :
                        #décompose la ligne de la séquence pour récupérer la 
                        #1ere position, la dernière et la séquence
                        sbjct = re.match(
                            '(Sbjct  )([0-9]+)[ ]+([A-Z-]+)[ ]+([0-9]+)', line)
                        start = sbjct.group(2)
                        stop = sbjct.group(4)
                        seq += sbjct.group(3)

                    #Si ce n'est pas la 1ere ligne de la sbjct et que c'est le 
                    #1er alignement
                    elif re_sbjct.match(line) and start != '' and nscore < 2 :
                        #décompose la ligne de la séquence pour récupérer la 
                        #la dernière position et la séquence
                        sbjct = re.match(
                            '(Sbjct  )([0-9]+)[ ]+([A-Z-]+)[ ]+([0-9]+)', line)
                        stop = sbjct.group(4)
                        seq += sbjct.group(3)  

                #Ajout dans le dictionnaire
                addInDict(dict_sbjct, ref.group(3), [start, stop, seq])                
                    
            #S'il ne trouve pas une ligne commençant par >         
            else :
                line = file1.readline()

    #Crée un fichier où tous les alignements récupérés sont ajoutés 
    fileOut = open(sys.argv[1]+'.faa', 'w')

    #Récupération de la sequence query
    fileQuery = open(sys.argv[1]+'.fasta','r')
    line1 = fileQuery.readline()
    line2 = fileQuery.readline()
    fileOut.write(line1 + line2)

    for key in dict_sbjct.keys() : 
        fileOut.write('>' + key + ' [' + dict_sbjct[key][0] + '...' + 
            dict_sbjct[key][1] + ']' + '\n' + dict_sbjct[key][2] + '\n')

    fileQuery.close()
    fileOut.close()

if __name__=='__main__':
    main()
