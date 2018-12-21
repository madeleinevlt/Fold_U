#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys


def AssignationSS(ss) :

    ''' Convertit les symboles des structures 
    secondaires de dssp pour qu'ils 
    correspondent à ceux de Psipred'''

    SSpsipred = ''
    if ss == 'H' or ss == 'I' :
        SSpsipred = 'H'
    elif ss == 'E' or ss == 'B' :
        SSpsipred = 'E'
    elif ss == 'G' :
        SSpsipred = 'G'
    else :
        SSpsipred = 'C'
    return SSpsipred


def main():

    #Ouverture du fichier dans lequel sont ajoutées l'assignation de structure 
    #secondaire et l'accessibilité au solvant pour la template concernée
    fileOut = open('dsspAll_out.txt', 'a')

    #liste récupérant la structure secondaire de la template
    list_ss = []

    #liste récupérant l'accessibilité au solvant 
    #pour chaque résidue de la template
    list_acc = []

    #fichier de sortie dssp : dssp<nom_query>.out
    #récupère juste le nom de la query
    template = sys.argv[1][4:-4]

    #Ecriture du nom de la template dans le fichier dsspAll_out.txt
    fileOut.write('>' + template + '\n') 

    #Ouverture du fichier de sortie dssp de la template
    with open(sys.argv[1], 'r') as dssp_file :

        line = dssp_file.readline()


        while not re.search('  #  RESIDUE', line) :
            #Lecture sans rien faire tant que l'on arrive pas aux lignes donnant 
            #la structure secondaire et l'accessibilité au solvant 
            line = dssp_file.readline()

        for line in dssp_file :
            if line[13]!= '!' :    #Pas d'AA détecté

                #Ajout de l'assignation et de l'accessibilité
                list_ss.append(AssignationSS(line[16]))
                list_acc.append(int(line[35:38]))
    
    for i in range(len(list_ss)) :
        #Sous le nom de la template -> ajout de la structure secondaire
        fileOut.write(list_ss[i])

    fileOut.write('\n')

    for j in range(len(list_acc)) :
        #Sous la structure secondaire -> ajout de l'accessibilité au solvant
        fileOut.write(str(list_acc[j]) + ' ')

    fileOut.write('\n')


    fileOut.close()

if __name__=='__main__':
    main()