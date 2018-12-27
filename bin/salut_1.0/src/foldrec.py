#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ce script cree un fichier .foldrec
Il part d'un fichier d'alignement des profils (query vs template)
Il lit les scores et trie les alignements selon ces scores (ordre croissant)
Il lance un psi-pred sur la query entiere qui ensuite sera decoupee selon l'alignement
et ajout les gaps
Il lit un fichier de sortie de dssp sur les templates :
- recupere la template correspondante
- decoupe le resultat selon l'alignement
- rajoute les gaps
La fonction cree deux fichiers :
- .foldrec1
- .foldrec2
Ces deux fichiers seront concatenes en un fichier .foldrec
"""

import sys
import re
import os
import numpy as np
import blosum as bl

###############################TRI DE L'ALIGNEMENT PROFIL PROFIL#####################################


def tri_insertion(dico_score, list_score):
    '''
    Effectue un tri par insertion en se basant sur les scores donnes en input
    '''
    for i in range(1,len(list_score)):
        en_cours = list_score[i]
        en_cours_info = dico_score["info"][i]
        en_cours_query = dico_score["query"][i]
        en_cours_template = dico_score["template"][i]
        en_cours_normscore = dico_score["normscore"][i]
        j = i
        #décalage des éléments du list_score
        while j>0 and list_score[j-1]>en_cours:
            list_score[j]=list_score[j-1]
            dico_score["info"][j]=dico_score["info"][j-1]
            dico_score["query"][j]=dico_score["query"][j-1]
            dico_score["template"][j]=dico_score["template"][j-1]
            dico_score["normscore"][j]=dico_score["normscore"][j-1]
            j = j-1
            #on insère l'élément à sa place
            list_score[j]=en_cours
            dico_score["info"][j]=en_cours_info
            dico_score["query"][j]=en_cours_query
            dico_score["template"][j]=en_cours_template
            dico_score["normscore"][j]=en_cours_normscore
    return(dico_score)


def tri_croissant(dico_score) :
    '''
    Trie le dictionnaire dans l'ordre croissant des scores
    '''
    list_score = [] #On va recuperer la liste des scores pour l'ordonner
    for info in dico_score["info"] :
        list_score.append(float(info.split(" | ")[2]))
    #On a la liste des scores, maintenant on va les trier
    dico_score = tri_insertion(dico_score, list_score)
    return(dico_score)


def lecture_scores(filescores):
    scores = []

    with open(filescores, 'r') as f :
        for line in f:
            l = line[:-1].split(' ')
            scores.append(float(l[1]))

    smin = min(scores)
    smax = max(scores)

    return(smin, smax)



def lecture_aligns(filealign, filescores):
    '''
    Lit un fichier score de sortie de comparaison profil-profil
    Cree un dictionnaire contenant :
    - les informations sur la comparaison
    - la sequence query alignee
    - la sequence template alignee
    '''
    #meanscores, stdscores = lecture_scores(filescores)
    smin, smax = lecture_scores(filescores)

    dico_score = {}
    dico_score["info"]=[]
    dico_score["query"]=[]
    dico_score["template"]=[]
    dico_score["normscore"]=[]
    score_file = open(filealign, "r")
    Q_T = 1 #flag pour savoir s'il s'agit de la query ou du template
    for line in score_file :
        if re.search("Query", line) :
            dico_score["info"].append(line[:-1]) 
            score = float(line.split(" | ")[2])
            #dico_score["normscore"].append((score - meanscores) / stdscores)
            dico_score["normscore"].append((score - smin) / (smax - smin))
        else :
            if Q_T == 1 : #il s'agit de la query
                dico_score["query"].append(line[:-1])
                Q_T = 2
            else : #il s'agit de la template
                dico_score["template"].append(line[:-1])
                Q_T = 1
    score_file.close()
    return(dico_score)


##################################         PSI-PRED QUERY           #################################


def horiz_lecture(filename):
    '''
    Lit un fichier horiz (sortie de psipred) et recupere la prediction et les scores
    '''
    file = open(filename, "r")
    struct = ""
    conf = ""
    for line in file :
        if re.search("Pred:", line) :
            struct = struct + line[6:-1]
        elif re.search("Conf:", line) :
            conf = conf + line[6:-1]

    file.close()
    return([struct, conf])


def secondary_structure_query(dico_score, query_pred) :
    '''
    Decoupe les prediction psi-pred et ajoute les gaps selon la query alignee
    On aura une liste de listes avec pour chaque liste la pred et le score de confiance
    '''
    query_prediction = []
    #On gere les query alignee une par une
    for i in range(len(dico_score["query"])) : #on parcourt toutes les seq query de la sortie profil/profil
        seq = dico_score["query"][i] #On s'interesse a la query i
        df = dico_score["info"][i].split(" | ")[5] #On recupere le numero des aa debut et fin
        debut = int(df.split("-")[0]) #numero du premier aa
        fin = int(df.split("-")[1]) #numero du dernier aa
        #On prend la partie qui nous interesse (entre debut et fin)
        pred = query_pred[0][(debut-1):fin]
        conf = query_pred[1][(debut-1):fin]
        #Maintenant on va rajouter les gaps
        for j in range(len(seq)) :
            if seq[j] == "-" : # s'il s'agit d'un gap
                pred = pred[:j]+"-"+pred[j:] #On ajoute le gap a la position j des structures
                conf = conf[:j]+"-"+conf[j:] #On ajoute un gap a la position j des scores
        #On a fini avec la seq i, on passe a la suivante
        query_prediction.append([pred, conf])
    return(query_prediction)


################################  secondary structure template  #####################################


def secondary_structure_template(dico_score) :
    '''
    Decoupe les assignations dssp et ajoute les gaps selon les templates alignees    
    On aura une liste de listes avec pour chaque liste la pred,
    le score de confiance (9) et l'accessibilite au solvent
    '''
    #Lecture fichier dssp
    dsspfile = open("bin/salut_1.0/data/dsspAll_out.txt","r")
    dssp = {}
    for line in dsspfile :
        if re.search(">", line) :
            nametemplate = line.split(">")[1][:-1]
            dssp[nametemplate] = [] #On initialise la cle du dico avec la template
        else :
            dssp[nametemplate].append(line[:-1])
    dsspfile.close()
    template_prediction = []
    #On gere les template une par une
    for i in range(len(dico_score["template"])) :#On parcourt toutes les seq template de la sortie profil/profil
        seq = dico_score["template"][i] #On s'interesse a la template i
        df = dico_score["info"][i].split(" | ")[6] #On recupere le numero des aa debut et fin
        debut = int(df.split("-")[0]) #numero du premier aa
        fin = int(df.split("-")[1]) #numero du dernier aa
        nametemplate = dico_score["info"][i].split(" | ")[1]
        #On prend la partie qui nous interesse (entre debut et fin)
        if nametemplate in dssp.keys() : #il y avait bien le pdb et l'assignation
            assign = dssp[nametemplate][0][(debut-1):fin]
            conf = "9" * len(assign) #Il s'agit d'une assignation, donc le score de confiance est a 9
            acc = dssp[nametemplate][1].split(" ")[(debut-1):fin] #accessibilite au solvent
        else :
            assign = "X" * len(seq)
            conf = "X" * len(seq)
            acc = ["X"] * len(seq)
        #Maintenant on va rajouter les gaps
        for j in range(len(seq)) :
            if seq[j] == "-" : # s'il s'agit d'un gap
                assign = assign[:j]+"-"+assign[j:] #On ajoute le gap a la position j des structures
                conf = conf[:j]+"-"+conf[j:] #On ajoute un gap a la position j des scores
                acc = acc[:j]+["-1"]+acc[j:] #On ajoute un gap a la position j des accessibilites
        #On a fini avec la seq i, on passe a la suivante
        template_prediction.append([assign, conf, acc])
    return(template_prediction)



#####################################    creation fichiers     ######################################

def calc_similarity(alignQ, alignT):
    """
    Calcule du pourcentage de similarité des 2 séquences alignés à partir de la BLOSUM62
    """
    n = len(alignT)
    cpt = 0
    sim = 0
    aa, blosum = bl.load_blosum("bin/salut_1.0/data/BLOSUM62.txt")

    for i in range(n):
        if alignQ[i] in aa and alignT[i] in aa:
            cpt += 1
            if blosum[alignQ[i]][alignT[i]] > 0:
                sim += 1

    return(sim * 100 / cpt)


def calcul_information(dico_score, query_len) :
    '''
    Calcule les informations suivantes :
    - query coverage
    - % identity
    - % similarité
    - % Gaps
    - Alignment length
    Les renvoie sous forme d'un dictionnaire
    '''
    dico_align = {} #Il s'agira d'un dico avec des listes dans l'ordre
    dico_align["query_coverage"] = []
    dico_align["identity"] = []
    dico_align["similarity"] = []
    dico_align["gaps"] = []
    dico_align["align_len"] = []
    
    #On parcourt les alignements du dico_score, il y a autant de seq templates que de scores
    for i in range(len(dico_score["template"])) :
        df = dico_score["info"][i].split(" | ")[5] #On recupere le numero des aa debut et fin
        debut = int(df.split("-")[0]) #numero du premier aa
        fin = int(df.split("-")[1]) #numero du dernier aa)
        query_coverage = float((fin - debut + 1) / query_len) * 100
        align_len = len(dico_score["query"][i])
        gap = 0
        identity = 0
        similarity = 0
        for j in range(len(dico_score["template"][i])) : #On parcourt tous les aa de la seq
            #On compte les gaps
            if dico_score["query"][i][j] == "-" : #il s'agit d'un gap
                gap = gap + 1
            #On compte les identites
            if dico_score["query"][i][j] == dico_score["template"][i][j] : #aa identique
                identity = identity + 1
        dico_align["query_coverage"].append(query_coverage)
        dico_align["gaps"].append((gap / align_len) * 100)
        dico_align["identity"].append((identity / align_len) * 100)
        dico_align["similarity"].append(calc_similarity(dico_score["query"][i], dico_score["template"][i]))
        dico_align["align_len"].append(align_len)
    return(dico_align)


def foldrec(name, query_len, dico_score, dico_align, query_prediction, template_prediction) :
    '''
    Ecrit un fichier foldrec en deux parties
    '''
    fold = open("data/queries/"+name+"/"+name+"part1.foldrec", "w")
    fold.write("*** HITS RANKED ***\n\n")
    fold.write("SEQUENCE QUERY FILE : {}, {} aa\n\n".format(name, query_len))
    fold.write("#| Score | Normalized score | Q. Length | T. Length | Q. begin-end |  T. begin-end | HITS\n")
    fold.write("---------------------------------------------------------\n")

    align = open("data/queries/"+name+"/"+name+"part2.foldrec", "w")
    align.write("\n\n\n*** ALIGNMENTS DETAILS ***\n\n")
    numero = 1
    for i in range((len(dico_score["info"])-1), -1, -1) : #On parcourt tous les alignements a l'envers
        score = float(dico_score["info"][i].split(" | ")[2])
        norm = dico_score["normscore"][i]
        Q_len = int(dico_score["info"][i].split(" | ")[3])
        T_len = int(dico_score["info"][i].split(" | ")[4])
        Q_BE = dico_score["info"][i].split(" | ")[5] #debut et fin de la query
        T_BE = dico_score["info"][i].split(" | ")[6] #debut et fin de la template
        fold.write("{0:>5}{1:>9}{2:>9}{3:>6}{4:>6}{5:>9}{6:>9}".format(numero, round(score, 3), round(norm, 3), Q_len, T_len, Q_BE, T_BE))
        hit = dico_score["info"][i].split(" | ")[1]
        scop = open("data/templates/{}/{}.scop_id".format(hit, hit), "r")
        hit = hit +" :  " + scop.readline() #il y a le \n dans la ligne
        fold.write("   {}".format(hit))

        align.write("No {}\n".format(numero))
        align.write("Alignment : {}, {} aa. vs  {}".format(name, query_len, hit))
        align.write("Score :{0:>9} | Normalized score :{1:>9} | Query coverage : {2:}% | ".format(round(score, 3), round(norm, 3), round(dico_align["query_coverage"][i],2)))
        align.write("Identity :{0:>9}% | Similarity :{1:>9}% | Gaps :{2:>9}%".format(round(dico_align["identity"][i],2), round(dico_align["similarity"][i],2), round(dico_align["gaps"][i],2)))
        align.write(" | Alignment length :{0:>6}\n\n".format(dico_align["align_len"][i]))
        #Alignement query
        align.write("Query{0:>8} {1:}{2:>6}\n".format(int(Q_BE.split("-")[0]), dico_score["query"][i], int(Q_BE.split("-")[1])))
        align.write("Query{0:>8} {1:}{2:>6}\n".format(int(Q_BE.split("-")[0]), query_prediction[i][0], int(Q_BE.split("-")[1])))
        align.write("Query{0:>8} {1:}{2:>6}\n\n".format(int(Q_BE.split("-")[0]), query_prediction[i][1], int(Q_BE.split("-")[1])))
        #Alignement template
        align.write("Template{0:>5} {1:}{2:>6}\n".format(int(T_BE.split("-")[0]), dico_score["template"][i], int(T_BE.split("-")[1])))
        align.write("Template{0:>5} {1:}{2:>6}\n".format(int(T_BE.split("-")[0]), template_prediction[i][0], int(T_BE.split("-")[1])))
        align.write("Template{0:>5} {1:}{2:>6}\n\n".format(int(T_BE.split("-")[0]), template_prediction[i][1], int(T_BE.split("-")[1])))
        align.write("{}\n\n".format(" ".join(template_prediction[i][2])))


        numero += 1
        #scop.close()
    align.close()
    fold.close()


def main() :
    namefile = sys.argv[1]
    filescores = namefile.split('.')[0] + "_scores.txt"
    fasta = sys.argv[2]
    name = sys.argv[3]

    #psi-pred sur la query
    query_pred = horiz_lecture("data/queries/"+name+"/"+name+".horiz")

    #Gerer l'information sur la longueur de la query
    f = open(fasta, "r") 
    f.readline()
    query_len = len(f.readline()) - 1
    f.close()

    #Recuperation des seq alignees dans profil/profil
    dico = (lecture_aligns(namefile, filescores))
    #Tri du dictionnaire selon l'ordre croissant des scores
    dico = (tri_croissant(dico))

    #Structures secondaires de la query
    query_prediction = secondary_structure_query(dico, query_pred)
    template_prediction = secondary_structure_template(dico)

    #Calcul des informations concernant l'alignement
    dico_align = calcul_information(dico, query_len)
    
    #Foldrec
    foldrec(name, query_len, dico, dico_align, query_prediction, template_prediction)



if __name__=='__main__':
    main()