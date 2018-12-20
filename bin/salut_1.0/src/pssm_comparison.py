#!/usr/bin/env python

import sys
import numpy as np
import blosum as bl
from random import randint


def load_pssm(namefile, aa) :
    """
    Lecture d'un fichier .aamtx et renvoie matrice PSSM
    """
    with open (namefile, 'r') as f :
        for line in f :
            if line[0] == '>' :
                name = line[1:-1]
            elif line[0] in aa :
                seq = line[:-1]
                pssm = []
            else :
                pssm.append([float(j) for j in line[:-1].split(' ') if j])

    return(name, list(seq), pssm)


def calc_score(vect1, vect2) :
    """
    Calcule du score d'alignement entre 2 aa pour une position donnée 
    """
    s = 0
    n1, n2 = len(vect1), len(vect2)
    for i in range(n1) :
        for j in range(n2) :
            s += vect1[i]*vect2[j]
    
    return(s)


def calc_gap_optimal(nq, nt, pssmQ, pssmT, blosum):
    """
    Fonction qui calcule et renvoie les pénalités d'ouverture (po) et d'extension de gap (pe) optimales
    """

    gaps = np.arange(0.1, 10, 0.1)
    A = []
    scores = []
    blsmean = blosum.stack().mean()
    blsstd = blosum.stack().std()

    for gap in gaps:
        A.append((gap - blsmean) / blsstd)

    # Scores - random PSSM
    for i in range(nq):
        for j in range(nt):
            ri, rj = randint(0, nq-1), randint(0, nt-1) #sélection des indices de la pssm de manière aléatoires
            scores.append(calc_score(pssmQ[ri], pssmT[rj]))

    smean = np.array(scores).mean()
    sstd = np.array(scores).std()
    po = []
    pe = []
    for i in range(len(A)):
        po.append((A[i] * sstd) + smean)
        pe.append(po[i] * 10 / 100)

    return(-np.array(po).mean(), -np.array(pe).mean())


def init_matInsert(i, j, po):
    """
    Initialisation matrice des ajouts de gaps (au niveau de la Query (Q)
    ou au niveau de la template (T)) entre aaQ (en position i) et aaT (en position j)
    """

    if i > 0:
        return -np.inf
    else:
        if j > 0 and i == 0:
            return po*j
        else:
            return 0


def init_matMatch(i, j):
    """
    Initialisation matrice des matchs (M) entre aaQ (en position i) et aaT (en position j)
    """
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return -np.inf
        else:
            return 0


def scores_propositions(i, j, M, Q, T, pssmQ, pssmT, p, type_score):
    """
    Calcul des propositions de scores d'alignement aaQ - aaT
    """
    propositions = []

    if type_score == 'M':
        score = calc_score(pssmQ[i-1], pssmT[j-1])
        propositions.append(score + M[i-1,j-1])
        propositions.append(Q[i,j])
        propositions.append(T[i,j])
    elif type_score == 'Q':
        propositions.append(M[i,j-1] + p)
        propositions.append(Q[i,j-1] + p)
        propositions.append(T[i,j-1] + p)
    elif type_score == 'T':
        propositions.append(M[i-1,j] + p)
        propositions.append(Q[i-1,j] + p)
        propositions.append(T[i-1,j] + p)

    return(propositions)


def crea_matrix(fileQ, fileT):
    """
    Création et initialisation des matrices des matchs (M) et des insertions de gaps (Q et T)
    (Méthode affine)
    """
    aa, blosum = bl.load_blosum("BLOSUM62.txt")
    nameQ, seqQ, pssmQ = load_pssm(fileQ, aa)
    nameT, seqT, pssmT = load_pssm(fileT, aa)

    if len(seqQ) != len(pssmQ):
        print(">>> ERREUR : La longueur de la séquence QUERY est différente du nombre de lignes lues dans la PSSM QUERY\n")
        sys.exit(1)
    elif len(seqT) != len(pssmT):
        print(">>> ERREUR : La longueur de la séquence TEMPLATE est différente du nombre de lignes lues dans la PSSM TEMPLATE\n")
        sys.exit(1)

    nq, nt = len(seqQ)+1, len(seqT)+1
    po, pe = calc_gap_optimal(len(seqQ), len(seqT), pssmQ, pssmT, blosum)

    M = np.zeros((nq, nt)) # Matchs
    Q = np.zeros((nq, nt)) # Ajout gap query
    T = np.zeros((nq, nt)) # Ajout gap template

    for i in range(nq) :
        for j in range(nt) :
            if i == 0 or j == 0 :
            # Initialisation des matrices d'alignements Match (M) / InsertGapQuery (Q) / InsertGapTemplate (T)
                M[i,j] = init_matMatch(i, j)
                Q[i,j] = init_matInsert(j, i, po)
                T[i,j] = init_matInsert(i, j, po)
            else :
            # Remplissage des matrices M, Q et T
                mat = [M[i-1,j-1], Q[i,j-1], T[i-1,j]]
                if mat.index(max(mat)) > 0 :
                    p = pe
                else:
                    p = po
                Q[i,j] = max(scores_propositions(i, j, M, Q, T, pssmQ, pssmT, p, 'Q'))
                T[i,j] = max(scores_propositions(i, j, M, Q, T, pssmQ, pssmT, p, 'T'))
                M[i,j] = max(scores_propositions(i, j, M, Q, T, pssmQ, pssmT, p, 'M'))

    return(nameQ, nameT, M, Q, T, seqQ, seqT, po, pe)


def align_matrix(M, Q, T, seqQ, seqT):
    """
    Lecture de l'alignement optimal / global (algorithme de Needleman et Wunsch)
    """
    alignQ, alignT = [], []
    i, j = len(seqQ), len(seqT)
    score = max(M[i,j], Q[i,j], T[i,j])

    while i > 0 and j > 0:
        mat = [M[i-1,j-1], Q[i,j-1], T[i-1,j]]
        # dans le cas où seqQ[i] et seqT[j] sont alignés
        if mat.index(max(mat)) == 0 and i > 0 and j > 0: # M
            i, j = i-1, j-1
            alignQ.append(seqQ[i])
            alignT.append(seqT[j])

        # cas où il y a ajout d'un gap au niveau de seqQ[i]
        elif mat.index(max(mat)) == 1 and j > 0: # Q
            i, j = i, j-1
            alignQ.append("-")
            alignT.append(seqT[j])

        # cas où il y a ajout d'un gap au niveau de seqT[j]
        elif mat.index(max(mat)) == 2 and i > 0: # T
            i, j = i-1, j
            alignQ.append(seqQ[i])
            alignT.append("-")

    alignQ.reverse()
    alignT.reverse()

    return("".join(alignQ), "".join(alignT), score)


def align_truncate(align, seq, deb, fin):
    """
    Fonction qui tronque un alignment (retire les gaps en amont + aval) et renvoie le résultat + l'intervalle des aa alignés
    """
    interv = []
    tmp = align.find('-')

    if tmp == -1: 
        return(align, [1, len(align)])
    else:
        interv.append(align.find(seq[deb]) +1)
        interv.append(align.rfind(seq[fin]) +1)

    trunc = align[interv[0]-1 : interv[1]]
 
    return(trunc, interv)


def crea_output(nameQ, nameT, truncQ, truncT, nq, nt, intervQ, intervT, scoreGlobal) :
    """
    Création et remplissage du fichier en output (query_name.aln)
    """
    fd = open("alignments/{}.aln".format(nameQ), "a+")
    fd.write("Query | {} | {:.3f} | {} | {} | {}-{} | {}-{}\n".format(nameT, scoreGlobal, nq, nt, intervQ[0], intervQ[1], intervT[0], intervT[1]))
    fd.write("{}\n{}\n".format(truncQ, truncT))
    fd.close()

    fscore = open("alignments/{}_scores.txt".format(nameQ), "a+")
    fscore.write("{} {:.3f}\n".format(nameT, scoreGlobal))
    fscore.close()



def main():
    # Fichiers input :
    if len(sys.argv) < 3 :
        print(">>> ERREUR : PSSM de la Query et/ou Template manquante(s)\n")
        sys.exit(1)
    else :
        fileQ = sys.argv[1]
        fileT = sys.argv[2]

    # Construction et lecture de l'alignement entre Query et Template
    nameQ, nameT, M, Q, T, seqQ, seqT, po, pe = crea_matrix(fileQ, fileT)
    print("Alignement : {} (Q) vs {} (T)".format(nameQ, nameT))
    print("Pénalités optimales : ouverture de gap = {:.2f}, extension de gap = {:.3f}".format(po, pe))

    alignQ, alignT, scoreGlobal = align_matrix(M, Q, T, seqQ, seqT)
    truncT, intervT = align_truncate(alignT, seqT, 0, -1)
    truncQ = alignQ[ intervT[0]-1 : intervT[1] ]
    intervQ = [ alignQ.find(seqQ[0])+1 ,  len(truncQ) - truncQ.count('-')]

    # Création de l'output
    crea_output(nameQ, nameT, truncQ, truncT, len(seqQ), len(seqT), intervQ, intervT, scoreGlobal)



if __name__=='__main__':
    main()
    print()