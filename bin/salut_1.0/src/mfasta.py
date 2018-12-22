#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

def crea_mfasta(prot, fd):
    fd.write("{:s}".format(prot[0]))
    for i in range(1, len(prot)):
        fd.write("{:s}".format(prot[i]))


def main():

    if len(sys.argv) != 2 :
        print(">>> ERREUR : nombre d'arguments - mfasta\n")
        sys.exit(1)
    else :
        query = sys.argv[1]
        file = "data/queries/"+query+"/"+query+".mfasta"

    # Stockage des données du fichier.mfasta
    with open (file, 'r') as f :
        mfasta = []
        tmp = []
        cpt = 0
        for line in f:

            if line[0] == '>' and tmp != []:
                mfasta.append(tmp)
                tmp = []
                cpt += 1

            if line == "> "+query+"\n":
                qcpt = cpt

            tmp.append(line)
        mfasta.append(tmp)

    # Réécriture dans un mfasta corrigé
    os.remove("data/queries/{}/{}.mfasta".format(query, query))
    fd = open("data/queries/{}/{}.mfasta".format(query, query), "a+")
    crea_mfasta(mfasta[qcpt], fd)

    for i in range(len(mfasta)):
        if i != qcpt:
            crea_mfasta(mfasta[i], fd)

    fd.close()


if __name__=='__main__':
    main()


