#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


def load_blosum(namefile) :
    """
    Chargement de la matrice BlOSUM62 dans une dataframe (pandas)
    """
    blosum = []
    with open(namefile, 'r') as f :
        for line in f:
            if line[0] != '#' :
                l = line.strip().split(' ')
                l = [i for i in l if i != '']
                if not l[1].isalpha() :
                    l = [float(i) for i in l[1:] if i]
                blosum.append(l)
    aa = blosum[0]
    blosum = blosum[1:]
    blosum = pd.DataFrame(np.matrix(blosum), columns=aa, index=aa)
    return(aa, blosum)