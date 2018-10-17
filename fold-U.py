#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS
import src.parsing as parse
import src.threading as threading
import numpy as np

METAFOLD_FILE = "data/METAFOLD.list"
FOLDREC_FILE = "data/foldrec/Cohesin.foldrec"
NB_TEMPLATES = 10
DIST_RANGE = [5, 10]

if __name__ == "__main__":
    metafold_dict = parse.metafold(METAFOLD_FILE)
    alignment_list = parse.foldrec(FOLDREC_FILE, NB_TEMPLATES, metafold_dict)
    for ali in alignment_list:
        print(len(ali.query_residues), len(ali.template.residues))
        matrix = threading.calc_dist_matrix(
            ali.query_residues, ali.template.residues, DIST_RANGE)
        #threading.display_matrix(matrix)
