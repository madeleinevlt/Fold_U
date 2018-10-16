#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS
import src.parsing as parse
import src.threading as threading

METAFOLD_FILE = "data/METAFOLD.list"
FOLDREC_FILE = "data/foldrec/Agglutinin.foldrec"
NB_TEMPLATES = 10
DIST_RANGE = [5, 10] 

if __name__ == "__main__":
    metafold_dict = parse.metafold(METAFOLD_FILE)
    alignment_list = parse.foldrec(FOLDREC_FILE, NB_TEMPLATES, metafold_dict)
    for ali in alignment_list:
        matrix = threading.calc_dist_matrix(ali.query_seq, ali.template.seq, DIST_RANGE)
    print(matrix)
