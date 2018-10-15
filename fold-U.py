#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS
import src.parsing as parse
import src.threading as threading

if __name__ == "__main__":

    nb_templates = 10
    dist_range = [5, 10]
    alignment_list = parse.foldrec(nb_templates)
    for ali in alignment_list:
        matrix = threading.calc_dist_matrix(ali.query_seq, ali.template.seq, dist_range)
    print(matrix)
