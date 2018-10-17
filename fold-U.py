#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        fold-U.py FILE [--nb_templates NUM]

    Options:
        -h, --help                   Show this
        -n NUM, --nb_templates NUM   First n templates to retrieve from
                                     the foldrec file [default:100]

"""

# IMPORTS
from docopt import docopt
import src.parsing as parse
import src.threading as threading
import numpy as np

METAFOLD_FILE = "data/METAFOLD.list"
DIST_RANGE = [5, 10]

if __name__ == "__main__":

    # Parse command line
    arguments = docopt(__doc__, version='fold-U 1.0')
    foldrec_file = arguments["FILE"]
    # First n templates
    if arguments["--nb_templates"]:
        nb_templates = int(arguments["--nb_templates"])
    else:
        nb_templates = 100

    metafold_dict = parse.metafold(METAFOLD_FILE)
    alignment_list = parse.foldrec(foldrec_file, nb_templates, metafold_dict)
    for ali in alignment_list:
        print(ali.template.pdb)
        matrix = threading.calc_dist_matrix(
            ali.query_residues, ali.template.residues, DIST_RANGE)
        #threading.display_matrix(matrix)
