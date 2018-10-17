#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        fold-U.py FILE [--nb_templates NUM] [--metafold METAFOLD] [--dope DOPE]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              The foldrec file [default:100]
        -m METAFOLD, --metafold METAFOLD      Path to the metafold.list file
        -d DOPE, --dope DOPE                  Path to the dope.par file
"""

# Local modules
import src.parsing as parse
import src.threading as threading

# Third-party modules
import numpy as np
from docopt import docopt


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

    if arguments["--metafold"]:
        METAFOLD = arguments["--metafold"]
    if arguments["--dope"]:
        DOPE = arguments["--dope"]

    # Parse Metafold file
    metafold_dict = parse.metafold(METAFOLD)
    # Parse Foldrec file
    alignment_list = parse.foldrec(foldrec_file, nb_templates, metafold_dict)
    # Parse DOPE file
    dope_df = parse.dope(DOPE)

    for ali in alignment_list:
        print(ali.template.pdb)
        matrix = threading.calc_dist_matrix(
            ali.query_residues, ali.template.residues, DIST_RANGE)
    # threading.display_matrix(matrix)
