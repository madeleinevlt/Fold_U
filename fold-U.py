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
from multiprocessing import Pool, cpu_count


DIST_RANGE = [5, 15]

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


    scores = {}

    for ali in alignment_list:
        #print(ali.template.pdb)
        matrix, dist_dict = threading.calc_dist_matrix(
            ali.query_residues[:30], ali.template.residues[:30], DIST_RANGE)
        #threading.display_matrix(matrix)

        energy_matrix = threading.convert_dist_to_energy(matrix, dist_dict, dope_df)
        scores[ali.template.name] = np.sum(energy_matrix)

    print(max(scores))
