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

# Third-party modules
from multiprocessing import Pool, cpu_count
from functools import partial
from operator import itemgetter
import numpy as np
from docopt import docopt

# Local modules
import src.parsing as parse
import src.threading as threading




DIST_RANGE = [5, 15]


def loop(alignement):
    """
        blabla

        Args:
            alignment:

        Returns:
            best template's name:
            best template's score:

    """
    # Calculate the distance matrix
    matrix, dist_dict = threading.calc_dist_matrix(
        alignement.query_residues, alignement.template.residues, DIST_RANGE)
    # Convert distances into energies based on DOPE energies
    energy_matrix = threading.convert_dist_to_energy(matrix, dist_dict, DOPE_DF)
    return alignement.template.name, np.nansum(energy_matrix)




if __name__ == "__main__":

    # Parse command line
    ARGUMENTS = docopt(__doc__, version='fold-U 1.0')
    FOLDREC_FILE = ARGUMENTS["FILE"]
    # First n templates
    if ARGUMENTS["--nb_templates"]:
        NB_TEMPLATES = int(ARGUMENTS["--nb_templates"])
    else:
        NB_TEMPLATES = 100
    # METAFOLD file
    if ARGUMENTS["--metafold"]:
        METAFOLD = ARGUMENTS["--metafold"]
    else:
        METAFOLD = "data/METAFOLD.list"
    # DOPE file
    if ARGUMENTS["--dope"]:
        DOPE = ARGUMENTS["--dope"]
    else:
        DOPE = "data/dope.par"

    # Parse Metafold file
    METAFOLD_DICT = parse.metafold(METAFOLD)
    # Parse Foldrec file
    ALIGNMENT_LIST = parse.foldrec(FOLDREC_FILE, NB_TEMPLATES, METAFOLD_DICT)
    # Parse DOPE file
    DOPE_DF = parse.dope(DOPE)


    # Parallelization of the main loop
    POOL = Pool(processes=cpu_count())
    # Necessary to pass arguments to parallelized function
    FUNC = partial(loop, DIST_RANGE, DOPE_DF)
    SCORES = POOL.imap(FUNC, ALIGNMENT_LIST)
    POOL.close()
    POOL.join()


    ########## RESULTS ##########


    BEST_TEMPLATE, BEST_SCORE = min(SCORES, key=itemgetter(1))
    print("\n\nBest template: {}\nBest energy: {:.2f}".format(BEST_TEMPLATE, BEST_SCORE))
