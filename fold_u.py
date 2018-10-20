#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        fold_u.py FILE [--nb_templates NUM] [--metafold METAFOLD] [--dope DOPE] [--scores SCORES]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              The foldrec file [default: 100]
        -m METAFOLD, --metafold METAFOLD      Path to the metafold.list file
                                              [default: data/METAFOLD.list]
        -d DOPE, --dope DOPE                  Path to the dope.par file
                                              [default: data/dope.par]
        -s SCORES, --scores SCORES            Path to the results file
                                              [default: res/threading_scores.out]
"""

# Third-party modules
from multiprocessing import Pool, cpu_count
from docopt import docopt
import numpy as np

# Local modules
import src.parsing as parse
import src.threading as threading
import src.writing as write




DIST_RANGE = [5, 15]



def process(ali):
    """
        Does the threading and gives the score for a given Alignment object.

        Args:
            ali (alignment object): An object of the class Alignment

        Returns:
            tupple: threading's result for the given alignment as a tuple
                    (template's score, template's name)

    """
    # Calculate the distance matrix
    query = ali.query_residues
    template = ali.template.residues
    energy_matrix = threading.calc_dist_convert_energy(query, template, DIST_RANGE, DOPE_DICT)
    return np.nansum(energy_matrix), ali.template.name


if __name__ == "__main__":


    ### Parse command line
    ######################

    ARGUMENTS = docopt(__doc__, version='fold-U 1.0')
    FOLDREC_FILE = ARGUMENTS["FILE"]
    # Process the first n templates only
    NB_TEMPLATES = int(ARGUMENTS["--nb_templates"])
    # METAFOLD file
    METAFOLD = ARGUMENTS["--metafold"]
    # DOPE file
    DOPE = ARGUMENTS["--dope"]
    # SCORES file
    SCORES_FILE = ARGUMENTS["--scores"]

    ### Parse data files
    ####################

    # Parse Metafold file
    METAFOLD_DICT = parse.metafold(METAFOLD)
    # Parse Foldrec file
    ALIGNMENT_LIST = parse.foldrec(FOLDREC_FILE, NB_TEMPLATES, METAFOLD_DICT)
    # Parse DOPE file
    DOPE_DICT = parse.dope(DOPE)


    ### Main calculations
    #####################

    # Parallelization of the main loop: threading calculations
    POOL = Pool(processes=cpu_count())
    # Necessary to pass arguments to parallelized function
    SCORES = POOL.imap(process, ALIGNMENT_LIST)
    POOL.close()
    POOL.join()

    ### Results
    ###########
    write.scores(SCORES_FILE, sorted(SCORES))
