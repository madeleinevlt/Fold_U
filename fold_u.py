#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        fold_u.py FILE [--nb_templates NUM] [--nb_pdb NB_PDB] [--output PATH]
                       [--metafold METAFOLD] [--dope DOPE]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              the foldrec file [default: 100]
        -p NB_PDB, --nb_pdb NB_PDB            Number of pdb to create
                                              [default: 10]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: res/threading]
        -m METAFOLD, --metafold METAFOLD      Path to the metafold.list file
                                              [default: data/METAFOLD.list]
        -d DOPE, --dope DOPE                  Path to the dope.par file
                                              [default: data/dope.par]
"""

# Third-party modules
from multiprocessing import Pool, cpu_count
from docopt import docopt
from schema import Schema, And, Use, SchemaError
import numpy as np

# Local modules
import src.parsing as parsing
import src.score as score


DIST_RANGE = [5, 15]
GAP_PENALTY = 1


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.

        Args:
            void

        Returns:
            void
    """
    schema = Schema({
        'FILE': Use(open, error='FILE should be readable'),
        '--metafold': Use(open, error='METAFOLD should be readable'),
        '--dope': Use(open, error='DOPE should be readable'),
        '--nb_templates': And(Use(int), lambda n: 1 <= n <= 413,\
                                error='--nb_templates=NUM should be integer 1 <= N <= 413'),
        '--nb_pdb': And(Use(int), lambda n: 1 <= n <= 413,\
                                error='--nb_pdb=NB_PDB should be integer 1 <= N <= 413'),
        # -The output PATH is created (if not exists) at the end of the program
        # so we skip the check.
        object: object})
    try:
        schema.validate(ARGUMENTS)
    except SchemaError as err:
        exit(err)


def process(ali):
    """
        Does the threading and gives the score for a given Alignment object.

        Args:
            ali (alignment object): An object of the class Alignment

        Returns:
            tupple: threading's result for the given alignment as a tuple
                    (template's score, template's name)

    """
    # Calculate the energy matrix
    energy_matrix = ali.calculate_energy(DIST_RANGE, GAP_PENALTY, DOPE_DICT)
    # for i, res in enumerate(energy_matrix):
    #     for j, res2 in enumerate(energy_matrix):
    #         print("{:^5.1f}".format(energy_matrix[i, j]),end="")
    #     print("")
    return np.nansum(energy_matrix), ali.template.name


if __name__ == "__main__":


    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='fold-U 1.1')

    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    FOLDREC_FILE = ARGUMENTS["FILE"]
    # Process the first n templates only
    NB_TEMPLATES = int(ARGUMENTS["--nb_templates"])
    # Create the first n pdb templates
    NB_PDB = int(ARGUMENTS["--nb_pdb"])
    # METAFOLD file
    METAFOLD = ARGUMENTS["--metafold"]
    # DOPE file
    DOPE = ARGUMENTS["--dope"]
    # SCORES file
    OUTPUT = ARGUMENTS["--output"]

    ### Parse data files
    ####################

    # Parse Metafold file
    METAFOLD_DICT = parsing.parse_metafold(METAFOLD)
    # Parse Foldrec file
    ALIGNMENT_DICT = parsing.parse_foldrec(FOLDREC_FILE, NB_TEMPLATES, METAFOLD_DICT)
    DOPE_DICT = parsing.parse_dope(DOPE)


    ### Main calculations
    #####################

    # Parallelization of the main loop: threading calculations
    POOL = Pool(processes=cpu_count())
    # Necessary to pass ARGUMENTS to parallelized function
    THREADING_SCORE = score.Score(POOL.imap(process, ALIGNMENT_DICT.values()))
    POOL.close()
    POOL.join()

    ### Results : Score and PDB files
    #################################
    THREADING_SCORE.write_score(OUTPUT, NB_PDB, ALIGNMENT_DICT)
