#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        fold_u.py FOLDREC_FILE [--nb_templates NUM] [--nb_pdb NUM] [--output PATH]
                               [--metafold FILE] [--dope FILE] [--blosum FILE]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              the foldrec file [default: 100]
        -p NUM, --nb_pdb NUM                  Number of pdb to create
                                              [default: 10]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: res/threading]
        -m FILE, --metafold FILE              Path to the metafold.list file
                                              [default: data/METAFOLD.list]
        -d FILE, --dope FILE                  Path to the dope.par file
                                              [default: data/dope.par]
        -b FILE, --blosum FILE                  Path to the dope.par file
                                              [default: data/BLOSUM62.txt]
"""

# Third-party modules
from multiprocessing import Pool, cpu_count
import numpy as np
from tqdm import tqdm
from docopt import docopt
from schema import Schema, And, Use, SchemaError

# Local modules
import src.parsing as parsing
from src.score import Score
from src.scoring import *

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
        'FOLDREC_FILE': Use(open, error='FOLDREC_FILE should be readable'),
        '--metafold': Use(open, error='METAFOLD_FILE should be readable'),
        '--dope': Use(open, error='DOPE_FILE should be readable'),
        '--nb_templates': And(Use(int), lambda n: 1 <= n <= 413,\
                                error='--nb_templates=NUM should be integer 1 <= N <= 413'),
        '--nb_pdb': And(Use(int), lambda n: 1 <= n <= 413,\
                                error='--nb_pdb=NUM should be integer 1 <= N <= 413'),
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
    # Return the sum of energy matrix with numpy's "Nan" interpreted as zeros with
    # alignment's informations
    return np.nansum(energy_matrix), ali.num, ali.template.name


if __name__ == "__main__":


    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='fold-U 1.1')

    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    FOLDREC_FILE = ARGUMENTS["FOLDREC_FILE"]
    # Process the first n templates only
    NB_TEMPLATES = int(ARGUMENTS["--nb_templates"])
    # Create the first n pdb templates
    NB_PDB = int(ARGUMENTS["--nb_pdb"])
    # METAFOLD file
    METAFOLD_FILE = ARGUMENTS["--metafold"]
    # DOPE file
    DOPE_FILE = ARGUMENTS["--dope"]
    # DOPE file
    BLOSUM_FILE = ARGUMENTS["--blosum"]
    # SCORES file
    OUTPUT_PATH = ARGUMENTS["--output"]

    ### Parse data files
    ####################

    # Parse Metafold file
    METAFOLD_DICT = parsing.parse_metafold(METAFOLD_FILE)
    # Parse Foldrec file
    ALIGNMENT_DICT = parsing.parse_foldrec(FOLDREC_FILE, NB_TEMPLATES, METAFOLD_DICT)
    # Parse DOPE file
    DOPE_DICT = parsing.parse_dope(DOPE_FILE)
    # Parse BLOSSUM file
    BLOSSUM_DICT = parsing.parse_blosum(BLOSUM_FILE)


    ### TEST SCORES ###
    for ali in ALIGNMENT_DICT.values():
        ss_score(ali)
        blossum_score(ali, BLOSSUM_DICT)
        break
    ###################


    ### Main calculations
    #####################

    # Parallelization of the main loop: threading calculations
    POOL = Pool(processes=cpu_count())
    # tqdm module enables an ETA progress bar of alignments
    # imap_unordered can smooth things out by yielding faster-calculated values
    # ahead of slower-calculated values.
    print("\nProcessing threading on templates ...\n\n")
    RESULTS = [res_ali for res_ali in tqdm(POOL.imap_unordered(process,\
                ALIGNMENT_DICT.values()), total=len(ALIGNMENT_DICT.values()))]
    POOL.close()
    POOL.join()

    # Calculate the threading score of all alignments
    THREADING_SCORE = Score(RESULTS)

    ### Results : Score and PDB files
    #################################
    THREADING_SCORE.write_score(OUTPUT_PATH, NB_PDB, ALIGNMENT_DICT)
    print("\nThe program ended successfully !\nThe results are stored in " + OUTPUT_PATH)
