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
from docopt import docopt

# Local modules
import src.parsing as parse
import src.threading as threading


DIST_RANGE = [5, 10]

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

    for ali in ALIGNMENT_LIST:
        print(ali.template.pdb)
        matrix = threading.calc_dist_matrix(
            ali.query_residues, ali.template.residues, DIST_RANGE)
    # threading.display_matrix(matrix)
