#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Usage:
        ./fold_u QUERY_SEQUENCE_PATH   [--nb_pdb NUM] [--output PATH] [--dssp PATH] [--cpu NUM]
                                       [--metafold FILE] [--dope FILE] [--benchmark FILE]

    Arguments:
        QUERY_SEQUENCE_PATH                   Path to the QUERY sequence to predicted

    Options:
        -h, --help                            Show this
        -p NUM, --nb_pdb NUM                  Number of pdb to create
                                              [default: 10]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: ./results]
        -a PATH, --dssp PATH                  Path to the dssp software
                                              binary [default: /usr/bin/mkdssp]
        -c NUM, --cpu NUM                     Number of cpus to use for parallelisation. By default
                                              using all available (0).
                                              [default: 0]
        -m FILE, --metafold FILE              Path to the metafold.list file
                                              [default: data/metafold.list]
        -d FILE, --dope FILE                  Path to the dope.par file
                                              [default: data/dope.par]
        -b FILE, --benchmark FILE             Path to the benchmark.list file
                                              [default: data/benchmark.list]
"""



#First Step : Generating files on HOMESTRAD database (independant Task)
#subprocess.call("./src/partie_amont/salut_1.0/salut1.sh", shell=True)

#2nd Step :


#3rd
