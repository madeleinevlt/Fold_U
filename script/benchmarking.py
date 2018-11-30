#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# """
# .. module:: benchmarking
#    :synopsis: This module runs the program "fold_u" on all of the benchmarks. It then generates
#               plots to visualize the contribution of each and every score to the re-ranking of
#               models/templates.
#               Three plots are generated, one for each structure type of benchmark
#               "Fold", "Superfamily" and "Family".
# """
"""
    Usage:
        ./script/benchmarking.py [--nb_templates NUM] [--output PATH] [--dssp PATH] [--sscore SCORE]
                                 [--cpu NUM]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates with the best
                                              score [default: 100]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and plot)
                                              [default: ./results/top_n]
        -d PATH, --dssp PATH                  Path to the dssp software
                                              binary [default: /usr/local/bin/mkdssp]
        -s SCORE, --sscore SCORE              Score for which you wish to see the statistics:
                                              "alignment", "threading", "modeller",
                                              "secondary_structure", "solvent_access" or all of them
                                              at once: "sum_scores".
                                              [default: sum_scores]
        -c NUM, --cpu NUM                     Number of cpus to use for parallelisation
                                              [default: 2]
"""

# Third-party modules
import os
import subprocess
from multiprocessing import cpu_count
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler
from docopt import docopt
from schema import Schema, And, Use, SchemaError


def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '--nb_templates': And(Use(int), lambda n: 1 <= n <= 100,
                              error='--nb_templates=NUM should be integer 1 <= N <= 100'),
        '--dssp': Use(open, error='dssp/mkdssp should be readable'),
        '--sscore': And(Use(str), lambda s: s in ["alignment", "threading", "modeller",
                                                  "secondary_structure", "solvent_access",
                                                  "sum_scores"],
                        error='SCORES should be an existing score'),
        '--cpu': And(Use(int), lambda n: 1 <= n <= cpu_count(),
                     error='--cpus=NUM should be integer 1 <= N <= ' + str(cpu_count())),
        # The output PATH is created (if not exists) at the end of the program so we skip the check.
        object: object})
    try:
        schema.validate(ARGUMENTS)
    except SchemaError as err:
        exit(err)


def plot_benchmark(output_path, struct, scores, rank, benchmarking_scores, selected_score):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            output_path (str): The path to store png file.
            struct (str): One of the three following : "Family", "Superfamily",
                          "Fold".
            scores (list): A list of score name.
            rank (list): A list of rank from 1 to N.
            benchmarking_scores (dict): a dictionnary containing benchmarking
                                        data for each type of score for
                                        different types of structures.

    """
    os.makedirs(output_path, exist_ok=True)
    plt.figure(num=struct)  # Window's name
    # Plot all scores
    if selected_score == "sum_scores":
        ali_struct = benchmarking_scores[scores[0]][struct].values
        thr_struct = benchmarking_scores[scores[1]][struct].values
        mod_struct = benchmarking_scores[scores[2]][struct].values
        ss_struct = benchmarking_scores[scores[3]][struct].values
        acc_struct = benchmarking_scores[scores[4]][struct].values
        sum_struct = benchmarking_scores[scores[5]][struct].values

        plt.plot(rank, ali_struct, "b", label=scores[0])
        plt.plot(rank, thr_struct, "#ffa201", label=scores[1])
        plt.plot(rank, mod_struct, "#EE82EE", label=scores[2])
        plt.plot(rank, ss_struct, "#00B200", label=scores[3])
        plt.plot(rank, acc_struct, "#7a9a91", label=scores[4])
        plt.plot(rank, sum_struct, "r", label=scores[5])
        plt.plot([0, len(ali_struct)], [0, max(ali_struct)], "k", label="random")
        plt.title("Global scores comparison using " + struct + " benchmarks")
        plt.ylabel("Benchmark")
        plt.xlabel("rank")
        plt.legend(loc="lower right")
        plt.savefig(output_path + "/" + struct + "_plot.png")
    # Plot scores individually
    else:
        score_struct = benchmarking_scores[selected_score][struct].values
        plt.plot(rank, score_struct, "b", label=selected_score)
        plt.plot([0, len(score_struct)], [0, max(score_struct)], "k", label="random")
        plt.title(selected_score + " score using " + struct + " benchmarks")
        plt.ylabel("Benchmark")
        plt.xlabel("rank")
        plt.legend(loc="lower right")
        plt.savefig(output_path + "/" + selected_score + "_" + struct + "_plot.png")
    plt.show()



def top_n(structures, scores, top_n, benchmarking_scores):
    """
        Show statistics based on the benchmark.list files separately for each fold-type: "Fold",
        "Family", "Superfamily".
        Represent the strength/weaknesses of the different scores independantly and/or combined.

        Args:
            structures (list): List containing fold types: "Family", "Superfamily", "Fold"
            scores (str): the score you want some stats on
            top_n (str): a maximum rank number
            benchmarking_scores (dict): a dictionnary containing benchmarking
                                        data for each type of score for
                                        different types of structures.

        Returns:
            a str "top_n_results" table summarizing the top_n results

    """
    rank = {}
    max_rank = {}
    for struct in structures:
        rank[struct] = benchmarking_scores[scores][struct][top_n-1]
        max_rank[struct] = max(benchmarking_scores[scores][struct])
    line1 = "top{0}\t{1}/{2}\t\t{3}/{4}\t\t{5}/{6}\n".format(top_n,
                                                             rank["Family"],
                                                             max_rank["Family"],
                                                             rank["Superfamily"],
                                                             max_rank["Superfamily"],
                                                             rank["Fold"],
                                                             max_rank["Fold"])
    line2 = "\t{0:>5.2f} % {1:>13.2f} %{2:>14.2f} %"\
        .format((rank["Family"]/max_rank["Family"])*100,
                (rank["Superfamily"]/max_rank["Superfamily"])*100,
                (rank["Fold"]/max_rank["Fold"])*100)
    top_n_results = line1 + line2
    return top_n_results


if __name__ == "__main__":
    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__, version='fold_u 1.2')
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    # Process the first n templates only
    NB_TEMPLATES = int(ARGUMENTS["--nb_templates"])
    # OUTPUT file
    OUTPUT_PATH = ARGUMENTS["--output"]
    # DSSP path
    DSSP_PATH = ARGUMENTS["--dssp"]
    # Number of cpus for parallelisation
    NB_PROC = ARGUMENTS["--cpu"]
    # Selected score you want to have info about
    SELECTED_SCORE = ARGUMENTS["--sscore"]
    # The 3 different structures from benchmark
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    # all the possible scores useful for plots
    SCORES = ["alignment", "threading", "modeller", "secondary_structure", "solvent_access", "sum_scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each structure (= 3 columns)
    BENCHMARKING_SCORES = {}
    for score in SCORES:
        BENCHMARKING_SCORES[score] = pd.DataFrame(np.zeros((399, 3)), columns=STRUCTURES)
    # For each query,
    ALL_FOLDRECS = os.listdir("data/foldrec")
    print("Processing all benchmarks ...\n")
    for ind, query in enumerate(ALL_FOLDRECS, 1):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not already generated
        if not os.path.isfile("results/" + query + "/scores.csv"):
            print("\nProcessing query {} / {} : {}\n".format(ind, len(ALL_FOLDRECS), query))
            p = subprocess.Popen(["./fold_u", "data/foldrec/" + query + ".foldrec",
                                  "-o", "results/" + query, "--dssp", DSSP_PATH,
                                  "--cpu", NB_PROC, "--nb_templates", "100"],
                                 stdout=subprocess.PIPE).communicate()[0]
            rows, columns = os.popen('stty size', 'r').read().split()
            print("\n" + "-"*int(columns))
        # Score results are stored in a pandas DataFrame
        query_scores = pd.read_csv("results/" + query + "/scores.csv", index_col=0)
        for score in SCORES:
            # The DataFrame is sorted by the current score
            query_score = query_scores.sort_values(by=score, ascending=False)
            # Initialization of the dictionary of counts
            structures_count = {}
            for structure in STRUCTURES:
                structures_count[structure] = 0
            # Cumulative sum of benchmark according to the structure
            for i, struct_type in enumerate(query_score["benchmark"]):
                for structure in STRUCTURES:
                    if struct_type == structure:
                        structures_count[structure] += 1
                    BENCHMARKING_SCORES[score][structure][i+1] += structures_count[structure]

    RANK = [i for i in BENCHMARKING_SCORES[SCORES[0]][STRUCTURES[0]].index]

    # plot settings
    COLORS = cycler('color', ['#EE6666', '#3388BB', '#9988DD',
                              '#EECC55', '#88BB44', '#FFBBBB'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
           axisbelow=True, grid=True, prop_cycle=COLORS)
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='gray')
    plt.rc('ytick', direction='out', color='gray')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=1.5)

    #OUTPUT_PATH = "results/plot/"
    for structure in STRUCTURES:
        plot_benchmark(OUTPUT_PATH, structure, SCORES, RANK, BENCHMARKING_SCORES, SELECTED_SCORE)
    print("\nThe plots are stored in " + OUTPUT_PATH + "\n")
    #OUTPUT_PATH = "results/top_N/"
    n_top_n = [5]
    if NB_TEMPLATES <= 5:
        n_top_n = [5]
    if NB_TEMPLATES >= 10:
        n_top_n = [5, 10]
    if NB_TEMPLATES >= 50:
        n_top_n = [5, 10, 50]
    if NB_TEMPLATES >= 100:
        n_top_n = [5, 10, 50, 100]
    print("Table summarizing the top {} results.\n".format(n_top_n))
    print("\tFamily\t\tSuperfamily\tFold\n")
    for top in n_top_n:
        print(top_n(STRUCTURES, SELECTED_SCORE, top, BENCHMARKING_SCORES))
        print("\t----------------------------------------")
    # with  open(output_path + "top_N_stats.txt", "w") as fileout:
        # fileout.write(top_n_results)
