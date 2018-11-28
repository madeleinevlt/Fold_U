#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. module:: benchmarking
   :synopsis: This module runs the program "fold_u" on all of the benchmarks. It then generates
              plots to visualize the contribution of each and every score to the re-ranking of
              models/templates.
              Three plots are generated, one for each structure type of benchmark
              "Fold", "Superfamily" and "Family".
"""

# Third-party modules
import os
import subprocess
from multiprocessing import cpu_count
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler


def plot_benchmark(output_path, struct, scores, rank, benchmarking_scores):
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

    ali_struct = benchmarking_scores[scores[0]][struct].values
    thr_struct = benchmarking_scores[scores[1]][struct].values
    mod_struct = benchmarking_scores[scores[2]][struct].values
    ss_struct = benchmarking_scores[scores[3]][struct].values
    sum_struct = benchmarking_scores[scores[4]][struct].values

    plt.figure(num=struct)  # Window's name
    plt.plot(rank, ali_struct, "b", label=scores[0])
    plt.plot(rank, thr_struct, "#ffa201", label=scores[1])
    plt.plot(rank, mod_struct, "#EE82EE", label=scores[2])
    plt.plot(rank, ss_struct, "#00B200", label=scores[3])
    plt.plot(rank, sum_struct, "r", label=scores[4])
    plt.plot([0, len(ali_struct)], [0, max(ali_struct)], "k", label="random")
    plt.title("Global scores comparison using " + struct + " benchmarks")
    plt.ylabel("benchmark")
    plt.xlabel("rank")
    plt.legend(loc="lower right")
    plt.savefig(output_path + struct + "_plot.png")
    plt.show()
    print("The plot for '" + struct + "' is stored in " + output_path)

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
        rank[struct] = benchmarking_scores[scores][struct][top_n]
        max_rank[struct] = max(benchmarking_scores[scores][struct])
    line1 =  "top{0}\t{1}/{2}\t\t{3}/{4}\t\t{5}/{6}\n".format(top_n,
                                                        rank["Family"],
                                                        max_rank["Family"],
                                                        rank["Superfamily"],
                                                        max_rank["Superfamily"],
                                                        rank["Fold"],
                                                        max_rank["Fold"])
    line2 =  "\t{0:.2f} %\t\t{1:.2f} %\t\t{2:.2f} %"\
                .format((rank["Family"]/max_rank["Family"])*100,
                        (rank["Superfamily"]/max_rank["Superfamily"])*100,
                        (rank["Fold"]/max_rank["Fold"])*100)
    top_n_results = line1 + line2
    return top_n_results


if __name__ == "__main__":
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    SCORES = ["alignment", "threading", "modeller", "secondary_structure", "sum scores"]
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
                                  "-o", "results/" + query, "--cpu", str(cpu_count())],
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

    OUTPUT_PATH = "results/plot/"
    for structure in STRUCTURES:
        plot_benchmark(OUTPUT_PATH, structure, SCORES, RANK, BENCHMARKING_SCORES)
    print("\nThe plots are stored in " + OUTPUT_PATH + "\n")
    OUTPUT_PATH = "results/top_N/"
    print("\tFamily\t\tSuperfamily\tFold\n")
    for i in [5, 10, 50, 100]:
        print(top_n(STRUCTURES, "sum scores", i, BENCHMARKING_SCORES))
        print("\t----------------------------------------")
    # with  open(output_path + "top_N_stats.txt", "w") as fileout:
        # fileout.write(top_n_results)
