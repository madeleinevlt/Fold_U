#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: benchmarking
   :synopsis: This module implements the functions necessary for the benchmark
                of our program.
"""
# Third-party modules
import os
import subprocess
from multiprocessing import cpu_count
import numpy as np
import pandas as pd
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

    ali_structure = benchmarking_scores[scores[0]][struct].values
    thr_structure = benchmarking_scores[scores[1]][struct].values
    sum_structure = benchmarking_scores[scores[2]][struct].values

    plt.figure(num=struct) # Window's name
    plt.plot(rank, ali_structure, "b", label=scores[0])
    plt.plot(rank, thr_structure, "#ffa201", label=scores[1])
    plt.plot(rank, sum_structure, "r", label=scores[2])
    plt.plot([0, len(ali_structure)], [0, max(ali_structure)], "k", label="random")
    plt.title("Global scores comparison using " + struct + " benchmarks")
    plt.ylabel("benchmark")
    plt.xlabel("rank")
    plt.legend(loc="lower right")
    plt.savefig(output_path + struct + "_plot.png")
    plt.show()

def top_n(structures, scores, topn, benchmarking_scores):
    """
        Create a top_N_stats.txt file showing statistics based on the
        benchmark.list files

        Args:
            output_path (str): The path to store stat file.
            score (str): the score you want some stats on
            n (str): a maximum rank number

        Return:
            a str "top_n_results" table summarizing the topN results

    """

    rank = {}
    max_rank = {}
    for struct in structures:
        rank[struct] = benchmarking_scores[scores][struct][topn]
        max_rank[struct] = max(benchmarking_scores[scores][struct])
    line1 = "\tFamily\t\tSuperfamily\tFold\n"
    line2 =  "top{0}\t{1}/{2}\t\t{3}/{4}\t\t{5}/{6}\n".format(topn,
                                                        rank["Family"],
                                                        max_rank["Family"],
                                                        rank["Superfamily"],
                                                        max_rank["Superfamily"],
                                                        rank["Fold"],
                                                        max_rank["Fold"])
    line3 =  "    %\t{0:.2f}\t\t{1:.2f}\t\t{2:.2f}"\
                .format(rank["Family"]/max_rank["Family"],
                        rank["Superfamily"]/max_rank["Superfamily"],
                        rank["Fold"]/max_rank["Fold"])
    top_n_results = line1 + line2 + line3
    return top_n_results


if __name__ == "__main__":
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    SCORES = ["alignment", "threading", "sum scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each
    # structure (= 3 columns)
    BENCHMARKING_SCORES = {}
    for score in SCORES:
        BENCHMARKING_SCORES[score] = pd.DataFrame(np.zeros((406, 3)),
                                                  columns=STRUCTURES)
    # For each query,
    ALL_FOLDRECS = os.listdir("data/foldrec")
    for ind, query in enumerate(ALL_FOLDRECS, 1):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not
        if not os.path.isfile("results/" + query + "/scores.csv"):
            print("\nProcessing query {} / {} : {}\n"\
                    .format(ind, len(ALL_FOLDRECS), query))
            p = subprocess.Popen(["./fold_u",
                                  "data/foldrec/" + query + ".foldrec",
                                  "-o", "results/" + query,
                                  "--cpu", str(cpu_count())],\
                                  stdout=subprocess.PIPE).communicate()[0]
            for i in p.decode("UTF-8").split("\n")[1:]:
                print(i)
            print("--------------------------------------------------------------------")
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
    for i in [10, 50, 100]:
        print(top_n(STRUCTURES, "sum scores", i, BENCHMARKING_SCORES))
        print("\t----------------------------------------")
    #with  open(output_path + "top_N_stats.txt", "w") as fileout:
        #fileout.write(top_n_results)
