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


def plot_structure_benchmark(output_path, benchmarking_scores, struct, scores, rank):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            output_path (str): The path to store png file.
            structure (str): One of the three following : "Family", "Superfamily", "Fold".
            scores (list): A list of score name.
            rank (list): A list of rank from 1 to N.
    """
    os.makedirs(output_path, exist_ok=True)

    ali_struct = benchmarking_scores[scores[0]][struct].values
    thr_struct = benchmarking_scores[scores[1]][struct].values
    ss_struct = benchmarking_scores[scores[2]][struct].values
    sum_struct = benchmarking_scores[scores[3]][struct].values

    plt.figure(num=struct)  # Window's name
    plt.plot(rank, ali_struct, "b", label=scores[0])
    plt.plot(rank, thr_struct, "#ffa201", label=scores[1])
    plt.plot(rank, ss_struct, "#00B200", label=scores[2])
    plt.plot(rank, sum_struct, "r", label=scores[3])
    plt.plot([0, len(ali_struct)], [0, max(ali_struct)], "k", label="random")
    plt.title("Global scores comparison using " + struct + " benchmarks")
    plt.ylabel("benchmark")
    plt.xlabel("rank")
    plt.legend(loc="lower right")
    plt.savefig(output_path + struct + "_plot.png")
    plt.show()
    print("The plot for '" + struct + "' is stored in " + output_path)

def top_N(output_path, BENCHMARKING_SCORES, score, n):
    """
        Create a top_N_stats.txt file showing statistics based on the benchmark.list files

        Args:
            output_path (str): The path to store stat file.
            score (str): the score you want some stats on
            n (str): a maximum rank number

    """
    os.makedirs(output_path, exist_ok=True)

    structures = ["Family", "Superfamily", "Fold"]
    rank = {}
    max_rank = {}
    for structure in structures:
        rank[structure] = BENCHMARKING_SCORES[score][structure][n]
        max_rank[structure] = max(BENCHMARKING_SCORES[score][structure])
    str = "          Family      Superfamily      Fold \ntop {0}    {1}/{2}       {3}/{4}      {5}/{6} \n%           {7:.2f}          {8:.2f}          {9:.2f}".format(n, rank["Family"], max_rank["Family"],rank["Superfamily"],
        max_rank["Superfamily"], rank["Fold"], max_rank["Fold"],
        rank["Family"]/max_rank["Family"],rank["Superfamily"]/max_rank["Superfamily"],
        rank["Fold"]/max_rank["Fold"])
    with  open(output_path + "top_N_stats.txt", "w") as fileout:
        fileout.write(str)


if __name__ == "__main__":
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    SCORES = ["alignment", "threading", "secondary_structure", "sum scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each structure (= 3 columns)
    BENCHMARKING_SCORES = {}
    for score in SCORES:
        BENCHMARKING_SCORES[score] = pd.DataFrame(np.zeros((405, 3)), columns=STRUCTURES)
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
    COLORS = cycler('color',
                    ['#EE6666', '#3388BB', '#9988DD',
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
        plot_structure_benchmark(OUTPUT_PATH, BENCHMARKING_SCORES, structure, SCORES, RANK)
    OUTPUT_PATH = "results/top_N/"
    top_N(OUTPUT_PATH, BENCHMARKING_SCORES, "sum scores", 10)
