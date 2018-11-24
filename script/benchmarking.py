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


def plot_structure_benchmark(benchmarking_scores, struct, scores, rank):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            struct (str) : one of the three following : "Family", "Superfamily", "Fold"
            scores (list):
            rank (list):
    """
    output_path = "results/plot/"
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


if __name__ == "__main__":
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    SCORES = ["alignment", "threading", "secondary_structure", "sum scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each structure (= 3 columns)
    BENCHMARKING_SCORES = {}
    for score in SCORES:
        BENCHMARKING_SCORES[score] = pd.DataFrame(np.zeros((413, 3)), columns=STRUCTURES)
    # For each query,
    ALL_FOLDRECS = os.listdir("data/foldrec")
    print("Processing all benchmarks ...\n")
    for ind, query in enumerate(ALL_FOLDRECS, 1):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not
        if not os.path.isfile("results/" + query + "/scores.csv"):
            print("\nProcessing query {} / {} : {}\n".format(ind, len(ALL_FOLDRECS), query))
            p = subprocess.Popen(["./fold_u", "data/foldrec/" + query + ".foldrec",
                                  "-o", "results/" + query, "--cpu", str(cpu_count())],
                                 stdout=subprocess.PIPE).communicate()[0]
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

    for structure in STRUCTURES:
        plot_structure_benchmark(BENCHMARKING_SCORES, structure, SCORES, RANK)
