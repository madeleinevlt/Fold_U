#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Third-party modules
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler


def plot_benchmark(structure, scores, rank):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            structure (str) : one of the three following : "Family", "Superfamily", "Fold"
            scores (list):
            rank (list):
    """
    OUTPUT_PATH = "results/plot/"
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    ali_structure = benchmarking_scores[scores[0]][structure].values
    thr_structure = benchmarking_scores[scores[1]][structure].values
    ss_structure = benchmarking_scores[scores[2]][structure].values
    sum_structure = benchmarking_scores[scores[3]][structure].values

    plt.figure(num=structure) # Window's name
    plt.plot(rank, ali_structure, "b", label=scores[0])
    plt.plot(rank, thr_structure, "#ffa201", label=scores[1])
    plt.plot(rank, ss_structure, "#00B200", label=scores[1])
    plt.plot(rank, sum_structure, "r", label=scores[2])
    plt.plot([0,len(ali_structure)],[0,max(ali_structure)], "k", label="random")
    plt.title("Global scores comparison using " + structure + " benchmarks")
    plt.ylabel("benchmark")
    plt.xlabel("rank")
    plt.legend(loc="lower right")
    plt.savefig(OUTPUT_PATH + structure + "_plot.png")
    plt.show()
    print("\nThe plots are stored in " + OUTPUT_PATH)


if __name__ == "__main__":
    structures = ["Family", "Superfamily", "Fold"]
    scores = ["alignment", "threading", "secondary_structure", "sum scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each structure (= 3 columns)
    benchmarking_scores = {}
    for score in scores:
        benchmarking_scores[score] = pd.DataFrame(np.zeros((413,3), dtype=np.int64), columns=structures)
    # For each query,
    for query in os.listdir("data/foldrec"):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not
        if not os.path.isfile("results/" + query + "/scores.csv"):
            subprocess.Popen(["./fold_u", "data/foldrec/" + query + ".foldrec", "-o", "results/" + query, "--cpu", "8"], stdout=subprocess.PIPE).communicate()
        # Score results are stored in a pandas DataFrame
        query_scores = pd.read_csv("results/" + query + "/scores.csv", index_col=0)
        for score in scores:
            # The DataFrame is sorted by the current score
            query_score = query_scores.sort_values(by=score, ascending=False)
            # Initialization of the dictionary of counts
            structures_count = {}
            for structure in structures:
                structures_count[structure] = 0
            # Cumulative sum of benchmark according to the structure
            for i, type in enumerate(query_score["benchmark"]):
                for structure in structures:
                    if type == structure:
                        structures_count[structure] += 1
                    benchmarking_scores[score][structure][i+1] += structures_count[structure]

    rank = [i for i in benchmarking_scores[scores[0]][structures[0]].index]

    # plot settings
    colors = cycler('color',
                    ['#EE6666', '#3388BB', '#9988DD',
                     '#EECC55', '#88BB44', '#FFBBBB'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
           axisbelow=True, grid=True, prop_cycle=colors)
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='gray')
    plt.rc('ytick', direction='out', color='gray')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=1.5)

    for structure in structures:
        plot_benchmark(structure, scores, rank)
