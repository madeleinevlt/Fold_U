#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Third-party modules
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler


def plot_benchmark(output_path, structure, scores, rank):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            output_path (str): The path to store png file.
            structure (str): One of the three following : "Family", "Superfamily", "Fold".
            scores (list): A list of score name.
            rank (list): A list of rank from 1 to N.
    """
    os.makedirs(output_path, exist_ok=True)

    ali_structure = benchmarking_scores[scores[0]][structure].values
    thr_structure = benchmarking_scores[scores[1]][structure].values
    sum_structure = benchmarking_scores[scores[2]][structure].values

    plt.figure(num=structure) # Window's name
    plt.plot(rank, ali_structure, "b", label=scores[0])
    plt.plot(rank, thr_structure, "#ffa201", label=scores[1])
    plt.plot(rank, sum_structure, "r", label=scores[2])
    plt.plot([0,len(ali_structure)],[0,max(ali_structure)], "k", label="random")
    plt.title("Global scores comparison using " + structure + " benchmarks")
    plt.ylabel("benchmark")
    plt.xlabel("rank")
    plt.legend(loc="lower right")
    plt.savefig(output_path + structure + "_plot.png")
    plt.show()

def top_N(output_path, score, n):
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
        rank[structure] = benchmarking_scores[score][structure][n]
        max_rank[structure] = max(benchmarking_scores[score][structure])
    str = "          Family      Superfamily      Fold \ntop {0}    {1}/{2}       {3}/{4}      {5}/{6} \n%           {7:.2f}          {8:.2f}          {9:.2f}".format(n, rank["Family"], max_rank["Family"],rank["Superfamily"],
        max_rank["Superfamily"], rank["Fold"], max_rank["Fold"],
        rank["Family"]/max_rank["Family"],rank["Superfamily"]/max_rank["Superfamily"],
        rank["Fold"]/max_rank["Fold"])
    with  open(output_path + "top_N_stats.txt", "w") as fileout:
        fileout.write(str)


if __name__ == "__main__":
    structures = ["Family", "Superfamily", "Fold"]
    scores = ["alignment", "threading", "sum scores"]
    # A dictionary of pandas DataFrames is created for each score
    # Each DataFrame will contain the cumulative sum of benchmarks for each structure (= 3 columns)
    benchmarking_scores = {}
    for score in scores:
        benchmarking_scores[score] = pd.DataFrame(np.zeros((406,3), dtype=np.float64), columns=structures)
    # For each query,
    all_foldrecs = os.listdir("data/foldrec")
    for ind, query in enumerate(all_foldrecs, 1):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not
        if not os.path.isfile("results/" + query + "/scores.csv"):
            print("\nProcessing query {} / {} : {}\n".format(ind, len(all_foldrecs), query))
            p = subprocess.Popen(["./fold_u", "data/foldrec/" + query + ".foldrec", "-o", "results/" + query, "--cpu", "4"], stdout=subprocess.PIPE).communicate()[0]
            for i in p.decode("UTF-8").split("\n")[1:]:
                print(i)
            print("--------------------------------------------------------------------")
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
    colors = cycler('color', ['#EE6666', '#3388BB', '#9988DD',
                              '#EECC55', '#88BB44', '#FFBBBB'])
    plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
           axisbelow=True, grid=True, prop_cycle=colors)
    plt.rc('grid', color='w', linestyle='solid')
    plt.rc('xtick', direction='out', color='gray')
    plt.rc('ytick', direction='out', color='gray')
    plt.rc('patch', edgecolor='#E6E6E6')
    plt.rc('lines', linewidth=1.5)

    output_path = "results/plot/"
    for structure in structures:
        plot_benchmark(output_path, structure, scores, rank)
    print("\nThe plots are stored in " + output_path + "\n")
    output_path = "results/top_N/"
    top_N(output_path, "sum scores", 10)
