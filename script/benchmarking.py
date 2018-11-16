#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Third-party modules
import os
from subprocess import call
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

structures = ["Family", "Superfamily", "Fold"]
scores = ["alignment", "threading", "sum scores"]
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
        call([" ./fold_u data/foldrec/" + query + ".foldrec -o results/" + query], shell=True)
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


def plot_benchmark(benchmark_type):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            benchmark_type (str) : one of the three following : "Family", "Superfamily", "Fold"
    """
    ali_structure = benchmarking_scores["alignment"][benchmark_type]  
    thr_structure = benchmarking_scores["threading"][benchmark_type]  
    #blo_structure = benchmarking_scores["blosum"][benchmark_type]  #score pretty bad so not useful but kept for the report
    sum_structure = benchmarking_scores["sum scores"][benchmark_type]  

    len_matrix = len(benchmarking_scores["sum scores"][benchmark_type].index) 
    rank = [ i for i in benchmarking_scores["sum scores"]["Family"].index]

    plt.plot(rank, ali_structure.values, "b", label="alignment")
    plt.plot(rank, thr_structure.values,"#ffa201", label="threading")
    #plt.plot(rank, blo_structure.values, "#33b74b", label="blosum")
    plt.plot(rank, sum_structure.values, "r", label="sum scores")
    plt.title("Global scores comparison using "+ benchmark_type +" Benchmarks")
    plt.ylabel('Benchmarks')
    plt.xlabel('rank')
    plt.legend(loc='lower right')
    plt.axis([0, 413, 0, max(ali_structure+1)])
    plt.show()

