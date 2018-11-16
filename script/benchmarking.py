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
