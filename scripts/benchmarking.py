#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    This module runs the program "fold_u" on all of the queries IF their results do not already
    exist respectively (checks if scores.csv is generated in results/foldrec_name/scores.csv).
    It then generates a plot to visualize the contribution of each and every score to the
    re-ranking of models/templates. You can choose to see the statistics for one particular score,
    or all scores combined (summed and normalized), or all the scores at the same time.
    A table is also printed in the terminal for the TOP N statistics, presenting the number and
    percentage of benchmarks for the TOP N found.

    Usage:
        ./script/benchmarking.py UNIREF_DB [--selected_score SCORE] [--cpu NUM] [--output PATH]

    Arguments:
        UNIREF_DB                             Path to Uniref database.

    Options:
        -h, --help                            Show this
        -s SCORE, --selected_score SCORE      Score for which you wish to see the statistics:
                                              "alignment", "threading", "modeller",
                                              "secondary_structure", "solvent_access"
                                              or "sum_scores",
                                              or all of them at once: "all" [default: all]
        -c NUM, --cpu NUM                     Number of cpus to use for parallelisation. By default
                                              using all available (0).
                                              [default: 0]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and plot)
                                              [default: ./results/plots]
"""

# Third-party modules
import os
import subprocess
from multiprocessing import cpu_count
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cycler
from docopt import docopt
from schema import Schema, And, Use, SchemaError

# plot settings
COLORS = cycler('color', ['#EE6666', '#3388BB', '#9988DD', '#EECC55', '#88BB44', '#FFBBBB'])
plt.rc('axes', facecolor='#E6E6E6', edgecolor='none',
       axisbelow=True, grid=True, prop_cycle=COLORS)
plt.rc('grid', color='w', linestyle='solid')
plt.rc('xtick', direction='out', color='gray')
plt.rc('ytick', direction='out', color='gray')
plt.rc('patch', edgecolor='#E6E6E6')
plt.rc('lines', linewidth=1.5)

def check_args():
    """
        Checks and validates the types of inputs parsed by docopt from command line.
    """
    schema = Schema({
        '--selected_score': And(Use(str), lambda s: s in ["alignment", "threading",
                                                          "modeller", "secondary_structure",
                                                          "solvent_access", "sum_scores", "all"],
                                error='SCORES should be an existing score'),
        '--cpu': And(Use(int), lambda n: 0 <= n <= cpu_count(),
                     error='--cpus=NUM should be integer 1 <= N <= ' + str(cpu_count())),
        # The output PATH is created (if not exists) at the end of the program so we skip the check.
        object: object})
    try:
        schema.validate(ARGUMENTS)
    except SchemaError as err:
        exit(err)

def run_command(command):
    """
        Run the FOLD_U program. It keeps on reading the stdout, checks for the return code and
        displays the output messages and the progress bar in real time.

        Args:
            command (str): The command line to run
    """
    # Run the command
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    # This part is to print the running program messages
    while True:
        line = process.stdout.readline()
        if not line:
            break
        else:
            print(line.decode('utf-8').strip())

def create_benchmarking_scores_dict(uniref, scores, structures, nb_proc):
    """
        Create a dictionary of scores with key = a score, value = a pandas dataframe
        which contains the cumulative sum of benchmark for each benchmark type and for
        all benchmarks (=4 columns).

        Args:
            uniref (str): Path to the uniref90 database
            scores (list): A list of score name
            structures (list): List containing fold types: "Family", "Superfamily", "Fold"
            nb_proc (int): Number of processors.

        Returns:
            dict: The dictionary of scores with key = a score, value = a pandas dataframe.
    """
    benchmarking_scores = {}
    # We search the minimum rank for the plot
    min_rank = 900
    for score in scores:
        benchmarking_scores[score] = pd.DataFrame(np.zeros((405, 3), dtype=int), columns=structures)
    # For each query,
    all_foldrecs = os.listdir("data/queries")
    print("\n\nProcessing all benchmarks ...\n")
    for ind, query in enumerate(all_foldrecs, 1):
        query = query.split(".")[0]
        # The Fold_U program is run on the current query if results are not already generated
        if not os.path.isfile("results/" + query + "/scores.csv"):
            print("\nProcessing query {} / {} : {}\n".format(ind, len(all_foldrecs), query))
            run_command(["./fold_u", "data/queries/" + query + "/"+ query + ".fasta", uniref, "-o",
                         "results/" + query, "-c", str(nb_proc)])
            _, columns = os.popen('stty size', 'r').read().split()
            print("\n" + "-" * int(columns))
        # Score results are stored in a pandas DataFrame
        query_scores = pd.read_csv("results/" + query + "/scores.csv", index_col=0)
        if len(query_scores) < min_rank:
            min_rank = len(query_scores)
        for score in scores:
            # The DataFrame is sorted by the current score
            query_score = query_scores.sort_values(by=score, ascending=False)
            # Initialization of the dictionary of counts
            structures_count = {}
            for structure in structures:
                structures_count[structure] = 0
            # Cumulative sum of benchmark according to the structure
            for i, struct_type in enumerate(query_score["benchmark"]):
                for structure in structures:
                    if struct_type == structure:
                        structures_count[structure] += 1
                    benchmarking_scores[score][structure][i+1] += structures_count[structure]
            benchmarking_scores[score]["total"] = benchmarking_scores[score]\
                    .apply(lambda row: row.Family+ row.Superfamily + row.Fold, axis=1)
    print("done.\n")
    return (benchmarking_scores, min_rank)

def plot_benchmark(output_path, min_rank, scores, benchmarking_scores, selected_score):
    """
        Create one plot for one benchmark type for all the foldrec files.

        Args:
            output_path (str): The path to store png file.
            scores (list): A list of score name.
            min_rank (int): The last value of the rank.
            benchmarking_scores (dict): A dictionnary containing benchmarking
                                        data for each type of score for
                                        different types of structures.
            selected_score (str): All scores or selected score for the plot.

    """
    rank = [i for i in range(0, min_rank)]
    os.makedirs(output_path, exist_ok=True)
    # Plot creation
    plt.figure(num="Enrichment")  # Window's name
    plt.ylabel("Benchmark")
    plt.xlabel("rank")

    # Plot all scores
    if selected_score == "all":
        ali_struct = benchmarking_scores[scores[0]]["total"].values[:min_rank]
        thr_struct = benchmarking_scores[scores[1]]["total"].values[:min_rank]
        mod_struct = benchmarking_scores[scores[2]]["total"].values[:min_rank]
        ss_struct = benchmarking_scores[scores[3]]["total"].values[:min_rank]
        acc_struct = benchmarking_scores[scores[4]]["total"].values[:min_rank]
        co_ev_struct = benchmarking_scores[scores[5]]["total"].values[:min_rank]
        sum_struct = benchmarking_scores[scores[6]]["total"].values[:min_rank]

        plt.plot(rank, ali_struct, "b", label=scores[0])
        plt.plot(rank, thr_struct, "#ffa201", label=scores[1])
        plt.plot(rank, mod_struct, "#EE82EE", label=scores[2])
        plt.plot(rank, ss_struct, "#00B200", label=scores[3])
        plt.plot(rank, acc_struct, "#7a9a91", label=scores[4])
        plt.plot(rank, co_ev_struct, "#660033", label=scores[5])
        plt.plot(rank, sum_struct, "r", label=scores[6])
        plt.plot([0, len(ali_struct)], [0, max(ali_struct)], "k", label="random")
        plt.legend(loc="lower right")
        plt.title("Enrichment plot")
        plt.savefig(output_path + "/" + "all_scores_plot.png")
    # Plot scores individually
    else:
        score_struct = benchmarking_scores[selected_score]["total"].values[:min_rank]
        plt.plot(rank, score_struct, "b", label=selected_score)
        plt.plot([0, len(score_struct)], [0, max(score_struct)], "k", label="random")
        plt.title("Enrichment plot of " + selected_score + " score")
        plt.legend(loc="lower right")
        plt.savefig(output_path + "/" + selected_score + "_plot.png")

def print_top_n(selected_score, top, benchmarking_scores):
    """
        Show statistics based on the benchmark.list files separately for each fold-type: "Fold",
        "Family", "Superfamily".
        Represent the strength/weaknesses of the different scores independantly and/or combined.

        Args:
            selected_score (str): The score you want some stats on.
            top (str): a maximum rank number.
            benchmarking_scores (dict): a dictionnary containing benchmarking data for each type
                                        of score for different types of structures.

        Returns:
            a str "top_results" table summarizing the top results.

    """
    rank = {}
    max_rank = {}
    if selected_score == "all":
        selected_score = "sum_scores"
    for struct in benchmarking_scores[selected_score].columns.values:
        rank[struct] = benchmarking_scores[selected_score][struct][top-1]
        max_rank[struct] = max(benchmarking_scores[selected_score][struct])
    line1 = "top {}\t\t{}/{}\t\t{}/{}\t\t{}/{}\t\t{}/{}\n"\
             .format(top, rank["Family"], max_rank["Family"], rank["Superfamily"],
                     max_rank["Superfamily"], rank["Fold"], max_rank["Fold"], rank["total"],
                     max_rank["total"])
    line2 = "\t\t{:<5.1f}%\t\t{:<5.1f}%\t\t{:<5.1f}%\t\t{:<5.1f}%"\
        .format((rank["Family"]/max_rank["Family"])*100,
                (rank["Superfamily"]/max_rank["Superfamily"])*100,
                (rank["Fold"]/max_rank["Fold"])*100,
                (rank["total"]/max_rank["total"])*100)
    top_results = line1 + line2
    return top_results

def print_table(selected_score, benchmarking_scores):
    """
        Print a table of top N.

        Args:
            structures (list): List containing fold types: "Family", "Superfamily", "Fold".
            selected_score (str): The score you want some stats on.
            benchmarking_scores (dict): a dictionnary containing benchmarking data for each type
                                        of score for different types of structures.
    """
    top_n_list = [5, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300, 350]
    print("\nTable summarizing the top N results.\n")
    print("\t\tFamily\t\tSuperfamily\tFold\t\tTotal\n")
    for top_n in top_n_list:
        print(print_top_n(selected_score, top_n, benchmarking_scores))
        print("-----------------------------------------------------------------------")

if __name__ == "__main__":
    START_TIME = datetime.now()
    ### Parse command line
    ######################
    ARGUMENTS = docopt(__doc__)
    # Check the types and ranges of the command line arguments parsed by docopt
    check_args()

    # Uniref database
    UNIREF = ARGUMENTS["UNIREF_DB"]
    # OUTPUT file
    OUTPUT_PATH = ARGUMENTS["--output"]
    # Number of cpus for parallelisation
    NB_PROC = cpu_count() if int(ARGUMENTS["--cpu"]) == 0 else int(ARGUMENTS["--cpu"])
    # Selected score you want to have info about
    SELECTED_SCORE = ARGUMENTS["--selected_score"]
    # The 3 different structures from benchmark
    STRUCTURES = ["Family", "Superfamily", "Fold"]
    # all the possible scores useful for plots
    SCORES = ["alignment", "threading", "modeller", "secondary_structure", "solvent_access",
              "co_evolution", "sum_scores"]

    (BENCHMARKING_SCORES, MIN_RANK) = create_benchmarking_scores_dict(UNIREF, SCORES, STRUCTURES,
                                                                      NB_PROC)

    print_table(SELECTED_SCORE, BENCHMARKING_SCORES)
    plot_benchmark(OUTPUT_PATH, MIN_RANK, SCORES, BENCHMARKING_SCORES, SELECTED_SCORE)
    print("\nTotal runtime: {} seconds".format(str(datetime.now() - START_TIME)))
    plt.show()
    print("\n\nThe plot is stored in " + OUTPUT_PATH + "\n")
