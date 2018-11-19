"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
from conkit.applications import CCMpredCommandline
import subprocess
import re
import numpy as np
from Bio.SeqUtils import seq1

# Local modules
from src.residue import Residue
from src.alignment import Alignment
from src.query import Query
from src.template import Template


def parse_alignement(aln_file):
    """
       Extract and calculate the tops positions of residues based on
       co-evolution score

        Args:
            aln_file (str): multiple alignement file in .aln

        Returns:
            dict: A dictionary with key = ranking of tops positions and value = tuple of index of positions ~ [i,j]
    """
    top_position_dict = {}

    #Parsing aln file & Reindexing
    with open(aln_file, "r") as alnfile:
        query = alnfile.readline().split('\n')[0]

    #Reindexing
    list_index_pos_nongaps = [i for i, e in enumerate(query) if e != "-"] #save "non_gap" positions index
    #Calculs of co-ev score & generation of top_ouputs matrix
    ccmpred_cline = CCMpredCommandline(
        cmd ='./bin/CCMpred/bin/ccmpred', alnfile= aln_file, matfile= "contact.mat"
    )
    ccmpred_cline()
    #default_value for the tops ~ 30
    subprocess.call(
        "./bin/CCMpred/scripts/top_couplings.py contact.mat > top_output.mat", shell=True
    )

    #further ameliorations
    #directly recuperate values of top_output mat for the next operations in stdout
    #in order to avoid to reopen the file top_output
    #check subprocess.call options or subprocess.Popen

    #open top_output.mat
    with open("top_output.mat", "r") as tp_file:
        liste_data = tp_file.read().split('\n')
    #delete header and file_end character ""
        liste_data = liste_data[1:len(liste_data)-1]
        #verif si les positions tops ne correspondent pas des gaps_positions
        for i, value in enumerate(liste_data):
            value = value.split("\t")
            if int(value[0]) in list_index_pos_nongaps and int(value[1]) in list_index_pos_nongaps:
                pos_i = list_index_pos_nongaps.index(int(value[0]))
                pos_j = list_index_pos_nongaps.index(int(value[1]))
                top_position_dict[i] = (pos_i, pos_j)
    return top_position_dict


def parse_metafold(metafold_file):
    """
        Extracts the name of the metafold as a key and the associated pdb as
        a value in a dictionary

        Args:
            metafold_file (str): A METAFOLD.list file containing for a given line a
                           template name and a pdb file name associated

        Returns:
            dict: A dictionary with key = template name and value = pdb file
    """
    metafold_dict = {}  # Initalization of the dictionary
    with open(metafold_file, "r") as file:
        for line in file:
            metafold_dict[line.split()[0]] = line.split()[1]
    return metafold_dict


def parse_dope(dope_file):
    """
        Extracts 30 dope energy values for the 20*20 residus-CA pairs and stores it
        in a dictionary with key = res_1-res_2 and value = an array of the 30 dope scores.

        Args:
            dope_file (str): The file dope.par containing energy values for each
                       amino acid pair.

        Returns:
            dict: A dictionary with key = res_1-res_2 and value = an array of
            the 30 dope energy values.
    """
    dope_dict = {}
    with open(dope_file, "r") as file:
        for line in file:
            # get the line with C-alpha for both amino acids
            if line[4:6] == "CA" and line[11:13] == "CA":
                res_1 = seq1(line[0:3])
                res_2 = seq1(line[7:10])
                dope_dict[res_1+res_2] = np.fromstring(line[14:-1], dtype=float, sep=" ")
    return dope_dict


def parse_benchmark(benchmark_file, foldrec_file, alignment_dict):
    """
        Extract the line of the benchmark file containing the benchmarks of
        the studied query. Each benchmark template have an assigned structure
        corresponding to the degree of similarity with the query
        ("Family", "Superfamily" or "Fold"). For a given benchmark template
        object, the structure is stored in the benchmark item.

        Args:
            benchmark_file (str): The path to the benchmark file.
            foldrec_file (str): The path to the foldrec file.
            alignment_dict (dictionary): A dictionary with key = template name
                                         and value = an Alignment object.
    """
    # The name of the query is retrieved from the foldrec file pathway
    query_reg = re.compile("^.*\\/(\\S*)\\.")
    query_name = re.search(query_reg, foldrec_file)
    with open(benchmark_file, "r") as file:
        for line in file:
            # Only the line containing the query name is extracted
            if query_name.group(1) in line:
                # List containing for each benchmark the template name and
                # the associated fold_type
                template_fold_type = line[0:-1].split(": ")[1].split(", ")
                for i in template_fold_type:
                    # The fold_type of the benchmark template is stored in
                    # the benchmark item
                    template = i.split(" ")[0]
                    fold_type = i.split(" ")[1]
                    if template in alignment_dict:
                        alignment_dict[template].template.set_benchmark(fold_type)
                # No need to continue because all benchmarks are stored
                break


def parse_foldrec(foldrec_file, nb_templates, metafold_dict):
    """
        Extracts the score, the template name, the query and template sequences for
        the first n alignments from a file containing N profil-profil alignments and
        creates a list of Alignment objects. It also gives the template's pdb file name
        and gets all the coordinates of the CA atoms in the template's Residue list.

        Args:
            foldrec_file (str): The file containing N profil-profil alignments and their
                                corresponding scores.
            nb_templates (int): Number of alignments to retrieve from the file and chosen
                                by the user.
            metafold_dict (dictionary): A dictionary with key = template name
                                        and value = pdb file.

        Returns:
            dict: A dictionary with key = template name and value = an Alignment object.
    """
    # Regex :
    num_reg = re.compile("^No\\s*([0-9]+)")
    template_name_reg = re.compile("Alignment :.*vs\\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\\s+([-0-9\\.]+)")
    query_seq_reg = re.compile("^Query\\s*([0-9]+)\\s*([A-Z0-9-]+)\\s*([0-9]+)")
    template_seq_reg = re.compile("^Template\\s*[0-9]+\\s*([A-Z-]+)")

    alignment_dict = {}
    count_templates = 0

    with open(foldrec_file, "r") as file:
        query_reg_count = 0
        template_reg_count = 0
        for line in file:
            # The loop is break when the required nb of templates is reached :
            if count_templates == nb_templates:
                break
            # Search a regex for the current line :
            num_found = re.search(num_reg, line)
            template_name_found = re.search(template_name_reg, line)
            score_found = re.search(score_reg, line)
            query_seq_found = re.search(query_seq_reg, line)
            template_seq_found = re.search(template_seq_reg, line)
            # A num is found :
            if num_found:
                num = int(num_found.group(1))
            # A template name is found :
            if template_name_found:
                template_name = template_name_found.group(1)
            # A score is found :
            if score_found:
                score = float(score_found.group(1))
            # A query sequence is found :
            if query_seq_found:
                if query_reg_count == 0:
                    query_first = int(query_seq_found.group(1))
                    query_seq = [Residue(name) for name in list(query_seq_found.group(2))]
                    query_last = int(query_seq_found.group(3))
                    query_reg_count += 1
                elif query_reg_count == 1:
                    for ind, sec_struct in enumerate(list(query_seq_found.group(2))):
                        query_seq[ind].secondary_struct = sec_struct
                    query_reg_count += 1
                elif query_reg_count == 2:
                    for ind, ss_conf in enumerate(list(query_seq_found.group(2))):
                        if ss_conf != "-":
                            ss_conf = int(ss_conf)
                        query_seq[ind].ss_confidence = ss_conf
                    query_reg_count = 0
            # A template sequence is founds :
            if template_seq_found:
                if template_reg_count == 0:
                    template_seq = [Residue(name) for name in list(template_seq_found.group(1))]
                    template_reg_count += 1
                elif template_reg_count == 1:
                    for ind, sec_struct in enumerate(list(template_seq_found.group(1))):
                        template_seq[ind].secondary_struct = sec_struct
                    template_reg_count = 0
                    # Add a new alignment object in the list :
                    ali = Alignment(num, score,
                                    Query(query_seq, query_first, query_last),
                                    Template(template_name, template_seq))
                    ali.template.set_pdb_name(metafold_dict)
                    ali.template.parse_pdb()
                    alignment_dict[template_name] = ali
                    count_templates += 1
    return alignment_dict
