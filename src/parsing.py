"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
import re
import numpy as np
from Bio.SeqUtils import seq1

# Local modules
from src.residue import Residue
from src.alignment import Alignment


def parse_metafold(metafold_file):
    """
        Extracts the name of the metafold as a key and the associated pdb as
        a value in a dictionary

        Args:
            metafold_file: A METAFOLD.list file containing for a given line a
                           template name and a pdb file name associated

        Returns:
            dictionary: A dictionary with key = template name and value = pdb file
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
            dope file: The file dope.par containing energy values for each
                       amino acid pair.

        Returns:
            dictionary: A dictionary with key = res_1-res_2 and
                        value = an array of the 30 dope energy values.
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


def parse_foldrec(foldrec_file, nb_templates, metafold_dict):
    """
        Extracts the score, the template name, the query and template sequences for
        the first n alignments from a file containing N profil-profil alignments and
        creates a list of Alignment objects. It also gives the template's pdb file name
        and gets all the coordinates of the CA atoms in the template's Residue list.

        Args:
            foldrec_file (file): The file containing N profil-profil alignments and their
                                 corresponding scores.
            nb_templates (int): Number of alignments to retrieve from the file and chosen
                                by the user.
            metafold_dict (dictionary): A dictionary with key = template name and
                                        value = pdb file.

        Returns:
            Dictionary of Alignments: A dictionary with key = template name and
                                      value = an Alignment object.
    """

    # Regex :
    template_name_reg = re.compile("Alignment :.*vs\\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\\s+([-0-9\\.]+)")
    query_seq_reg = re.compile("^Query\\s*[0-9]+\\s*([A-Z-]+)")
    template_seq_reg = re.compile("^Template\\s*[0-9]+\\s*([A-Z-]+)")

    alignment_dict = {}
    count_templates = 0

    with open(foldrec_file, "r") as file:
        prev_line = file.readline()
        for line in file:
            # The loop is break when the required nb of templates is reached :
            if count_templates == nb_templates:
                break
            # Search a regex for the current line :
            template_name_found = re.search(template_name_reg, line)
            score_found = re.search(score_reg, line)
            query_seq_found = re.search(query_seq_reg, line)
            template_seq_found = re.search(template_seq_reg, line)
            # A template name is found :
            if template_name_found:
                template_name = template_name_found.group(1)
            # A score is found :
            if score_found:
                score = float(score_found.group(1))
            # A query sequence is found :
            if query_seq_found and prev_line == '\n':
                query_seq = [Residue(name) for name in list(query_seq_found.group(1))]
            # A template sequence is founds :
            if template_seq_found and prev_line == '\n':
                template_seq = [Residue(name) for name in list(template_seq_found.group(1))]
                # Add a new alignment object in the list :
                ali = Alignment(score, query_seq, template_name, template_seq)
                ali.template.set_pdb_name(metafold_dict)
                ali.template.parse_pdb()
                alignment_dict[template_name] = ali
                count_templates += 1
            prev_line = line
    return alignment_dict
