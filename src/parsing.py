"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
import re
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1

# Local modules
from src.classes import Residue
from src.classes import Alignment



def metafold(metafold_file):
    """
        Extracts the name of the metafold as a key and the associated pdb as a value in a dictionary

        Args:
            file: A METAFOLD.list file containing for a given line a template name and a pdb file
            name associated

        Returns:
            dictionary: A dictionary with key = template name and value = pdb file
    """
    metafold_dict = {}  # Initalization of the dictionary
    with open(metafold_file, "r") as file:
        for line in file:
            metafold_dict[line.split()[0]] = line.split()[1]
    return metafold_dict


def dope(dope_file):
    """
        Extracts 10 energy values for the 20*20 residus-CA pairs and stores it
        in double indexed and labelled pandas DataFrame.

        Args:
            dope_file: The file dope_file.par containing energy values for each
            amino acid pair.

        Returns:
            DataFrame : amino acids indexed matrix containing lists of energies

     """
    # set up aa liste for rownames & colnames of the DataFrame
    aa1 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',\
    'S', 'T', 'V', 'W', 'Y']

    # set up matrix object 20*20
    dope_df = pd.DataFrame(index=aa1, columns=aa1, dtype=object)

    with open(dope_file, 'r') as file:
        for line in file:
            if line[4:6] == 'CA' and line[11:13] == 'CA':
                # get the line with C-alpha for both amino acids
                res_1 = protein_letters_3to1[line[0:3].title()]
                res_2 = protein_letters_3to1[line[7:10].title()]
                energy_res1_res2 = list(map(float, line[14:-1].split(" ")))
                dope_df[res_1][res_2] = energy_res1_res2
    return dope_df


def foldrec(foldrec_file, nb_templates, metafold_dict):
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
            alignment_list (list of Alignment): A list of Alignment objects.
    """

    # Regex :
    template_name_reg = re.compile("Alignment :.*vs\\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\\s+([-0-9\\.]+)")
    query_seq_reg = re.compile("^Query\\s*[0-9]+\\s*([A-Z-]+)")
    template_seq_reg = re.compile("^Template\\s*[0-9]+\\s*([A-Z-]+)")

    alignment_list = []
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
                ali.template.get_pdb(metafold_dict)
                ali.template.get_all_ca_coords()
                alignment_list.append(ali)
                count_templates = count_templates + 1
            prev_line = line
    return alignment_list
