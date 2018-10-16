"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# IMPORTS
import re
from Bio.PDB import PDBParser
import src.classes as cl

def metafold(metafold_file):
    """
        Extracts the name of the metafold as a key and the associated pdb as a value in a dictionary

        Args:
            file: A METAFOLD.list file containing for a given line a template name and a pdb file
            name associated

        Returns:
            dictionary: A dictionary with key = template name and value = pdb file
    """
    metafold_dict = {} # Initalization of the dictionary
    with open(metafold_file, "r") as file:
        for line in file:
            metafold_dict[line.split()[0]] = line.split()[1]
    return metafold_dict

def get_ca_coords(alignment):
    """
        Takes a pdb file and creates a list of Residue objects.
        Each Residue object contains the number the name and the CA coordinates.

        Args:
            alignment: An Alignment object

        Returns:
            void:
    """
    print(alignment.score)
    pdb = PDBParser(PERMISSIVE=2, QUIET=True) # QUIET = True : Warnings issued are suppressed

    structure = pdb.get_structure(alignment.template.pdb, 'data/pdb/'+alignment.template.pdb)

    for atom in structure.get_atoms():
        if atom.name == "CA":
            for res in alignment.template.residues:
                if res != "-":
                    res.ca_coords = atom.get_vector()


def foldrec(foldrec_file, nb_templates, metafold_dict):
    """
        Extracts the score, the template name, the query sequence and the query template
        for each alignment and creates a list of alignment objects.

        Args:
            nb_template: Number of templates chosen

        Returns:
            list of Alignment: A list of alignment objects
    """

    # Regex :
    template_name_reg = re.compile("Alignment :.*vs\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\s+([-0-9\.]+)")
    query_seq_reg = re.compile("^Query\s*[0-9]+\s*([A-Z-]+)")
    template_seq_reg = re.compile("^Template\s*[0-9]+\s*([A-Z-]+)")

    alignment_list = []

    with open(foldrec_file, "r") as file:
        prev_line = file.readline()
        for line in file:
            # The loop is break when the required nb of templates is reached :
            # if count == nb_templates:
            #     break
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
                query_seq = [cl.Residue(name) for name in list(query_seq_found.group(1))]
            # A template sequence is founds :
            if template_seq_found and prev_line == '\n':
                template_seq = [cl.Residue(name) for name in list(template_seq_found.group(1))]
                # Add a new alignment object in the list :
                ali = cl.Alignment(score, query_seq, template_name, template_seq)
                ali.template.pdb = metafold_dict[template_name]
                alignment_list.append(ali)
                get_ca_coords(ali)
            prev_line = line
    return alignment_list
