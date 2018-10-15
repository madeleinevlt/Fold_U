"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# IMPORTS
import re
import src.classes as cl
from Bio.PDB import * # Biopython

def metafold(metafold_file):
    """
        Extracts the name of the metafold as a key and the associated pdb as a value in a dictionary

        Args:
            metafold_file: The file METAFOLD.list containing this information

        Returns:
            metafold_dict: The dictionary with key = metafold and value = pdb file
    """
    metafold_dict = {}
    with open(metafold_file,"r") as f:
        for line in f:
            metafold_dict[line.split()[0]] = line.split()[1]
    return metafold_dict

def get_CA_coords(alignment_list,count):
    """
        Takes a pdb file and creates a list of Residue objects.
        Each Residue object contains the number the name and the CA coordinates.

        Args:
            pdb: The name of a pdb file

        Returns:
            residues_list: A list of Residue objects
    """
    p = PDBParser(QUIET=True) # QUIET = True : Warnings issued are suppressed
    pdb = p.get_structure(alignment_list[count].template.pdb, 'data/pdb/'+alignment_list[count].template.pdb)
    res_num = 0
    for atom in pdb.get_atoms():
        if atom.name == "CA":
            while alignment_list[count].template.seq[res_num].name == '-':
                res_num = res_num + 1
            alignment_list[count].template.seq[res_num].CA_coords = atom.get_vector()
            res_num = res_num + 1

def foldrec(nb_templates):
    """
        Extracts the score, the template name, the query sequence and the query template for each alignment
        and creates a list of alignment objects.

        Args:
            nb_template: Number of templates chosen

        Returns:
            alignment_list : A list of alignment objects
    """
    metafold_dict = metafold("data/METAFOLD.list")

    # Regex :
    template_name_reg = re.compile("Alignment :.*vs\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\s+([-0-9\.]+)")
    query_seq_reg = re.compile("^Query\s*[0-9]+\s*([A-Z-]+)")
    template_seq_reg = re.compile("^Template\s*[0-9]+\s*([A-Z-]+)")

    alignment_list = []
    count = 0

    with open("data/Agglutinin.foldrec","r") as f:
        prev_line = f.readline()
        for line in f:
            # The loop is break when the required nb of templates is reached :
            if (count == nb_templates):
                break
            # Search a regex for the current line :
            template_name_found = re.search(template_name_reg,line)
            score_found = re.search(score_reg,line)
            query_seq_found = re.search(query_seq_reg,line)
            template_seq_found = re.search(template_seq_reg,line)
            # A template name is found :
            if (template_name_found):
                template_name = template_name_found.group(1)
            # A score is found :
            if (score_found):
                score = float(score_found.group(1))
            # A query sequence is found :
            if (query_seq_found and prev_line == '\n'):
                query_seq = [cl.Residue(name) for name in list(query_seq_found.group(1))]
            # A template sequence is founds :
            if (template_seq_found and prev_line == '\n'):
                template_seq = [cl.Residue(name) for name in list(template_seq_found.group(1))]
                # Add a new alignment object in the list :
                alignment_list.append(cl.Alignment(score, query_seq, template_name, template_seq))
                alignment_list[count].template.pdb = metafold_dict[template_name]
                get_CA_coords(alignment_list,count)
                count = count + 1
            prev_line = line
    return(alignment_list)
