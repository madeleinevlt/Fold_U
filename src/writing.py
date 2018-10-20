"""
    .. module:: writing
      :synopsis: This module implements all the functions to create result files.
"""

# Third-party modules
import os
from Bio.SeqUtils import seq3


def scores_and_pdb(res_path, score_list, nb_pdb, alignment_dict):
    """
        Creates a scores.out file containing for each line a score and the template's
        name associated. Creates the pdb files of the n best templates (best scores).

        Args:
            res_path (str): The path of the directory where to stock the created files.
            score_list (list): A list of tupples (score, template's name)
            nb_pdb (int): Number of pdb to create using the n first templates.
            alignment_dict (dictionary): A dictionary containing Alignment objects.

        Returns:
            void
    """
    os.makedirs(res_path, exist_ok=True)
    os.makedirs(res_path+"/pdb", exist_ok=True)
    with open(res_path+"/scores.out", "w") as file:
        template_num = 1
        for score, name in score_list:
            # Write a line in the score.out file containing the
            # score and the name of the current ranked template
            file.write("{:<5d}{}\n".format(int(score), name))
            # Only nb_pdb pdb files are created
            if (template_num <= nb_pdb):
                write_pdb(res_path, name, template_num, alignment_dict)
                template_num = template_num + 1
    print("\nThe program ended successfully !\nThe results are stored in " + res_path)


def write_pdb(res_path, template_name, template_num, alignment_dict):
    """
        Write a pdb file by threading the query sequence on the template CA coordinates.

        Args:
            res_path (str): The path of the directory where to stock the created files.
            template_name (str): Name of the template used.
            template_num (int): Rank of the template used.
            alignment_dict (dictionary): A dictionary containing Alignment objects.

        Returns:
            void
    """
    with open(res_path+"/pdb/template_"+str(template_num)+".pdb", "w") as file:
        query = alignment_dict[template_name].query_residues
        template = alignment_dict[template_name].template.residues
        # Extra informations on the template used to generate the pdb file
        file.write("REMARK Threading of query sequence on the Template {:d} ({:s}).\n".format(template_num, template_name))
        for res_num, res_q in enumerate(query):
            res_t = template[res_num]
            if (res_q == "-" or res_t == "-"):
                continue
            # An "ATOM" line of the created pdb file 
            file.write("{:6s}{:5d} {:^4s} {:3s}{:>2s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
            .format("ATOM", res_num+1, "CA", seq3(res_q.name).upper(), "A", res_num+1, res_t.ca_coords[0], res_t.ca_coords[1], res_t.ca_coords[2], 1.00, 0, "C"))
        # The two last lines of the created pdb file ("END" and "TER" lines)
        file.write("{:6s}{:5d}{:>9s}{:>2s}{:4d}\nEND\n".format("TER", res_num+1, seq3(res_q.name).upper(), "A", res_num+1))
