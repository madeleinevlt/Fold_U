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
    os.makedirs(res_path+"/pdb", exist_ok=True)
    with open(res_path+"/scores.out", "w") as file:
        template_num = 1
        for score, name in score_list:
            # Write a line in the score.out file containing the
            # score and the name of the current ranked template
            file.write("{:<8d}{}\n".format(int(score), name))
            # Only nb_pdb pdb files are created
            if template_num <= nb_pdb:
                pdb_filename = res_path+"/pdb/template_"+str(template_num)+".pdb"
                write_pdb(pdb_filename, alignment_dict[name])
                template_num += 1
    print("\nThe program ended successfully !\nThe results are stored in " + res_path)


def write_pdb(pdb_filename, ali):
    """
        Write a pdb file by threading the query sequence on the template CA coordinates.

        Args:
            pdb_filename (str): Name of the pdb file to create.
            ali (object): Alignment object used.

        Returns:
            void
    """
    with open(pdb_filename, "w") as file:
        # Extra informations on the template used to generate the pdb file
        file.write("REMARK Threading of query sequence on the {:s} template.\n"\
            .format(ali.template.name))
        for res_num, res_q in enumerate(ali.query_residues):
            res_t = ali.template.residues[res_num]
            if res_num > len(ali.template.residues) or res_num == len(ali.template.residues) - 1:
                break
            if res_q.name == "-" or res_t.name == "-":
                continue
            # An "ATOM" line of the created pdb file
            file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
                .format("ATOM", res_num+1, "CA", seq3(res_q.name).upper(), "A", res_num+1,\
                res_t.ca_coords[0], res_t.ca_coords[1], res_t.ca_coords[2], 1.00, 0, "C"))
        # The two last lines of the created pdb file ("END" and "TER" lines)
        file.write("{:6s}{:5d}{:>9s}{:>2s}{:4d}\nEND\n"\
            .format("TER", res_num+1, seq3(res_q.name).upper(), "A", res_num+1))
