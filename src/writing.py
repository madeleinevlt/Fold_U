"""
    .. module:: writing
      :synopsis: This module implements all the functions to create result files.
"""

# Third-party modules
import os

def write_scores(res_path, score_list, nb_pdb, alignment_dict):
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
                pdb_filename = res_path + "/pdb/top_" + str(template_num) + ".pdb"
                alignment_dict[name].write_pdb(pdb_filename)
                template_num += 1
    print("\nThe program ended successfully !\nThe results are stored in " + res_path)
