"""
    .. module:: Score
      :synopsis: This module implements the Score class.
"""

# Third-party modules
import os


class Score:
    """
    .. class:: Score

      This class groups informations about a score.

    Attributes:
        iterator (iterator): An iterator of the generated scores.
    """

    def __init__(self, iterator):
        self.iterator = iterator

    def write_score(self, res_path, nb_pdb, alignment_dict):
        """
            Creates a scores.out file containing for each line a score and the template's
            name associated. Creates the pdb files of the n best templates (best scores).

            Args:
                res_path (str): The path of the directory where to stock the created files.
                nb_pdb (int): Number of pdb to create using the n first templates.
                alignment_dict (dictionary): A dictionary containing Alignment objects.
        """
        os.makedirs(res_path+"/pdb", exist_ok=True)
        with open(res_path+"/scores.out", "w") as file:
            template_num = 1
            file.write("{:<18s} {:<11s} {:<9s} {:<9s} {:<9s}\n"\
                .format("template_name", "benchmark", "alignment", "threading", "blosum"))
            for _, ali_score, thr_score, blosum_score, name, benchmark in sorted(self.iterator):
                # Write a line in the score.out file containing the
                # score and the name of the current ranked template
                file.write("{:<18s} {:<11s} {:<9d} {:<9d} {:<9d}\n"\
                    .format(name, benchmark, int(ali_score), int(thr_score), int(blosum_score)))
                # Only nb_pdb pdb files are created
                if template_num <= nb_pdb:
                    pdb_filename = res_path + "/pdb/top_" + str(template_num) + ".pdb"
                    alignment_dict[name].write_pdb(pdb_filename)
                    template_num += 1
