"""
    .. module:: Score
      :synopsis: This module implements the Score class.
"""

# Third-party modules
import os
import pandas as pd


def normalize_score(score_type):
    """
        Normalization of a score using the min-max scaling method (values between 0 and 1).

        Args:
            score_type (Pandas Series): A score

        Returns:
            Pandas Series: The score normalized
    """
    try:
        normalized = (score_type - min(score_type)) / (max(score_type) - min(score_type))
    except:
        normalized = 0
    return normalized
    # Use sklearn module instead
    # https://web.archive.org/web/20160520170701/http://chrisalbon.com:80/python/pandas_normalize_column.html


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
            Creates a scores.csv file containing for each line a template's name and its different
            scores. Creates the pdb files of the n best sum scores.

            Args:
                res_path (str): The path of the directory where to stock the created files.
                nb_pdb (int): Number of pdb to create using the n first templates.
                alignment_dict (dictionary): A dictionary containing Alignment objects.
        """
        os.makedirs(res_path+"/pdb", exist_ok=True)

        # A dataframe is created with pandas and elements of the iterator are stored
        scores_df = pd.DataFrame(columns=['benchmark', 'alignment', 'threading', 'modeller'])
        for _, ali_score, thr_score, modeller_score, _, name, benchmark in sorted(self.iterator):
            scores_df.loc[name] = [benchmark, ali_score, thr_score, modeller_score]

        # The first row is removed because it corresponds to the query
        scores_df = scores_df.drop(scores_df.index[0])
        # Normalization of the scores
        for index in ['alignment', 'threading', 'modeller']:
            scores_df[index] = normalize_score(scores_df[index])
        # Sum of the different scores and normalization
        scores_df['sum scores'] = normalize_score(scores_df['alignment']
                                                  + scores_df['threading']
                                                  + scores_df['modeller'])
        # Sort of the templates according to the sum score
        scores_df = scores_df.sort_values(by="sum scores", ascending=False)
        # A csv file containing the normalized scores is created
        scores_df.to_csv(res_path+"/scores.csv")

        # Only nb_pdb pdb files are created
        for i in range(nb_pdb):
            pdb_filename = res_path + "/pdb/top_" + str(i+1) + ".pdb"
            alignment_dict[scores_df.index[i]].write_pdb(pdb_filename)
