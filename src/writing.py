"""
    .. module:: writing
      :synopsis: This module implements all the functions to create result files.
"""




def scores(filename, score_list):
    """
        Creates a file containing for each line a score and the template's
        name associated.

        Args:
            filename (str): The path of the file the create
            score_list (list): A list of tupples (score, template's name)

        Returns:
            void
    """
    with open(filename, "w") as file:
        for score, name in score_list:
            file.write("{}\t{}\n".format(round(score, 1), name))
