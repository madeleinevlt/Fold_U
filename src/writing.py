"""
    .. module:: writing
      :synopsis: This module implements all the functions to create result files.
"""




def scores(filename, score_list):
    """
        Creates a file containing for each line a score and the template's
        name associated.

        Args:
            filename (str): The path of the file to create
            score_list (list): A list of tupples (score, template's name)

        Returns:
            void
    """
    print("\nThe program ended successfully !\nThe results are stored in " + filename)
    with open(filename, "w") as file:
        for score, name in score_list:
            file.write("{:<5d}{}\n".format(int(score), name))
