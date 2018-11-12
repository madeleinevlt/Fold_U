"""
.. module:: Query
   :synopsis: This module implements the Query class.
"""


class Query:
    """
    .. class:: Query

      This class groups informations about a query sequence.

    Attributes:
        residues (list of Residue objects): Template's sequence of residues as list of Residues
                                            objects
        first (int): First residue of the query sequence.
        last (int): Last residue of the query sequence.
    """

    def __init__(self, residues, first, last):
        self.residues = residues
        self.first = first
        self.last = last

    def get_size(self):
        """
            Get the size of the query sequence.

            Returns:
                int: The size of the query sequence.

        """
        return len(self.residues)
