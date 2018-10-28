"""
.. module:: Query
   :synopsis: This module implements the Query class.
"""


class Query:
    """
    .. class:: Query

      This class groups informations about a query sequence.

    Attributes:
        residues: Template's sequence of residues as list of Residues objects
        first (int): First residue of the query sequence.
        last (int): Last residue of the query sequence.
    """

    def __init__(self, residues, first, last):
        self.residues = residues
        self.first = first
        self.last = last


    def display(self):
        return "".join(str(res.name) for res in self.residues)
