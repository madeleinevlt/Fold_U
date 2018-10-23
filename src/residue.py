"""
.. module:: Residue
   :synopsis: This module implements the Residue class.
"""

# Third-party modules
import numpy as np

# Local modules
from src.atom import Atom


class Residue:
    """
    .. class:: Residue

      This class groups informations about a residue.

    Attributes:
        name: Name of the residue (1 letter code)
        ca_coords: 3D coordinates of the residue
    """

    def __init__(self, name):
        self.name = name
        self.ca_atom = Atom("CA")
        self.cb_atom = Atom("CB")
        self.n_atom = Atom("N")

    def calculate_distance(self, residue):
        """
            Calculate Euclidian distance between two residues with THE most efficient method.
            Formula: distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)

            Args:
                residue (object): An object of the Residue class.

            Returns:
                dist (float): The calculated distance.
        """
        a_min_b = self.ca_atom.coords - residue.ca_atom.coords
        dist = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
        return dist
