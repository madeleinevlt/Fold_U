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
        name (str): Name of the residue (1 letter code)
        ca_coords (Numpy array): 3D coordinates of the residue
    """

    def __init__(self, name):
        self.name = name
        self.ca_atom = Atom("CA")
        self.c_atom = Atom("C")
        self.n_atom = Atom("N")
        self.ss = None
        self.ss_confidence = None

    def calculate_distance(self, residue):
        """
            Calculate Euclidian distance between two residues with THE most efficient method.

            :math: `distance = \\sqrt{(x_a-x_b)^2+(y_a-y_b)^2+(z_a-z_b)^2}`

            Args:
                residue (object): An object of the Residue class.

            Returns:
                float: The calculated distance.
        """
        a_min_b = self.ca_atom.coords - residue.ca_atom.coords
        dist = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
        return dist
