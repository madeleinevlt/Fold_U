"""
.. module:: Residue
   :synopsis: This module implements the Residue class.
"""


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
        self.ca_coords = None

    def __str__(self):
        if self.ca_coords is None:
            return "<" + self.name + " | " + "empty coordinates>"
        coords = [str(coord) for coord in self.ca_coords]
        return "<" + self.name + " | " + str(coords) + ">"

    def set_ca_coords(self, coords):
        """
            Sets the coordinates of the Residue.

            Args:
                coords (np array): x,y,z coordinates.

            Returns:
                void
        """
        self.ca_coords = coords
