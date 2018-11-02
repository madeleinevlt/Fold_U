"""
.. module:: Atom
   :synopsis: This module implements the Atom class.
"""


class Atom:
    """
    .. class:: Atom

      This class groups informations about an atom.

    Attributes:
        name: Name of the atom
        coords: 3D coordinates of the atom
    """

    def __init__(self, name):
        self.name = name
        self.coords = None

    def __str__(self):
        if self.coords is None:
            return "<" + self.name + " | " + "empty coordinates>"
        coords = [str(coord) for coord in self.coords]
        return "<" + self.name + " | " + str(coords) + ">"

    def set_coords(self, coords):
        """
            Sets the coordinates of the atom.

            Args:
                coords (np array): x,y,z coordinates.

            Returns:
                void
        """
        self.coords = coords
