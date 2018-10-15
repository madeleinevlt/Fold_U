"""
.. module:: classes
   :synopsis: This module implements Alignment, Template and Residue classes
"""


class Alignment:
    """
    .. class:: Alignment
      This class groups informations about an alignment.

    Attributes:
        score: Alignment score
        query_seq: Query sequence
        template: Instance of a Template object.
        .. code-block:: python
            Template(template_name, template_seq)
    """

    def __init__(self, score, query_seq, template_name, template_seq):
        """The constructor of an instance of the Alignment class."""
        self.score = score
        self.query_seq = query_seq
        self.template = Template(template_name, template_seq)

class Template:
    """
    .. class:: Template
      This class groups informations about a template sequence/structure.

    Attributes:
        name: Name of the template
        seq: Template's sequence of residues as list of Residues objects
        pdb: PDB code of the template
    """

    def __init__(self, name, seq):
        """The constructor of an instance of the Template class."""
        self.name = name
        self.seq = seq
        self.pdb = None

class Residue:
    """
    .. class:: Residue
      This class groups informations about a residue.

    Attributes:
        name: Name of the residue (3 letters code)
        CA_coords: 3D coordinates of the residues
    """

    def __init__(self, name):
        """The constructor of an instance of the Residue class."""
        self.name = name
        self.CA_coords = None
