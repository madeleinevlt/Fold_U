"""
.. module:: classes
   :synopsis: This module implements Alignment, Template and Residue classes
"""

class Alignment:
    """Class that groups informations about an alignment."""
    def __init__(self, score, query_seq, template_name, template_seq):
        """The constructor of an instance of the Alignment class."""
        self.score = score
        self.query_seq = query_seq
        self.template = Template(template_name, template_seq)

class Query:
    """Class that groups informations about a template."""
    def __init__(self, name, seq):
        """The constructor of an instance of the Template class."""



class Template:
    """Class that groups informations about a template."""
    def __init__(self, name, seq):
        """The constructor of an instance of the Template class."""
        self.name = name
        self.seq = seq
        self.pdb = None

class Residue:
    """Class that groups informations about a residue."""
    def __init__(self, res_name):
        """The constructor of an instance of the Residue class."""
        self.res_name = res_name
        self.CA_coords = None
