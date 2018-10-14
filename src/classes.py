"""
.. module:: classes
   :synopsis: This module implements Alignment, Template and Residue classes
"""

class Alignment:
    """Class that groups informations about an alignment."""
    def __init__(self, score, querySeq, templateName, templateSeq):
        """The constructor of an instance of the Alignment class."""
        self.score = score
        self.querySeq = querySeq
        self.template = Template(templateName, templateSeq)

class Template:
    """Class that groups informations about a template."""
    def __init__(self, name, seq):
        """The constructor of an instance of the Template class."""
        self.name = name
        self.seq = seq
        self.pdbID = None


class Residue:
    """Class that groups informations about a residue."""
    def __init__(self, chain, resNum):
        """The constructor of an instance of the Residue class."""
        self.resNum = resNum
        self.resName = chain[self.resNum].get_resname()
        self.CA_coords = chain[self.resNum]['CA'].get_vector()
