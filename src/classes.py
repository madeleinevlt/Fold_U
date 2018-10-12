class Alignment:
	"""Class that groups informations about an alignment."""
    def __init__(self, score, idSCOP, query, template):
    	"""The constructor of an instance of the Alignment class."""
    	self.score =score
    	self.idSCOP = idSCOP
        self.query = query
        self.template = template

class Residue:
	"""Class that groups informations about a residue."""
    def __init__(self, chain, resNum):
    	 """The constructor of an instance of the Residue class."""
        self.resNum = resNum
        self.resName = chain[resNum].get_resname()
        self.CA = chain[self.resNum]['CA'].get_vector()