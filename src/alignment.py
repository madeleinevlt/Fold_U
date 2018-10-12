class Alignment:
    """Class that groups informations about an alignment."""
    def __init__(self, score, templateName):
        """The constructor of an instance of the Alignment class."""
        self.score = score
        self.templateName = templateName
        #self.query = query
        #self.template = template