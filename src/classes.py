class Alignment:
    def __init__(self, chain, resNum):

class Residue:
    def __init__(self, chain, resNum):
        self.chain = chain.id
        self.resNum = resNum
        self.resName = chain[resNum].get_resname()
