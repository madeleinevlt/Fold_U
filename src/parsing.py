from Bio.PDB import *
import re
# import alignment as ali

def parsing_foldrec():
    score_templateName = re.compile("^\s+[0-9]+\s+([0-9]+\.[0-9]+)(\s|[0-9]|\.|E|-|\+)+([A-Za-z].*):")
    aligmentList = []
    with open("data/Agglutinin.foldrec","r") as f:
        for line in f:
            for i in re.finditer(score_templateName,line):
                a = Alignment(i.group(1),i.group(3).strip(),template,cible)
                aligmentList.append(a)

if __name__ == "__main__":
    parsing_foldrec()

"""
    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(opt.input,opt.input+".H")

    resList = [] # List of Residue instances

    for chain in pdb.get_chains():
        # First residue of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))

        # Last residue of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
            last = res.get_id()[1]

        for resNum in range(first,last+1):
            r = classes.Residue(chain, index, chain.id, resNum)
            resList.append(r)
"""