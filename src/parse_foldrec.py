from Bio.PDB import *
import re

def parsing_foldrec(foldrec,top):
    pattn = re.compile('Query')
    for q in re.finditer(pattn,foldrec):
        print(q)

        #"^\s+[0-9]+\s+([0-9]+\.[0-9]+)(\s|[0-9]|\.|E|-|\+)+([A-Za-z].*):"



if __name__ == "__main__":
    with open("data/Agglutinin.foldrec","r") as f:
        foldrec = f.read()


'''

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

'''