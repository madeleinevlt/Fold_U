import re
import src.alignment as ali

def parsing_foldrec(nb_template):
    """
        Parsing of a foldrec file

        Args: nb_template : Number of template chosen

        Returns: List of alignment objects

    """
    score_template_name = re.compile("^\s+[0-9]+\s+([\-0-9]+\.[0-9]+)(\s|[0-9]|\.|E|-|\+)+([A-Za-z].*):")
    query = re.compile("^Query\s*[0-9]+\s*([A-Z-]+)")
    template = re.compile("^Template\s*[0-9]+\s*([A-Z-]+)")

    alignment_list = []
    flag = False
    cpt = 0

    with open("data/Agglutinin.foldrec","r") as f:
        prev_line = f.readline()
        for line in f:
            if (line == "*** ALIGNMENTS DETAILS ***\n"):
                flag = True
            # Parsing of all scores and template names
            if (flag == False):
                reg1 = re.search(score_template_name,line)
                if (reg1):
                    a = ali.Alignment(float(reg1.group(1)),reg1.group(3).strip())
                    alignment_list.append(a)
            # Parsing of all template and query sequences
            else:
                reg2 = re.search(query,line)
                reg3 = re.search(template,line)
                if (reg2):
                    if (prev_line == '\n'):
                        alignment_list[cpt].query = reg2.group(1)
                if (reg3):
                    if (prev_line == '\n'):
                        alignment_list[cpt].template = reg3.group(1)
                        cpt = cpt + 1
            prev_line = line
    return(alignment_list)


def get_pdb(pdb):
    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(pdb,pdb)

    residues_list = [] # List of Residue instances

    for chain in pdb.get_chains():
        # First residue of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))

        # Last residue of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
            last = res.get_id()[1]

        for resNum in range(first,last+1):
            r = Residue(chain, resNum)
            residues_list.append(r)
    return(residues_list)

if __name__ == "__main__":
    nb_template = 10
    parsing_foldrec(nb_template)
