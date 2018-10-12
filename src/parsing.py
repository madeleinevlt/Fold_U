import re
import alignment as ali

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
                        alignment_list[cpt].query = reg3.group(1)
                        cpt = cpt + 1      
            prev_line = line
    return(alignment_list)


if __name__ == "__main__":
    nb_template = 10
    parsing_foldrec(nb_template)
