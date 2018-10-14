"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# IMPORTS
import re
import src.classes as classe


def metafold(metafold_file):
    """
        Extract the name of the metafold as a key and the associated pdb as a value in a dictionary

        Args:
            metafold_file: the file METAFOLD.list containing this information

        Returns:
            metafold_dict: the dictionary with key = metafold and value = pdb file
    """
    metafold_dict = {}
    with open(metafold_file,"r") as f:
        for line in f:
            metafold_dict[line.split()[0]] = line.split()[1]
    return(metafold_dict)


def foldrec(nb_templates, metafold_dict):
    """
        Extract the score, the template name, the query sequence and the query template for each alignment
        and create a list of alignment objects.
            
        Args:
            nb_template: Number of templates chosen

        Returns:
            alignment_list : A list of alignment objects
    """
    # Regex :
    templateName_reg = re.compile("Alignment :.*vs\s+([A-Za-z0-9-_]+)")
    score_reg = re.compile("^Score :\s+([-0-9\.]+)")
    querySeq_reg = re.compile("^Query\s*[0-9]+\s*([A-Z-]+)")
    templateSeq_reg = re.compile("^Template\s*[0-9]+\s*([A-Z-]+)")

    alignment_list = []
    count = 0

    with open("data/Agglutinin.foldrec","r") as f:
        prev_line = f.readline()
        for line in f:
            # The loop is break when the required nb of templates is reached :
            if (count == nb_templates):
                break
            # Search a regex for the current line :
            templateName_found = re.search(templateName_reg,line)
            score_found = re.search(score_reg,line)
            querySeq_found = re.search(querySeq_reg,line)
            templateSeq_found = re.search(templateSeq_reg,line)
            # A template name is found :
            if (templateName_found):
                templateName = templateName_found.group(1)
            # A score is found :
            if (score_found):
                score = float(score_found.group(1))
            # A query sequence is found :
            if (querySeq_found and prev_line == '\n'):
                querySeq = querySeq_found.group(1)
            # A template sequence is founds :
            if (templateSeq_found and prev_line == '\n'):
                templateSeq = templateSeq_found.group(1)
                # Add a new alignment object in the list :
                alignment_list.append(classe.Alignment(score, querySeq, templateName, templateSeq))
                alignment_list[count].template.pdbID = metafold_dict[templateName]
                count = count + 1
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