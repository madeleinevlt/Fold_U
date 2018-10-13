def metafold(metafold_file):
    """
        Extract the name of the meatafold as a key and the associated pdb as a value in a dictionary

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