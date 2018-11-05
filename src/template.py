"""
.. module:: Template
   :synopsis: This module implements the Template class.
"""

# Third-party modules
import numpy as np


class Template:
    """
    .. class:: Template

      This class groups informations about a template sequence/structure.

    Attributes:
        name (str): Name of the template
        residues (list of Residue object): Template's sequence of residues as list of Residues objects
        pdb (str): PDB filename of the template
    """

    def __init__(self, name, residues):
        self.name = name
        self.residues = residues
        self.pdb = None
        self.benchmark = ""

    def set_pdb_name(self, metafold_dict):
        """
            Get the PDB file name of the current template from the template's name.

            Args:
                metafold_dict: A dictionary with key = template name and value = pdb file
        """
        self.pdb = metafold_dict[self.name]

    def parse_pdb(self):
        """
            Parse the pdb file and set the CA coordinates.
        """
        count_res = 0
        nb_atoms = 0
        with open("data/pdb/" + self.name + "/" + self.pdb, 'r') as file:
            for line in file:
                line_type = line[0:6].strip()
                name_at = line[12:16].strip()
                if line_type == "ATOM" and (name_at == "N" or name_at == "CA" or name_at == "C"):
                    x_coord = float(line[30:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:54].strip())
                    if count_res <= len(self.residues):
                        # Skip gaps in the template
                        while self.residues[count_res].name == "-":
                            count_res += 1
                        if line_type == "ATOM" and name_at == "N":
                            self.residues[count_res].ca_atom\
                            .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        elif line_type == "ATOM" and name_at == "CA":
                            self.residues[count_res].c_atom\
                            .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        elif line_type == "ATOM" and name_at == "C":
                            self.residues[count_res].n_atom\
                            .set_coords(np.array([x_coord, y_coord, z_coord]))
                            nb_atoms += 1
                        if nb_atoms == 3:
                            count_res += 1
                            nb_atoms = 0

    def get_fasta_file(self):
        """
        Retrieve the full FASTA amino acid sequence (gapless) of the template,
        and write it to an alignment file.

        Returns:
            str: Path to the FASTA file containing the sequence of the template.
        """
        url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=" + self.pdb + "&compressionType=uncompressed"
        file_name = None
        try:
            file_name = wget.download(url, bar=None)
        except Exception as e:
            print("Errors encountered while downloading FASTA file of "
                  + self.pdb + ".atm template.")
            print("Initial error: " + str(e))
        return file_name

    def add_gaps_in_template_sequence(self):
        """
            Sets the benchmark attribute of the template.
            The benchmark attribute represents the fold family type of the
            template: "Family", "Superfamily" or "Fold", that is to say how similar
            the template is from the query. This is necessary to be
            able to benchmark the new scoring functions.

            Args:
                fold_type (str): "Family", "Superfamily" or "Fold"

            Returns:
                void
        """
        self.benchmark = fold_type
