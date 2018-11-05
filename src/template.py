"""
.. module:: Template
   :synopsis: This module implements the Template class.
"""

# Third-party modules
import numpy as np
import os
import wget
from Bio import SeqIO
from modeller import *


class Template:
    """
    .. class:: Template

      This class groups informations about a template sequence/structure.

    Attributes:
        name: Name of the template
        residues: Template's sequence of residues as list of Residues objects
        pdb: PDB filename of the template
    """

    def __init__(self, name, residues):
        self.name = name            # ex: Agglutinin
        self.residues = residues
        self.pdb = None             # ex: 1jlxa1
        self.first = None

    def display(self):
        return "".join(str(res.name) for res in self.residues)

    def set_pdb_name(self, metafold_dict):
        """
            Get the PDB name of the current template from the template's name.

            Args:
                metafold_dict: A dictionary with key = template name and value = pdb file

            Returns:
                void
        """
        self.pdb = metafold_dict[self.name].split(".")[0]

    def parse_pdb(self):
        """
            Parse the pdb file and set the CA coordinates.

            Args:
                void

            Returns:
                void
        """
        count_res = 0
        nb_atoms = 0
        flag = True
        with open("data/pdb/" + self.name + "/" + self.pdb + ".atm", 'r') as file:
            for line in file:
                line_type = line[0:6].strip()
                name_at = line[12:16].strip()
                if line_type == "ATOM" and (name_at == "N" or name_at == "CA" or name_at == "C"):
                    # In some PDBs the residues ids do not start at 1
                    # so we remember the first id number
                    if flag:
                        self.first = int(line[22:26].strip())
                        flag = False
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

        Args:
            void

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
            This function replaces the missing residues of the template's
            sequence

            Returns:
                str: Path to the new (gapless) PDB file of the template.
        """
        # Fetch FASTA file from RCSB PDB database
        fasta_file = self.get_fasta_file()
        # Parse the file to get sequence
        sequence = list(SeqIO.parse(fasta_file, "fasta"))[0].seq
        # Add gaps "-" to the initial sequence inplace of the missing residues
        new_seq = ""
        self_ind = 0
        for ind, res in enumerate(sequence):
            if self_ind < len(self.residues) and res == self.residues[self_ind].name:
                new_seq += res
                self_ind += 1
            else:
                new_seq += "-"
        return new_seq


def write_new_alignment(self, res_path):
    """
        This function launches MODELLER to rebuild the PDB structures
        containing missing residues. Some structures have unresolved parts
        which are consequently not integrated in the PDB file because no
        coordinates could be determined for the faulty residues.
        We can use MODELLER to "fill in" these missing residues by treating
        the original structure (without the missing residues) as a template,
        and building a comparative model using the full sequence.

        Returns:
            str: Path to the new (gapless) PDB file of the template.
    """
    new_template_seq = self.add_gaps_in_template_sequence()
    ali_out_dir = res_path + "modeller/alignments/"
    with open(ali_out_dir + self.template.name+".ali", "w") as ali_out:
        template_len = len([res for res in self.template.residues if res.name != '-'])
        ali_out.write(">P1;" + self.template.pdb)
        ali_out.write("\nstructure:" + self.template.pdb
                      + ":" + str(self.template.first) + ":@:" + str(template_len) + ":@::::\n")
        ali_out.write(self.template.display() + "*")
        ali_out.write("\n\n>P1;"+self.template.pdb+"_fill")
        ali_out.write("\nsequence:::::::::\n")
        ali_out.write(new_template_seq + "*")


def build_new_pdb_with_modeller(self, res_path):
    """
        This function launches MODELLER to rebuild the PDB structures
        containing missing residues. Some structures have unresolved parts
        which are consequently not integrated in the PDB file because no
        coordinates could be determined for the faulty residues.
        We can use MODELLER to "fill in" these missing residues by treating
        the original structure (without the missing residues) as a template,
        and building a comparative model using the full sequence.

        Returns:
            str: Path to the new (gapless) PDB file of the template.
    """
