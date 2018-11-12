"""
.. module:: Template
   :synopsis: This module implements the Template class.
"""

# Local modules
from src.residue import Residue

# Third-party modules
import numpy as np
import os
import wget
from Bio import SeqIO
import modeller as m
import modeller.automodel as am


class Template:
    """
    .. class:: Template

      This class groups informations about a template sequence/structure.

    Attributes:
        name (str): Name of the template
        residues (list of Residue object): Template's sequence of residues as list of Residues
                                           objects
        benchmark (int): Fold family type of the template: 3: Family, 2: Superfamily,
                         1: Fold and 0: None
                         It tells how similar the template is from the query structure.
                         This is necessary to be able to benchmark the new scoring functions.
        pdb (str): PDB filename of the template
    """

    def __init__(self, name, residues):
        self.name = name            # ex: Agglutinin
        self.residues = residues
        self.benchmark = 0
        self.pdb = None

    def set_benchmark(self, fold_type):
        """
            Sets the benchmark attribute of the template.

            Args:
                fold_type (int): 3: Family, 2: Superfamily, 1: Fold and 0: None
        """
        self.benchmark = fold_type


    def set_pdb_name(self, metafold_dict):
        """
            Set the pdb file name of the template from the template's name.

            Args:
                metafold_dict: A dictionary with key = template name and value = pdb file
        """
        self.pdb = metafold_dict[self.name].split(".")[0]

    def parse_pdb(self, pdb_path):
        """
            Parse the pdb file and set the CA coordinates.
        """
        count_res = 0
        nb_atoms = 0
        flag = True
        with open(pdb_path + self.name + "/" + self.pdb + ".atm", 'r') as file:
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
                        if count_res == len(self.residues):
                            break
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
        if len(self.pdb) == 6:
            chain = self.pdb[-2].upper()
            name = self.pdb[:-2].upper()
        elif len(self.pdb) == 5:
            chain = self.pdb[-1].upper()
            name = self.pdb[:-1].upper()
        else:
            chain = "A"
            name = self.pdb.upper()
        url = "https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=fastachain&compression=NO&structureId=" + name + "&chainId=" + chain
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
                (str, str): Tuple: The full FASTA sequence and the new sequence
                with gaps instead of missing residues.
        """
        # Fetch FASTA file from RCSB PDB database
        fasta_file = self.get_fasta_file()
        # Parse the file to get sequence
        sequence = list(SeqIO.parse(fasta_file, "fasta"))[0].seq
        new_seq = ""
        self_ind = 0
        # Add gaps "-" to the initial sequence inplace of the missing residues
        for ind, res in enumerate(sequence):
            if self_ind < len(self.residues) and res == self.residues[self_ind].name:
                new_seq += res
                self_ind += 1
            else:
                new_seq += "-"
        # Assign the new residues
        new_residues = [Residue(name) for name in new_seq]
        self.residues = new_residues

        return str(sequence), new_seq

    def write_template_alignment(self):
        """
            This function creates a PIR formatted file for MODELLER, containing
            the alignment between the sequence of the template containing
            missing residues and the new sequence of the template in which
            holes were replaced by gaps "-".

            Returns:
                str: Path to the new (gapless) PDB file of the template.
        """
        fasta_seq, new_template_seq = self.add_gaps_in_template_sequence()
        with open("alignments/" + self.name+".ali", "w") as ali_out:
            ali_out.write(">P1;" + self.pdb)
            ali_out.write("\nstructure:" + self.pdb
                          + ":" + str(self.first)
                          + ":@:" + str(len(new_template_seq)) + ":@::::\n")
            ali_out.write(new_template_seq + "*")
            ali_out.write("\n>P1;"+self.pdb+"_fill")
            ali_out.write("\nsequence:::::::::\n")
            ali_out.write(fasta_seq + "*")

    def build_new_pdb_with_modeller(self, dir_to_atm, ali_dir):
        """
            This function launches MODELLER to rebuild the PDB structures
            containing missing residues. Some structures have unresolved parts
            which are consequently not integrated in the PDB file because no
            coordinates could be determined for the faulty residues.
            We can use MODELLER to "fill in" these missing residues by treating
            the original structure (without the missing residues) as a template,
            and building a comparative model using the full sequence.

            Args:
                dir_to_atm (str): Path to template's .atm file
                ali_dir (str): Path to the alignment files.

            Returns:
                str: Path to the new (gapless) PDB file of the template.
        """
        m.log.verbose()
        m.env = m.environ()

        # directories for input atom files
        m.env.io.atom_files_directory = [dir_to_atm]

        # class MyModel(automodel):
        #     def select_atoms(self):
        #         return selection(self.residue_range('133', '135'),
        #                          self.residue_range('217', '230'))

        a_model = am.automodel(m.env,
                               # alignment filename
                               alnfile=ali_dir + self.name + '.ali',
                               # codes of the templates
                               knowns=self.pdb,
                               # code of the target
                               sequence=self.pdb + "_fill")
        a_model.starting_model = 1
        a_model.ending_model = 1

        a_model.make()
        new_model_pdb = a_model.outputs[0]["name"]
        os.rename(new_model_pdb, dir_to_atm + "/" + self.pdb + ".atm")
