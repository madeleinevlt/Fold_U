"""
.. module:: Alignment
   :synopsis: This module implements the Alignment class.
"""

# Third-party modules
import numpy as np
from Bio.SubsMat import MatrixInfo
from Bio.SeqUtils import seq3
from Bio.PDB import PDBParser
from Bio.PDB import NACCESS as nac
import modeller as m
import modeller.automodel as am
from modeller.automodel import assess
import contextlib
import os
import re
import pathlib
import logging


logging.basicConfig(filename="run_warnings.log", level=logging.WARNING)


def process(dist_range, gap_penalty, dope_dict, output_path, naccess_bin_path, ali):
    """
        Generates the threading and the blosum scores for a given Alignment object.

        Returns:
            tuple: (Sum of the different scores, Number of the Alignment,
            Template's name, Template's benchmark)

    """
    #print(ali.template.name, ali.template.pdb)
    # Calculate the threading score of all alignments and find the initial templates
    threading_score = ali.calculate_threading_score(dist_range, gap_penalty, dope_dict)
    # Calculate the modeller score of all alignments
    modeller_score = ali.calculate_modeller_score(output_path)
    # Calculate secondary structure score
    ss_score = ali.calculate_ss_score()
    # Calculate the BLOSUM score of all alignments
    #blosum_score = ali.calculate_blosum_score()
    blosum_score = 0
    access_score = ali.calculate_access_score(naccess_bin_path, 20)
    return ali.num, ali.score, threading_score, modeller_score, ss_score, access_score, blosum_score,\
        ali.template.name, ali.template.benchmark


def clean_modeller_outputs():
    """
    Removes files matching the pattern in the directory given in argument.

    Args:
        dir (str): Path to directory to clean
        pattern (str): Regex expression for files to remove
    """
    modeller_dir = "."
    pattern = "^[^\\.].*\\.(?!pdb).*$"
    for file in os.listdir(modeller_dir):
        if re.search(pattern, file):
            os.remove(os.path.join(modeller_dir, file))


class Alignment:
    """
    .. class:: Alignment

      This class groups informations about an alignment.

    Attributes:
        num (int): Number of the alignment
        score (float): Score of the alignment
        query (Query object): Instance of a Query object as
                              ``Query(query_residues, query_first, query_last)``
        template (Template object): Instance of a Template object as
                                    ``Template(template_name, template_residues)``
    """

    def __init__(self, num, score, query, template):
        self.num = num
        self.score = score
        self.query = query
        self.template = template

    def calculate_threading_score(self, dist_range, gap_penalty, dope_dict):
        """
            Calculate the threading score of the query on the template sequence.
            1) For each pair residues of the query sequence the distance between them is
            calculated using the coordinates of the template sequence.The distance calculation
            is optimized.
            2) The calculated distance is then convert into a dope energy value.
            3) All elements of the energy matrix generated are finally sum.

            Args:
                dist_range (list of int): Range of distances in angströms. Distances
                                          within this range only are taken into account
                gap_penalty (int): Gap penalty
                dope_dict (dictionary): A dictionary with key = res_1-res_2 and
                                        value = an array of 30 dope energy values.

            Returns:
                float: The threading score calculated. Sum of the energy matrix.
        """
        query = self.query.residues
        template = self.template.residues

        query_size = self.query.get_size()
        # This sets a numpy matrix of shape query * query which will contain all
        # the energies corresponding to distances between all pairs of residues
        # between the query and itself: the coordinates are the ones from the template.
        energy = np.empty((query_size, query_size), dtype=object)
        # Filling the matrix afterwards with "NaN" is faster
        energy.fill(np.nan)

        for i, row_res in enumerate(query):
            # The gap was already treated
            if i <= (query_size - 3) and energy[i, (i + 2)] == gap_penalty:
                continue
            # There is a gap in the query or the template
            elif row_res.name == "-" or template[i].name == "-":
                # The whole line is set with gap penalty value
                energy[i, (i + 2):] = gap_penalty
                continue
            for j in range(i + 2, query_size - 1):
                col_res = query[j]
                # The gap was already treated
                if energy[i, j] == gap_penalty:
                    continue
                # There is a gap in the query or the template
                elif col_res.name == "-" or template[j].name == "-":
                    # The whole column is set with gap penalty value
                    energy[:(j - 1), j] = gap_penalty
                    continue
                else:
                    # Calculate to distance between two residues
                    dist = template[i].calculate_distance(template[j])
                    # Keep distances only in a defined range because we don't want to
                    # take into account directly bonded residues (dist < ~5 A) and too far residues
                    if dist_range[0] <= dist <= dist_range[1]:
                        # DOPE energy values spread between 0.25 and 15 by 0.5 intervals
                        # So 30 intervals and max value = 15
                        interval_index = round(int((dist * 30) / 15))
                        energy[i, j] = dope_dict[row_res.name + col_res.name][interval_index]
        # Return the sum of energy matrix with numpy's "Nan" interpreted as zeros
        return -np.nansum(energy)

    def calculate_blosum_score(self):
        """
            Calculate a blosum score using the blosum62 matrix. This matrix contains
            substitution scores for each amino acid pair. A positive score is given to the more
            likely substitutions while a negative score is given to the less likely substitutions.

            Returns:
                int: The blosum score calculated.
        """
        # dictionary containing substitution scores of the blosum62 matrix
        blosum62 = MatrixInfo.blosum62
        score = 0
        for ind, res_q in enumerate(self.query.residues):
            res_t = self.template.residues[ind]
            if res_q.name == "-" or res_t.name == "-":
                continue
            # The blosum62 matrix is symetric, thus we need to test
            # if (res_1, res_2) key is in the dictionary.
            if (res_q.name, res_t.name) in blosum62:
                score += blosum62[(res_q.name, res_t.name)]
            # If not, the corresponding key is (res_2, res_1)
            else:
                score += blosum62[(res_t.name, res_q.name)]
        return score

    def calculate_ss_score(self):
        """
        This function calculates a score based on the secondary structure prediction of the query.
        We consider only PSIPRED predictions with a confidence >= 7 on the scale 0-9.
        The reason is that we want secondary structure predictions that have at least a Q3 accuracy
        of 80% (cf. doi:[10.1186/1471-2164-11-S4-S4])
        We use the average three-state prediction accuracy (Q3) to measure the accuracy of the
        secondary structure prediction of PSIPRED.
        The formula is: score = (N - total_incorrect) / N

        With N = length of the gapless query and total_incorrect = all the incorrectly predicted
        secondary structures with a confidence score < 7 (True negatives)

        Returns:
            float: Q3, the secondary structure score, as the proportion of well predicted secondary
            structure predictions
        """
        total_incorrect = 0
        score = 0
        query_ind = 0
        templ_ind = 0
        ind = 0
        query_len = self.query.get_size()
        while ind < query_len:
            # Skip gaps in secondary structure predictions
            while query_ind < query_len and self.query.residues[query_ind].secondary_struct == "-":
                query_ind += 1
            while (templ_ind < query_len
                    and self.template.residues[templ_ind].secondary_struct == "-"):
                templ_ind += 1
            # Count incorrect secondary structure predictions of query according to the template
            if (query_ind < query_len and templ_ind < query_len
                and self.query.residues[query_ind].ss_confidence < 7
                    and self.query.residues[query_ind].secondary_struct
                    != self.template.residues[templ_ind].secondary_struct):
                total_incorrect += 1
            ind += 1
            query_ind += 1
            templ_ind += 1
        try:
            # Calculate Q3
            gapless_query_len = len([res for res in self.query.residues if res.name != "-"])
            score = (gapless_query_len - total_incorrect) / gapless_query_len
        except ZeroDivisionError as err:
            print(str(err), "\n\nError ss_score: the query seems to be of size null")
        return score

    def write_pdb(self, pdb_path):
        """
            Write a pdb file by threading the query sequence on the template CA coordinates.

            Args:
                pdb_path (str): Path of the pdb file to create.
        """
        with open(pdb_path, "w") as file:
            # Extra informations on the template used to generate the pdb file
            file.write("REMARK Threading of query sequence on the {:s} template #{:d}.\n"
                       .format(self.template.name, self.num))
            ind = 0
            count_atom = 1
            for count_res in range(self.query.first, self.query.last+1):
                res_t = self.template.residues[ind]
                res_q = self.query.residues[ind]
                if res_q.name == "-" or res_t.name == "-":
                    ind += 1
                    continue
                # # N "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"
                           .format("ATOM", count_atom, "N", seq3(res_q.name).upper(), "A",
                                   count_res,
                                   res_t.n_atom.coords[0],
                                   res_t.n_atom.coords[1],
                                   res_t.n_atom.coords[2],
                                   1.00, 0, "N"))
                count_atom += 1
                # CA "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"
                           .format("ATOM", count_atom, "CA", seq3(res_q.name).upper(), "A",
                                   count_res,
                                   res_t.ca_atom.coords[0],
                                   res_t.ca_atom.coords[1],
                                   res_t.ca_atom.coords[2],
                                   1.00, 0, "C"))
                count_atom += 1
                # C "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"
                           .format("ATOM", count_atom, "C", seq3(res_q.name).upper(), "A",
                                   count_res,
                                   res_t.c_atom.coords[0],
                                   res_t.c_atom.coords[1],
                                   res_t.c_atom.coords[2],
                                   1.00, 0, "C"))
                count_atom += 1
                ind += 1
            # The two last lines of the created pdb file ("END" and "TER" lines)
            file.write("END\n")

    def calculate_access_score(self, naccess_bin_path, threshold):
        '''
            Calculate the accessibility score between the predicted model
            and the template pdb structure.

            Args:
                naccess_bin (String): Path to the naccess binary
                threshold cutoff (int): Cutoff for relative accessibility residue values

            Returns:
                float: The Accessibility score calculated
        '''
        # Path to the predicted model PDB
        predicted_model_pdb = "data/pdb/"+self.template.name+"/"+self.template.pdb+".atm"
        initial_template = "data/pdb/"+self.template.name+"/"+self.template.pdb.split("_")[0]+".atm"
        #generate PDB_model
        structure_predicted_model = PDBParser(QUIET=True).get_structure(
            "predicted_model", predicted_model_pdb)[0]
        structure_template_model = PDBParser(QUIET=True).get_structure(
            "template_model", initial_template)[0]
        #generate accessibilities files
        rsa_data_predicted_model, asa_data_predicted_model = nac.run_naccess(
            structure_predicted_model, predicted_model_pdb, naccess= naccess_bin_path)
        rsa_data_template_model, asa_data_template_model = nac.run_naccess(
            structure_template_model, initial_template, naccess= naccess_bin_path)
        #Parse the naccess output .rsa file to retrieve the
        #relative  % of solvant accessible area for each CA
        predicted_model_rsa = nac.process_rsa_data(rsa_data_predicted_model)
        template_model_rsa = nac.process_rsa_data(rsa_data_template_model)
        #Keep only residues with a relative accessibilities  threshold
        predicted_model_accessible_residues = self.keep_accessible_residues(predicted_model_rsa, threshold)
        template_model_accessible_residues = self.keep_accessible_residues(template_model_rsa, threshold)
        #Get the common accesible residues
        common_residues_number = len(set(predicted_model_accessible_residues).intersection(template_model_accessible_residues))
        #normalization
        return(common_residues_number/len(predicted_model_accessible_residues))

    def keep_accessible_residues(self, naccess_rsa, threshold):
        """
            From the output of naccess we keep only accessible residues which have a
            all_atoms_rel value > 30 (arbitrary threshold)

            Args:
                naccess_rsa (dict): A dictionnary containing the output of naccess's calculations
                threshold (int): all_atom_rel value

            Returns:
                dict: Keys are the residue ids and as value their solvant accessible area
        """
        accessible_residues_dict = {}
        for (chain_id, res_id), data_dict in naccess_rsa.items():
            for key, val in data_dict.items():
                if key == "all_atoms_rel" and val >= threshold:
                    accessible_residues_dict[res_id[1]] = val
        return accessible_residues_dict

    def write_alignment_for_modeller(self, ali_path):
        """
            Modeller requires a specifically formatted alignment file ".ali",
            following the PIR type format.
            This function writes this .ali file for the current alignment between
            the query and the template.

            Example of an alignment.ali file for Modeller:
            >P1;5fd1
            structureX:5fd1:1:A:106:A:ferredoxin:Azotobacter vinelandii:1.90: 0.19
            AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
            EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*

            >P1;1fdx
            sequence:1fdx:1::54::ferredoxin:Peptococcus aerogenes:2.00:-1.00
            AYVINDSC--IACGACKPECPVNIIQGS--IYAIDADSCIDCGSCASVCPVGAPNPED-----------------
            -------------------------------*
        """
        with open(ali_path + self.template.name+".ali", "w") as ali_out:
            template_len = len([res for res in self.template.residues if res.name != '-'])
            ali_out.write(">P1;" + self.template.pdb)
            if self.template.first == 1:
                ali_out.write("\nstructure:" + self.template.pdb
                              + ":" + str(self.template.first) + ":@:"
                              + str(template_len)
                              + ":@::::\n")
            else:
                ali_out.write("\nstructure:" + self.template.pdb
                              + ":" + str(self.template.first) + ":@:"
                              + str(template_len + self.template.first - 1)
                              + ":@::::\n")
            ali_out.write(self.template.display() + "*")
            ali_out.write("\n>P1;query_" + self.template.name)
            ali_out.write("\nsequence:::::::::\n")
            ali_out.write(self.query.display(modeller=True) + "*")

    def calculate_modeller_score(self, res_path):
        """
            * This function constructs a single comparative model for the query
            sequence from the known template structure, using alignment.ali,
            a PIR format alignment of query and template. The final model is
            written into the PDB file.
            * This function also returns the DOPE assessed score of the model
            generated by MODELLER.
            DOPE is the most reliable score at separating native-like models
            from decoys (lower, i.e, more negative, DOPE scores tend to
            correlate with more native-like models).

            Args:
                res_path (str): Path to the results folder

            Returns:
                float: The DOPE score of the model generated by MODELLER.
        """
        root_dir = os.getcwd()
        modeller_out_dir = res_path + "/modeller/"
        ali_dir = "alignments/"
        pathlib.Path(modeller_out_dir + ali_dir).mkdir(parents=True, exist_ok=True)
        # MODELLER generates the result files in his current directory, so we must
        # go to the results directory and come back to root dir afterwards.
        os.chdir(modeller_out_dir)
        path_to_atm = root_dir + "/data/pdb/" + self.template.name
        # We reindex all the PDB files to avoid any problem with modeller
        self.template.reindex_pdb(1, path_to_atm, True)
        # Parse the new PDB to get new residues and their coordinates generated by MODELLER
        self.template.parse_pdb(path_to_atm + "/" + self.template.pdb + ".atm")
        # Insert the gaps of the new template at the same positions in the query
        #self.query.add_gaps_in_query(self.template)
        self.write_alignment_for_modeller("./alignments/")
        # request no verbose output (still gives few that will be silenced afterwards)
        #m.log.none()
        # Redirect Modeller's verbose into nothingness, nil, chaos and abysses.
        with contextlib.redirect_stdout(None):
            # create a new MODELLER environment to build this model in
            m.env = m.environ()

            # directories for input atom files
            m.env.io.atom_files_directory = [path_to_atm]
            a_model = am.automodel(m.env,
                                   # alignment filename
                                   alnfile=ali_dir + self.template.name + '.ali',
                                   # codes of the templates
                                   knowns=self.template.pdb,
                                   # code of the target
                                   sequence='query_' + self.template.name,
                                   # DOPEHR is very similar to DOPEHR but is obtained at
                                   # Higher Resolution (using a bin size of 0.125Å
                                   # rather than 0.5Å).
                                   assess_methods=assess.DOPEHR)
            a_model.very_fast()
            # index of the first and last model (determines how many models to calculate)
            a_model.starting_model = 1
            a_model.ending_model = 1
            modeller_dope_score = 0
            try:
                a_model.make()
            except m.ModellerError as err:
                logging.warning("Modeller error with " + self.template.name + " | "
                            + self.template.pdb + "\n", str(err))
                pass
        new_model_pdb = a_model.outputs[0]["name"]
        modeller_dope_score = a_model.outputs[0]["DOPE-HR score"]
        os.rename(new_model_pdb, path_to_atm + "/" + self.template.pdb + "_mod.atm")
        # Parse the new models (PDB) generated by MODELLER to get the residues and coordinates
        # of residues
        self.template.parse_pdb(path_to_atm + "/" + self.template.pdb + "_mod.atm")
        self.template.pdb = self.template.pdb + "_mod"
        # remove useless files generated by MODELLER
        clean_modeller_outputs()
        # Go back to root directory
        os.chdir(root_dir)
        return modeller_dope_score
