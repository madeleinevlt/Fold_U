"""
.. module:: Alignment
   :synopsis: This module implements the Alignment class.
"""

# Third-party modules
import sys
sys.path.append("../bin/CCMpred/scripts/")
from top_couplings import get_top_pairs

import numpy as np
from conkit.applications import CCMpredCommandline
from Bio.SubsMat import MatrixInfo
from Bio.SeqUtils import seq3


class Alignment:
    """
    .. class:: Alignment

      This class groups informations about an alignment.

    Attributes:
        num (int): Number of the alignment
        score: Score of the alignment
        query: Instance of a Query object as
               ``Query(query_residues, query_first, query_last)``
        template: Instance of a Template object as
                  ``Template(template_name, template_residues)``
    """

    def __init__(self, num, score, query, template):
        self.num = num
        self.score = score
        self.query = query
        self.template = template

    #@fn_timer
    def calculate_threading_score(self, dist_range, gap_penalty, dope_dict):
        """
            Calculate the threading score of the query on the template sequence.
            1) For each pair residues of the query sequence the distance between them is
            calculated using the coordinates of the template sequence.The distance calculation
            is optimized.
            2) The calculated distance is then convert into a dope energy value.
            3) All elements of the energy matrix generated are finally sum.

            Args:
                dist_range (list of int): Range of distances in angstroms. Distances
                                          within this range only are taken into account
                gap_penalty (int): Gap penalty
                dope_dict (dictionary): A dictionary with key = res_1-res_2 and
                                        value = an array of 30 dope energy values.

            Returns:
                float: The threading score calculated. Sum of the energy matrix.
        """
        query = self.query.residues
        template = self.template.residues

        query_size = len(query)
        # This sets a numpy matrix of shape query * query which will contain all
        # the energies corresponding to distances between all pairs of residues
        # between the query and itself: the coordinates are the ones from the template.
        energy = np.empty((query_size, query_size), dtype=object)
        # Filling the matrix afterwards with "NaN" is faster
        energy.fill(np.nan)

        for i, row_res in enumerate(query):
            # The gap was already treated
            if i <= (query_size - 3) and energy[i, (i+2)] == gap_penalty:
                continue
            # There is a gap in the query or the template
            elif row_res.name == "-" or template[i].name == "-":
                # The whole line is set with gap penalty value
                energy[i, (i+2):] = gap_penalty
                continue
            for j in range(i+2, query_size-1):
                col_res = query[j]
                # The gap was already treated
                if energy[i, j] == gap_penalty:
                    continue
                # There is a gap in the query or the template
                elif col_res.name == "-" or template[j].name == "-":
                    # The whole column is set with gap penalty value
                    energy[:(j-1), j] = gap_penalty
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
                        energy[i, j] = dope_dict[row_res.name+col_res.name][interval_index]
        # Return the sum of energy matrix with numpy's "Nan" interpreted as zeros
        return np.nansum(energy)

    def calculate_physics_score(self):
        """
            Calculate a physics-based score using the BLOSSUM62 matrix. This matrix contains
            substitution scores for each amino acid pair. A positive score is given to the more
            likely substitutions while a negative score is given to the less likely substitutions.

            Args:
                void

            Returns:
                int: The physics-based score calculated.
        """
        # dictionary containing substitution scores of the BLOSSUM62 matrix
        blosum62 = MatrixInfo.blosum62
        score = 0
        for ind, res_q in enumerate(self.query.residues):
            res_t = self.template.residues[ind]
            if res_q.name == "-" or res_t.name == "-":
                continue
            # The BLOSSUM62 matrix is symetric, thus we need to test
            # if (res_1, res_2) key is in the dictionary.
            if (res_q.name, res_t.name) in blosum62:
                score += blosum62[(res_q.name, res_t.name)]
            # If not, the corresponding key is (res_2, res_1)
            else:
                score += blosum62[(res_t.name, res_q.name)]
        return -score

    def write_pdb(self, pdb_path):
        """
            Write a pdb file by threading the query sequence on the template CA coordinates.

            Args:
                pdb_path (str): Path of the pdb file to create.

            Returns:
                void
        """
        with open(pdb_path, "w") as file:
            # Extra informations on the template used to generate the pdb file
            file.write("REMARK Threading of query sequence on the {:s} template #{:d}.\n"\
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
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
                    .format("ATOM", count_atom, "N", seq3(res_q.name).upper(), "A", count_res,\
                    res_t.n_atom.coords[0], res_t.n_atom.coords[1], res_t.n_atom.coords[2], 1.00, 0, "N"))
                count_atom += 1
                # CA "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
                    .format("ATOM", count_atom, "CA", seq3(res_q.name).upper(), "A", count_res,\
                    res_t.ca_atom.coords[0], res_t.ca_atom.coords[1], res_t.ca_atom.coords[2], 1.00, 0, "C"))
                count_atom += 1
                # C "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
                    .format("ATOM", count_atom, "C", seq3(res_q.name).upper(), "A", count_res,\
                    res_t.c_atom.coords[0], res_t.c_atom.coords[1], res_t.c_atom.coords[2], 1.00, 0, "C"))
                count_atom += 1
                ind += 1
            # The two last lines of the created pdb file ("END" and "TER" lines)
            file.write("END\n")


def calculate_co_evolution_score(aln_file,distance_matrix):
    """
        Calculates co-evolution score of the query based on the MSA alignment
        Co-evolution score measures co-occurence of a pair of amino acid in
        ortholog sequences. Two amino acid have co-evoluated if the occurence
        of one of this amino never occur whithout the other.

        Args:
            aln_file:clustal file in .aln
            distance_matrix: distance matrix of the query against itself
        Returns:
            score_co-evolution
    """
    #Parsing aln file in order to get query + gaps
    with open(aln_file,'r') as file:
        data = file.read()
        query = data.split('\n', 1)[0]
        #repport index of gaps
    i_gap =[i for i, e in enumerate(query) if e == "-"]
    #Predict contacts
    ccmpred_cline = CCMpredCommandline(
        cmd ='bin/CCMpred/bin/ccmpred', alnfile= aln_file, matfile= "contact.mat"
    )
    ccmpred_cline()
    #Extract 30 top coupling
    subprocess.call(
        "./bin/CCMpred/scripts/top_couplings.py contact.mat > top_output.mat",
        shell=True
    )
    #parsing top_coupling
    dic_top_score= {}
    with open("top_output.mat", "r") as file:
        index = 1
        file.readline()
        for line in file:
            if float(line[0:2]) not in i_gap and float(line[3:5]) not in i_gap:
                dic_top_score[index] = [float(line[0:2]), float(line[3:5]),
                float(line[6:19])]
                index += 1
    #re-indexing
    for i, val in dic_top_score.items():
        for pos_gap in range(len(i_gap)-1):
            if val[0] > i_gap[pos_gap] and val[0] < i_gap[pos_gap+1]:
                val[0] = val[0]+1
            if val[1] > i_gap[pos_gap] and val[1] < i_gap[pos_gap+1]:
                val[1] = val[1]+1
    #compare top 30/ distance_matrix
    for i, val in dic_top_score.items():
        if distance_matrix[val[0],val[1]] < 7:
            score_co-evo = score_co-evo +1
        else:
            score_co-evo = score_co-evo -1
    return score_co-evo
