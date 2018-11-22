"""
.. module:: Alignment
   :synopsis: This module implements the Alignment class.
"""

# Third-party modules
import sys
from CCMpred.scripts.top_couplings import get_top_pairs
import numpy as np
from Bio.SubsMat import MatrixInfo
from Bio.SeqUtils import seq3
import pandas as pd


def process(dist_range, gap_penality, dope_dict, top_couplings_dict, query_index, ali):
    """
        Generates the threading and the blosum scores for a given Alignment object.

        Returns:
            tuple: (Sum of the different scores, Number of the Alignment,
            Template's name, Template's benchmark)

    """
    # Calculate the threading score of all alignments
    print(ali.template.name)
    threading_score = ali.calculate_threading_score(dist_range, gap_penality, dope_dict)
    blosum_score = ali.calculate_blosum_score()
    distance_matrix = ali.calculate_distance(query_index)
    print(distance_matrix)
    # for i in distance_matrix.shape[0]:
    #     for j in distance_matrix.shape[1]:
    #         print(distance_matrix[i, j], end=" ")
    ccmpred_score = ali.calculate_contact_score(top_couplings_dict, distance_matrix)
    return ali.num, ali.score, threading_score, blosum_score, ccmpred_score,\
           ali.template.name, ali.template.benchmark


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
        self.aln = None

    def calculate_distance(self, query_index):
        """
        Calculate distance
        """
        distance = np.empty((query_index[1], query_index[1]), dtype=object)
        k = 0
        while k < self.query.first-1:
            distance[k, (k+1):] = "x"
            k += 1
        print("###############################")
        query_size = self.query.get_size()
        query_len = len([res for res in self.query.residues if res != "-"])
        template_len = len([res for res in self.template.residues if res != "-"])
        print(query_len, template_len)
        print("\n\n")

        for i, row_res in enumerate(query_residues):
            # There is a gap in the query or the template
            if row_res.name == "-" or template_residues[i].name == "-":
                # The whole line is set with gap penalty value
                distance[i+k, (i+k+1):] = "x"
                continue
            for j in range(i + 1, query_size - 1):
                col_res = query_residues[j]
                # There is a gap in the query or the template
                if col_res.name == "-" or template_residues[j].name == "-":
                    # The whole column is set with gap penalty value
                    distance[:(j+k-1), j+k] = "x"
                    continue
                else:
                    # Calculate to distance between two residues
                    # print(i+k, j+k)
                    distance[i+k, j+k] = template_residues[i].calculate_distance(template_residues[j])
                    # Keep distances only in a defined range because we don't want to
                    # take into account directly bonded residues (dist < ~5 A) and too far residues
        k = i
        while query_index[1] > k > self.query.last-1:
            distance[k, (k+1):] = "x"
            k += 1
        return distance


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
                           .format("ATOM", count_atom, "N", seq3(res_q.name).upper(), "A",\
                                   count_res,
                                   res_t.n_atom.coords[0],\
                                   res_t.n_atom.coords[1],\
                                   res_t.n_atom.coords[2],\
                                   1.00, 0, "N"))
                count_atom += 1
                # CA "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"
                           .format("ATOM", count_atom, "CA", seq3(res_q.name).upper(), "A",\
                                   count_res,
                                   res_t.ca_atom.coords[0],\
                                   res_t.ca_atom.coords[1],\
                                   res_t.ca_atom.coords[2],\
                                   1.00, 0, "C"))
                count_atom += 1
                # C "ATOM" line
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"
                           .format("ATOM", count_atom, "C", seq3(res_q.name).upper(), "A",\
                                   count_res,
                                   res_t.c_atom.coords[0],\
                                   res_t.c_atom.coords[1],\
                                   res_t.c_atom.coords[2],\
                                   1.00, 0, "C"))
                count_atom += 1
                ind += 1
            # The two last lines of the created pdb file ("END" and "TER" lines)
            file.write("END\n")


    def calculate_contact_score(self, top_couplings_dict, distance_matrix):
        """
            Compare top 30 contacts in the query with corresponding calculated
            distances

            Args:
                top_couplings_dict: top ranking couplings indexes in the query
                distance_matrix: distance matrix of the query against itself
            Returns:
                contact_score
        """
        TP = 0
        for top_position in top_couplings_dict.values():
            #as matrix is triange, get matrix [i,j]
            if distance_matrix[top_position[0], top_position[1]] != None:
                if distance_matrix[top_position[0], top_position[1]] < 8:
                    TP += 1
            elif distance_matrix[top_position[1], top_position[0]] != None:
                if distance_matrix[top_position[1], top_position[0]] < 8:
                    TP += 1
            contact_score = TP/len(top_couplings_dict)
        return contact_score
