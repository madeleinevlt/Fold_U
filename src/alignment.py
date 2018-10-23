"""
.. module:: Alignment
   :synopsis: This module implements the Alignment class.
"""

# Third-party modules
import numpy as np
from Bio.SeqUtils import seq3

# Local modules
from src.template import Template
from src.query import Query


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

    def __init__(self, num, score, query_residues, query_first, query_last, template_name, template_residues):
        self.num = num
        self.score = score
        self.query = Query(query_residues, query_first, query_last)
        self.template = Template(template_name, template_residues)

    #@fn_timer
    def calculate_energy(self, dist_range, gap_penalty, dope_dict):
        """
            Calculate the matrix of distances of pair residues of the query sequence,
            using the coordinates of the template sequence: threading of the query
            on the template sequence. The distance calculation is optimized.
            The calculated distance is then convert into a dope energy value.

            Args:
                dist_range (list of int): Range of distances in angströms. Distances
                                          within this range only are taken into account
                gap_penalty (int): Gap penalty
                dope_dict (dictionary): A dictionary with key = res_1-res_2 and
                                        value = an array of 30 dope energy values.

            Returns:
                numpy 2D array: 2D matrix containing energy between pairs of residues
                                of the query sequence threaded on the template.
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
                        interval_index = round(int((dist * 30) / 15))
                        energy[i, j] = dope_dict[row_res.name+col_res.name][interval_index]
        return energy

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
            res_num = -1
            #for res_num, res_t in enumerate(self.template.residues):
            for ind in range(self.query.first, self.query.last+1):
                res_num += 1
                res_t = self.template.residues[res_num]
                res_q = self.query.residues[res_num]
                if res_q.name == "-" or res_t.name == "-":
                    continue
                # An "ATOM" line of the created pdb file
                file.write("{:6s}{:5d} {:^4s} {:>3s}{:>2s}{:4d}{:>12.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>12s}\n"\
                    .format("ATOM", ind, "CA", seq3(res_q.name).upper(), "A", ind,\
                    res_t.ca_coords[0], res_t.ca_coords[1], res_t.ca_coords[2], 1.00, 0, "C"))
            # The two last lines of the created pdb file ("END" and "TER" lines)
            file.write("END\n")
