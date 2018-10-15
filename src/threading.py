"""
    .. module:: threading
      :synopsis: This module implements all the functions to make the threading of
                    the query sequence on the templates.
"""


import numpy as np





def calc_dist_matrix(query, template, dist_range, gap_penalty):
    """
        Calculate the matrix of distances between the residues of the query sequence.
        Using the coordinates of the template sequence: threading of the query
        on the template sequence. The distance calculation is optimized.

        Args:
            query: List of residues of the query sequence
            template: List of residues of the template sequence
            dist_range: Range of distances in angströms.
                        Distances within this range only are taken into account
            gap_penalty: Penalization of gaps (integer)

        Returns:
            2D numpy matrix: The distance matrix between pairs of residues of the
            query sequence after being threaded on the template sequence.
    """

    size = len(query)
    matrix = np.empty((size, size), np.float)
    # use the query only for indexes
    for i, _ in enumerate(query):
        for j, _ in enumerate(query):
            # Penalize gaps
            if template[i] == "X" or template[j] == "X":
                matrix[i, j] = gap_penalty
            # One of the most efficient method to calculate the distances
            # https://stackoverflow.com/a/47775357/6401758
            # distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
            dist = np.linalg.norm(template[i].CA_coords - template[j].CA_coords)
            # Keep distances only in a defined range because we don't want to
            # take into account directly bonded residues (dist < ~5) and too far residues
            if dist_range[0] <= dist <= dist_range[1]:
                matrix[i, j] = dist
    return matrix
