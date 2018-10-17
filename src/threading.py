"""
    .. module:: threading
      :synopsis: This module implements all the functions to make the threading of
                 the query sequence on the templates.
"""

# IMPORTS
import numpy as np

def display_matrix(matrix):
    """
        Display the distance matrix in the terminal.

        Args:
            matrix (numpy ndarray): Matrix of shape query_length * query_length
                                    containing distances and gaps represented as
                                    "*".

        Returns:
            void
    """
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    for i in range(0, rows):
        for j in range(0, cols):
            if (matrix[i][j] == "*"):
                print("{:^4s}".format(matrix[i][j]), end="")
            else:
                print("{:4.1f}".format(matrix[i][j]), end="")
        print("\n")
    print("\n\n")


def calc_dist_matrix(query, template, dist_range):
    """
        Calculate the matrix of distances of pair residues of the query sequence,
        using the coordinates of the template sequence: threading of the query
        on the template sequence. The distance calculation is optimized.

        Args:
            query (list of Residue): List of residues of the query sequence
            template (list of Residue): List of residues of the template sequence
            dist_range (list of int): Range of distances in angströms.
                        Distances within this range only are taken into account

        Returns:
            2D numpy matrix: The distance matrix between pairs of residues of the
            query sequence after being threaded on the template sequence.
    """

    query_size = len(query)
    # This matrix holds distances and gaps ("*") for all pairs of residues of
    # the query sequence after being threaded on the template sequence
    matrix = np.full((query_size, query_size), fill_value=np.inf, dtype=object)
    # The distance matrix is symmetric so we make the calculations only for
    # the upper right triangular matrix. This saves have computation time.
    for i in range(len(query)):
        # We do not calculate distances between bonded residues and between
        # themselves so the 2nd loop starts at i+2
        for j in range(i+2, len(query)):
            row_res = query[i]
            col_res = query[j]
            # A gap represented by "-" = no distance calculation.
            # The whole line / row in the matrix will necessarily be "*"
            # This saves computation time
            if row_res.name == "-" or template[i].name == "-":
                matrix[i, j:] = "*"
                break
            elif col_res.name == "-" or template[j].name == "-":
                matrix[:i, j] = "*"
                break
            # One of the most efficient method to calculate the distances.
            # https://stackoverflow.com/a/47775357/6401758
            # distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
            #print(template[i].name,template[i].ca_coords,template[j].ca_coords,template[j].name)
            dist = np.linalg.norm(template[i].ca_coords - template[j].ca_coords)
            # Keep distances only in a defined range because we don't want to
            # take into account directly bonded residues (dist < ~5 A) and too far residues
            if dist_range[0] <= dist <= dist_range[1]:
                matrix[i, j] = dist
    return matrix
