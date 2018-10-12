#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def calc_dist_matrix(query, template, dist_range, gap_penalty):
    """
        Calculate the matrix of distances between the residues of the query sequence.
        Using the coordinates of the template sequence: threading of the query
        on the template sequence.
        The distance calculation is optimized.

        Args:
            query: List of residues of the query sequence
            template: List of residues of the template sequence
            dist_range: Range of distances in angströms. Distances within this
                        range only are taken into account
            gap_penalty: Penalization of gaps (integer)
        Returns:
            matrix: The distance matrix between pairs of residues of the query
                     sequence after being threaded on the template sequence.
    """

    size = len(query)
    matrix = np.empty((size, size), np.float)
    for i, res_row in enumerate(query):
        for j, res_col in enumerate(query):
            # Penalize gaps
            if template[i] == "X" or template[j] == "X":
                matrix[i, j] = gap_penalty
            # distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
            # The method here is a more efficient way of calculating the
            # distance then the numpy function np.linalg.norm(A-B)
            # https://stackoverflow.com/a/47775357/6401758
            a_min_b = template[i].coords - template[j].coords
            dist = np.sqrt(np.einsum('ij,ij->i', a_min_b, a_min_b))
            # Keep distances only in a defined range
            if dist_range[0] <= dist <= dist_range[1]:
                matrix[i, j] = dist
    return matrix






if __name__ == '__main__':

    dist_matrix = calc_dist_matrix(query, template, dist_range)
