"""
    .. module:: threading
      :synopsis: This module implements all the functions to make the threading of
                 the query sequence on the templates.
"""

# Third-party modules
import numpy as np
import time
from functools import wraps


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print("Total time running %s: %s seconds" %
               (function.__name__, str(t1-t0)))
        return result
    return function_timer

#@fn_timer
def calc_dist_convert_energy(query, template, dist_range, dope):
    """
        Calculate the matrix of distances of pair residues of the query sequence,
        using the coordinates of the template sequence: threading of the query
        on the template sequence. The distance calculation is optimized.
        The calculated distance is then convert into a dope energy value.

        Args:
            query (list of Residue): List of residues of the query sequence
            template (list of Residue): List of residues of the template sequence
            dist_range (list of int): Range of distances in angströms. Distances
                                      within this range only are taken into account
            dope (dictionary): A dictionary with key = res_1-res_2 and value = an array of
                               30 dope energy values.

        Returns:
            numpy 2D array: 2D matrix containing energy between pairs of residues
                            of the query sequence threaded on the template.
    """

    query_size = len(query)
    energy = np.empty((query_size, query_size), dtype=object)
    energy.fill(np.nan)

    for i, row_res in enumerate(query):
        if row_res.name == "-" or template[i].name == "-":
            energy[i, (i+2):] = 10
            break
        for j in range(i+2, query_size):
            col_res = query[j]
            if col_res.name == "-" or template[j].name == "-":
                energy[:i, j] = 10
                break
            else:
                # THE most efficient method to calculate Euclidian distances.
                # Formula: distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
                a_min_b = template[j].ca_coords - template[i].ca_coords
                dist = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
                # Keep distances only in a defined range because we don't want to
                # take into account directly bonded residues (dist < ~5 A) and too far residues
                if dist_range[0] <= dist <= dist_range[1]:
                    interval_index = round(int((dist * 30) / 15))
                    energy[i, j] = dope[row_res.name+col_res.name][interval_index]
    return energy
