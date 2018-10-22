"""
    .. module:: threading
      :synopsis: This module implements all the functions to make the threading of
                 the query sequence on the templates.
"""

# Third-party modules
from functools import wraps
import time
import numpy as np



def fn_timer(function):
    """This function is a wrapper to benchmarck a function
    when trying to optimize it. It returns the total running time of the function

        Args:
            function: The function to benchmarck

        Returns:
            float: The time the function took to run

    """
    @wraps(function)
    def function_timer(*args, **kwargs):
        """It calculates the time during the run of a function"""
        t_0 = time.time()
        result = function(*args, **kwargs)
        t_1 = time.time()
        print("Total time running {:s}: {:s} seconds"\
            .format(function.__name__, str(t_1-t_0)))
        return result
    return function_timer



#@fn_timer
def calc_energy(ali, dist_range, gap_penality, dope):
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
            gap_penality (int): Gap penality
            dope (dictionary): A dictionary with key = res_1-res_2 and value = an array of
                               30 dope energy values.

        Returns:
            numpy 2D array: 2D matrix containing energy between pairs of residues
                            of the query sequence threaded on the template.
    """
    query = ali.query_residues
    template = ali.template.residues

    query_size = len(query)
    # This sets a numpy matrix of shape query * query which will contain all
    # the energies corresponding to distances between all pairs of residues
    # between the query and itself: the coordinates are the ones from the template.
    energy = np.empty((query_size, query_size), dtype=object)
    # Filling the matrix afterwards with "NaN" is faster
    energy.fill(np.nan)

    for i, row_res in enumerate(query):
        # The gap was already treated
        if i <= (query_size - 3) and energy[i, (i+2)] == gap_penality:
            continue
        # There is a gap in the query or the template
        elif row_res.name == "-" or template[i].name == "-":
            # The whole line is set with gap penalty value
            energy[i, (i+2):] = gap_penality
            continue
        for j in range(i+2, query_size-1):
            col_res = query[j]
            # The gap was already treated
            if energy[i, j] == gap_penality:
                continue
            # There is a gap in the query or the template
            elif col_res.name == "-" or template[j].name == "-":
                # The whole column is set with gap penalty value
                energy[:(j-1), j] = gap_penality
                continue
            else:
                # THE most efficient method to calculate Euclidian distances.
                # Formula: distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)
                #print(template[i].name, template[j].name)
                a_min_b = template[i].ca_coords - template[j].ca_coords
                dist = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
                # Keep distances only in a defined range because we don't want to
                # take into account directly bonded residues (dist < ~5 A) and too far residues
                if dist_range[0] <= dist <= dist_range[1]:
                    interval_index = round(int((dist * 30) / 15))
                    energy[i, j] = dope[row_res.name+col_res.name][interval_index]
    return energy
