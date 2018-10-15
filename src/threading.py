"""
    .. module:: threading
      :synopsis: This module implements all the functions to make the threading of
                    the query sequence on the templates.
"""


import numpy as np


def get_pdb(pdb):
    p = PDBParser(QUIET=True) # QUIET=T : Warnings issued are suppressed
    pdb = p.get_structure(pdb,pdb)

    resList = [] # List of Residue instances

    for chain in pdb.get_chains():
        # First residue of the current chain
        first = next(res.id[1] for res in chain.get_residues() if (res.id[1] < 99999))

        # Last residue of the current chain
        for res in chain.get_residues():
            if (res.get_id()[0] != ' '):
                break
            last = res.get_id()[1]

        for resNum in range(first,last+1):
            r = Residue(chain, resNum)
            resList.append(r)
    return resList


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
    for i, res_row in enumerate(query):
        for j, res_col in enumerate(query):
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






if __name__ == '__main__':

    query = "AGLPVIMCLKSNNHQKYLRYQSDNIQQYGLLQFSADKILDPLAQFEVEPSKTYDGLVHIKSRYTNKYLVRWSPNHYWITASANEPDENKSNWACTLFKPLYVEEGNMKKVRLLHVQLGHYTQNYTVGGSFVSYLFAESSQIDTGSKDVFHVID"
    template = get_pdb("../data/1jlxa1.atm")
    dist_range = [5, 10]
    gap_penalty = -3
