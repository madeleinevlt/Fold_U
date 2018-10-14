
# coding: utf-8

# In[ ]:


import numpy as np

def parsing_dope(dope_file):
    """
        Args:
            dope file.par
	    example 1 line
            "aa1 atom aa2 atom 10 values of energy for each amino acids pair at distance intervalles [0.25,15] with a step of 0.5"
            ...
        Returns:
            matrix of lists for each amino acids pair
	    1 list = 10 values of energy for each amino acids pair at distance intervalles [0.25,15] with a step of 0.5
     """
    #set up aa dictionnary	
    d = {'CYS': 0, 'ASP': 1, 'SER': 2, 'GLN': 3, 'LYS': 4,
        'ILE': 5, 'PRO': 6, 'THR': 7, 'PHE': 8, 'ASN': 9,
        'GLY': 10, 'HIS': 11, 'LEU': 12, 'ARG': 13, 'TRP': 14,
        'ALA': 15, 'VAL': 16, 'GLU': 17, 'TYR': 18, 'MET': 19}

    #set up matrix object
    mat = np.empty((20,20,30), dtype = object)

    with open(dope_file,'r') as dope_file:
        for line in dope_file:
            if line[4:6] == 'CA' and line[11:13]== 'CA':
                    # get the line with C-alpha for both amino acids
                    res_1 = line[0:3]
                    res_2 =line[7:10]
                    energy_res1_res2 = list(map(float, line[14:-1].split(" ")))
                    mat[d[res_1]][d[res_2]] = energy_res1_res2

