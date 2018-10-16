
def dope(dope_file):
   """
       Extract 10 energy values for the 20*20 residus-CA pairs and store it in indexed  double matrix
       
       Args:
           dope_file: the file dope_file.par containing energy values for each amino acid pair
       Returns:
           mat : amino acids indexed matrix containing lists of energies

    """
   #set up aa liste for rownames & colnames
   aa3 = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
          'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

   #set up matrix object 20*20 
   mat = pd.DataFrame(index = aa3, columns = aa3, dtype = object)

   with open(dope_file,'r') as dope_file:
       for line in dope_file:
           if line[4:6] == 'CA' and line[11:13]== 'CA':
                   # get the line with C-alpha for both amino acids
                   res_1 = line[0:3]
                   res_2 =line[7:10]
                   energy_res1_res2 = list(map(float, line[14:-1].split(" ")))
                   mat[res_1][res_2] = energy_res1_res2

   return mat