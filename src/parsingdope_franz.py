#!/usr/bin/env python
#! -*- coding: utf-8 -*-
#Python 3

#written by  AKE Franz-Arnold
#aerod7710@gmail.com
#~ October 2018 ~ 00h22
#Parsing a dop.file


#Loading packages
import pandas as pd

def get_dopetab(Path_file):
    '''
    Function for getting from a dope.par file an dope table
    return a DataFrame [20*20] of residues
    in which values are a long string of energy_value per residues
    with sep =  "**"

    input = Path of dope.par file
    output = Dope_table
    '''
    #Reading data
    dope_data = pd.read_csv(Path_file, header = None, sep = " ")

    #Deleting rows that are not CA
    dope_data = dope_data.loc[dope_data[1] == "CA"]
    dope_data = dope_data.loc[dope_data[3] == "CA"]

    #drop some columns (of CAs)
    dope_data.drop([1,3], axis = 1, inplace = True)

    #reindexation
    dope_data.reset_index(drop = True, inplace = True)
    dope_data.columns = range(dope_data.shape[1])

    #convert num_values of dope_data into string type
    dope_data_nums = dope_data.iloc[:,2:].astype(str)

    #Creating & Filling dope_tab
    dope_tab = pd.DataFrame(index = dope_data.loc[:19,1],
                            columns = dope_data.loc[:19,1]).astype(str)
    #...transforming values of energy into list of string !
    list_test = dope_data_nums.iloc[:,:].values.tolist()
    list_test = list(map('**'.join, list_test)) #separator of values ~ "**"

    #filling dope_tab
    start = 0
    end = dope_tab.shape[0]
    for i in range(dope_tab.shape[0]):
        dope_tab.iloc[:dope_tab.shape[0],i] = list_test[start:end]
        start = end
        end += dope_tab.shape[0]

    #returning dope_table
    return dope_tab
