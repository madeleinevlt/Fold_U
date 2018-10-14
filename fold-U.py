#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# IMPORTS
import src.parsing as parse

if __name__ == "__main__":
    
    metafold_dict = parse.metafold("data/METAFOLD.list")
    nb_templates = 10
    alignment_list = parse.foldrec(nb_templates, metafold_dict)
    for i in alignment_list:
        print(i.template.name)
        print(i.template.pdbID)
        print('\n')