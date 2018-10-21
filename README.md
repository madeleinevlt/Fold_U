# Fold U: A Protein Structure Prediction Program
![Fold U version](https://img.shields.io/badge/Fold%20U-1.1-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
![Sphinx version](https://img.shields.io/badge/sphinx%20build-1.7.4-brightgreen.svg)

This program is the second step (downstream) of a protein structure prediction project. This step consists of threading a query sequence on different given models.


Our project is part of the Meet-U 2018-2019 competition.
Meet-U is a collaborative pedagogical and research initiative between several Universities of Paris area. The course is intended to Master students (2nd year) in Bioinformatics. For more details refer to [http://www.meet-u.org/](http://www.meet-u.org/).

## Installation

### Clone the repository
```shell
git clone https://github.com/meetU-MasterStudents/2018---2019-Equipe-4.git
cd 2018---2019-Equipe-4
```

### Requirements
Install the few required packages / modules:
```shell
pip install -r src/requirements.txt
```

## Run the program

### Run the test
```shell
./fold_u.py data/foldrec/Agglutinin.foldrec
```
### Get help
```
./fold_u.py -h

    Usage:
        fold_u.py FILE [--nb_templates NUM] [--nb_pdb NB_PDB] [--metafold METAFOLD] [--dope DOPE] [--output PATH]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              The foldrec file [default: 100]
        -p NB_PDB, --nb_pdb NB_PDB            Number of pdb to create
                                              [default: 10]
        -m METAFOLD, --metafold METAFOLD      Path to the metafold.list file
                                              [default: data/METAFOLD.list]
        -d DOPE, --dope DOPE                  Path to the dope.par file
                                              [default: data/dope.par]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: res/threading]

```

## Documentation

The documentation of the project was generated with Sphinx and Read The Docs theme.
You can access it by opening the `docs/build/html/index.html` file
with your favorite web browser.

```shell
firefox docs/build/html/index.html
```

## Authors
- [Gabriel Cretin](https://github.com/gabrielctn)
- [Hélène Kabbech](https://github.com/kabhel)
- [Tom Gutman](https://github.com/tomgutman)
- [Flora Mikaeloff](https://github.com/FloraMika)
- [Franz-Arnold Ake](https://github.com/franzx5) 

## License

This project is licensed under the MIT License.
