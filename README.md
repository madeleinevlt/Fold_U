![Fold-U release](https://img.shields.io/badge/fold--u-v1.2-blue.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
![Sphinx version](https://img.shields.io/badge/sphinx%20build-1.7.4-brightgreen.svg)  

<br>  

# Fold U: A Protein Structure Prediction Program

<p align="center">
  <img width="400" src="img/logo_foldu.png" alt="logo_foldu"/>
</p>

This program is the second step (downstream) of a protein structure prediction project. This step consists of threading a query sequence on different given templates.


Our project is part of the Meet-U 2018-2019 competition.
Meet-U is a collaborative pedagogical and research initiative between several Universities of Paris area. The course is intended to Master students (2nd year) in Bioinformatics. For more details, please refer to [http://www.meet-u.org/](http://www.meet-u.org/).

## Installation

### Clone the repository
```
git clone https://github.com/meetU-MasterStudents/Fold_U.git
cd Fold_U
```

### Requirements
Install the few required packages / modules:
```
pip install -r src/requirements.txt
```

## Run the program
fold_u takes in input N profil-profil alignments and their corresponding scores (foldrec file).

### Run the test
```
./fold_u.py data/foldrec/Agglutinin.foldrec
```
### Get help
```
./fold_u -h

    Usage:
        ./fold_u FOLDREC_FILE [--nb_templates NUM] [--nb_pdb NUM] [--output PATH]
                               [--metafold FILE] [--dope FILE] [--benchmark FILE]

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              the foldrec file [default: 100]
        -p NUM, --nb_pdb NUM                  Number of pdb to create
                                              [default: 10]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: res/threading]
        -m FILE, --metafold FILE              Path to the metafold.list file
                                              [default: data/metafold.list]
        -d FILE, --dope FILE                  Path to the dope.par file
                                              [default: data/dope.par]
        -b FILE, --benchmark FILE             Path to the benchmark.list file
                                              [default: data/benchmark.list]
```

## Documentation

The documentation of our program is generated with Sphinx and Read The Docs theme and is accessible at the following address:
[https://meetu-masterstudents.github.io/Fold_U/](https://meetu-masterstudents.github.io/Fold_U/)

## Authors

We are master students in bioinformatics at Paris Diderot University.
- [Gabriel Cretin](https://github.com/gabrielctn)
- [Hélène Kabbech](https://github.com/kabhel)
- [Tom Gutman](https://github.com/tomgutman)
- [Flora Mikaeloff](https://github.com/FloraMika)
- [Franz-Arnold Ake](https://github.com/franzx5)

## Acknowledgment

Thanks to [Maïté Cretin](https://www.linkedin.com/in/maitewho/) for the nice logo.

## License

This project is licensed under the MIT License.
