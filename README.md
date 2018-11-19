[![Fold-U release](https://img.shields.io/badge/fold--u-v1.2-blue.svg)](https://github.com/meetU-MasterStudents/Fold_U/releases/tag/v1.2)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python version](https://img.shields.io/badge/python-3-brightgreen.svg)
[![Documentation Status](https://readthedocs.org/projects/fold-u/badge/?version=latest)](https://fold-u.readthedocs.io/en/latest/?badge=latest)

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
pip install -r requirements.txt
```

## Run the program
fold_u takes in input **N profil-profil alignments and their corresponding score** (foldrec file) and returns a `score.csv` file and the **top N pdb structures**.

### Toy example

The `scores.csv` and the **top 10 pdb structures** of the His_biosynth query sequence are stored in `results/His_biosynth` folder.
The alignment, threading and blosum scores are normalized using the **min-max scaling method** (values between 0 and 1). The last score represents the sum of these 3 scores. It has also been normalized.
```
./fold_u data/foldrec/His_biosynth.foldrec -o results/His_biosynth
```

### Generation of all the results and plot of the different scores
The `benchmark_rank.png` generated plot represents the cumulative sum of benchmarks encountered along the ranking (from rank 1 to rank 412). A cross means a family or superfamily type banchmark.
```
./script/benchmarking
```

<p align="center">
  <img width="450" src="img/His_biosynth_benchmark_rank.png" alt="benchmark_rank.png"/>
</p>

### Get help
```
./fold_u -h

    Usage:
        ./fold_u FOLDREC_FILE [--nb_templates NUM] [--nb_pdb NUM] [--output PATH]
                              [--metafold FILE] [--dope FILE] [--benchmark FILE] [--cpu NUM]

    Arguments:
        FOLDREC_FILE                          N profile * profile alignment and
                                              their corresponding score

    Options:
        -h, --help                            Show this
        -n NUM, --nb_templates NUM            First n templates to retrieve from
                                              the foldrec file [default: 413]
        -p NUM, --nb_pdb NUM                  Number of pdb to create
                                              [default: 10]
        -o PATH, --output PATH                Path to the directory containing
                                              the result files (scores and pdb)
                                              [default: ./results]
        -m FILE, --metafold FILE              Path to the metafold.list file
                                              [default: data/metafold.list]
        -d FILE, --dope FILE                  Path to the dope.par file
                                              [default: data/dope.par]
        -b FILE, --benchmark FILE             Path to the benchmark.list file
                                              [default: data/benchmark.list]
        -c NUM, --cpu NUM                     Number of cpus to use for parallelisation
                                              [default: 2]
```

## Documentation

The documentation of our program is generated with Sphinx and and built on [Read The Docs](https://fold-u.readthedocs.io/en/latest).

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
