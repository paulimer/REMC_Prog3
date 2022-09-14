# REMC Reimplementation - Short Project
This is the repo for the short project in Programming 3 & Project Management course of Paris-Cité University.

It aims to reimplement the replica exchange Monte Carlo algorithm for protein described in [this article](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-342) in 2D dimensions.

## Installation
### With Docker
-  First install Docker : <https://docs.docker.com/get-docker/>
-  Then run :

``` sh
docker pull paulimer/remc
docker run --rm -it paulimer/remc
```
-  You can exit the container with `exit` or Ctrl-D
### With conda
-  Install Miniconda : <https://docs.conda.io/en/latest/miniconda.html>
-  Clone the repository
``` sh
git clone git@github.com:paulimer/REMC_Prog3.git
```
-  Move to the root of the directory

``` sh
cd REMC_Prog3
```
-  Create the conda environment

``` sh
conda env create -f conda_env.yml
```
-  Activate the conda environment

``` sh
conda activate remc
```

-  Make the main file executable

``` sh
chmod u+x ./main.py
```

`

## Usage
-  you can run the REMC folding algorithm with insulin as a demo protein

``` sh
./main.py --file ./data/insulin.fasta
```
-  You can also customize the different parameters.

| type of argument   | argument              | meaning                                    | default |
|:-------------------|:----------------------|:-------------------------------------------|--------:|
| positional         | nb_replicas           | the number of replicas                     | 5       |
|                    | steps                 | the maximum number of exchange steps       | 20      |
|                    | local_steps           | the maximum number of monte carlo search   | 50      |
|                    | t_min                 | the minimum temperature of the replicas    | 160     |
|                    | t_max                 | the maximum temperature of the replicas    | 220     |
| optional arguments | -h, --help            | prints a help message                      |         |
|                    | --file [FILE]         | read the protein sequence from this file   |         |
|                    | --sequence [SEQUENCE] | reads the sequence from the given SEQUENCE |         |
|                    | --energy              | the energy to achieve                      |         |

One of `--sequence` or `--file` is required
-  You can also use you own sequences, with the `--sequence` optional argument or the `--file` one. You cannot use your own fasta file with the docker container
