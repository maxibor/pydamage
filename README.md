# ![pydamage logo](docs/img/logo.png)

[![pydamage CI](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions) [![Documentation Status](https://readthedocs.org/projects/pydamage/badge/?version=latest)](https://pydamage.readthedocs.io/en/latest/?badge=latest)

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
It uses a process akin to a likelihood ratio test to attempt to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Install

Pydamage is not yet on *pypi* nor *conda*, but you can already install it using pip, provided that you have access to this repository.

### Install dependencies in conda environment

```bash
git clone git@github.com:maxibor/pydamage.git
cd pydamage
conda env create -f environment.yml
conda activate pydamage
```

### Install pydamage

- from source

```bash
python setup.py install
```

- from Github using pip

```bash
pip install git+ssh://git@github.com/maxibor/pydamage.git
```

## Documentation

[pydamage.readthedocs.io](https://pydamage.readthedocs.io)

## Help

```bash
$ pydamage --help
Usage: pydamage [OPTIONS] BAM

  PyDamage: Damage parameter estimation for ancient DNA
  Author: Maxime Borry
  Contact: <borry[at]shh.mpg.de>
  Homepage & Documentation: github.com/maxibor/pydamage

  BAM: path to BAM/SAM/CRAM alignment file

Options:
  --version              Show the version and exit.
  -w, --wlen INTEGER     Window length for damage modeling  [default: 20]
  -p, --process INTEGER  Number of processes  [default: 2]
  -m, --mini INTEGER     Minimum reads aligned to consider reference
                         [default: 2000]

  -c, --cov FLOAT        Minimum coverage to consider reference  [default:
                         0.5]

  -s, --show_al          Show alignments representations
  -pl, --plot            Make the damage plots
  --verbose              Verbose mode
  -o, --outdir PATH      Output directory  [default: pydamage_results]
  --force                Force overwriting of results directory
  --help                 Show this message and exit.
```

> **pydamage logic: `n_reads >=minimum reads OR coverage >= minimum coverage`**
