# ![pydamage logo](docs/img/logo.png)

[![pydamage CI](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions)

Pydamage is a Python software to automate the process of identifying assembled contigs with characteristic patterns of ancient DNA.
It uses a process akin to a likelihood ratio test to attempt to discriminate between truly ancient, and modern contigs likely originating from sample contamination.

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

## Help

```bash
$ pydamage --help
Usage: pydamage [OPTIONS] BAM

Options:
  --version              Show the version and exit.
  -w, --wlen INTEGER     Window length from beginning of read  [default: 20]
  -p, --process INTEGER  Number of processes  [default: 2]
  -m, --mini INTEGER     Minimum reads required to be aligned to a reference to estimate damage
                         [default: 2000]
  -c, --cov FLOAT        Minimum coverage to consider reference to estimate damage [default:
                         0.5]
  -s, --show_al          Show alignments representations
  -pl, --plot            Make the damage plots
  --verbose              Verbose mode
  -o, --outdir PATH      Output directory  [default: pydamage_results]
  --help                 Show this message and exit.
```

> **pydamage logic: `n_reads >=minimum reads OR coverage >= minimum coverage`**
