# Introduction

[![](../img/logo.png)](https://github.com/maxibor/pydamage)

[![](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions) [![](https://readthedocs.org/projects/pydamage/badge/?version=latest)](https://pydamage.readthedocs.io/en/latest/?badge=latest)

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

## CLI help

Command line interface help message

```bash
pydamage --help
```

See also [CLI](CLI)
