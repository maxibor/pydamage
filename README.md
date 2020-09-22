[![pydamage logo](https://github.com/maxibor/pydamage/raw/master/docs/img/logo.png)](https://github.com/maxibor/pydamage)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/maxibor/pydamage?include_prereleases&label=version)](https://github.com/maxibor/pydamage/releases) [![pydamage CI](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions) [![Documentation Status](https://readthedocs.org/projects/pydamage/badge/?version=latest)](https://pydamage.readthedocs.io/en/latest/?badge=latest) [![PyPI](https://img.shields.io/pypi/v/pydamage)](https://pypi.org/project/pydamage/) [![Conda](https://img.shields.io/conda/v/maxibor/pydamage)](https://anaconda.org/maxibor/pydamage)

# PyDamage

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
After modelling the ancient DNA damage using the C to T transitions, Pydamage uses a likelihood ratio test to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Installation

### With [conda](https://docs.conda.io/en/latest/) (recommended)

```bash
conda install -c conda-forge -c bioconda -c maxibor pydamage
```

### With pip

```bash
pip install pydamage
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/maxibor/pydamage.git@dev
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:maxibor/pydamage.git
cd pydamage
git checkout dev
conda env create -f environment.yml
conda activate pydamage
pip install -e .
```


## Quick start

```bash
pydamage aligned.bam reference.fa
```

## CLI help

Command line interface help message

```bash
pydamage --help
```

## Documentation

[pydamage.readthedocs.io](https://pydamage.readthedocs.io)