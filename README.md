# [![pydamage logo](https://github.com/maxibor/pydamage/raw/master/docs/img/logo.png)](https://github.com/maxibor/pydamage)

[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/maxibor/pydamage?include_prereleases&label=version)](https://github.com/maxibor/pydamage/releases) [![pydamage CI](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions) [![Documentation Status](https://readthedocs.org/projects/pydamage/badge/?version=latest)](https://pydamage.readthedocs.io/en/latest/?badge=latest) [![PyPI](https://img.shields.io/pypi/v/pydamage)](https://pypi.org/project/pydamage/) [![Conda](https://img.shields.io/conda/v/maxibor/pydamage)](https://anaconda.org/maxibor/pydamage)

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
After modelling the ancient DNA damage using the C to T transitions, Pydamage uses a likelihood ratio test to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Install

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
pip install git+ssh://git@github.com/maxibor/pydamage.git
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:maxibor/pydamage.git
cd pydamage
conda env create -f environment.yml
conda activate pydamage
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
  -w, --wlen INTEGER     Window length for damage modeling  [default: 35]
  -p, --process INTEGER  Number of processes  [default: 2]
  -m, --mini INTEGER     Minimum reads aligned to consider reference
                         [default: 1000]

  -c, --cov FLOAT        Minimum coverage to consider reference  [default:8]

  -s, --show_al          Show alignments representations
  -pl, --plot            Make the damage plots
  --verbose              Verbose mode
  -o, --outdir PATH      Output directory  [default: pydamage_results]
  --force                Force overwriting of results directory
  --help                 Show this message and exit.
```

> **pydamage logic: `n_reads >=minimum reads OR coverage >= minimum coverage`**
