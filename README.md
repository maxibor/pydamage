<p align="center">
    <a href="https://github.com/maxibor/pydamage"><img src="https://github.com/maxibor/pydamage/raw/master/docs/img/logo.png" alt="PyDamage"></a>
</p>

<p align="center">
    <a href="https://github.com/maxibor/pydamage/releases"><img src="https://img.shields.io/github/v/release/maxibor/pydamage?include_prereleases&label=version"/></a>
    <a href="https://github.com/maxibor/pydamage/actions"><img src="https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg"/></a>
    <a href="https://pydamage.readthedocs.io"><img src="https://readthedocs.org/projects/pydamage/badge/?version=latest"/></a>
    <a href="https://pypi.org/project/pydamage/"><img src="https://img.shields.io/badge/install%20with-pip-blue"/></a>
    <a href="https://anaconda.org/bioconda/pydamage"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat"/></a>
</p>

# PyDamage

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
After modelling the ancient DNA damage using the C to T transitions, Pydamage uses a likelihood ratio test to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Installation

### With [conda](https://docs.conda.io/en/latest/) (recommended)

```bash
conda install -c bioconda pydamage
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
pydamage --outdir result_directory analyze aligned.bam
```

> Note that if you specify `--outdir`, it has to be before the PyDamage subcommand, example: `pydamage --outdir test filter pydamage_results.csv`

## CLI help

Command line interface help message

```bash
pydamage --help
```

## Documentation

[pydamage.readthedocs.io](https://pydamage.readthedocs.io)

## Cite

PyDamage has been published in PeerJ: [10.7717/peerj.11845](https://doi.org/10.7717/peerj.11845)

```
@article{borry_pydamage_2021,
    author = {Borry, Maxime and Hübner, Alexander and Rohrlach, Adam B. and Warinner, Christina},
    doi = {10.7717/peerj.11845},
    issn = {2167-8359},
    journal = {PeerJ},
    language = {en},
    month = {July},
    note = {Publisher: PeerJ Inc.},
    pages = {e11845},
    shorttitle = {PyDamage},
    title = {PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly},
    url = {https://peerj.com/articles/11845},
    urldate = {2021-07-27},
    volume = {9},
    year = {2021}
}

```