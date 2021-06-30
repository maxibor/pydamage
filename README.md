<p align="center">
    <a href="https://github.com/maxibor/pydamage"><img src="https://github.com/maxibor/pydamage/raw/master/docs/img/logo.png" alt="PyDamage"></a>
</p>

<p align="center">
    <a href="https://github.com/maxibor/pydamage/releases"><img src="https://img.shields.io/github/v/release/maxibor/pydamage?include_prereleases&label=version"/></a>
    <a href="https://github.com/maxibor/pydamage/actions"><img src="https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg"/></a>
    <a href="https://readthedocs.org/projects/pydamage/badge/?version=latest"><img src="https://readthedocs.org/projects/pydamage/badge/?version=latest"/></a>
    <a href="https://pypi.org/project/pydamage/"><img src="https://img.shields.io/pypi/v/pydamage"/></a>
    <a href="https://anaconda.org/maxibor/pydamage"><img src="https://img.shields.io/conda/v/maxibor/pydamage"/></a>
</p>

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

```
@article{Borry2021_pydamage,
    author = {Borry, Maxime and Huebner, Alexander and Rohrlach, Adam B and Warinner, Christina G},
    doi = {10.1101/2021.03.24.436838},
    elocation-id = {2021.03.24.436838},
    eprint = {https://www.biorxiv.org/content/early/2021/03/24/2021.03.24.436838.full.pdf},
    journal = {bioRxiv},
    publisher = {Cold Spring Harbor Laboratory},
    title = {PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly},
    url = {https://www.biorxiv.org/content/early/2021/03/24/2021.03.24.436838},
    year = {2021}
}
```