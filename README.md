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

PyDamage is a fast Python tool to automate ancient DNA (aDNA) damage identification and estimation for assembled contigs in *de novo* metagenomic assembly.

After modeling ancient DNA damage using C to T / G to A transitions, PyDamage uses a likelihood ratio test and GLM accuracy prediction to discriminate between truly ancient and modern contigs originating from sample contamination.

---

## Installation

### With [`uv`](https://docs.astral.sh/uv/) (recommended)

```bash
uv pip install pydamage
# or run directly without manual install:
uv run --with pydamage pydamage --help
```

### With [bioconda](https://anaconda.org/bioconda/pydamage)

```bash
conda install -c bioconda pydamage
```

### With `pip`

```bash
pip install pydamage
```

---

## Development Setup

Clone the repository and synchronize the virtual environment using `uv`:

```bash
git clone https://github.com/maxibor/pydamage.git
cd pydamage

# Install all dependencies and development/profiling tools
uv sync --all-groups

# Run the test suite
uv run pytest
```

---

## Quick Start

```bash
# Run analysis on a BAM/SAM/CRAM file
pydamage --outdir result_directory analyze aligned.bam

# Run multi-threaded analysis with plotting
pydamage --outdir result_directory analyze -p 8 --plot aligned.bam

# Filter results by accuracy and q-value
pydamage filter result_directory/pydamage_results.csv

# Generate bin damage plot
pydamage binplot result_directory/pydamage_results.csv bin.fasta
```

---

## Documentation

Full documentation is available at [pydamage.readthedocs.io](https://pydamage.readthedocs.io).

---

## Cite

PyDamage is published in PeerJ: [10.7717/peerj.11845](https://doi.org/10.7717/peerj.11845)

```bibtex
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