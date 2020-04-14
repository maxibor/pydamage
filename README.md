[![](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions)

<img src="docs/img/logo.png" alt="pydamage logo" width="200"/>

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

```
python setup.py install
```

- from Github using pip

```
pip install git+ssh://git@github.com/maxibor/pydamage.git
```

## Help

```bash
$ pydamage --help
Usage: pydamage [OPTIONS] BAM

Options:
  --version              Show the version and exit.
  -w, --wlen INTEGER     Window length from beginning of read  [default: 30]
  -p, --process INTEGER  Number of processes/CPUs to use  [default: 2]
  -m, --mini INTEGER     Minimum reads required to be aligned to a reference
                         to estimate damage [default: 2000]
  -s, --show_al          Show alignments representations
  --verbose              Verbose mode
  -o, --output PATH      Output file  [default: ./contigs.csv]
  --help                 Show this message and exit.
```
## FAQ

### Estimation fails but I can see a sort of damage-pattern!

The modelled curve is constructed based on the classic damage pattern of a smooth but consistent decrease in frequencies of C to T. Other curves, such as decreased frequencies on the first base due to lower effiency of adapter ligation e.g. [Seguin-Orlando et al 2013](https://doi.org/10.1371/journal.pone.0078575), may not work well.
