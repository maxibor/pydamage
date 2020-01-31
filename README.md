[![](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions)

<img src="docs/img/logo.png" alt="pydamage logo" width="200"/>

Pydamage, is a Python software to automate the process of contig damage identification and estimation. 
It uses a process akin to a likelihood ratio test to attempt to discriminate between truly ancient, and modern contigs originating from sample contamination.

## Install

Pydamage is not yet on *pip* or *conda*

```bash
git clone git@github.com:maxibor/pydamage.git
cd pydamage
conda create -f environment.yml
conda activate pydamage
python setup.py install
```

## Help

```bash
$ pydamage --help
Usage: pydamage [OPTIONS] BAM

Options:
  --version              Show the version and exit.
  -w, --wlen INTEGER     Window length from beginning of read  [default: 30]
  -p, --process INTEGER  Number of processes  [default: 2]
  -m, --mini INTEGER     Minimum reads aligned to consider reference
                         [default: 2000]
  -s, --show_al          Show alignments representations
  --verbose              Verbose mode
  -o, --output PATH      Output file  [default: ./contigs.csv]
  --help                 Show this message and exit.
```

