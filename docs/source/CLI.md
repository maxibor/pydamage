# CLI

Command Line Interface

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