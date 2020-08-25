# Introduction

[![](../img/logo.png)](https://github.com/maxibor/pydamage)

[![](https://github.com/maxibor/pydamage/workflows/pydamage_ci/badge.svg)](https://github.com/maxibor/pydamage/actions) [![](https://readthedocs.org/projects/pydamage/badge/?version=latest)](https://pydamage.readthedocs.io/en/latest/?badge=latest)

Pydamage, is a Python software to automate the process of contig damage identification and estimation.
It uses a process akin to a likelihood ratio test to attempt to discriminate between truly ancient, and modern contigs originating from sample contamination.

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

## CLI help

Command line interface help message

```bash
pydamage --help
```

See also [CLI](CLI)

## Output

Pydamage generates both a tabular and a visual output.

The tabular output is a comma-separated file with the following information for all analysed reference genomes:

  * `reference`: name of the reference genome/ contig
  * `null_model_p0`: probability of the null model
  * `null_model_p0_stdev`: standard error of the probability of the null model
  * `damage_model_p`: probability of the damage model
  * `damage_model_p_stdev`: standard error of the probability of the damage model
  * `damage_model_pmin`: ?
  * `damage_model_pmin_stdev`: ?
  * `damage_model_pmax`: ?
  * `damage_model_pmax_stdev`: ?
  * `pvalue`: p-value calculated from the likelihood-ratio test-statistic using a $\chi^2$ distribution
  * `qvalue`: p-value corrected for multiple testing using Benjamini-Hochberg procedure
  * `RMSE`: residual mean standard error of the model fit of the damage model
  * `nb_reads_aligned`: number of aligned reads
  * `coverage`: average coverage along the reference genome

Finally, the remaining columns indicate the frequency of C to T and G to A transitions at every read position observed in the sequencing data.

The visual output are PNG files, one per reference contig. In it, we plot the frequency of observed C to T and G to A transitions at the 5' end of the sequencing data and overlay it with the fitted models for both the damage and the null model including the 95% confidence intervals. Furthermore, we provide a "residuals versus fits" plot in order to allow the user to evaluate the fit of the model. Finally, the plot contains informtion on the average coverage along the reference genome and the p-value calculated from the likelihood-ratio test-statistic using a $\chi^2$ distribution.
