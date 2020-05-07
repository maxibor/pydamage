#!/usr/bin/env python3

import click
from pydamage.main import analyze
from pydamage import __version__

@click.command()
@click.version_option(__version__)
@click.argument('bam', type=click.Path(exists=True))
@click.option('-w',
              '--wlen',
              default=20,
              type=int,
              show_default=True,
              help='Window length for damage modeling')
@click.option('-p',
              '--process',
              default=2,
              type=int,
              show_default=True,
              help='Number of processes')
@click.option('-m',
              '--mini',
              default=1000,
              type=int,
              show_default=True,
              help='Minimum reads aligned to consider reference')
@click.option('-c',
              '--cov',
              default=8,
              type=float,
              show_default=True,
              help='Minimum coverage to consider reference')
@click.option('-s',
              '--show_al',
              is_flag=True,
              help='Show alignments representations')
@click.option('-pl',
              '--plot',
              is_flag=True,
              help='Make the damage plots')
@click.option('--verbose', is_flag=True, help='Verbose mode')
@click.option('-o',
              '--outdir',
              type=click.Path(writable=True, dir_okay=True),
              default="pydamage_results",
              show_default=True,
              help="Output directory")
@click.option('--force', 
              is_flag=True, 
              help='Force overwriting of results directory')
    
def cli(no_args_is_help=True, **kwargs):
    """\b
    PyDamage: Damage parameter estimation for ancient DNA
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/pydamage

    BAM: path to BAM/SAM/CRAM alignment file
    """
    analyze(**kwargs)


if __name__ == "__main__":
    cli()
