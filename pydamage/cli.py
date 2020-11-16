#!/usr/bin/env python3

import click
from pydamage.main import analyze, analyze_group
from pydamage import __version__


@click.command()
@click.version_option(__version__)
@click.argument("bam", type=click.Path(exists=True))
@click.option(
    "-w",
    "--wlen",
    default=35,
    type=int,
    show_default=True,
    help="Window length (in bp) for damage modeling",
)
@click.option(
    "-p",
    "--process",
    default=2,
    type=int,
    show_default=True,
    help="Number of processes for parallel computing",
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(writable=True, dir_okay=True),
    default="pydamage_results",
    show_default=True,
    help="Output directory",
)
@click.option("-s", "--show_al",
              is_flag=True,
              help="Display alignments representations")
@click.option("-pl", "--plot",
              is_flag=True,
              help="Write damage fitting plots to disk")
@click.option("-vv", "--verbose", is_flag=True, help="Verbose mode")
@click.option("-f", "--force",
              is_flag=True,
              help="Force overwriting of results directory")
@click.option("-g", "--group",
              is_flag=True,
              help="Use entire BAM file as single reference for analyis "
              "(ignore reference headers)")
def cli(no_args_is_help=True, **kwargs):
    """\b
    PyDamage: Damage parameter estimation for ancient DNA
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/pydamage

    BAM: path to BAM/SAM/CRAM alignment file. MD tags need to be set.
    """
    analyze(**kwargs)


if __name__ == "__main__":
    cli()
