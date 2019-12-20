#!/usr/bin/env python3

import click
from .main import main


@click.command()
@click.argument('bam', type=click.Path(exists=True))
@click.option('-w',
              '--wlen',
              default=30,
              type=int,
              show_default=True,
              help='Window length from beginning of read')
@click.option('-p',
              '--process',
              default=2,
              type=int,
              show_default=True,
              help='Number of processes')
@click.option('-m',
              '--mini',
              default=2000,
              type=int,
              show_default=True,
              help='Minimum reads aligned to consider reference')
@click.option('-s',
              '--show_al',
              is_flag=True,
              help='Show alignments representations')
@click.option('-o',
              '--output',
              type=click.Path(writable=True),
              default="./contigs.csv",
              show_default=True,
              help="Output file")
def cli(no_args_is_help=True, **kwargs):
    main(**kwargs)


if __name__ == "__main__":
    cli()
