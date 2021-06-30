#!/usr/bin/env python3

import click
from pydamage.main import pydamage_analyze
from pydamage.citation import get_citation
from pydamage.filter import apply_filter
from pydamage import __version__
from collections import OrderedDict


class NaturalOrderGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(NaturalOrderGroup, self).__init__(name=name, commands=None, **attrs)
        if commands is None:
            commands = OrderedDict()
        elif not isinstance(commands, OrderedDict):
            commands = OrderedDict(commands)
        self.commands = commands

    def list_commands(self, ctx):
        return self.commands.keys()


@click.group(cls=NaturalOrderGroup)
@click.version_option(__version__)
@click.pass_context
@click.option(
    "-o",
    "--outdir",
    type=click.Path(writable=True, dir_okay=True),
    default="pydamage_results",
    show_default=True,
    help="Output directory",
)
def cli(ctx, outdir):
    """\b
    PyDamage: Damage parameter estimation for ancient DNA
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/pydamage

    """

    ctx.ensure_object(dict)

    ctx.obj["outdir"] = outdir
    pass


@cli.command()
@click.pass_context
@click.argument("bam", type=click.Path(exists=True))
@click.option(
    "-w",
    "--wlen",
    default=20,
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
    "-s", "--show_al", is_flag=True, help="Display alignments representations"
)
@click.option("-pl", "--plot", is_flag=True, help="Write damage fitting plots to disk")
@click.option("-vv", "--verbose", is_flag=True, help="Verbose mode")
@click.option(
    "-f", "--force", is_flag=True, help="Force overwriting of results directory"
)
@click.option(
    "-g",
    "--group",
    is_flag=True,
    help="Use entire BAM file as single reference for analyis "
    "(ignore reference headers)",
)
def analyze(ctx, no_args_is_help=True, **kwargs):
    """\b
    Run the PyDamage analysis

    BAM: path to BAM/SAM/CRAM alignment file. MD tags need to be set.
    """
    pydamage_analyze(**kwargs, **ctx.obj)


@cli.command()
@click.pass_context
@click.argument("csv", type=click.Path(exists=True))
@click.option(
    "-t",
    "--threshold",
    default=0.5,
    type=float,
    show_default=True,
    help="Predicted accuracy filtering threshold. Set to 0 for finding threshold with kneed method",
)
def filter(ctx, no_args_is_help=True, **kwargs):
    """\b
    Filter PyDamage results on predicted accuracy and qvalue thresholds.

    CSV: path to PyDamage result file
    """

    apply_filter(**kwargs, **ctx.obj)


@cli.command()
def cite():
    """Get pydamage citation in bibtex format"""
    get_citation()


if __name__ == "__main__":
    cli()
