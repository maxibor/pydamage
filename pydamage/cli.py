#!/usr/bin/env python3

import click


@click.command()
@click.argument('bam', type=click.Path(exists=True))
def main(bam):
    print(bam)


if __name__ == "__main__":
    main()
