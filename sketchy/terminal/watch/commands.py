import click
import shutil

from pathlib import Path

from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--dir', '-d', help='Directory with live FASTQ files from basecalling.'
)
@click.option(
    '--taxid', '-t', help='Taxonomic identifiers (from available indices) '
                          'to pre-filter for, otherwise use all available.'
)
@click.option(
    '--index', help='Show a list of genotype indices available in Sketchy.'
)
def watch(dir, taxid, index):

    """ Watch live run and pre-filter reads by taxonomic identifier """

    if index:
        show_indices()



