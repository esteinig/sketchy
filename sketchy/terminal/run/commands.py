import os
import click

from sketchy.sketchy import SketchyWrapper
from pathlib import Path

TEMPLATES = ['kpneumoniae', 'saureus', 'mtuberculosis']


@click.command()
@click.option(
    '--fastq',
    '-f',
    default=None,
    type=Path,
    required=True,
    help='Path to input Fastq containing basecalled nanopore reads'
)
@click.option(
    '--sketch',
    '-s',
    required=True,
    type=Path,
    help='Path to reference sketch local files or species template'
)
@click.option(
    '--ranks',
    '-r',
    default=20,
    type=int,
    help='Output highest ranking sum of shared hashes [20]'
)
@click.option(
    '--outdir',
    '-o',
    default='sketchy_out',
    type=Path,
    help='Output directory for data and plots [sketchy_out]'
)
@click.option(
    '--prefix',
    '-p',
    default='sketchy',
    type=str,
    help='Prefix for output files [sketchy].'
)
@click.option(
    '--limit',
    '-l',
    default=None,
    type=int,
    help='Maximum number of reads to predict on [all]'
)
@click.option(
    '--palette',
    '-c',
    default='YlGnBu',
    type=str,
    help='Brewer color palette for plots [YlGnBu]'
)
@click.option(
    '-b',
    '--stable',
    type=int,
    default=1000,
    help='Stability parameter to compute stable breakpoints, in reads [1000]'
)
@click.option(
    '-t',
    '--threads',
    type=int,
    default=4,
    help='Threads for sketch queries in Mash [4]'
)
@click.option(
    '-h',
    '--home',
    type=Path,
    default=Path.home() / '.sketchy',
    help='Sketchy path to reference sketch home directory. '
         'Can be set via environmental variable: SKETCHY_PATH'
)
@click.option(
    '-q',
    '--quiet',
    is_flag=True,
)
def run(
    fastq,
    sketch,
    ranks,
    outdir,
    prefix,
    limit,
    palette,
    stable,
    threads,
    home,
    quiet
):

    """ Sketchy main access wrapper for Rust pipeline """

    # TODO: Needs checks for user input etc.

    try:
        sketchy_path = os.environ['SKETCHY_PATH']
    except KeyError:
        sketchy_path = home

    try:
        sk = str(sketch.name)
        temp = sk.split('_')[0]
        if sk in TEMPLATES:
            sketch_file = sketchy_path / Path(f'{temp}/{sk}_15_1000')
        elif temp in TEMPLATES:
            sketch_file = sketchy_path / Path(f'{temp}/{sk}')
        else:
            sketch_file = sketch
    except IndexError:
        sketch_file = sketch

    sw = SketchyWrapper(
        fastx=fastq,
        sketch=sketch_file,
        prefix=prefix,
        outdir=outdir,
        verbose=not quiet
    )

    sw.run(
        ranks=ranks,
        limit=limit,
        stable=stable,
        threads=threads,
        palette=palette
    )
