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
    default=10,
    type=int,
    help='Output highest ranking sum of shared hashes [10]'
)
@click.option(
    '--outdir',
    '-o',
    default='sketchy',
    type=Path,
    help='Output directory for data and plots [sketchy]'
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
    default=1000,
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
    type=float,
    default=1000.,
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
    help="Run without logging output or progress bar."
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

    if limit == -1:
        limit = None

    if 0. < stable <= 1.:
        if limit is None:
            stable = 1000
        else:
            stable = int(limit*stable)
    else:
        stable = int(stable)

    try:
        sketch_name = str(sketch.name)
        sketch_prefix = sketch_name.split('_')[0]
        if sketch_name in TEMPLATES:  # default sketch
            sketch_file = sketchy_path / Path(f'{sketch_name}_15_1000')
        elif sketch_prefix in TEMPLATES:  # prefixed sketches
            sketch_file = sketchy_path / Path(f'{sketch_name}')
        else:
            sketch_file = sketch  # sketch path if no prefix found
    except IndexError:
        sketch_file = sketch  # sketch path if no prefix found

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
        palette=palette,
        image_format='pdf'
    )
