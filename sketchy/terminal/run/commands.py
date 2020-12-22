import os
import click

from sketchy.sketchy import SketchyWrapper
from pathlib import Path


@click.command()
@click.option(
    '--fastx',
    '-f',
    default=None,
    type=Path,
    required=True,
    help='Path to input Fastx containing basecalled nanopore reads'
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
    '--sketchy',
    type=Path,
    default=Path.home() / '.sketchy',
    help='Sketchy path to reference sketch home directory. '
         'Can be set via environmental variable: SKETCHY_PATH'
)
@click.option(
    '--quiet',
    is_flag=True,
    help="Run without logging output or progress bar."
)
@click.option(
    '--no_plot',
    is_flag=True,
    help="Do not plot the results; output only table [false]"
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help="Use this Matplotlib backend [default]"
)
@click.option(
    '--image_format',
    type=str,
    default="png",
    help="Use this image format [png]"
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
    sketchy,
    no_plot,
    mpl_backend,
    image_format,
    quiet
):

    """ Sketchy wrapper for Rust streaming algorithm on completed read sets"""

    try:
        sketchy_path = Path(os.environ['SKETCHY'])
    except KeyError:
        sketchy_path = sketchy

    if limit == -1:
        limit = None

    if 0. < stable <= 1.:
        if limit is None:
            stable = 1000
        else:
            stable = int(limit*stable)
    else:
        stable = int(stable)

    if not (sketchy_path / sketch).exists():
        raise ValueError(f"Could not find database sketch: {sketch} at {sketchy_path}")
    else:
        sketch_file

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
        plot=not no_plot,
        palette=palette,
        image_format=image_format,
        mpl_backend=mpl_backend
    )
