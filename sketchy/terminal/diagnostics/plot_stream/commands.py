import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--sssh',
    '-s',
    type=Path,
    required=True,
    help='Path to sum of ranked sums shared hashes data file from evaluation',
)
@click.option(
    '--stable',
    '-st',
    type=int,
    required=False,
    default=100,
    help='Stability parameter passed to: sketchy stream to compute stable breakpoint for each feature [none]'
)
@click.option(
    '--color',
    '-c',
    type=str,
    default='YlGnBu',
    help='Color palette for output plots [YlGnBu]'
)
@click.option(
    '--max_ranks',
    '-m',
    type=int,
    default=5,
    help='In the SSSH plots, show max feature values / prediction ranks [5]'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose logging output [false]'
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help='Matplotlib backend [default]'
)
def plot_stream(sssh, stable, max_ranks, color, mpl_backend, verbose):

    """ Diagnostic plots for output from stream client """

    sd = SketchyDiagnostics()

    sd.process_sssh(sssh_file=sssh, stable=stable, max_ranks=max_ranks)

