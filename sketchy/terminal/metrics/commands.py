import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--match_data',
    '-m',
    required=True,
    type=Path,
    help='Path to match data from task'
)
def metrics(
    match_data
):

    """ Compute {accuracy, precision, recall, F1} for matched data """

    sd = SketchyDiagnostics(outdir=None, verbose=True, mpl_backend=None)
    sd.get_metrics(data=match_data)