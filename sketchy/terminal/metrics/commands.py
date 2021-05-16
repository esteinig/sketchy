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
@click.option(
    '--average',
    '-a',
    required=False,
    type=str,
    default="micro",
    help='Average method for precision and recall of multilabel features'
)
def metrics(
    match_data, average
):

    """ Compute {accuracy, precision, recall, F1} for matched data """

    sd = SketchyDiagnostics(outdir=None, verbose=True, mpl_backend=None)
    sd.get_metrics(data=match_data, multi_average=average)