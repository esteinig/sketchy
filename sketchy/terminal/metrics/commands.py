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
    '--outfile',
    '-o',
    required=False,
    type=Path,
    default=Path("metrics.tsv"),
    help='Path to metric output file [metrics.tsv]'
)
@click.option(
    '--average',
    '-a',
    required=False,
    type=str,
    default="micro",
    help='Average method for precision and recall of multilabel features'
)
@click.option(
    '--force_db',
    '-f',
    required=False,
    type=str,
    default="",
    help='Force database [none]'
)
@click.option(
    '--read_levels',
    '-r',
    required=False,
    type=str,
    default="200,1000",
    help='Force database [none]'
)
def metrics(
    match_data, average, force_db, outfile, read_levels
):

    """ Compute {accuracy, precision, recall, F1} for matched data """
    read_levels = [int(n.strip()) for n in read_levels.split(',')]
    sd = SketchyDiagnostics(outdir=None, verbose=True, mpl_backend=None)
    sd.get_metrics(data=match_data, outfile=outfile, multi_average=average, force_db=force_db, read_levels=read_levels)