import click
import pandas

from sketchy.sketchy import SketchyEvaluation
from pathlib import Path


@click.command()
@click.option(
    '--ssh', '-s', type=Path, required=True,
    help='Path to sum of shared hashes data file from prediction',
)
@click.option(
    '--data', '-d', type=Path, required=True,
    help='Path to sketch data with genotype feature columns',
)
@click.option(
    '--feature', '-f', default=None,  type=str,
    help='Select subset of feature columns for evaluation; '
         'comma-delimited string [all]'
)
@click.option(
    '--output', '-o', default='sketchy.png', type=str,
    help='Output file for feature plots; image format by extension [sketchy.png]'
)
@click.option(
    '--color', '-c', default='YlGnBu', type=str,
    help='Color scheme for plots [YlGnBu]'
)
@click.option(
    '--stable', default=500,  type=int,
    help='Number of reads to define a stable breakpoint [500]'
)
@click.option(
    '--top', default=5,  type=int,
    help='Select top most frequent features values for plots [5]'
)
@click.option(
    '--limit', default=None, type=int,
    help='Limit the number of reads to show in plots [all]'
)
def evaluate(
    ssh, data, feature, stable, top, color, output, limit
):
    """ Evaluate sum of shared hashes with Sketchy """

    if feature is not None:
        feature = feature.split(',')

    se = SketchyEvaluation(
        ssh_data=ssh, feature_data=data, top=top, limit=limit
    )
    se.evaluate(
        features=feature, stable=stable, color=color, fout=output
    )
