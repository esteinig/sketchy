import click
from numpy import where

from pathlib import Path
from sketchy.sketchy import Evaluation


@click.command()
@click.option(
    '--index',
    '-i',
    type=Path,
    required=True,
    help='Path to genotype feature index used in: sketchy evaluate',
)
@click.option(
    '--key',
    '-k',
    type=Path,
    required=True,
    help='Path to genotype feature index key for translation from numeric.',
)
@click.option(
    '--sssh',
    '-ss',
    type=Path,
    required=True,
    help='Path to sum of ranked sums shared hashes data file from evaluation',
)
@click.option(
    '--ssh',
    '-s',
    type=Path,
    required=False,
    default=None,
    help='Path to sum of shared hashes data file from prediction',
)
@click.option(
    '--stable',
    '-st',
    type=int,
    required=False,
    default=None,
    help='Stability parameter passed to: sketchy evaluate',
)
@click.option(
    '--palette',
    '-c',
    type=str,
    default='YlGnBu',
    help='Color palette for output plots [YlGnBu]'
)
@click.option(
    '--prefix',
    '-p',
    type=str,
    default='sketchy',
    help='Output prefix for all files [sketchy]'
)
@click.option(
    '--format',
    '-f',
    type=str,
    default='png',
    help='Output image format [png]'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose logging output [false]'
)
def plot(sssh, ssh, index, key, stable, color, prefix, format, verbose):

    """ Plot output from Sketchy Rust pipeline """

    eve = Evaluation(
        sssh=sssh,
        index=index,
        key=key,
        stable=stable,
        ssh=ssh,
        verbose=verbose
    )

    eve.plot_feature_evaluations(
        plot_file=Path(prefix + f'.{format}'),
        break_file=Path(prefix + f'.data.tsv'),
        color=color,
        stable_point=False,
        break_point=True
    )
