import click
import pandas
from dateutil import parser as dp

from pathlib import Path
from sketchy.utils import SketchySimulator, MutuallyExclusiveOption


@click.command()
@click.option(
    "--fastx",
    "-f",
    type=Path,
    help="Path to Fast{a,q} input file used in evaluation",
    default=None,
    required=True,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["index"]
)
@click.option(
    "--index",
    "-i",
    type=Path,
    help="Path to input file read index from `sketchy utils fx-sort`",
    default=None,
    required=False,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["fastx"]
)
@click.option(
    "--evaluation",
    "-e",
    type=Path,
    help="Path to evaluation file from `sketchy evaluate`",
    default=None,
    required=True,
)
def fx_time(fastx, evaluation, index):

    """ Experimental: compute time to prediction from evaluation output and input reads """

    sim = SketchySimulator(
        fastx=fastx, fastx_index=index
    )

    if index:
        fx = sim.fastx_index
    else:
        fx = sim.get_run_index()

    feature_predictions = pandas.read_csv(
        evaluation, sep='\t', names=['feature', 'prediction', 'first', 'stable']
    )

    for i, row in feature_predictions.iterrows():
        first_delta, stable_delta = compute_feature_time_delta(fx, row)
        print(
            f"{row['feature']}\t{row['prediction']}\t{row['first']}\t"
            f"{first_delta}\t{row['stable']}\t{stable_delta}"
        )


def compute_feature_time_delta(fx, row):

    try:
        first = int(
            row['first']
        )-1
        first_feature = fx.iloc[first]
    except (TypeError, KeyError, ValueError):
        first_feature = None

    try:
        stable = int(
            row['stable']
        ) - 1
        stable_feature = fx.iloc[stable]
    except (TypeError, KeyError, ValueError):
        stable_feature = None

    first_read = fx.iloc[0]

    if first_feature is not None:
        first_feature_delta = \
            dp.parse(first_feature.start_time) - dp.parse(first_read.start_time)
    else:
        first_feature_delta = None

    if stable_feature is not None:
        stable_feature_delta = \
            dp.parse(stable_feature.start_time) - dp.parse(first_read.start_time)
    else:
        stable_feature_delta = None

    return first_feature_delta, stable_feature_delta
