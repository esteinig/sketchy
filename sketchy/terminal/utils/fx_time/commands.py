import click
import pandas
from dateutil import parser as dp

from pathlib import Path
from sketchy.utils import SketchySimulator, MutuallyExclusiveOption


@click.command()
@click.option(
    "--fastq",
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
    help="Path to evaluation file containing predictions (data.tsv)",
    default=None,
    required=False,
)
@click.option(
    "--prefix",
    "-p",
    type=str,
    help="Output prefix for time data: {prefix}.time.tsv [skecthy]",
    default='sketchy',
    required=False,
)
@click.option(
    "--delta",
    "-d",
    type=str,
    help="Compute time delta between 'first' read or start time of run "
         "!! GMT !! in format: '20/11/20 16:20:00' [first]",
    default=None,
    required=False,
)
def fx_time(fastq, evaluation, index, prefix, delta):

    """ Compute time of prediction from reads and evaluations """

    sim = SketchySimulator(
        fastx=fastq, fastx_index=index
    )

    if index:
        fx = sim.fastx_index
    else:
        fx = sim.get_run_index()

    fx.sort_index().to_csv(f'{prefix}.time.tsv', sep='\t', index_label='read')

    if evaluation is not None:
        feature_predictions = pandas.read_csv(
            evaluation, sep='\t', header=0
        )

        dates = [
            compute_feature_time_delta(fx, row, delta=delta)
            for _, row in feature_predictions.iterrows()
        ]

        feature_predictions['time'] = dates

        feature_predictions.to_csv(
            f'{prefix}.time.data.tsv', sep='\t', index=False
        )


def compute_feature_time_delta(fx, row, delta: str = None):

    try:
        stable = int(
            row['stability']
        )
        if stable == -1:
            stable_feature = None
        else:
            stable_feature = fx.iloc[stable]
    except (TypeError, KeyError, ValueError):
        stable_feature = None

    if stable_feature is not None:
        if delta:
            if delta == 'first':
                first_read = fx.iloc[0]
                stable_feature_date = dp.parse(stable_feature.start_time) - \
                    dp.parse(first_read.start_time)
            else:
                read_time = dp.parse(stable_feature.start_time).replace(tzinfo=None)
                start_time = dp.parse(delta, dayfirst=True).replace(tzinfo=None)
                stable_feature_date = read_time - start_time
        else:
            stable_feature_data = dp.parse(stable_feature.start_time)
            stable_feature_date = stable_feature_data.strftime("%d-%m-%Y %H:%M:%S")
    else:
        stable_feature_date = None

    return stable_feature_date
