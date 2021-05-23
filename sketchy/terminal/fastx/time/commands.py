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
def time(fastq, evaluation, index, prefix, delta):

    """ Compute time of prediction from reads and evaluations """

    sim = SketchySimulator(
        fastx=fastq, fastx_index=index
    )

    if index:
        fx = sim.fastx_index
    else:
        fx = sim.get_run_index()

    fx.index.name = 'read'
    fx.sort_index().to_csv(f'{prefix}.time.tsv', sep='\t', index_label='read')


def compute_time_delta(fx, read: int = 100, delta: str = None):

    read_data = fx[fx['read'] == read]

    if delta == 'first':
        first_read = fx.iloc[0]
        readd = dp.parse(read_data.start_time) - \
            dp.parse(first_read.start_time)
    elif delta == 'null':  # nextflow setting
        readd = dp.parse(read_data.start_time)
        readd = readd.strftime("%d-%m-%Y %H:%M:%S")
    else:
        read_time = dp.parse(read_data.start_time).replace(tzinfo=None)
        start_time = dp.parse(delta, dayfirst=True).replace(tzinfo=None)
        readd = read_time - start_time


    return readd
