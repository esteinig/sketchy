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
    "--read",
    "-r",
    type=int,
    help="Print time delta at read [100]",
    default=100,
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
    default='first',
    required=False,
)
def time(fastq, read, index, prefix, delta):

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

    compute_time_delta(fx, read=read, delta=delta)

def compute_time_delta(fx, read: int = 100, delta: str = None):

    read_time = fx.loc[fx['read'] == read, 'start_time']
    print(read_time)
    if delta == 'first':
        first_read = fx.iloc[0]
        readd = dp.parse(read_time) - \
            dp.parse(first_read.start_time)
    elif delta == 'null':  # nextflow setting
        readd = dp.parse(read_time)
        readd = readd.strftime("%d-%m-%Y %H:%M:%S")
    else:
        read_time = dp.parse(read_time).replace(tzinfo=None)
        start_time = dp.parse(delta, dayfirst=True).replace(tzinfo=None)
        readd = read_time - start_time


    return readd
