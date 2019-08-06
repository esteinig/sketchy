import click
import pandas

from pathlib import Path
from sketchy.utils import query_fastq, construct_fastq


@click.command()
@click.option(
    '--fpath', '-f', type=Path, help='Path to Fastq file to select from.'
)
@click.option(
    '--output', '-o', type=Path, help='Output to Fastq file.'
)
@click.option(
    '--sample', '-s', type=int or float, default=0.1,
    help='Number of reads to sample or fraction of reads to sample.'
)
@click.option(
    '--replace', '-r', is_flag=True,
    help='Sample with replacement.'
)
def fq_sample(fpath, output, sample, replace):
    """ Sample reads from fastq with or without replacement (sketchy.nf) """

    df, records = query_fastq(fpath, full=True)

    if isinstance(sample, int):
        df = df.sample(n=sample, replace=replace)
    else:
        df = df.sample(frac=sample, replace=replace)

    construct_fastq(output, headers=df['header'], records=records)







