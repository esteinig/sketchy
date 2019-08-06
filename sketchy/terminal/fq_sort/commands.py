import click

from pathlib import Path
from sketchy.utils import query_fastq, construct_fastq


@click.command()
@click.option(
    '--fastq', '-f', help='Input fastq file to sort.', type=Path
)
@click.option(
    '--output', '-o', help='Fastq sorted by increasing datetime.', type=Path
)
@click.option(
    '--shuffle', '-s', is_flag=True, help='Shuffle start datetimes in fastq.'
)
@click.option(
    '--start',  is_flag=True, help='Print start datetime of first read in fastq.'
)
def fq_sort(fastq, output, shuffle, start):
    """ Sort or shuffle reads by start date time in header. """

    df, records = query_fastq(fpath=fastq, full=True)

    if start:
        print(str(df.loc[0, 'date']))
        exit(0)

    if shuffle:
        df = df.sample(frac=1)

    headers = df['header'].tolist()

    construct_fastq(
        output=output, headers=headers, records=records
    )
