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
@click.option(
    '--outlist', '-l', help='List output of read headers, lengths and date times', type=Path
)
def fq_sort(fastq, output, shuffle, start, outlist):
    """ Sort or shuffle reads by start date time in header. """

    if output:
        df, records = query_fastq(fpath=fastq, full=True)
    else:
        df = query_fastq(fpath=fastq)

    if start:
        print(
            str(df.loc[0, 'date'])
        )
        exit(0)

    if shuffle:
        df = df.sample(frac=1)

    headers = df['header'].tolist()

    if output:
        construct_fastq(
            output=output, headers=headers, records=records
        )

    if outlist:
        df.to_csv(outlist, sep='\t', index=False)