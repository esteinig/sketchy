import click

from pathlib import Path
from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--fastq', '-f', help='Input fastq file to sort.', type=Path
)
@click.option(
    '--output', '-o', help='Fastq sorted by increasing datetime.', type=Path
)
@click.option(
    '--shuffle', '-s', is_flag=True, help='Shuffle start datetimes or  in fastq.'
)
@click.option(
    '--start',  is_flag=True, help='Print start datetime of first read in fastq.'
)
@click.option(
    '--sorted_out',  default=None, type=Path,
    help='Sorted tab-delimited data frame of read IDs and start time.'
)
def sort(fastq, output, shuffle, start, sorted_out):
    """ Sort basecalled reads by start time in header """

    sketchy = Sketchy()
    sketchy.sort_fastq(
        file=fastq, fastq=output, shuffle=shuffle,
        first_read=start, outfile=sorted_out
    )
