import click

from sketchy.sketchy import Sketchy

@click.command()
@click.option(
    '--fastq', '-f', help='Input FASTQ file to slice and predict'
)
@click.option(
    '--output', '-o', help='Sorted by start datetime FASTQ output file'
)
@click.option(
    '--shuffle', '-s', is_flag=True, help='Shuffle start datetimes in FASTQ'
)
def sort(fastq, output, shuffle):
    """ Sort a read file (.fq) by read start times from read header """

    sketchy = Sketchy()
    sketchy.sort_fastq(
        file=fastq, fastq=output, shuffle=shuffle
    )
