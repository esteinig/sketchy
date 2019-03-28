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
    '--shuffle', '-s', is_flag=True, help='Shuffle start datetimes or  in FASTQ'
)
@click.option(
    '--bootstraps', '-b', default=None, type=int, help='Number of bootstrap samples.'
)
@click.option(
    '--replacement', '-r', is_flag=True, help='Replacement sampling for bootstraps.'
)
@click.option(
    '--prefix', '-p', default='boot_', type=int, help='Prefix for bootstrap output files.'
)
@click.option(
    '--sample_reads', default=100, type=int, help='Sample size of reads for each bootstrap replicate.'
)
def sort(fastq, output, shuffle, bootstraps, replacement, prefix, sample_reads):
    """ Sort basecalled reads by start time (Albacore) """

    sketchy = Sketchy()
    sketchy.sort_fastq(
        file=fastq, fastq=output, shuffle=shuffle, nbootstrap=bootstraps,
        replacement=replacement, prefix=prefix, sample_size=sample_reads
    )
