import click
import pandas

from pathlib import Path
from sketchy.utils import filter_fastq

@click.command()
@click.option(
    '--fpath', '-f', type=Path, help='Path to Fastq file to select from.'
)
@click.option(
    '--output', '-o', type=Path, help='Output to Fastq file.'
)
@click.option(
    '--ids', '-i', type=Path,
    help='Path to file containing the read '
         'IDs to select in second column.'
)
@click.option(
    '--column', '-c', default=1,
    help='Column index that contains the IDs.'
)
@click.option(
    '--sep', '-s', default='\t',
    help='File separator to read columns.'
)
def fq_filter(fpath, ids, output, column, sep):
    """ Filter reads by headers; used for species grep (sketchy.nf) """

    ids_df = pandas.read_csv(ids, sep=sep)
    read_ids = [str(_) for _ in ids_df.iloc[:, column].tolist()]

    filter_fastq(fpath=fpath, output=output, records=read_ids)










