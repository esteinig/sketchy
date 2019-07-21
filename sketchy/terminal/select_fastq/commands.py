import click
import pandas

from pathlib import Path
from sketchy.minhash import MashSketch


@click.command()
@click.option(
    '--fpath', '-f', type=Path, help='Path to Fastq file to select from.'
)
@click.option(
    '--output', '-o', type=Path, help='Output to Fastq file.'
)
@click.option(
    '--ids', '-i', type=Path,
    help='Path to tab-delimited file containing the read '
         'IDs to select in second column.'
)
def select_fastq(fpath, ids, output):
    """ Select reads by ID from output of Kraken2 grep for species """

    import sys

    ids_df = pandas.read_csv(ids, sep='\t')

    read_ids = [str(_) for _ in ids_df.iloc[:, 1].tolist()]

    def process(lines=None):
        ks = ['name', 'sequence', 'optional', 'quality']
        return {k: v for k, v in zip(ks, lines)}

    n = 4
    with output.open('w') as outfq:
        with fpath.open() as fh:
            lines = []
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == n:
                    record = process(lines)
                    read_id = str(
                        record['name'].split()[0].replace('@', '')
                    )
                    if read_id in read_ids:
                        outfq.write(
                            f"{record['name']}\n{record['sequence']}\n"
                            f"{record['optional']}\n{record['quality']}\n"
                        )
                    lines = []








