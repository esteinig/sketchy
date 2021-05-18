import click
import pandas

from pysam import FastxFile
from pathlib import Path
from sketchy.utils import get_output_handle


@click.command()
@click.option(
    "--fasta",
    "-f",
    type=Path,
    help="Path to fasta directory to bootstrap sample",
    required=True,
)
@click.option(
    "--fasta",
    "-f",
    type=Path,
    help="Path to fasta directory to bootstrap sample",
    required=True,
)
def filter(fpath, ids, output, column, sep):
    """ Filter reads by external file of read header names """

    ids_df = pandas.read_csv(ids, sep=sep, header=None)

    read_ids = set(
        [str(read_id) for read_id in ids_df.iloc[:, column].tolist()]
    )

    with FastxFile(fpath) as fin, \
            get_output_handle(output) as fout:
        for read in fin:
            if read.name in read_ids:
                fout.write(str(read)+"\n")
