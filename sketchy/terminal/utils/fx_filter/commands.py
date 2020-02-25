import click
import pandas

from pysam import FastxFile
from pathlib import Path
from sketchy.utils import get_output_handle


@click.command()
@click.option(
    "--fpath",
    "-f",
    type=Path,
    help="Path to Fast{a,q} input file.",
    required=True,
)
@click.option(
    "--output",
    "-o",
    type=str,
    help="Output to Fast{a,q} file. Default stdout [-]",
    default="-",
)
@click.option(
    "--ids",
    "-i",
    type=Path,
    help=(
        "Path to file containing the read IDs to get from Fast{a,q}. "
    ),
    required=True,
)
@click.option(
    "--column",
    "-c",
    default=1,
    help="Column index that contains the IDs (0-based). [1]",
    type=int,
)
@click.option(
    "--sep",
    "-s",
    default="\t",
    help="File separator to read columns. ['\\t']",
    type=str,
)
def fx_filter(fpath, ids, output, column, sep):
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
