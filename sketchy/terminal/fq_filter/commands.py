import click
import pandas
import pysam
import sys

from pathlib import Path
from sketchy.utils import filter_fastq


@click.command()
@click.option(
    "--fpath",
    "-f",
    type=Path,
    help="Path to Fast{a,q} file to select from.",
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
        "Path to file containing the read IDs to get from fast{a,q}. "
        "A header line is expected."
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
def fq_filter(fpath: Path, ids: Path, output: str, column: int, sep: str):
    """ Filter reads by headers; used for species grep (sketchy.nf) """
    if output == "-":
        output = sys.stdout
    else:
        p = Path(output)
        if not p.parent.is_dir():
            raise NotADirectoryError(
                "Directory specified for output file does not exist: {}".format(
                    p.parent
                )
            )
        output = p.open("w")

    ids_df = pandas.read_csv(ids, sep=sep)
    read_ids = set([str(read_id) for read_id in ids_df.iloc[:, column].tolist()])

    with pysam.FastxFile(fpath) as fastq_in:
        filter_fastq(fastq_in=fastq_in, fastq_out=output, records=read_ids)
    output.close()
