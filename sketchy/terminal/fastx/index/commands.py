import click

from pathlib import Path
from sketchy.utils import SketchySimulator


@click.command()
@click.option(
    "--fastx",
    "-f",
    type=Path,
    help="Path to Fast{a,q} input file to index reads from headers",
    required=True,
)
@click.option(
    "--sort_by",
    "-s",
    type=str,
    help="Sort read index by column: runid, sampleid, barcode, name, start_time [start_time]",
    default="start_time",
)
@click.option(
    "--output",
    "-o",
    type=Path,
    default="read_index.tsv",
    help="Output sorted read index to tab-delimited file [read_index.tsv]",
)
def index(fastx, output, sort_by):

    """ Create read index from information in read headers (Guppy) """

    sim = SketchySimulator(fastx=fastx)
    sim.get_run_index(fout=output, sort_by=sort_by)
