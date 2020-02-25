import click

from pathlib import Path
from sketchy.utils import SketchySimulator


@click.command()
@click.option(
    '--fastx', '-f', default=None, type=Path, required=True,
    help='Path to Fast{a,q} file to simulate run [required]'
)
@click.option(
    '--index', '-i', default=None, type=Path, required=False,
    help='Path to Fast{a,q} index file from previous simulation [none]'
)
@click.option(
    '--outdir', '-o', default='run_sim', type=Path, required=False,
    help='Output directory for run simulation [run_sim]'
)
@click.option(
    '--reads_per_file', '-r', default=100, type=int, required=False,
    help='Number of reads per Fast{a,q} to simulate live basecalling [4000]'
)
@click.option(
    '--barcodes', '-b', default=None, type=str, required=False,
    help='Barcode integers, comma-separated e.g. 1,2,3,4 [None]'
)
@click.option(
    '--speedup', '-s', default=1.0, type=float, required=False,
    help='Speed up the simulation by this factor [1.0]'
)
@click.option(
    '--start_time_regex', default=None, type=str, required=False,
    help='Read start time regex to parse from read headers.'
)
@click.option(
    '--barcode_regex', default=None, type=str, required=False,
    help='Barcode regex to parse from read headers.'
)
@click.option(
    '--test', is_flag=None, help='Run simple test mode; see docs [false]'
)
def simulate(
    fastx, index, outdir, reads_per_file, barcodes, speedup, start_time_regex, barcode_regex, test
):

    """ Experimental: simulate a live sequencing run with optional barcodes """

    # Remember to check if there are multiple runs in same file and
    # concat them with a standard gap between them (rather than real time)

    sim = SketchySimulator(
        fastx=fastx,
        fastx_index=index
    )
