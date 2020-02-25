import click
import logging
import shutil

from pathlib import Path
from sketchy.utils import SketchySimulator, create_fastx_index, get_output_handle, PoreLogger


@click.command()
@click.option(
    "--fastx",
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
def fx_sort(fastx, output):

    """ Sort reads by 'start_time' in basecalled header from Guppy"""

    logger = PoreLogger(level=logging.INFO).logger

    logger.info('Creating index for random access to reads')
    fx, build_read = create_fastx_index(fastx=fastx)

    logger.info('Creating start time index for all reads')
    sim = SketchySimulator(fastx=fastx)
    run_index = sim.get_run_index()

    logger.info(f'Writing sorted reads by start time date')
    with get_output_handle(fpath=output) as fout:
        for i, row in run_index.iterrows():
            read_str = build_read(
                fx[row['name']], comment=sim.create_header_comment(row)
            )
            fout.write(read_str + '\n')

    logger.info(f'Completed')

    Path(str(fastx) + '.fxi').unlink()
