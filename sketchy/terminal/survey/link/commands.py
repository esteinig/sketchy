import click
import os
import pandas
import logging
import random

from pathlib import Path
from sketchy.utils import PoreLogger


@click.command()
@click.option(
    '--iid', '-i', type=Path,
    help='Path to isolate ID file from Pathfinder Survey'
)
@click.option(
    '--directory', '-d', type=Path,
    help='Path to directory from which to extract files'
)
@click.option(
    '--column', '-c', default=None, type=str,
    help='Use a header and column name to parse isolate IDs [None]'
)
@click.option(
    '--extension', '-e', type=str, default='.fasta',
    help='File extension to link isolate IDs to file [.fasta]'
)
@click.option(
    '--symlink', '-s', is_flag=True,
    help='Symlink the detected files to the output directory.'
)
@click.option(
    '--outdir', '-o', type=Path, default="sketchy_link",
    help='If symlink, output directory for symbolic links to files [sketchy_link]'
)
@click.option(
    '--bootstrap', '-b', type=int, default=None,
    help='Bootstrap sample isolate IDs'
)
def link(iid: Path, directory, column, extension, outdir, symlink, bootstrap):

    """ Link ID file to FASTA, e.g. from filtered Pathfinder Survey """

    if column is None:
        iids = pandas.read_csv(
            iid, index_col=0, header=0, dtype=str
        )
    else:
        iids = pandas.read_csv(
            iid, usecols=[column], header=0, sep='\t'
        )
        iids.rename(columns={column: 'iid'}, inplace=True)

    if symlink:
        outdir.mkdir(exist_ok=True, parents=True)

    log = PoreLogger(level=logging.INFO).logger

    if bootstrap:
        isolates = random.choices(iids.iid, k=bootstrap)
    else:
        isolates = iids.iid

    for i in isolates:
        iid_path = directory / str(i + extension)
        if iid_path.exists():
            if symlink:
                sym_path = (outdir / str(i + extension)).absolute()
                log.info(f'Symlink: {iid_path.absolute()} to {sym_path.absolute()}')
                os.symlink(
                    str(iid_path.absolute()), str(sym_path.absolute())
                )
            else:
                print(f"{iid_path}")
        else:
            log.debug(f'Could not find: {iid_path}')


