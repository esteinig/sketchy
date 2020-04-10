import click
import json

from pathlib import Path
from sketchy.sketchy import LineageIndex


@click.command()
@click.option(
    '--index', '-i', type=Path, required=True,
    help='Path to feature index input file'
)
@click.option(
    '--drop', '-d', type=str, required=False, default=None,
    help='Comma separated string of columns to drop'
)
@click.option(
    '--prefix', '-p', type=Path, required=False, default="index",
    help='Prefix for prepared feature index output files'
)
def prepare(index, drop, prefix):

    """ Prepare a feature index file for evaluation in Rust """

    idx = LineageIndex(index_file=index)

    idx.write(file=Path(f"{prefix}.reference.tsv"), idx=False, header=True)

    if drop is not None:
        if ',' in drop:
            drop = drop.split(',')
        else:
            drop = [drop]

        idx.index = idx.index.drop(columns=drop)

    _, index_key = idx.prepare_columns(integers=True)

    idx.write(file=Path(f"{prefix}.tsv"), idx=False, header=False)

    with Path(f"{prefix}.json").open('w') as fout:
        json.dump(index_key, fout, sort_keys=False)
