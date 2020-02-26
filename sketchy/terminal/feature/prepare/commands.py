import click
import json

from pathlib import Path
from sketchy.sketchy import LineageIndex


@click.command()
@click.option(
    '--index', '-i', type=Path, required=True,
    help='Path to feature index input file.'
)
@click.option(
    '--drop', '-d', type=str, required=False, default=None,
    help='Comma separated string of column names to drop.'
)
@click.option(
    '--output', '-o', type=Path, required=False, default="index.prepped.tsv",
    help='Path to dropped feature index output file.'
)
def prepare(index, drop, output):

    """ Prepare a feature index file for evaluation in Rust """

    idx = LineageIndex(index_file=index)

    idx.write(file=Path(output.stem + ".sorted.tsv"), idx=True, header=True)

    if drop is not None:
        if ',' in drop:
            drop = drop.split(',')
        else:
            drop = [drop]

        idx.index = idx.index.drop(columns=drop)

    _, index_key = idx.prepare_columns(integers=True)

    idx.write(file=output, idx=False, header=False)

    with Path(output.stem + ".json").open('w') as fout:
        json.dump(index_key, fout, sort_keys=False)
