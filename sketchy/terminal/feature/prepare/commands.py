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
    '--output', '-o', type=Path, required=False, default="index.prepped.tsv",
    help='Path to dropped feature index output file.'
)
@click.option(
    '--exclude_numeric', '-e', is_flag=True,
    help='Treat integer columns as categorical.'
)
def prepare(index, exclude_numeric, output):

    """ Prepare a feature index file for evaluation in Rust """

    idx = LineageIndex(index_file=index)

    idx.write(file=Path(output.stem + ".sorted.tsv"), idx=True, header=True)

    _, index_key = idx.prepare_columns(integers=not exclude_numeric)

    idx.write(file=output, idx=False, header=False)

    with Path(output.stem + ".json").open('w') as fout:
        json.dump(index_key, fout, sort_keys=False)
