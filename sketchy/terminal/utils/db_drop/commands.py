import click
import os
import pandas
import logging

from pathlib import Path
from sketchy.sketchy import LineageIndex


@click.command()
@click.option(
    '--index', '-i', type=Path, required=True,
    help='Path to create index input file.'
)
@click.option(
    '--output', '-o', type=Path, required=False, default="index.dropped.tsv",
    help='Path to dropped create index output file [index.dropped.tsv]'
)
@click.option(
    '--columns', '-c', type=str, default=None,
    help='Comma-delimited string of columns to db_drop or "clean" [clean]'
)
def db_drop(index, columns, output):

    """ Drop columns from the create index file """

    li = LineageIndex(index_file=index)

    if isinstance(columns, str):
        columns = columns.split(',')

    df = li.drop_columns(columns=columns)

    df.to_csv(output, sep='\t', index=True, header=True)


