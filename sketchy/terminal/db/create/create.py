import click
import json

from pathlib import Path
from sketchy.sketchy import SketchyDatabase


@click.command()
@click.option(
    '--sketch', '-s', type=Path, required=True,
    help='Path to reference sketch file'
)
@click.option(
    '--genotypes', '-g', type=Path, required=True,
    help='Path to sketch genotypes file'
)
@click.option(
    '--outdir', '-o', type=Path, required=False, default="sketchy-db",
    help='Name and prefix of database files [db]'
)
@click.option(
    '--id_column', '-i', type=str, required=False, default="uuid",
    help='Column name containing genome identifiers matching file stems used in sketch construction [uuid]'
)
@click.option(
    '--drop', '-d', type=str, required=False, default=None,
    help='Comma separated string of column names to drop [none]'
)
@click.option(
    '--numeric', '-n', is_flag=True,
    help='Include columns with numeric formats [false]'
)
def create(sketch, genotypes, outdir, id_column, drop, numeric):

    """ Create a reference database for Sketchy """

    db = SketchyDatabase(
        sketch=sketch, genotypes=genotypes
    )

    db.create_database(id_column=id_column, outdir=outdir, drop=drop, numeric=False)


