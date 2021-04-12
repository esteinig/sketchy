import click

from pathlib import Path
from sketchy.sketchy import SketchyDatabase


@click.command()
@click.option(
    '--outdir', '-o', type=Path, required=False, default="sketchy-db",
    help='Name and prefix of database files [db]'
)
@click.option(
    '--full', '-f', is_flag=True,
    help='Include columns with numeric formats [false]'
)
def download(outdir, full):

    """ Download default sketches """

    db = SketchyDatabase(
        sketch_file=None, genotype_file=None
    )

    db.download_databases(outdir=outdir, full=full)


