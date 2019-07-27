import click

from pathlib import Path
from sketchy.cloud.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '--p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [ ~/.sketchy ]'
)
def db_list(path):
    """ List currently available database sketches """

    gcs = GoogleCloudSketch(sketch_path=path / 'db')
    gcs.list_dbs()

