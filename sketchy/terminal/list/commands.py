import os
import click

from pathlib import Path
from sketchy.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '-p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [~/.sketchy]'
)
def list(path):

    """ List currently available database sketches in repository """

    try:
        sketchy_path = Path(
            os.environ['SKETCHY_PATH']
        )
    except KeyError:
        sketchy_path = path

    gcs = GoogleCloudSketch(sketch_path=sketchy_path)
    gcs.list_sketches()

