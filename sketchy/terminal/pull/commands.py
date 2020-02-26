import click
from pathlib import Path
from sketchy.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '--p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [~/.sketchy ]'
)
@click.option(
    '--full_sketch', '-f', is_flag=True,
    help='Pull the full default sketch collections [false]'
)
def pull(path, full_sketch):

    """ Pull default sketch collection into local storage """

    path.mkdir(exist_ok=True, parents=True)
    GoogleCloudSketch(sketch_path=path, full_sketch=full_sketch).pull()
